require 'rbbt'
require 'rbbt/util/misc'
require 'rbbt/workflow'

require 'rbbt/sources/organism'
require 'rbbt/sources/uniprot'
require 'structure/ssw'
require 'structure/interactome_3d'
require 'structure/pdb_helper'
require 'structure/uniprot'
require 'structure/appris'
require 'structure/neighbours'
require 'structure/COSMIC'

Workflow.require_workflow 'Genomics'
Workflow.require_workflow 'Translation'
Workflow.require_workflow 'PdbTools'

require 'rbbt/entity/mutated_isoform'

module Structure
  extend Workflow

  def self.alignment_map(alignment_source, alignment_target)
    map = {}

    offset_source, alignment_source = alignment_source.match(/^(_*)(.*)/).values_at( 1, 2)
    offset_target, alignment_target = alignment_target.match(/^(_*)(.*)/).values_at( 1, 2)
 
    gaps_source = 0 
    gaps_target = 0
    alignment_source.chars.zip(alignment_target.chars).each_with_index do |p,i|
      char_source, char_target = p
      gaps_source += 1 if char_source == '-'
      gaps_target += 1 if char_target == '-'
      map[i + offset_source.length - gaps_source] = i + offset_target.length - gaps_target if char_source == char_source  and char_source != "-"
    end
    
    map
  end

  def self.sequence_map(source_sequence, target_sequence)
    source_alignment, target_alignment = SmithWaterman.align(source_sequence, target_sequence)
    Structure.alignment_map(source_alignment, target_alignment)
  end

  def self.match_position(protein_position, protein_alignment, chain_alignment)
    if Array === protein_position
      alignment_map(protein_alignment, chain_alignment).chunked_values_at protein_position.collect{|p| p.to_i}
    else
      alignment_map(protein_alignment, chain_alignment)[protein_position]
    end
  end

  def self.mutated_isoforms_to_residue_list(mutated_isoforms)
    MutatedIsoform.setup(mutated_isoforms, "Hsa")
    residues = TSV.setup({}, :key_field => "Ensembl Protein ID", :fields => ["Residues"], :type => :flat, :cast => :to_i, :namespace => "Hsa")

    MutatedIsoform.setup(mutated_isoforms, "Hsa")
    mutated_isoforms.each do |mi|
      next unless mi.consequence == "MISS-SENSE"
      protein = mi.protein
      position = mi.position
      residues[protein] ||=[]
      residues[protein] << position
    end

    residues
  end

  input :sequence, :text, "Protein sequence"
  input :positions, :array, "Positions within protein sequence"
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  task :sequence_position_in_pdb => :yaml do |protein_sequence, protein_positions, pdb, pdbfile|
    atoms = PDBHelper.atoms(pdb, pdbfile)

    chains = {}
    atoms.split("\n").each do |line|
      chain = line[20..21].strip
      aapos = line[22..25].to_i
      aa    = line[17..19]

      next if aapos < 0

      chains[chain] ||= Array.new
      chains[chain][aapos] = aa
    end

    alignments = {}
    chains.each do |chain,chain_sequence|
      log pdb, "Pdb #{ pdb}, chain #{ chain }."

      chain_sequence = chain_sequence.collect{|aa| aa.nil? ? '?' : Misc::THREE_TO_ONE_AA_CODE[aa.downcase]} * ""

      chain_alignment, protein_alignment = SmithWaterman.align(chain_sequence, protein_sequence)

      alignments[chain] = Structure.match_position(protein_positions, protein_alignment, chain_alignment)
    end

    alignments.delete_if{|c,p| p.nil? or p.empty?}

    alignments
  end
  export_exec :sequence_position_in_pdb

  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  input :chain, :string, "PDB chain"
  input :positions, :array, "Position within PDB chain"
  input :sequence, :text, "Protein sequence"
  task :pdb_chain_position_in_sequence => :array do |pdb, pdbfile, chain, positions, protein_sequence|
    atoms = PDBHelper.atoms(pdb, pdbfile)

    chains = {}
    atoms.split("\n").each do |line|
      pdb_chain = line[20..21].strip
      aapos = line[22..25].to_i
      aa    = line[17..19]

      next if aapos < 0

      chains[pdb_chain] ||= Array.new
      chains[pdb_chain][aapos] = aa
    end

    chain_sequence = chains[chain].collect{|aa| aa.nil? ? '?' : Misc::THREE_TO_ONE_AA_CODE[aa.downcase]} * ""
    protein_alignment, chain_alignment = SmithWaterman.align(protein_sequence, chain_sequence)

    Structure.match_position(positions, chain_alignment, protein_alignment)
  end
  export_exec :pdb_chain_position_in_sequence


  dep :sequence_position_in_pdb
  input :sequence, :text, "Protein sequence"
  input :position, :integer, "Position inside sequence"
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  input :distance, :float, "Distance"
  task :amino_acid_neighbours_in_pdb => :array do |sequence,position,pdb,pdbfile, distance|
    alignments = step(:sequence_position_in_pdb).load.collect{|k,pos| [k,pos] * ":"}

    pdbfile ||= Open.read("http://www.pdb.org/pdb/files/#{ pdb }.pdb.gz")

    neighbours = PdbTools.job(:pdb_close_positions, name, :pdb => pdb, :pdbfile => pdbfile, :distance => distance).run

    positions_in_chains = neighbours.values_at(*alignments).flatten.uniq

    positions_in_chains.collect do |p| 
      chain, pos = p.split(":")
      Structure.job(:pdb_chain_position_in_sequence, name, :sequence => sequence, :pdb => pdb, :pdbfile => pdbfile, :chain =>  chain, :position => pos.to_i).run
    end.compact
  end
  export_asynchronous :amino_acid_neighbours_in_pdb

  input :sequence, :text, "Protein sequence"
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  task :alignment_map => :tsv do |protein_sequence, pdb, pdbfile|
    atoms = PDBHelper.atoms(pdb, pdbfile)

    chains = {}
    atoms.split("\n").each do |line|
      chain = line[20..21].strip
      aapos = line[22..25].to_i
      aa    = line[17..19]
      
      next if aapos < 0

      chains[chain] ||= Array.new
      chains[chain][aapos] = aa
    end

    result = TSV.setup({}, :key_field => "Sequence position", :fields => ["Chain:Position in PDB"], :type => :flat)
    chains.each do |chain,chain_sequence|
      chain_sequence = chain_sequence.collect{|aa| aa.nil? ? '?' : Misc::THREE_TO_ONE_AA_CODE[aa.downcase]} * ""

      chain_alignment, protein_alignment = SmithWaterman.align(chain_sequence, protein_sequence)

      map = Structure.alignment_map(protein_alignment, chain_alignment)
      map.each do |seq_pos, chain_pos|
        if result[seq_pos].nil?
          result[seq_pos] = [[chain, chain_pos] * ":"]
        else
          result[seq_pos] << [chain, chain_pos] * ":"
        end
      end
    end

    result
  end
  export_exec :alignment_map

  input :uniprot, :string, "UniPro/SwissProt Accession"
  task :i3d_protein_pdbs => :array do |uniprot|
    begin
      CMD.cmd("cat '#{Interactome3d.proteins_tsv.produce.find}' | grep '#{ uniprot}'").read.split("\n").collect{|l| l.split("\t").last.split("|")}.flatten
    rescue
      Log.debug("Error in grep: #{$!.message}")
      []
    end
  end

  input :uniprot, :string, "UniPro/SwissProt Accession"
  task :i3d_interaction_pdbs => :array do |uniprot|
    begin
      CMD.cmd("cat '#{Interactome3d.interactions_tsv.produce.find}' | grep '#{ uniprot}'").read.split("\n").collect{|l| l.split("\t").last.split("|")}.flatten
    rescue
      Log.debug("Error in grep: #{$!.message}")
      []
    end.select{|file| file.index uniprot}
  end

  input :residues, :tsv, "Proteins and their affected residues", nil
  task :annotate_residues_UNIPROT => :tsv do |residues|
    tsv = TSV.setup({}, :key_field => "Ensembl Protein ID", :fields => ["Residue", "UniProt Features", "UniProt Feature locations", "UniProt Feature Descriptions"], :type => :double)

    iso2uni = Organism.protein_identifiers("Hsa").index :target => "UniProt/SwissProt Accession", :persist => true
    iso2sequence = Organism.protein_sequence("Hsa").tsv :type => :single, :persist => true


    residues.each do |isoform, list|
      uniprot = iso2uni[isoform]
      next if uniprot.nil?

      features = Structure.corrected_uniprot_features(uniprot, iso2sequence[isoform])
      overlapping = [[],[],[],[]]
      list.flatten.each do |position|
        position = position.to_i
        features.select{|info|
          case info[:type]
          when "VAR_SEQ", "CONFLICT", "CHAIN", "UNSURE"
            false
          when "DISULFID", "CROSSLNK", "VARIANT"
            info[:start] == position or info[:end] == position
          else
            info[:start] <= position and info[:end] >= position
          end
        }.each{|info|
          overlapping[0] << position
          overlapping[1] << info[:type]
          overlapping[2] << [info[:start], info[:end]] * ":"
          overlapping[3] << (info[:description] || "").strip.sub(/\.$/,'')
        }
      end

      tsv[isoform] = overlapping
    end

    tsv
  end
  export_asynchronous :annotate_residues_UNIPROT

  input :residues, :tsv, "Proteins and their affected residues", nil
  task :annotate_residues_Appris => :tsv do |residues|
    tsv = TSV.setup({}, :key_field => "Ensembl Protein ID", :fields => ["Residue", "Appris Features", "Appris Feature locations", "Appris Feature Descriptions"], :type => :double)

    iso2sequence = Organism.protein_sequence("Hsa").tsv :type => :single, :persist => true

    residues.each do |isoform, list|

      features = Structure.appris_features(isoform)

      overlapping = [[],[],[],[]]
      list.flatten.each do |position|
        position = position.to_i
        features.select{|info|
          info[:start] <= position and info[:end] >= position
        }.each{|info|
          overlapping[0] << position
          overlapping[1] << info[:type]
          overlapping[2] << [info[:start], info[:end]] * ":"
          overlapping[3] << (info[:description] || "").strip.sub(/\.$/,'')
        }
      end

      tsv[isoform] = overlapping
    end

    tsv
  end
  export_asynchronous :annotate_residues_Appris

  input :residues, :tsv, "Proteins and their affected residues", nil
  task :annotate_variants_COSMIC => :tsv do |residues|

    cosmic_residue_mutations = Structure.COSMIC_residues
    
    isoform_matched_variants = {}
    residues.each do |protein, positions|
      positions.flatten.each do |position|
        matching_mutations = cosmic_residue_mutations[[protein, position]*":"]
        next if matching_mutations.nil? or matching_mutations.empty?
        isoform_matched_variants[protein] ||= []
        isoform_matched_variants[protein].concat matching_mutations
      end
    end

    cosmic_mutation_annotations = Structure.COSMIC_mutation_annotations

    res = TSV.setup({}, :key_field => "Ensembl Protein ID", :fields => ["Genomic Mutation"].concat(cosmic_mutation_annotations.fields), :type => :double)

    isoform_matched_variants.each do |protein, mutations|
      values = []
      mutations.each do |mutation|
        values << [mutation].concat(cosmic_mutation_annotations[mutation])
      end

      res[protein] = Misc.zip_fields(values)
    end

    res
  end
  export_asynchronous :annotate_variants_COSMIC

  input :residues, :tsv, "Proteins and their affected residues", nil
  task :annotate_variants_UNIPROT => :tsv do |residues|

    uniprot_residue_mutations = Structure.UniProt_residues
    
    isoform_matched_variants = {}
    residues.each do |protein, positions|
      positions.flatten.each do |position|
        matching_mutations = uniprot_residue_mutations[[protein, position]*":"]
        next if matching_mutations.nil? or matching_mutations.empty?
        isoform_matched_variants[protein] ||= []
        isoform_matched_variants[protein].concat matching_mutations
      end
    end

    uniprot_mutation_annotations = Structure.UniProt_mutation_annotations

    res = TSV.setup({}, :key_field => "Ensembl Protein ID", :fields => uniprot_mutation_annotations.fields, :type => :double)

    isoform_matched_variants.each do |protein, mutations|
      values = []
      mutations.each do |mutation|
        values << uniprot_mutation_annotations[mutation]
      end

      res[protein] = Misc.zip_fields(values)
    end

    res
  end
  export_asynchronous :annotate_residues_UNIPROT

  input :mutated_isoforms, :array, "e.g. ENSP0000001:A12V", nil
  input :genomic_mutations, :array, "e.g. 1:173853127:T", nil
  input :organism, :string, "Organism code", "Hsa/jan2013"
  task :mutated_isoform_annotation => :string do |mis,muts,organism|
    if mis.nil? 
      Workflow.require_workflow "Sequence"
      raise ParameterException, "No mutated_isoforms or genomic_mutations specified" if muts.nil? 
      job = Sequence.job(:mutated_isoforms_for_genomic_mutations, name, :mutations => muts, :organism => organism)
      tsv = job.run
      mis = tsv.values.compact.flatten
    end

    residues = Structure.mutated_isoforms_to_residue_list(mis)

    Open.write(file(:uniprot), Structure.job(:annotate_residues_UNIPROT, name, :residues => residues).run.to_s)
    Open.write(file(:appris), Structure.job(:annotate_residues_Appris, name, :residues => residues).run.to_s)
    Open.write(file(:COSMIC), Structure.job(:annotate_variants_COSMIC, name, :residues => residues).run.to_s)
    Open.write(file(:uniprot_variants), Structure.job(:annotate_variants_UNIPROT, name, :residues => residues).run.to_s)

    "DONE"
  end
  export_asynchronous :mutated_isoform_annotation

  input :mutated_isoforms, :array, "e.g. ENSP0000001:A12V", nil
  task :mutated_isoform_neighbour_annotation => :string do |mis|
    residues = Structure.mutated_isoforms_to_residue_list(mis)

    neighbour_residues = {}
    residues.each do |protein, list|
      list = list.flatten
      neighbours = Structure.neighbours(protein, list.flatten)
      next if neighbours.nil? or neighbours.empty?

      neighbour_residues[protein] = neighbours
    end
    residues.annotate neighbour_residues
    neighbour_residues.type = :flat

    Open.write(file(:neighbour_uniprot), Structure.job(:annotate_residues_UNIPROT, name, :residues => neighbour_residues).clean.run.to_s)
    Open.write(file(:neighbour_appris), Structure.job(:annotate_residues_Appris, name, :residues => neighbour_residues).clean.run.to_s)
    Open.write(file(:neighbour_COSMIC), Structure.job(:annotate_variants_COSMIC, name, :residues => neighbour_residues).clean.run.to_s)

    "DONE"
  end
  export_asynchronous :mutated_isoform_neighbour_annotation


end

#require 'structure/cosmic_feature_analysis'


if defined? Entity and defined? MutatedIsoform and Entity === MutatedIsoform
  module MutatedIsoform
    property :pdbs_and_positions => :single do
      return [] if pdbs.nil?
      pdbs.collect do |pdb, info|
        [pdb, Structure.job(:sequence_position_in_pdb, "Protein: #{ self }", :sequence => protein.sequence, :organism => organism, :position => position, :pdb => pdb).run]
      end
    end
  end
end
