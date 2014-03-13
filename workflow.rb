require 'rbbt'
require 'rbbt/util/misc'
require 'rbbt/workflow'

require 'rbbt/sources/organism'
require 'rbbt/sources/uniprot'

require 'structure/alignment'
require 'structure/pdb_alignment'

require 'structure/interactome_3d'
require 'structure/pdb_helper'
require 'structure/uniprot'
require 'structure/appris'
require 'structure/neighbours'
require 'structure/COSMIC'

Workflow.require_workflow 'Genomics'
Workflow.require_workflow 'Translation'
Workflow.require_workflow 'PdbTools'
Workflow.require_workflow "Sequence"

require 'rbbt/entity/mutated_isoform'

module Structure
  extend Workflow

  def self.mutated_isoforms_to_residue_list(mutated_isoforms)
    log :mutated_isoform_to_residue_list, "Find residues affected in each isoform" do
      residues = {}

      Annotated.purge(mutated_isoforms).each do |mi|
        protein, _sep, change = mi.partition ":"
        if change.match(/([A-Z])(\d+)([A-Z])$/)
          next if $1 == $3
          position = $2
          residues[protein] ||=[]
          residues[protein] << position.to_i
        end
      end

      TSV.setup(residues, :key_field => "Ensembl Protein ID", :fields => ["Residues"], :type => :flat, :cast => :to_i, :namespace => "Hsa")
    end
  end

  #{{{ ALIGNMENTS

  # In structure/pdb_alignment
  desc <<-EOF
Translate the positions inside a given amino-acid sequence to positions in the sequence of a PDB by
aligning them

The pdb can be specified as a url or a PDB code, or a pdbfile can be provided directly.
  EOF
  input :sequence, :text, "Protein sequence"
  input :positions, :array, "Positions within protein sequence"
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  task :sequence_position_in_pdb => :yaml
  export_exec :sequence_position_in_pdb


  # In structure/pdb_alignment
  desc <<-EOF
Translate the positions of amino-acids in a particular chain of the provided PDB into
positions inside a given sequence

The pdb can be specified as a url or a PDB code, or a pdbfile can be provided directly.
  EOF
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  input :chain, :string, "PDB chain"
  input :positions, :array, "Position within PDB chain"
  input :sequence, :text, "Protein sequence"
  task :pdb_chain_position_in_sequence => :array 
  export_exec :pdb_chain_position_in_sequence


  desc <<-EOF
Find all pairs of residues in a PDB that fall within 'distance' of each other.

The pdb can be specified as a url or a PDB code, or a pdbfile can be provided directly.
  EOF
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  input :distance, :float, "Distance", 5
  task :neighbour_map => :tsv do |pdb, pdbfile, distance|
    tsv = PDBHelper.pdb_close_residues(distance, pdb, pdbfile)
    TSV.setup tsv, :key_field => "Residue", :fields => ["Neighbours"], :type => :flat
  end
  export_asynchronous :neighbour_map

  desc <<-EOF
Use a pdb to find the residues neighbouring, in three dimensional space, a particular residue in a given sequence

The pdb can be specified as a url or a PDB code, or a pdbfile can be provided directly.
  EOF
  input :sequence, :text, "Protein sequence"
  input :positions, :integer, "Positions inside sequence"
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  input :chain, :string, "Check only a particular chain", 7
  input :distance, :float, "Distance", 5
  task :neighbours_in_pdb => :yaml
  export_asynchronous :neighbours_in_pdb

  # In structure/pdb_alignment
  desc <<-EOF
Find the correspondance between sequence positions in a PDB and in a given sequence

The pdb can be specified as a url or a PDB code, or a pdbfile can be provided directly.
  EOF
  input :sequence, :text, "Protein sequence"
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  task :pdb_alignment_map => :tsv 
  export_exec :pdb_alignment_map


  #{{{ ANNOTATIONS


  desc <<-EOF
Given a set of proteins and resudies inside these proteins, finds the protein features that overlap, as annotated in UniProt.

The proteins and residues are specified as a TSV file, with Ensembl Protein IDs
as key, and residues as values. Types `flat` (separated by more tabs),
and `double` (separated by '|') are supported.
  EOF
  input :residues, :tsv, "Proteins and their affected residues", nil
  task :annotate_residues_UNIPROT => :tsv do |residues|
    raise ParameterException, "No residues provided" if residues.nil?
    tsv = TSV.setup({}, :key_field => "Ensembl Protein ID", :fields => ["Residue", "UniProt Features", "UniProt Feature locations", "UniProt Feature Descriptions"], :type => :double)

    iso2uni = Organism.protein_identifiers("Hsa").index :target => "UniProt/SwissProt Accession", :persist => true
    iso2sequence = Organism.protein_sequence("Hsa").tsv :type => :single, :persist => true

    missing = []
    residues.each do |isoform, list|
      list = Array === list ? list.flatten : [list]

      uniprot = iso2uni[isoform]
      if uniprot.nil?
        missing << isoform
        next
      end

      features = Structure.corrected_uniprot_features(uniprot, iso2sequence[isoform])
      overlapping = [[],[],[],[]]
      list.each do |position|
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

    if missing.any?
      Log.warn "Some isoforms failed to translate to UniProt: #{missing.length}"
      set_info(:missing, missing)
    end

    tsv
  end
  export_asynchronous :annotate_residues_UNIPROT

  desc <<-EOF
Given a set of proteins and resudies inside these proteins, finds the protein features that overlap, as annotated in Appris.

The proteins and residues are specified as a TSV file, with Ensembl Protein IDs
as key, and residues as values. Types `flat` (separated by more tabs),
and `double` (separated by '|') are supported.
  EOF
  input :residues, :tsv, "Proteins and their affected residues", nil
  task :annotate_residues_Appris => :tsv do |residues|
    tsv = TSV.setup({}, :key_field => "Ensembl Protein ID", :fields => ["Residue", "Appris Features", "Appris Feature locations", "Appris Feature Descriptions"], :type => :double)

    iso2sequence = Organism.protein_sequence("Hsa").tsv :type => :single, :persist => true

    residues.each do |isoform, list|
      list = Array === list ? list.flatten : [list]

      features = Structure.appris_features(isoform)

      overlapping = [[],[],[],[]]
      list.each do |position|
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

  desc <<-EOF
Given a set of proteins and resudies inside these proteins, finds the mutations registered in COSMIC that affect those residues, and provide
some annotations from the samples that contained them

The proteins and residues are specified as a TSV file, with Ensembl Protein IDs
as key, and residues as values. Types `flat` (separated by more tabs),
and `double` (separated by '|') are supported.

The result is the proteins along with the overlapping features and some information about them
  EOF
  input :residues, :tsv, "Proteins and their affected residues", nil
  task :annotate_variants_COSMIC => :tsv do |residues|

    cosmic_residue_mutations = Structure.COSMIC_residues

    isoform_matched_variants = {}
    residues.through do |protein,list|
      list = Array === list ? list.flatten : [list]

      list.each do |position|
        matching_mutations = cosmic_residue_mutations[[protein, position]*":"]
        next if matching_mutations.nil? or matching_mutations.empty?
        isoform_matched_variants[protein] ||= []
        isoform_matched_variants[protein].concat matching_mutations
      end
    end

    cosmic_mutation_annotations = Structure.COSMIC_mutation_annotations

    res = TSV.setup({}, :key_field => "Ensembl Protein ID", :fields => ["Residue", "Genomic Mutation"].concat(cosmic_mutation_annotations.fields), :type => :double)

    isoform_matched_variants.each do |protein, mutations|
      values = []
      mutations.each do |mutation|
        begin
          annotations = cosmic_mutation_annotations[mutation]
          raise "No annotations for #{ mutation } in #{cosmic_mutation_annotations.persistence_path}" if annotations.nil?
          aa_mutation = (annotations.first || []).compact.first
          raise "No AA mutation" if aa_mutation.nil?
          if aa_mutation.match(/(\d+)/)
            residue = $1
          else
            raise "AA mutation does not have a position: #{ aa_mutation }" if aa_mutation.nil?
          end
          values << [residue,mutation].concat(annotations || [])
        rescue
          Log.exception $!
          next
        end
      end

      res[protein] = Misc.zip_fields(values)
    end

    res
  end
  export_asynchronous :annotate_variants_COSMIC


  desc <<-EOF
Given a set of proteins and resudies inside these proteins, finds the known variants that overlap, as annotated in UniProt.

The proteins and residues are specified as a TSV file, with Ensembl Protein IDs
as key, and residues as values. Types `flat` (separated by more tabs),
and `double` (separated by '|') are supported.

The result is the proteins along with the overlapping features and some information about them
  EOF
  input :residues, :tsv, "Proteins and their affected residues", nil
  task :annotate_variants_UNIPROT => :tsv do |residues|

    uniprot_residue_mutations = Structure.UniProt_residues

    isoform_matched_variants = {}
    residues.each do |protein, list|
      list = Array === list ? list.flatten : [list]


      list.each do |position|
        key = [protein, position]*":"

        matching_mutations = uniprot_residue_mutations[key]
        next if matching_mutations.nil? or matching_mutations.empty?
        isoform_matched_variants[protein] ||= []
        isoform_matched_variants[protein].concat matching_mutations
      end
    end

    uniprot_mutation_annotations = Structure.UniProt_mutation_annotations

    res = TSV.setup({}, :key_field => "Ensembl Protein ID", :fields => ["Residue"] + uniprot_mutation_annotations.fields, :type => :double)

    isoform_matched_variants.each do |protein, mutations|
      values = []
      mutations.each do |mutation|
        begin
          annotations = uniprot_mutation_annotations[mutation]
          raise "No annotations for #{ mutation } in #{uniprot_mutation_annotations.persistence_path}" if annotations.nil?
          residue = annotations.first.scan(/\d+/)
          values << [residue].concat(annotations || [])
        rescue
          Log.exception $!
        end
      end

      res[protein] = Misc.zip_fields(values)
    end

    res
  end
  export_asynchronous :annotate_variants_UNIPROT


  #{{{ MUTATED ISOFORMS ANALYSIS

  helper :residue_neighbours do |residues|
    log :residue_neighbours, "Find neighbouring residues"
    all_neighbours = TSV.setup({}, :key_field => "Isoform:residue", :fields => ["Ensembl Protein ID", "Residue", "PDB", "Neighbours"], :type => :list)

    residues.with_monitor :desc => "Finding neighbours" do
      residues.through do |protein, list|
        list = list.flatten.compact.uniq
        neighbours = Structure.neighbours_i3d(protein, list)
        next if neighbours.nil? or neighbours.empty?
        if all_neighbours.empty?
          all_neighbours = neighbours
        else
          all_neighbours.merge! neighbours
        end
      end
    end

    all_neighbours
  end

  helper :mutated_isoforms do |mutations,organism|
    raise ParameterException, "No mutated_isoforms or genomic_mutations specified" if mutations.nil? 

    log :mutated_isoforms, "Extracting mutated_isoforms from genomic_mutations" do

      job = Sequence.job(:mutated_isoforms_for_genomic_mutations, name, :mutations => mutations, :organism => organism)

      mis = Set.new
      job.run(true).path.traverse do |m,_mis|
        mis.merge _mis
      end

      mis.to_a
    end
  end

  desc <<-EOF
Given a set of protein mutations, finds the affected protein features and variants in Appris, COSMIC, UniProt

The results for the different annotation types are save as job files
  EOF
  input :mutated_isoforms, :array, "e.g. ENSP0000001:A12V", nil
  input :genomic_mutations, :array, "e.g. 1:173853127:T", nil
  input :organism, :string, "Organism code", "Hsa"
  task :mutated_isoform_annotation => :string do |mis,muts,organism|
    mis = mutated_isoforms muts, organism if mis.nil?

    residues = Structure.mutated_isoforms_to_residue_list(mis)

    Open.write(file(:uniprot), Structure.job(:annotate_residues_UNIPROT, name, :residues => residues).run.to_s)
    Open.write(file(:appris), Structure.job(:annotate_residues_Appris, name, :residues => residues).run.to_s)
    Open.write(file(:COSMIC), Structure.job(:annotate_variants_COSMIC, name, :residues => residues).run.to_s)
    Open.write(file(:uniprot_variants), Structure.job(:annotate_variants_UNIPROT, name, :residues => residues).run.to_s)

    "DONE"
  end
  export_asynchronous :mutated_isoform_annotation


  desc <<-EOF
Given a set of protein mutations, finds the neighbouring protein features and variants in Appris, COSMIC, UniProt.

The results for the different annotation types are save as job files. This method
consider only the features of neighbouring residues, no the residue itself.
  EOF
  input :mutated_isoforms, :array, "e.g. ENSP0000001:A12V", nil
  input :genomic_mutations, :array, "e.g. 1:173853127:T", nil
  input :organism, :string, "Organism code", "Hsa"
  task :mutated_isoform_neighbour_annotation => :tsv do |mis,muts,organism|
    mis = mutated_isoforms muts, organism if mis.nil?

    residues = Structure.mutated_isoforms_to_residue_list(mis)

    neighbours = self.residue_neighbours residues

    neighbour_residues = {}

    log :neighbour_residues, "Sorting residues"
    neighbours.pthrough do |iso,values|
      protein, _pos = iso.split ":"
      positions = values.last.split ";"
      neighbour_residues[protein] ||= []
      neighbour_residues[protein].concat positions
      neighbour_residues[protein].uniq!
    end
    residues.annotate neighbour_residues

    log :annotation, "Annotating neighbours"
    Open.write(file(:neighbour_uniprot), Structure.job(:annotate_residues_UNIPROT, name, :residues => neighbour_residues).clean.run.to_s)
    Open.write(file(:neighbour_uniprot_variants), Structure.job(:annotate_variants_UNIPROT, name, :residues => neighbour_residues).clean.run.to_s)
    Open.write(file(:neighbour_appris), Structure.job(:annotate_residues_Appris, name, :residues => neighbour_residues).clean.run.to_s)
    Open.write(file(:neighbour_COSMIC), Structure.job(:annotate_variants_COSMIC, name, :residues => neighbour_residues).clean.run.to_s)

    neighbours
  end
  export_asynchronous :mutated_isoform_neighbour_annotation


  desc <<-EOF
Given a set of protein mutations, finds residues that affect protein interfaces, as represented
in PDB models of protein complexes in Interactome3d

The results for the different annotation types are save as job files. This method
consider only the features of neighbouring residues, no the residue itself.
  EOF
  input :mutated_isoforms, :array, "e.g. ENSP0000001:A12V", nil
  input :genomic_mutations, :array, "e.g. 1:173853127:T", nil
  input :organism, :string, "Organism code", "Hsa"
  task :mutated_isoform_interface_residues => :tsv do |mis,muts,organism|
    mis = mutated_isoforms muts, organism if mis.nil?

    residues = Structure.mutated_isoforms_to_residue_list(mis)

    neighbour_residues = {}

    residues.ppthrough_callback do |neig|
      next if neig.nil?
      neighbour_residues.merge!(neig)
    end

    residues.with_monitor :desc => "Processing residues" do
    residues.ppthrough(7) do |protein, list|
      list = list.flatten.compact.sort
      begin
        neighbours = Structure.interface_neighbours_i3d(protein, list)
        next if neighbours.nil? or neighbours.empty?
        neighbours
      rescue
        Log.exception $!
        next
      end
    end
    end

    TSV.setup(neighbour_residues, :key_field => "Isoform:residue:partner", :fields => ["Partner", "PDB", "Partner residues"], :type => :double)
  end
  export_asynchronous :mutated_isoform_interface_residues

  desc <<-EOF
Given a set of protein mutations, finds residues that affect protein interfaces, as represented
in PDB models of protein complexes in Interactome3d

The results for the different annotation types are save as job files. This method
consider only the features of neighbouring residues, no the residue itself.
  EOF
  input :mutated_isoforms, :array, "e.g. ENSP0000001:A12V", nil
  input :genomic_mutations, :array, "e.g. 1:173853127:T", nil
  input :organism, :string, "Organism code", "Hsa"
  task :mutated_isoform_interface_residues_old => :tsv do |mis,muts,organism|
    mis = mutated_isoforms muts, organism if mis.nil?

    residues = Structure.mutated_isoforms_to_residue_list(mis)

    neighbour_residues = TSV.setup({}, :key_field => "Isoform:residue:partner", :fields => ["Partner", "PDB", "Partner residues"], :type => :double)

    residues.with_monitor :desc => "Processing residues" do
    residues.pthrough do |protein, list|
      list = list.flatten.compact.sort
      neighbours = Structure.interface_neighbours_i3d(protein, list)
      next if neighbours.nil? or neighbours.empty?
      neighbour_residues.merge!(neighbours)
    end
    end
    neighbour_residues
  end
  export_asynchronous :mutated_isoform_interface_residues
end

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
