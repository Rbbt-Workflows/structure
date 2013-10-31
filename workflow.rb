require 'rbbt'
require 'rbbt/util/misc'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/sources/uniprot'
require 'ssw'
require 'interactome_3d'
require 'pdb_helper'

Workflow.require_workflow 'Genomics'
Workflow.require_workflow 'Translation'
Workflow.require_workflow 'PdbTools'
Workflow.require_workflow 'Appris'

module Structure
  extend Workflow

  def self.alignment_map(alignment1, alignment2)
    map = {}

    offset1, alignment1 = alignment1.match(/^(_*)(.*)/).values_at 1, 2
    offset2, alignment2 = alignment2.match(/^(_*)(.*)/).values_at 1, 2

    gaps1 = 0 
    gaps2 = 0
    alignment1.chars.zip(alignment2.chars).each_with_index do |p,i|
      char1, char2 = p
      gaps1 += 1 if char1 == '-'
      gaps2 += 1 if char2 == '-'
      map[i + offset1.length + gaps1] = i + offset2.length + gaps2 if char1 == char2  and char1 != "-"
    end
    
    map
  end

  def self.match_position(protein_position, protein_alignment, chain_alignment)
    if Array === protein_position
      alignment_map(protein_alignment, chain_alignment).chunked_values_at protein_position.collect{|p| p.to_i}
    else
      alignment_map(protein_alignment, chain_alignment)[protein_position]
    end
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
  task :pdb_chain_position_in_sequence => :array do |pdb, pdbfile, chain, chain_position, protein_sequence|
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

    protein_positions = protein_position.split(/,/).collect{|p| p.strip.to_i}
    Structure.match_position(chain_positions, chain_alignment, protein_alignment)
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

    neighbours = PdbTools.job(:pdb_close_positions, name, :pdb => pdbfile, :distance => distance).run

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
end

require 'cosmic_feature_analysis'

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
