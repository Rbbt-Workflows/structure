module Structure
  # In structure/pdb_alignment
  input :sequence, :text, "Protein sequence"
  input :positions, :array, "Positions within protein sequence"
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  task :sequence_position_in_pdb => :tsv
  export_exec :sequence_position_in_pdb

  # In structure/pdb_alignment
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  input :chain, :string, "PDB chain"
  input :positions, :array, "Position within PDB chain"
  input :sequence, :text, "Protein sequence"
  task :pdb_chain_position_in_sequence => :tsv 
  export_exec :pdb_chain_position_in_sequence

  input :distance, :float, "Distance", 5
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  task :neighbour_map => :tsv 
  export_asynchronous :neighbour_map

  input :sequence, :text, "Protein sequence"
  input :positions, :array, "Positions inside sequence"
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  input :chain, :string, "Check only a particular chain", nil
  input :distance, :float, "Distance", 5
  task :neighbours_in_pdb => :tsv
  export_asynchronous :neighbours_in_pdb

  # In structure/pdb_alignment
  input :sequence, :text, "Protein sequence"
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  task :pdb_alignment_map => :tsv 
  export_exec :pdb_alignment_map

end
