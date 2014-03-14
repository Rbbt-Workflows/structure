module Structure
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
  task :sequence_position_in_pdb => :tsv
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
  task :pdb_chain_position_in_sequence => :tsv 
  export_exec :pdb_chain_position_in_sequence


  desc <<-EOF
Find all pairs of residues in a PDB that fall within 'distance' of each other.

The pdb can be specified as a url or a PDB code, or a pdbfile can be provided directly.
  EOF
  input :distance, :float, "Distance", 5
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  task :neighbour_map => :tsv 
  export_asynchronous :neighbour_map

  desc <<-EOF
Use a pdb to find the residues neighbouring, in three dimensional space, a particular residue in a given sequence

The pdb can be specified as a url or a PDB code, or a pdbfile can be provided directly.
  EOF
  input :sequence, :text, "Protein sequence"
  input :positions, :array, "Positions inside sequence"
  input :pdb, :string, "Option 1: Name of pdb to align (from rcsb.org)", nil
  input :pdbfile, :text, "Option 2: Content of pdb to align", nil
  input :chain, :string, "Check only a particular chain", nil
  input :distance, :float, "Distance", 5
  task :neighbours_in_pdb => :tsv
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

end
