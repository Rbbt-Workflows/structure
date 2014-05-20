Functionalities regarding protein Structures

This workflow offers several functionalities to explore the consequence of
protein mutations. It reports features that overlap the mutations or that are
in close physical proximity. 

The features reported include features domains, variants, helices, ligand
binding residues, catalytic sites, transmembrane domains, InterPro domains, or
known somatic mutations in different types of cancer. This information is
extracted from resources such as UniProt, COSMIC, InterPro and Appris.  It can
also identify mutations affecting the interfaces of protein complexes.

This workflow makes use of PDB files to calculate residues in close proximity.
This information is used to find features close to the mutations, at a distance
of 5 angstroms, or mutations in residues close to residues in a complex partner,
at a distance of up to 8 angstroms. 

PDBs are extracted from Interactome3d, which has organized thousands of PDBs,
including both experimental structures and structure models, of individual
proteins, and of protein complexes.

Pairwise (Smith-Watterman) alignment is used to fix all inconsistencies between
protein sequences in PDBs, and sequence differences between proteins in Ensembl and
proteins in UniProt.

# Tasks

The annotation tasks take `residues` as input. These residues are given as a
flat TSV file of proteins and positions over them. The result is a TSV file
containing the proteins and features that overlap. The fields of the TSV
file depend on the features of database used.

When PDBs are required, the PDB code or a URL can be specified using the `pdb`
parameter. Alternatively, the content of a PDB file can be provided using the
`pdbfile` parameter.

The main tasks are: `annotate`, `annotate_neighbours`, and `interfaces`. 
Or alternatively for mutated isoforms: `annotate_mi`, `annotate_mi_neighbours`, and `mi_interfaces`. 

## annotate

Annotates genomic mutations based on the protein features that are overlapping
amino-acid changes

## annotate_neighbours

Annotates genomic mutations based on the protein features that are in close
physical proximity to amino-acid changes

At a distance of 5 angstroms 

## interfaces

Find variants that affect residues in protein-protein interaction surfaces

Residues at a distance of 8 angstroms of a residue from an interaction partner

## annotate_mi

Annotates mutated isoforms based on the protein features that are overlapping
amino-acid changes

## annotate_mi_neighbours

Annotates mutated isoforms based on the protein features that are in close
physical proximity to amino-acid changes

At a distance of 5 angstroms 

## mi_interfaces

Find mutated_isoforms with affected residues in protein-protein interaction sufaces

Residues at a distance of 8 angstroms of a residue from an interaction partner

## mi_neighbours

Finds residues physical proximity to amino-acid changes in mutated isoforms

At a distance of 5 angstroms 

## neighbour_map

For a given PDB, find all pairs of residues in a PDB that fall within a given
'distance' of each other. It uses PDBs from Interactome3d for individual
proteins.

## neighbours_in_pdb

Use a pdb to find the residues neighbouring, in three dimensional space, a
particular residue in a given sequence. 

## pdb_alignment_map

Find the correspondence between sequence positions in a PDB and in a given
sequence. PDB positions are reported as `chain:position`.

## pdb_chain_position_in_sequence

Translate the positions of amino-acids in a particular chain of the provided
PDB into positions inside a given sequence.

## sequence_position_in_pdb

Translate the positions inside a given amino-acid sequence to positions in the
sequence of a PDB by aligning them
