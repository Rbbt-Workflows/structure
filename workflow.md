Functionalities regarding protein Structures

This workflow offers several functionalities to explore the consequence of
protein mutations. It reports features that overlap the mutations or that are
in close physical proximity. 

The features reported include features domains, variants, helices, ligand
binding residues, catalitic sites, transmembrane domains, InterPro domains, or
known somatic mutations in different types of cancer. This information is
extracted from resources such as UniProt, COSMIC, InterPro and Appris.  It can
also identify mutations affecting the interfaces of protein complexes.

This workflow makes use of PDB files to calculate residues in close proximity.
This information is used to find features close to the mutations, at a distance
of 5 amstrons, or mutations in residues close to residues in a complex partner. 
PDBs are extracted from Interactome3d, which has organized thousands of PDBs,
including both experimental structures and structure models, both of individual
proteins, and of protein complexes.

Pairwise (Smith-Watterman) alignment is used to fix all inconsistencies between
protein sequences and sequence differences between proteins in Ensembl and
proteins in UniProt.

# Tasks

The annotation tasks take `residues` as input. These residues are given as a
flat TSV file of proteins and positions over them. The result is a TSV file
containing the proteins and features that overlap. The fields of the TSV
file depend on the features of database used.


## annotate_residues_Appris


## annotate_residues_UNIPROT

Given a set of proteins and resudies inside these proteins, finds the protein
features that overlap, as annotated in UniProt.

## annotate_variants_COSMIC

Given a set of proteins and resudies inside these proteins, finds the mutations
registered in COSMIC that affect those residues, and provide some annotations
from the samples that contained them

## annotate_variants_UNIPROT

Given a set of proteins and resudies inside these proteins, finds the known
variants that overlap, as annotated in UniProt.

## neighbour_map

Find all pairs of residues in a PDB that fall within 'distance' of each other.

## neighbours_in_pdb

Use a pdb to find the residues neighbouring, in three dimensional space, a
particular residue in a given sequence

## pdb_alignment_map

Find the correspondance between sequence positions in a PDB and in a given
sequence

## pdb_chain_position_in_sequence

Translate the positions of amino-acids in a particular chain of the provided
PDB into positions inside a given sequence

## residue_interfaces

Finds residues that lay over protein complex interfaces. It does so by checking
PDBs of complexes extracted from Interactome3D and looking for residues close 
to residues in the complex partner (8 amstrongs).

## sequence_position_in_pdb

Translate the positions inside a given amino-acid sequence to positions in the
sequence of a PDB by aligning them


