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
of 5 angstroms, or mutations in residues close to residues in a complex partner,
at a distance of up to 8 angstroms. 

PDBs are extracted from Interactome3d, which has organized thousands of PDBs,
including both experimental structures and structure models, of individual
proteins, and of protein complexes.

Pairwise (Smith-Watterman) alignment is used to fix all inconsistencies between
protein sequences and sequence differences between proteins in Ensembl and
proteins in UniProt.

# Tasks

The annotation tasks take `residues` as input. These residues are given as a
flat TSV file of proteins and positions over them. The result is a TSV file
containing the proteins and features that overlap. The fields of the TSV
file depend on the features of database used.

When PDBs are required, the PDB code or a URL can be specified using the `pdb`
parameter. Alternatively, the content of a PDB file can be provided using the
`pdbfile` parameter.

## neighbour_map

For a given PDB, find all pairs of residues in a PDB that fall within a given
'distance' of each other. It uses PDBs from Interactome3d for individual
proteins.


## neighbours_in_pdb

Use a pdb to find the residues neighbouring, in three dimensional space, a
particular residue in a given sequence. 

## pdb_alignment_map

Find the correspondance between sequence positions in a PDB and in a given
sequence. PDB positions are reported as `chain:position`.

## pdb_chain_position_in_sequence

Translate the positions of amino-acids in a particular chain of the provided
PDB into positions inside a given sequence.

## annotate_residues_Appris

Given a set of proteins and resudies inside these proteins, finds the protein
features that overlap, as annotated in Appris.

## annotate_residues_InterPro

Given a set of proteins and resudies inside these proteins, finds the protein
domains that overlap, as annotated in InterPro.

Sequence alignment is used to correct discrepancies between sequences in
UniProt and in Ensembl.

## annotate_residues_UniProt

Given a set of proteins and resudies inside these proteins, finds the protein
features that overlap, as annotated in UniProt.

Sequence alignment is used to correct discrepancies between sequences in
UniProt and in Ensembl.

## annotate_residues_COSMIC

Given a set of proteins and resudies inside these proteins, finds the mutations
registered in COSMIC that affect those residues, and provide some annotations
from the samples that contained them

## residue_interfaces

Finds residues that lay over protein complex interfaces. It does so by checking
PDBs of complexes extracted from Interactome3D and looking for residues in
close physical proximity  to residues in the complex partner protein (8
angstroms).

## sequence_position_in_pdb

Translate the positions inside a given amino-acid sequence to positions in the
sequence of a PDB by aligning them

## annotated_variants

Annotates variants given as mutated isoforms or genomic mutations. It
translates them to protein and residue positions and annotates them using the
different annotation methods.

## annotated_variant_neighbours

Annotates variants given as mutated isoforms or genomic mutations. It
translates them to protein and residue positions, uses them to find
neighbouring residues, and annotates them using the different annotation methods.

## variant_interfaces

Identifies variants, given as mutated isoforms or genomic mutations, that
occur over protein-protein binding interfaces. 
