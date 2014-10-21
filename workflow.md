Functionalities regarding protein Structures

This workflow offers several functionalities to explore the consequence of
protein mutations. It reports features that overlap the mutations, or that are
in close physical proximity. 

The features reported include protein domains, variants, helices, ligand
binding residues, catalytic sites, transmembrane domains, InterPro domains, and
known somatic mutations in different types of cancer. This information is
extracted from resources such as UniProt, COSMIC, InterPro and Appris. It can
also identify mutations affecting the interfaces of protein complexes.

This workflow makes use of PDB files to calculate residues in close proximity.
This information is used to find features close to the mutations, at a distance
of 5 angstroms, or mutations in residues close to residues in a complex partner,
at a distance of up to 8 angstroms. 

PDBs are extracted from Interactome3d, which organized thousands of PDBs,
for both experimental structures and structure models, of individual
proteins and protein complexes.

Pairwise (Smith-Watterman) alignment is used to fix all inconsistencies between
protein sequences in PDBs, Uniprot and Ensembl Protein ID.

### Reference:

Vazquez M, Valencia A, Pons T. (2014) Structure-PPi:_ a module for the
annotation of cancer-related single-nucleotide variants at protein-protein
interfaces_. *Bioinformatics* (submitted)

# Tasks

The annotation tasks take either genomic mutations or mutated isoforms. Mutated
isoforms must be in reference to Ensembl Protein IDs (e.g.
ENSP00000449454:E433D). When inputing genomic mutations, the consequence of
this mutations in terms of mutated isoforms is computed automatically. Genomic
mutations can be reported in Chromosome:Position:AlternativeAllele format (e.g.
1:19949995:T) or in VCF format. Normaly genomic mutations are given always with
respect to the Watson or forward strand. This is the default. Sometimes,
however, mutations are given with respect to the strand that encodes de gene;
If this is the case, specify `watson` to be `false`. The version of the methods
using mutated isoforms as input have the `_mi_` term in their name. The
`organism` input is used to specify the version of the genome and gene set to
use, leave `Hsa` for the most recent `H`omo `sa`piens data, or `Hsa/may2009`
for the latest gene set of hg18. The date codes correspond the Ensembl
archives, which are used to download and process all the information
consistently for each version.

When PDBs are required, the PDB code or a URL can be specified using the `pdb`
parameter. Alternatively, the content of a PDB file can be provided using the
`pdbfile` parameter.

The result of these tasks are reported in TSV tables. These tables often encode
complex relationships. For instance a genomic mutation might have as
consequence several mutated isoforms, with each protein mutation overlapping
different domains or features. This multiplicity is encoded in the TSV in a
standard way: using tabs to separate fields, vertical bars (`|`) to separate
diferent values for each field, and semicolons (`;`) to further separate
multiple options for each result.

The main tasks are: `annotate`, `annotate_neighbours`, and `interfaces`. 
Or alternatively for mutated isoforms: `annotate_mi`, `annotate_mi_neighbours`, 
and `mi_interfaces`. 

All the workflow tasks can be executed from the web interface, programatically,
through the REST interface using `wget` or `curl`, or using the `rbbt` command.
See [the rbbt documentation](http://mikisvaz.github.io/rbbt/) for more
information.

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

