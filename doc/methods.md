Methods
=======

## Residue annotation

The `annotate_residues_*` methods take a TSV file with proteins as keys (in
`Ensembl Protein ID` format), and lists of positions in the sequence as values.
It then traverses the file annotating those residues that fall over features
annotated in the corresponding databases: UniProt, Appris, and COSMIC.

UniProt features are extracted from their corresponding UniProt entries. `Ensembl
Protein ID` codes are translated to `UniProt/SwissProt Accession`, the
UniProt web site is queried for the entry, and the features are extracted with
their location inside the sequence. Since UniProt sequence may differ from the
sequence of the Ensembl isoform, an alignment map  is produced between the
sequences and the locations of the features are translated in to locations
in the Ensembl sequence of the isoform. The `organism` parameter is used
to select the right build to use in extracting the basic Ensembl information 
(identifier equivalences, and protein sequences).

COSMIC mutations are translated into protein mutations and indexed. Residues
are annotated with the mutations in COSMIC that affect that same residue, and
with the sample annotations associated to those variants in COSMIC (tumor site,
histology, etc). Since the proteins in COSMIC are listed in hg19, performing
the match at the protein residue level allows us to match it with variants
found in previous builds (e.g. hg18), as long as the isoforms codes and
sequences are the same. No steps have been taken to solve any possible
inconsistencies that could still remain.

Appris information includes domains, transmembrane helices, ligand biding
residues and catalytic sites.

## Variant annotations

Variants can be specified as `genomic_mutations` (e.g. 12:9265129:A) or
`mutated_isoforms` (e.g. ENSP00000323929:P193L). Mutated isoforms will be
derived from the genomic mutations using the `Sequence` workflow. Mutated
isoforms are then translated into TSV with proteins and residue positions,
which is then forwarded to the residue annotation methods. The results of
residue annotations are mapped back to the mutated isoforms, and genomic
mutations if provided.

## Variant neighbour annotations

Variants are translated into proteins and residues, as decribed above. Then
these residues are matched to neighbours using the `neighbour_map` task. 
Neighbouring residues are then subject of annotations as detailed above.
The distance for neighbours is set to 4 angstroms.

## Residue and variant interfaces

A `neighbour_map` of PDB files containing complexes extracted from
Interactome3d is used to identify residues in our protein of interest
that are in close proximity. The distance for interface residues
is set to 8 angstroms. Variants are translated into residues prior to 
annotation, as is described above.

## Other uses: homology, and paralogy

This workflow exports some low-level methods to map sequences to PDB. Since
these methods can take any input sequence, they can be used for homolog
and paralog sequences as well.


