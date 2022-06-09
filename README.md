# OrthologAnalysis
A flexible pipeline for deducing orthology relationships between proteins in genomic (proteomic) datasets.

Data that can be generated and considered includes pairwise (BLASTP/DIAMOND) bi-directional best hits, synteny (conserved gene order), and membership in TIGRfam equivalog families (based on embedded trusted cutoffs).

Input is protein fasta files (1 per organism) and gff3 files (required for synteny analysis).
The fasta identifiers must match the ID values in the gff3 file.

Output is a tab-delimited table with the following columns
1. mcl cluster id
2. bbh cluster id
3. synteny family id
4. TIGRfam id
5. genome count
6. proteins count
7. product name
8. gene symbol
9. role category
10. EC number
11. TC number
12. - n. one column per organism containing the accession(s) for putative ortholog family members

