# OrthologAnalysis
A flexible pipeline for deducing orthology relationships between proteins in genomic (proteomic) datasets.

Data that can be generated and considered includes pairwise (BLASTP/DIAMOND) bi-directional best hits, synteny (conserved gene order), and membership in TIGRfam equivalog families (based on embedded trusted cutoffs).

Prerequisites:
perl
Bio::FeatureIO (for parsing gff files)
DIAMOND (https://github.com/bbuchfink/diamond) or BLASTP (https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
mcl (https://micans.org/mcl/)
hmmer (http://hmmer.org/)

This process also uses an altered version of MultiParanoid (

Input is protein fasta files (1 per organism) and gff3 files (required for synteny analysis).
The fasta identifiers must match the ID values in the gff3 file.

The files should all be in one directory. Running the ortholog_pipeline.sh script should be sufficient to run the process.
The shell script can be altered to customize analysis.

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
12. . . n. one column per organism containing the accession(s) for putative ortholog family members

