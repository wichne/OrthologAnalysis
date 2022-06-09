#!/bin/bash
export ORTHO=/projects/perlscripts/OrthologAnalysis

function joinArray { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }

# USAGE: ortholog_pipeline.sh NAMEBASE_FOR_ORTHO_FAMILIES
# This script will then scrape .faa and .gff files from the current directory as the input for the process
# You should have a pep and gff file for each genome you want to analyze.

# The following programs must be in the executable path:
# makeblastdb and blastp [NCBI BLAST]
# diamond
# mcxload and mcl [MCL]
# multi_bbh.pl
# myMultiParanoid.pl
# synteny_analysis.pl
# build_ortholog_table.pl
# hmmsearch [HMMer]

# This is the prefix for the ortho family id
namebase=$1

# delete any existing old run files
rm nr.all.fasta nr.all.fasta.clstr all.fasta all.fasta.* all_v_all.tab mcl_in.abc all.tab all.mci *.ortho synteny_ortho.out all.equivalog.out all.equivalog.tblout multiparanoid.out out.all.mci.I40 $namebase.outfile geneFamily

# Get the number of data sets. 
n=$(ls *.faa | wc -l)

echo "n=" $n

# make the data set string.
echo "Organism identifiers are drawn from the file prefix of the .gff files."

# Catenate the protein fasta files into all.pep
echo "Fasta file is built from .faa files"
cd input
cat *.faa >> all.fasta
# cd-hit -i all.faa -o nr.all.faa -c 1

# Format all.pep for blast
echo "Protein all-v-all search"
diamond makedb --in all.fasta -d all
# makeblastdb -in nr.all.faa -dbtype prot -parse_seqids

# Perform all vs all search
diamond blastp -d all -q all.fasta -f 6 qseqid qlen qstart qend sseqid sseqid slen sstart send bitscore bitscore ppos pident evalue -o all_v_all.tab
# blastp -db nr.all.faa -query nr.all.faa -num_threads 30 -evalue 1e-10 | tab_blast.pl > all_v_all.tab

# Use mcl to cluster based on blast output
echo "Cluster using mcl"
cut -f 1,5,14 all_v_all.tab > mcl_in.abc
mcxload -abc mcl_in.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o all.mci -write-tab all.tab
mcl all.mci -I 4 -use-tab all.tab
cp out.all.mci.I40 geneFamily

# Build pairwise orthologs by bidirectional best hit method
echo "Build pairwaise bbh clusters with multi_bbh"
# v2 does gffs; the other does info
$ORTHO/multi_bbh.pl -i all_v_all.tab -p .

# Stitch the pairwise orthologs into ortholog groups
echo "Build ortholog clusters with myMultiParanoid"
$ORTHO/myMultiParanoid.pl --inpath .

# analyze the bbh output for synteny
echo "Perform synteny analysis"
# need the gff files for this step
m=$(ls *.gff | wc -l)
if [ $n != $m ]
then
    echo "Can't perform synteny analysis: faa files ($n) not equal to gff files ($m)"
    exit
fi
$ORTHO/synteny_analysis.pl

# Do an HMM search against the equivalog models OR parse existing results
echo "Perform equivalog TIGRfam search"
hmmsearch -o all.equivalog.out --noali --acc --tblout all.equivalog.tblout --cut_tc --cpu 24 /projects/db/hmms/TIGRFAM/equivalog_models.HMM all.fasta

# Use the output to build an ortholog table
#build_ortholog_table.pl -u nels329 -p RenMan -i multiparanoid_15720.out -P Halo -c geneFamily > Halo_ortho.tsv
#/home/bifx/inparanoid_4.1/build_ortholog_table.pl -i multiparanoid_13653.out -p RenMan -P OD1 -c out.all.mci.I40 > OD1.orthologs
echo "Build ortholog table"
$ORTHO/build_ortholog_table.pl -namebase "../$namebase" -p .
cd ..
