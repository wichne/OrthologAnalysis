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

tmp=${namebase}_tmp
if [ -e $tmp ]
then
    echo $tmp " already exists"
    # rm -rf $tmp
else
    echo "Make temporary folder $tmp"
    mkdir $tmp
fi

# delete any existing old run files
if [ -e output ]
then
    
   cd output
   #rm all_v_all.tab geneFamily multiparanoid.out synteny_ortho.out all.equivalog.tblout $namebase.outfile 
   cd ..
else
    mkdir output
fi
   
#######################################################
# check the input files
$ORTHO/check_input.pl --glob gff --faext faa
if [ $? -ne 0 ]
then
    echo "Problem with input files? $?"
    exit
fi

if [ -e $tmp/all.fasta ]
then
    echo "Found existing all.fasta. Continue? (y/n)"
    read fasta_continue
    if [ $fasta_continue != "y" ]
    then
	exit
    fi
else
    # Catenate the protein fasta files into all.pep
    echo "Fasta file is built from .faa files"
    cat input/*.faa >> $tmp/all.fasta
    # cd-hit -i all.faa -o nr.all.faa -c 1
fi

###############################################################
# start by running OrthoFinder
# this generates a LOT of files, but we should be able to get the BBH info from it.
if [ -e output/OrthoFinder ]
then
    echo "Found existing OrthoFinder run. Should I use it? (y/n) "
    read use_OF
    if [ $use_OF != "y" ]
    then
	echo "Running OrthoFinder..."
        /projects/bifx/OrthoFinder/orthofinder -f ./input -t 32
        mv ./input/OrthoFinder ./output/
    fi
else
    echo "Running OrthoFinder..."
    /projects/bifx/OrthoFinder/orthofinder -f ./input -t 32
    mv ./input/OrthoFinder ./output/
fi


##########################################################################
#### pair-wise all-vs-all
echo "Perform pair-wise all versus all search."

if [ -e "$tmp/all_v_all.tab" ]
then
    echo "Found existing all_v_all.tab. Should I use it? (y/n)"
    read use_pairwise
    if [ $use_pairwise != "y" ]
    then
        cd $tmp
    	# Perform all vs all search
	    if [ ! -e all.dmnd ]
        then
            # Format all.fasta for diamond
            diamond makedb --in all.fasta -d all
            # makeblastdb -in all.fasta -dbtype prot -parse_seqids
        fi
        echo "Performing diamond search"
	    diamond blastp -d all -q all.fasta -f 6 qseqid qlen qstart qend sseqid sseqid slen sstart send bitscore bitscore ppos pident evalue -o all_v_all.tab --masking seg
	    # blastp -db all.fasta -query all.fasta -num_threads 30 -evalue 1e-10 | tab_blast.pl > all_v_all.tab
        cd ..
    fi
else
    cd $tmp
    # Perform all vs all search
    echo "Performing diamond search"
    if [ ! -e all.dmnd ]
    then
        # Format all.fasta for diamond
        diamond makedb --in all.fasta -d all
        # makeblastdb -in all.fasta -dbtype prot -parse_seqids
    fi
    diamond blastp -d all -q all.fasta -f 6 qseqid qlen qstart qend sseqid sseqid slen sstart send bitscore bitscore ppos pident evalue -o all_v_all.tab --masking seg
    # blastp -db all.fasta -query all.fasta -num_threads 30 -evalue 1e-10 | tab_blast.pl > all_v_all.tab
    cd ..
fi

############################################################################

# Incorporate OrthoMCL into this script. It is already implemented in build_ortholog_table.pl (reading in groups.txt file)


# Use mcl to cluster based on blast output
echo "Cluster using mcl"
if [ -e $tmp/out.all.mci.I40 ]
then
    echo "Found an existing out.all.mci.I40. Should I use it? (y/n)"
    read use_mci
    if [ $use_mci != "y" ]
    then
        cd $tmp
    	cut -f 1,5,14 all_v_all.tab > mcl_in.abc
	    mcxload -abc mcl_in.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o all.mci -write-tab all.tab
	    mcl all.mci -I 4 -use-tab all.tab
        cd ..
    fi
else
    cd $tmp
    if [ -z all_v_all.tab ]
    then
        echo "all_v_all.tab is empty. Exiting."
        exit
    fi

    cut -f 1,5,14 all_v_all.tab > mcl_in.abc
    mcxload -abc mcl_in.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o all.mci -write-tab all.tab
    mcl all.mci -I 4 -use-tab all.tab
    cd ..
fi

# Build clusters using mmseqs2
echo "Build clusters using mmseqs2"
if [ -e $tmp/mmseqs2.out_cluster.tsv ]
then
    echo "Found existing mmseqs2.out_cluster.tsv. Should I use it? (y/n) "
    read use_clusters
    if [ $use_clusters != "y" ]
    then
        cd $tmp
        mmseqs easy-cluster all.fasta mmseqs2.out mmseqs_tmp/ --cluster-reassign
        cd ..
    fi
else
    cd $tmp
    mmseqs easy-cluster all.fasta mmseqs2.out mmseqs_tmp/ --cluster-reassign
    cd ..
fi

# Build pairwise orthologs by bidirectional best hit method
echo "Build pairwaise bbh clusters with multi_bbh"
if [ -e output/multiparanoid.out ]
then
    echo "Found existing multiparanoid.out. Should I use it? (y/n) "
    read use_clusters
    if [ $use_clusters != "y" ]
    then
	cd $tmp
	# v2 does gffs; the other does info
	$ORTHO/multi_bbh.pl -i all_v_all.tab -p ../input -e faa

	# Stitch the pairwise orthologs into ortholog groups
	echo "Build ortholog clusters with myMultiParanoid"
	$ORTHO/myMultiParanoid.pl --inpath BBH
	cd ..
    fi
else
    cd $tmp
    # v2 does gffs; the other does info
    $ORTHO/multi_bbh.pl -i all_v_all.tab -p ../input
    
    # Stitch the pairwise orthologs into ortholog groups
    echo "Build ortholog clusters with myMultiParanoid"
    $ORTHO/myMultiParanoid.pl --inpath BBH
    cd ..
fi

# analyze the bbh output for synteny
echo "Perform synteny analysis" ### !!!!!!! THIS IS NOT WORKING !!!!!!!!!!!!!!!!
if [ -e $tmp/synteny_ortho.out ]
then
    echo "Found an existing synteny_ortho.out. Should I use it? (y/n) "
    read do_synteny
    if [ $do_synteny != "y" ]
    then
        cd $tmp
    	# need the gff files for this step
        if [ -e BBH/multiparanoid.out ]
        then
            bbhfile=BBH/multiparanoid.out
        elif [ -e mmseqs2.out_cluster.tsv ]
        then 
            bbhfile=mmseqs2.out_cluster.tsv
        fi

        $ORTHO/synteny_analysis.pl -c $bbhfile
        cd ..
    fi
else
    cd $tmp
    # need the gff files for this step
    if [ -e BBH/multiparanoid.out ]
    then
        bbhfile=BBH/multiparanoid.out
    elif [ -e mmseqs2.out_cluster.tsv ]
    then 
        bbhfile=mmseqs2.out_cluster.tsv
    fi

    $ORTHO/synteny_analysis.pl -c $bbhfile
    cd ..
fi

# Do an HMM search against the equivalog models OR parse existing results
echo "Perform equivalog TIGRfam search"
if [ -e $tmp/all.equivalog.tblout ]
then
    echo "Found existing TIGRfam search results. Use them? (y/n) "
    read do_hmm
    if [ $do_hmm != "y" ]
    then
        cd $tmp
    	hmmsearch -o all.equivalog.out --noali --acc --tblout all.equivalog.tblout --cut_tc --cpu 24 /projects/db/hmms/TIGRFAM/equivalog_models.HMM all.fasta
        cd ..
    fi
else
    cd $tmp
    hmmsearch -o all.equivalog.out --noali --acc --tblout all.equivalog.tblout --cut_tc --cpu 24 /projects/db/hmms/TIGRFAM/equivalog_models.HMM all.fasta
    cd ..
fi
       
cp $tmp/all_v_all.tab output/
cp $tmp/out.all.mci.I40 output/geneFamily
cp $tmp/BBH/multiparanoid.out ../output/
cp $tmp/synteny_ortho.out output/
cp $tmp/all.equivalog.tblout output/
cp $tmp/mmseqs2.out_cluster.tsv output/mmseqs2.out
# OrthoFinder output dir is now moved to output when it's made (see above)
#cp input/OrthoFinder/.../Orthologues/Ortho... output/
 
# Use the output to build an ortholog table
#build_ortholog_table.pl -u nels329 -p RenMan -i multiparanoid_15720.out -P Halo -c geneFamily > Halo_ortho.tsv
#/home/bifx/inparanoid_4.1/build_ortholog_table.pl -i multiparanoid_13653.out -p RenMan -P OD1 -c out.all.mci.I40 > OD1.orthologs
echo "Build ortholog table"
$ORTHO/build_ortholog_table.pl --namebase "$namebase" > $namebase.ortho_out

#echo "Deleting temporary files in $tmp"
#rm -rf $tmp

# cat $namebase.ortho_out | perl -ne 'chomp; @f=split(/\t/, $_); @g=(); for (my $i=0; $i<@f; $i++) { if ($i < 13) { push @g, $f[$i] } else { my @e = split(/\s+/, $f[$i]); push @g, scalar(@e); }} print join("\t", @g) . "\n";' > $namebase.ortho_counts.csv
