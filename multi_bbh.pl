#!/usr/bin/perl
use strict;
use Getopt::Std;
use Bio::FeatureIO;
use lib $ENV{ORTHO};
use OrthologAnalysis;

my $usage = "
    # This program looks for bi-directional best hits between organism pairs.
    # Input is blastp output from an all_v_all search between multiple genomes.
    # The top non-self hit for each gene to each genome is stored in a hash, then we go 
    # through each protein from one of the genomes and check for reciprocal
    # best hits. Cutoffs are >30% id over > 70% of the query length
    # for the first iteration. (The length cutoff only has to be unidirectional
    # to accommodate fusion proteins.) The avg and stdev for %id (AAI) are
    # calculated from the initial results. In the second iteration, the cutoff 
    # is set to avg-2*stdev, constraining the 'best' hits to reasonable scores.
    # Clear?

  USAGE: multi_bbh.pl -i inputfile -p results_dir [-e fasta_file_ext]

    -i - blast_tab.pl output from an all vs all blastp search.
    -p - file path to directory containing input files and all_v_all results
    -e - file extension used for input fasta files (default 'faa')
    ";

my %arg;
&getopts('i:p:e:h', \%arg);

if ($arg{'h'}) { die $usage }
my $file = $arg{'i'} or die "Need to provide input btab file with -i\n";
my $path = $arg{'p'};
my $ext = $arg{'e'} ? $arg{'e'} : 'faa';

my $VERBOSE = 1;

our $length_cutoff = 0.7;
our $id_cutoff = 30;

# still have the problem of assocating proteins with genomes, since it's not explicit in the
# pair-wise output. In this instance, we will use the fasta files, since they were the
# original input for the searches.
my ($PtoDB, $DBtoP) = &parse_protein_fasta_files($path, $ext);

my $HITS = &read_file($file, $PtoDB); # populates %HITS

our %BBH;
my $per_id_ref = &look_for_bbh(keys %$DBtoP); # populates %BBH

my @DB = sort keys %$per_id_ref;

# Roll through the pair-wise comparisons
for (my $i=0; $i<@DB; $i++) {
    my $db1 = $DB[$i];
    print STDERR "Doing $db1...\n";
    for (my $j=$i+1; $j<@DB; $j++) {
	my $db2 = $DB[$j];
	print STDERR "...with $db2\n";
	if (! defined $per_id_ref->{$db1}->{$db2}) { next } # If we don't have scores for these two dbs

	# Calculate the AAI for the two genomes so we can set a threshhold.
	my $index = scalar(@{$per_id_ref->{$db1}->{$db2}}); # how many bbhs are there?
	my $sum = 0;
	foreach (@{$per_id_ref->{$db1}->{$db2}}) { $sum += $_ } # sum the AAIDs
	my $avg = $sum/$index; # calculate the average

	# Calculate the stdev for bbh AAI
	my $x = 0;
	foreach (@{$per_id_ref->{$db1}->{$db2}}) { $x += ($_ - $avg)**2 }
	my $stdev = sqrt($x/$index);

	# The threshhold is 2 stdev from the average.
	my $keep_threshold = sprintf "%.2f\n",($avg-2*$stdev);
	if ($keep_threshold < $id_cutoff) { $keep_threshold = $id_cutoff }
	
	if ($VERBOSE) {
	    print "$db1 - $db2\n";
	    print "index =    $index\n";
	    printf "Avg per_id = %.2f\n",$avg;
	    printf "St dev     = %.2f\n",$stdev;
	    print "Keep threshold  = $keep_threshold\n";
	}
	
	open OUT, ">$path/${db1}_BBH_${db2}";	

	# we want to print out because this is the multiparanoid input
	# OID     bits    'species'               'IPscore'  gene        boostrap
	# 1       1302    dros_sh_dfl_12.pep      1.000   S0_Dshi_0002    100%
	# 1       1302    bin05.pep       1.000   bin05_02438     100%

	# some of it we'll make up (like IPscore and bootstrap)

	my $index = 0;
	foreach my $p (sort {$a cmp $b} keys %BBH) {
	    if ($PtoDB->{$p} ne $db1) { next }
	    if (defined $BBH{$p}->{$db2}
		&& $HITS->{$p}->{$db2}->{'per_id'} >= $keep_threshold
		) {
		$index++;

		# multiParanoid input
		print OUT "$index\t$HITS->{$p}->{$db2}->{bits}\t$db1\t$HITS->{$p}->{$db2}->{score}\t$p\t100%\n";
		print OUT "$index\t$HITS->{$BBH{$p}->{$db2}}->{$db1}->{bits}\t$db2\t$HITS->{$BBH{$p}->{$db2}}->{$db1}->{score}\t$BBH{$p}->{$db2}\t100%\n";

	    }
	}
	close OUT;
    }
}

exit();

sub look_for_bbh {
    my @DB = @_;
    my %per_id_ref;
    #my %BBH;
    
    # Go through each protein
    foreach my $p (keys %$HITS) {
	my $db = $PtoDB->{$p};
	if (!$db) { die "Where did this protein come from? It's not in any of the parsed fasta files: $p\n";}

	# go through each db
	foreach my $qdb (@DB) {
	    # skip if there's no hits
	    if (! defined $HITS->{$p}->{$qdb}->{'tophit'}) { next }

	    my $q = $HITS->{$p}->{$qdb}->{'tophit'}; # ie top hit from other genome
	    # Check for unidirectional hits. 
	    if (! defined $HITS->{$q}->{$db}) {   # perhaps because length cutoff wasnt' met for a fusion protein...
		# Not sure why this code is here. Do we really want to include 
		# unidirectional hits?
		# print OUT "$p\t$q\t$HITS{$p}->{$qdb}->{per_id}\t!!! no reciprocal hit\n";
		# $BBH{$p} = {$q};
	    }
	    # This is the magic - when two proteins are reciprocal best matches
	    elsif ($HITS->{$q}->{$db}->{'tophit'} eq $p) {
		if ($BBH{$p}->{$qdb} eq $q && $BBH{$q}->{$db} eq $p) { next } # we already found this one
		# print "$p\t$q\t$HITS{$p}->{$qdb}->{bits}\n";
		# Add per_ids to array for AAAI calculation
		push @{$per_id_ref{$db}->{$qdb}}, ($HITS->{$p}->{$qdb}->{'per_id'}, $HITS->{$q}->{$db}->{'per_id'});
		push @{$per_id_ref{$qdb}->{$db}}, ($HITS->{$p}->{$qdb}->{'per_id'}, $HITS->{$q}->{$db}->{'per_id'});
		$BBH{$p}->{$qdb} = $q;
		$BBH{$q}->{$db} = $p;
	    }
	}
    }
    return (\%per_id_ref);
}

sub read_file {
    my $file = shift;
    my $PtoDB = shift;
    
    print STDERR "Reading file $file\n";
    open my $IN, $file or die "Can't open file $file : $!\n";

    my $header;
    my $last;
    while (my $line = <$IN>) {
	my ($qryacc, $refacc, $bits, $per_id, $qry_len, $ref_len);
	# Need the header stuff for MSOAR, if we're using that, which we aren't.
	if ($line =~ /^\#/) {
	    $header .= $line;
	    next;
	}
	chomp $line;
	my @f = split/\t/, $line;
	if ($line && @f == 14) {
	    # tab_blast data
	    $qryacc = id_from_accession($f[0]);
	    $refacc = id_from_accession($f[4]);
	    $bits = $f[10];
	    $per_id = $f[12];
	    $qry_len = $f[1];
	    $ref_len = $f[6];
	} elsif ($line && @f == 12) {
	    # mtab data
	    $qryacc = id_from_accession($f[0]);
	    $refacc = id_from_accession($f[1]);
	    $bits = $f[11] * 1;
	    $per_id = $f[2];
	    $qry_len = 1;
	    $ref_len = 1;
	} else {
	    die "Don't recognize format of tabulated blast output\n";
	}

	if (! defined $PtoDB->{$qryacc}) {
	    die "Hey, fool. The accession $f[0] ($qryacc) isn't in the fasta file. Huh?\n";
	}
	
	if (! defined $PtoDB->{$refacc}) {
	    die "Hey, fool. The accession $f[4] ($refacc) isn't in the fasta file. Huh?\n";
	}
	
	# self hit
	if ($qryacc eq $refacc) {
	    $HITS->{$qryacc}->{self_hit} = $bits;
	    next;
	}
	# paralogous hit (in same genome)
	elsif ($PtoDB->{$qryacc} eq $PtoDB->{$refacc}) { next }
	
	# This is the meat right here. Comparison against length and id cutoffs
	# populates the Data table.
	if (!defined $HITS->{$qryacc}->{$PtoDB->{$refacc}} &&  # this is the top hit to this genome
	    $per_id >= $id_cutoff &&
	    $ref_len/$qry_len >= $length_cutoff) {
	    # why do I calculate this score and then not use it? I guess this gets piped to paranoid, and it doesn't like non-1 values?
	    my $score = $HITS->{$qryacc}->{self_hit} ? $bits/$HITS->{$qryacc}->{self_hit} : 1;
	    $HITS->{$qryacc}->{$PtoDB->{$refacc}} = {'header' => $header,
							 'tophit' => $refacc,
							 'per_id' => $per_id,
							 'bits'   => $bits,
							 'score'  => 1,
							 'blast_data' => $line};
	}
	$header = "";
    }
    close $IN;
    return $HITS;
}

