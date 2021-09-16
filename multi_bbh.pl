#!/usr/bin/perl
use strict;
use Getopt::Std;
use Bio::FeatureIO;
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

USAGE: multi_bbh.pl -i inputfile -n number_of_genomes \@info_files

-i - blast_tab.pl output from an all vs all blastp search.

";

my %arg;
&getopts('i:n:h', \%arg);

if ($arg{'h'}) { die $usage }
my $file = $arg{'i'} or die "Need to provide input btab file with -i\n";

my $VERBOSE = 1;

our $length_cutoff = 0.7;
our $id_cutoff = 30;

# still have the problem of assocating proteins with genomes, since it's not explicit in the
# pair-wise output. In this instance, we will use the fasta files, since they were the
# original input for the searches.
print STDERR "Parsing accessions into genomes....\n";
our %HITS; # Lookup table: prot acc -> genome; also holds blast hit info
our @DB;
foreach my $fa_file (glob("*.faa"), glob("*.fasta"), glob("*.pep")) {
    if ($fa_file =~ /^all[\.\_]/) {
		print STDERR "\t$fa_file looks like the agglomeration file, so skip it...\n";
		next;
    }

	# generate the db name
    my $db = $fa_file;
    $db =~ s/\.(pep|faa|fasta)$//;
    push @DB, $db;

    open(my $fa, $fa_file) or die "Can't open $fa_file for read: $!\n";
    while (my $h = <$fa>) {
		if ($h =~ />(\S+)/) {
			# my @acc = split("|", $1);
			# my $thisAcc = @acc > 1 ? $acc[1] : $acc[0];
			# $HITS{$thisAcc}->{'genome'} = $db;
			$HITS{$1}->{'genome'} = $db;
		}
	}
}
print STDERR "Done.\n";

&read_file($file); # populates %HITS

our %BBH;
my $per_id_ref = &look_for_bbh(@DB); # populates %BBH

my @DB = sort keys %$per_id_ref;

# Roll through the pair-wise comparisons
for (my $i=0; $i<@DB; $i++) {
    my $db1 = $DB[$i];
    print STDERR "Doing $db1...\n";
    for (my $j=$i; $j<@DB; $j++) {
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
		
		open OUT, ">${db1}.${db2}.ortho";	

		# we want to print out because this is the multiparanoid input
		# OID     bits    'species'               'IPscore'  gene        boostrap
		# 1       1302    dros_sh_dfl_12.pep      1.000   S0_Dshi_0002    100%
		# 1       1302    bin05.pep       1.000   bin05_02438     100%

		# some of it we'll make up (like IPscore and bootstrap)

		my $index = 0;
		foreach my $p (sort {$a cmp $b} keys %BBH) {
			if ($HITS{$p}->{'genome'} ne $db1) { next }
			if (defined $BBH{$p}->{$db2}
			     && $HITS{$p}->{$db2}->{'per_id'} >= $keep_threshold
			    ) {
				$index++;

				# multiParanoid input
				print OUT "$index\t$HITS{$p}->{$db2}->{bits}\t$db1\t$HITS{$p}->{$db2}->{score}\t$p\t100%\n";
				print OUT "$index\t$HITS{$BBH{$p}->{$db2}}->{$db1}->{bits}\t$db2\t$HITS{$BBH{$p}->{$db2}}->{$db1}->{score}\t$BBH{$p}->{$db2}\t100%\n";

			}
		}
		close OUT;
    }
}

exit();

sub look_for_bbh {
    my @DB = @_;
    my %per_id_ref;

    # Go through each protein
    foreach my $p (keys %HITS) {
		my $db = $HITS{$p}->{'genome'};

		# go through each db
		foreach my $qdb (@DB) {
			# skip if there's no hits
			if (! defined $HITS{$p}->{$qdb}->{'tophit'}) { next }

			my $q = $HITS{$p}->{$qdb}->{'tophit'}; # ie top hit from other genome
			# Check for unidirectional hits. 
			if (! defined $HITS{$q}->{$db}) {   # perhaps because length cutoff wasnt' met for a fusion protein...
				# Not sure why this code is here. Do we really want to include 
				# unidirectional hits?
				# print OUT "$p\t$q\t$HITS{$p}->{$qdb}->{per_id}\t!!! no reciprocal hit\n";
				# $BBH{$p} = {$q};
			}
			# This is the magic - when two proteins are reciprocal best matches
			elsif ($HITS{$q}->{$db}->{'tophit'} eq $p) {
		        # print "$p\t$q\t$HITS{$p}->{$qdb}->{bits}\n";
				# Add per_ids to array for AAAI calculation
				push @{$per_id_ref{$db}->{$qdb}}, ($HITS{$p}->{$qdb}->{'per_id'}, $HITS{$q}->{$db}->{'per_id'});
				push @{$per_id_ref{$qdb}->{$db}}, ($HITS{$p}->{$qdb}->{'per_id'}, $HITS{$q}->{$db}->{'per_id'});
				$BBH{$p}->{$qdb} = $q;
				$BBH{$q}->{$db} = $p;
			}
		}
    }
    return (\%per_id_ref);
}

sub read_file {
    my $file = shift;
    print STDERR "Reading file $file\n";
    open my $IN, $file or die "Can't open file $file : $!\n";

    my $header;
    my $last;
    while (my $line = <$IN>) {
		# Need the header stuff for MSOAR
		if ($line =~ /^\#/) {
			$header .= $line;
			next;
		}
		chomp $line;
		my @f = split/\t/, $line;
		if ($line && @f == 14) { 
			# tab_blast data

			if ($f[0] eq $f[4]) {
				$HITS{$f[0]}->{self_hit} = $f[10];
				next;
			} # self hit
			# determine if these are from the same genome to ignore paralog hits
			if (! (defined $HITS{$f[0]}->{'genome'})) {
				if ($f[0] =~ /.*\|(.*)/) {
					if (defined $HITS{$1}->{'genome'}) { $f[0] = $1 }
				} else {
					die "Hey, fool. The accessions in your fasta file don't match the gff file. '$f[0]'\n" . join("' '", keys %HITS) ."\n";
				}
			}
			if (! (defined $HITS{$f[4]}->{'genome'})) {
				if ($f[4] =~ /.*\|(.*)/) {
					if (defined $HITS{$1}->{'genome'}) { $f[4] = $1 }
				} else {
					die "Hey, fool. The accessions in your fasta file don't match the .info file. $f[4]\n" . join(" ", keys %HITS) ."\n";
				}
	    	}
		   
			if ($HITS{$f[0]}->{'genome'} eq $HITS{$f[4]}->{'genome'}) { next }

			# This is the meat right here. Comparison against length and id cutoffs
			# populates the Data table.
			if (!defined $HITS{$f[0]}->{$HITS{$f[4]}->{'genome'}} &&  # this is the top hit to this genome
  			     $f[12] >= $id_cutoff &&
			     $f[6]/$f[1] >= $length_cutoff) {
				my $score = $HITS{$f[0]}->{self_hit} ? $f[10]/$HITS{$f[0]}->{self_hit} : 1;
				$HITS{$f[0]}->{$HITS{$f[4]}->{'genome'}} = {'header' => $header,
										'tophit' => $f[4],
										'per_id' => $f[12],
										'bits'   => $f[10],
										'score'  => 1,
										'blast_data' => $line};
			}
			$header = "";
		} elsif ($line && @f == 12) {
			# mtab data
			if ($f[0] eq $f[1]) { 
				$HITS{$f[0]}->{self_hit} = $f[11] * 1;
				next;
			} # self hit
			if ($HITS{$f[0]}->{'genome'} eq $HITS{$f[1]}->{'genome'}) { next }

			# for this data, we don't have the query and reference seq lengths,
			# so we'll use a bit score cutoff instead
			if (!defined $HITS{$f[0]}->{$HITS{$f[1]}->{'genome'}} &&  # this is the top hit to this genome
			     $f[2] >= $id_cutoff &&
			     $f[11] > $HITS{$f[0]}->{self_hit}/10) {
				my $score = $HITS{$f[0]}->{self_hit} ? $f[11]/$HITS{$f[0]}->{self_hit} : 1;
				$HITS{$f[0]}->{$HITS{$f[1]}->{'genome'}} = {'header' => $header,
										'tophit' => $f[1],
										'per_id' => $f[2],
										'bits'   => $f[11],
										'score'  => 1,
										'blast_data' => $line};
			}
			$header = "";
		} else {
			die "Don't recognize format of tabulated blast output\n";
		}
    }
    close $IN;
}

