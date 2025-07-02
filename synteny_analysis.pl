#!/usr/bin/perl
use Getopt::Std;
use lib $ENV{ORTHO};
use OrthologAnalysis;
use strict;
use Cwd;

my %arg;
&getopts('c:',\%arg);
my $bbhfile = $arg{c} ? $arg{c} : "multiparanoid.out";
if (!-e $bbhfile) { die "$bbhfile can't be found in " . cwd() . "\n"; }
my $outfile = "synteny_ortho.out";
open(my $OUT, ">$outfile") or die "Can't open $outfile for write: $!\n";

# check each gene in a gff file for synteny family
# For each gene in a gff file, see if it is a member of an MCL
# If it's not, then we move on to the next and report no synteny
# If it is, compare the target gene to all other genes in the cluster,
# looking for primarily BBH, but secondarily MCL membership of upstream
# and down stream genes.
# 

# holds the db, sequence, position for each gene
# ie $LOOKUP{gene} = {set => db,
#                     seq => seqid,
#                     idx => position (array index in GENES)}
my %LOOKUP;
# holds genes in positional order by sequence (in case multiple) and db
# ie $GENES->{db}->{seqid} = [gene1, gene2, gene3,...]
my $GENES = {};
my @db; 
for my $gfile (glob("../input/*.gff")) {
    my $db = &process_gff_file($gfile, $GENES, \%LOOKUP);
	push @db, $db;
}
# Alternative ways to get same info
# get sequences and genes from the sets and make ordered structures
# 
#    foreach my $set ($set1, $set2) {
#	if (-e "$set.info" && ! defined $INFO{$set}) {
#	    &process_info_file($set, $GENES, \%LOOKUP);
#	} elsif (-e "$+set.gtf" && ! defined $INFO{$set}) {
#	    &process_gtf_file($set, $GENES, \%LOOKUP);
#	} elsif (-e "$set.gff" && ! defined $INFO{$set}) {
#	    &process_gff_file($set, $GENES, \%LOOKUP);
#	} else {
#	    die "Need a .gff, .gtf or .info file for every genome: need one for $set1 or $set2";
#	}
#    }

# read in the multiparanoid BBH cluster output
my %GENECLUST;
my %CLUSTER;
open(my $clusterIN, $bbhfile) or die "Can't open $bbhfile for read: $!.\n";
while (my $line = <$clusterIN>) {
    if ($line =~ /^clusterID/) { next } # this is the header line
    chomp $line;
	my ($clust_id,
			$genome,
			$protein,
			$is_seed,
			$confidence,
			$all_species,
			$tree_conflict);

	my @clust = split(/\t/, $line);
	if (@clust == 2) { # this is mmseqs2 output
		($clust_id,
		 $protein) = @clust;
	} elsif (@clust == 7) { 	# this is multiparanoid output
	    ($clust_id,
			$genome,
			$protein,
			$is_seed,
			$confidence,
			$all_species,
			$tree_conflict) = @clust;
	}
    $GENECLUST{$protein} = $clust_id;
	push @{$CLUSTER{$clust_id}}, $protein;
}

my %SYN;
my %SYNGROUP;
my $syngroupcount;
my %GROUP;
my %INFO;

# read in a BBH file
# with format
# id, bits?, genome, (1)?, protein_id, confidence(100%)?
while (my $file = glob("BBH/*_BBH_*")) {
	print STDERR "Examining $file\n";
    my %BBH;
    # get the set names
    my ($set1, $set2);
    if ($file =~ /BBH\/(.+)\_BBH\_(.+)/) {
		$set1 = $1;
		$set2 = $2;
    } else { die "Set names can't be determined from filename $file. Needs to be set1_BBH_set2\n";}

    # why am I doing this if I have the BBH clusters in GENECLUST?
	# because geneclust isn't necessarily 1:1
    # read the bbhs into a structure
    open (my $bbh_in, $file) or die "Can't open $file: $!\n";
    while (my $line1 = <$bbh_in>) {
		my @f1 = split(/\s+/, $line1);
		my @acc1 = split(/\|/, $f1[4]);
		# many of the accessions have 'db|acc' structure, so make sure we get an acc by 
		my $acc1 = @acc1 > 1 ? $acc1[1] : $acc1[0];

		# lines come in pairs, so read in second line
		my $line2 = <$bbh_in>;
		my @f2 = split(/\s+/, $line2);
		my @acc2 = split(/\|/, $f2[4]);
		my $acc2 = @acc2 > 1 ? $acc2[1] : $acc2[0];

		# make sure the lines pair appropriately by the bbh id
		if ($f1[0] != $f2[0]) { die "Problem with bbh file $file:\n$line1$line2"; }
		# I guess we won't worry for now about identical accs between genomes
		$BBH{$set1}->{$acc1} = $acc2;
		$BBH{$set2}->{$acc2} = $acc1;
    }
    
    # g and h are the two genes in each BBH pair (only for these two sets)
    while (my ($g, $h) = each %{$BBH{$set1}}) {
		# see if the upstream and downstream genes are bbhs with 
		# the bbh's upstream and downstream genes
		# First, get the 2 upstream and 2 downstream genes for both g and h
		if (! $LOOKUP{$g}->{'set'}) { warn "!No set for $g\n";}
		if (! $LOOKUP{$h}->{'set'}) { warn "No set for $g\n";}
		my $usg1 = $GENES->{$LOOKUP{$g}->{'set'}}->{$LOOKUP{$g}->{'seq'}}->[$LOOKUP{$g}->{'idx'}+1];
		my $dsg1 = $GENES->{$LOOKUP{$g}->{'set'}}->{$LOOKUP{$g}->{'seq'}}->[$LOOKUP{$g}->{'idx'}-1];
		my $ush1 = $GENES->{$LOOKUP{$h}->{'set'}}->{$LOOKUP{$h}->{'seq'}}->[$LOOKUP{$h}->{'idx'}+1];
		my $dsh1 = $GENES->{$LOOKUP{$h}->{'set'}}->{$LOOKUP{$h}->{'seq'}}->[$LOOKUP{$h}->{'idx'}-1];
		my $usg2 = $GENES->{$LOOKUP{$g}->{'set'}}->{$LOOKUP{$g}->{'seq'}}->[$LOOKUP{$g}->{'idx'}+2];
		my $dsg2 = $GENES->{$LOOKUP{$g}->{'set'}}->{$LOOKUP{$g}->{'seq'}}->[$LOOKUP{$g}->{'idx'}-2];
		my $ush2 = $GENES->{$LOOKUP{$h}->{'set'}}->{$LOOKUP{$h}->{'seq'}}->[$LOOKUP{$h}->{'idx'}+2];
		my $dsh2 = $GENES->{$LOOKUP{$h}->{'set'}}->{$LOOKUP{$h}->{'seq'}}->[$LOOKUP{$h}->{'idx'}-2];
		#	print "$usg\t$g\t$dsg  <=>  $ush\t$h\t$dsh\n$BBH{$usg}\t$BBH{$g}\t$BBH{$dsg}  <=>  $BBH{$ush}\t$BBH{$h}\t$BBH{$dsh}\n\n" if (($BBH{$usg} && ($BBH{$usg} eq $ush || $BBH{$usg} eq $dsh)) || ($BBH{$dsg} && ($BBH{$dsg} eq $ush || $BBH{$dsg} eq $dsh)));

		# check to see if these have already been put in a synteny group
		# and if they are the same. If so, we are done.
		# If not, apparently merging didn't work so well
		if (defined $SYNGROUP{$g} && defined $SYNGROUP{$h}){
			if ($SYNGROUP{$g} == $SYNGROUP{$h}) {
				# Hooray! These have already been done.
			} else {
				# Hrm.
				# print STDERR "Merging groups $SYNGROUP{$g} and $SYNGROUP{$h}\n";
				# my $deletegroup = $SYNGROUP{$h};
				# foreach my $i (keys %{$GROUP{$deletegroup}}) {
				#     $SYNGROUP{$i} = $SYNGROUP{$g};
				#     $GROUP{$SYNGROUP{$g}}->{$i} = 1;
				# }
				# delete $GROUP{$deletegroup};
			}
		} else {
			# Let's see if these guys are syntenic.
			# we will score all combinations, to account for (single) indels and reverse orientation
			my $syn_score = 0;
			foreach my $sg ($usg1, $dsg1, $usg2, $dsg2) {
				foreach my $sh ($ush1, $dsh1, $ush2, $dsh2) {
					# two points for a match to the paired BBH
					$syn_score += 2 if ($BBH{$set1}->{$sg} && ($BBH{$set1}->{$sg} eq $sh));
					# a bonus point if the BBH is a part of a global cluster.
#					$syn_score += 1 if ($GENECLUST{$set1}->{$sg} && ($GENECLUST{$set1}->{$sg} eq $GENECLUST{$set2}->{$sh}));
					$syn_score += 1 if ($GENECLUST{$sg} && ($GENECLUST{$sg} eq $GENECLUST{$sh}));
				}
			}

			if ($syn_score > 1) {
				# so now, if one is part of a group, we add the other to the group
				# otherwise, we start a new syngroup
				if (defined $SYNGROUP{$g}) {
					$SYNGROUP{$h} = $SYNGROUP{$g};
					$GROUP{$SYNGROUP{$g}}->{$set1} = $g;
					$GROUP{$SYNGROUP{$h}}->{$set2} = $h;
				} elsif (defined $SYNGROUP{$h}) {
					$SYNGROUP{$g} = $SYNGROUP{$h};
					$GROUP{$SYNGROUP{$h}}->{$set2} = $h;
					$GROUP{$SYNGROUP{$g}}->{$set1} = $g;
				} else {
					$syngroupcount++;
					$SYNGROUP{$g} = $syngroupcount;
					$GROUP{$SYNGROUP{$g}}->{$set1} = $g;
					$SYNGROUP{$h} = $SYNGROUP{$g};
					$GROUP{$SYNGROUP{$h}}->{$set2} = $h;
				}
			}
		}
    }
}

# set for print to output file
select $OUT;

# print a header row
print "Id\t" . join("\t", sort keys %$GENES) . "\n";

foreach my $i (sort {$a<=>$b} keys %GROUP) {
    print "$i";
    foreach my $g(sort keys %$GENES) {
		print "\t" . $GROUP{$i}->{$g};
    }
    #    print "$g\t";
    #    &explode_print(\%SYN, $g);
    print "\n";
}

exit();

sub explode_print {
    my $href = shift;
    my $g = shift;
    foreach my $h (@{$href->{$g}->{2}}) {
	print "$h\t";
	&explode_print($href, $h);
    }
}

sub process_info_file {
    my $set = shift;
    my $GENES = shift;
    my $LOOKUP = shift;
    
    my $file = $set . ".info";
    open(my $info, $file) or die "Can't open $file for read: $!\n";
    my %D;
    while (my $line = <$info>) {
	chomp; $line;
	my ($racc, $lo, $dir, $seq_acc) = split/\s+/, $line;
	my @acc = split/\|/, $racc;
	my $acc = @acc > 1 ? $acc[1] : $acc[0];
	$D{$acc} = { 'posn' => $lo,
			 'dir' => $dir,
			 'seq' => $seq_acc};
    }
    foreach my $racc (sort {$D{$a}->{seq} cmp $D{$b}->{seq} || $D{$a}->{posn} <=> $D{$b}->{posn}} keys %D) {
	my @acc = split/\|/, $racc;
	my $acc = @acc > 1 ? $acc[1] : $acc[0];
	push @{$GENES->{$set}->{$D{$acc}->{seq}}}, $acc;
	$LOOKUP->{$acc} = {'set' => $set,
			       'seq' => $D{$acc}->{seq},
			       'idx' => $#{$GENES->{$set}->{$D{$acc}->{seq}}}};
    }
}

sub process_gtf_file {
    my $set = shift;
    my $GENES = shift;
    my $LOOKUP = shift;

    my $file = $set . ".gtf";
    open(my $gtf, $set . ".gtf") or die "Can't open $file for read: $!\n";
    my %D;
    while (my $line = <$gtf>) {
	chomp $line;
	my ($seq_acc,
	    $source,
	    $type,
	    $start,
	    $end,
	    $dot,
	    $strand,
	    $frame,
	    $notes) = split(/\t/, $line);
	if ($type ne "CDS") { next }
	my @n = split(/;/, $notes);
	my $racc; # for 'raw accession'
	if ($n[0] =~ /ID=(.*)/) { $racc = $1 }
	else { warn "Wait, why is the ID the first field in '$notes'?\n" }

	# clean up the accession
	my @acc = split(/\|/, $racc);
	my $acc = @acc > 1 ? $acc[1] : $acc[0];
	$D{$acc} = { 'posn' => $start,
			 'dir' => $strand,
			 'seq' => $seq_acc};
    }
    foreach my $acc (sort {$D{$a}->{seq} cmp $D{$b}->{seq} || $D{$a}->{posn} <=> $D{$b}->{posn}} keys %D) {
	# build an array of the accessions in order, stored by seq_acc and genome
	push @{$GENES->{$set}->{$D{$acc}->{seq}}}, $acc;
	# and a lookup of the genome, seq_acc, and position index by accession
	$LOOKUP->{$acc} = {'set' => $set,
			       'seq' => $D{$acc}->{seq},
			       'idx' => $#{$GENES->{$set}->{$D{$acc}->{seq}}}};
    }
    
}

sub process_gff_file {
    # this goes through the gff file and builds a structure
    # by db and contig that holds all gene accs in order of position
    # And a structure allowing lookup of db and contig by acc
    
    my $file = shift;
    my $GENES = shift;
    my $LOOKUP = shift;
    print STDERR "Processing $file\n";

    $file =~ /(.*\/)?(.*)\.gff$/;
    my $set = $2; # = db

    if (! $set) { die "Can't get db/set from $file\n"; }
    
    open(my $gff, $file) or die "Can't open $file for read: $!\n";
    my %D;
    while (my $line = <$gff>) {
	chomp $line;
	my ($seq_acc,
	    $source,
	    $type,
	    $start,
	    $end,
	    $dot,
	    $strand,
	    $frame,
	    $notes) = split(/\t/, $line);
	if ($type ne "CDS") { next }
	my @n = split(/;/, $notes);

	# Enforce that fasta acc has to match ID field
	my $racc; # for 'raw accession'
	if ($n[0] =~ /ID=(.*)/) { $racc = $1; }
	else { die "What is up with your gff file notes field? Where is 'ID'? '$notes'\n"; }
	
	if (! $racc) { die "Couldn't find an accession in $notes\n"; }

	# clean up the accession
	my $acc = id_from_accession($racc);
	$D{$acc} = { 'posn' => $start,
			 'dir' => $strand,
			 'seq' => $seq_acc};
    }
    foreach my $acc (sort {$D{$a}->{seq} cmp $D{$b}->{seq} || $D{$a}->{posn} <=> $D{$b}->{posn}} keys %D) {
		# build an array of the accessions in order, stored by seq_acc and genome
		push @{$GENES->{$set}->{$D{$acc}->{seq}}}, $acc;
		# and a lookup of the genome, seq_acc, and position index by accession
		$LOOKUP->{$acc} = {'set' => $set,
					'seq' => $D{$acc}->{seq},
			        'idx' => $#{$GENES->{$set}->{$D{$acc}->{seq}}}};
    }
    return ($set)
}
