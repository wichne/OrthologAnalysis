#!/usr/bin/perl
#This is a script that makes output from the intermediate results files of the pipeline
use strict;
use Getopt::Long;
use lib $ENV{ORTHO};
use OrthologAnalysis;

$| = 1;

our (@dbs, $prefix, %DONE);
my $bbhfile = "multiparanoid.out";
my $cluster_file = "geneFamily";
my $synfile = "synteny_ortho.out";
my $equiv_file = "all.equivalog.tblout";
my $nameglob = "gff";
my $path = "input";
my $faext = "faa";

&GetOptions(
    "namebase=s" => \$prefix,
    "bbh=s" => \$bbhfile,
    "syn=s" => \$synfile,
    "mcl=s" => \$cluster_file,
    "hmm=s" => \$equiv_file,
    "glob=s" => \$nameglob,
    "path=s" => \$path,
    "faext=s" => \$faext,
    );

for my $f (glob("$path/*.$nameglob")) { $f =~ /$path\/(.*)\.$nameglob$/; push @dbs, $1 }
if (! @dbs) { die "Need to have a $nameglob file for each genome\n" }
$prefix or die "Please provide an orthlog family prefix with -n\n";
$bbhfile or die "Please provide multiparanoid.out file with -bbh";
#$synfile or die "Please provide synortho.out with -syn";
$cluster_file or die "Please provide out.mcl.all.I40 with -cluster";

# first we need to link all the protein identifiers with their databases
# go through each dataset and build PtoDB
my ($PtoDB, $DBtoP, $PtoANN) = &parse_gff_files($path, $nameglob);

print STDERR "Parsing BBH from '$path/$bbhfile'...\n";
my ($PtoB, $BtoP) = &read_BBH("$path/$bbhfile");
my ($PtoS, $StoP);
if ($synfile) {
    print STDERR "Parsing synteny from '$path/$synfile'...\n";
    ($PtoS, $StoP) = &read_SYN("$path/$synfile");
}
print STDERR "Parsing MCL clusters from '$path/$cluster_file'...\n";
my ($PtoM, $MtoP) = &read_MCL("$path/$cluster_file");
my ($PtoH, $HtoP);
if ($equiv_file) {
    print STDERR "Parsing equivalogs from '$path/$equiv_file'...\n";
    ($PtoH, $HtoP) = &read_HMM("$path/$equiv_file");
}

# Instantiate output file
open(my $OUT, ">$prefix.ortho_out") or die "Can't open $prefix.ortho_out for write: $!\n";
print $OUT "MCLFamID\tBBHFamID\tSynFamID\tHMMFam\tGenomeCount\tGeneCount\tProduct\tGene\tRole\tEC\tTC";
foreach my $d(@dbs) { print $OUT "\t$d" }
print $OUT "\n";

# Instantiate LOG file for testing
open (LOG, ">$prefix.$$.log");

# go through each dataset
foreach my $db (@dbs) {
    if ($db eq "all") { next }
    print STDERR "Finding orthologs for $db...";
    my $protCount = 0;

    my $filepath;
    if (-e "$path/$db.$nameglob") {
		$filepath = "$path/$db.$nameglob";
    }

    # to optimize analysis, we will first go through all the genes in a genome and count up how much
    # orthology evidence each gene has. 
    my ($with_MCL, $with_HMM, $with_BBH, $with_SYN);
    my %ORTHEV;
    foreach my $acc (sort {$a cmp $b} @{$DBtoP->{$db}}) {

		$ORTHEV{$acc} = 0;
		$protCount++;
		next if ($DONE{$acc});

		# count up the ortholog evidences applying weights that will sum to unique values
		if (defined ($PtoM->{$acc})) { $ORTHEV{$acc} += 1; $with_MCL++; }
		if (defined ($PtoB->{$acc})) { $ORTHEV{$acc} += 2; $with_BBH++; }
		if (defined ($PtoH->{$acc})) { $ORTHEV{$acc} += 4; $with_HMM++; }
		if (defined ($PtoS->{$acc})) { $ORTHEV{$acc} += 8; $with_SYN++; }

    }
    print STDERR "$protCount proteins found\nwith MCL: $with_MCL\nwith BBH: $with_BBH\nwith HMM: $with_HMM\nwith SYN: $with_SYN\n";

	# Now we will go through all the genes in order of decreasing
    # evidence, marking 'DONE' genes that get assigned to families
    for my $acc (sort {$ORTHEV{$b} <=> $ORTHEV{$a} || $a cmp $b } keys %ORTHEV) {
		# pull the annotation
		if (! $acc) { die "Why is there no accession?"}
		if ($DONE{$acc}) { next } # need this here, too to avoid double-printing
		
		my $orth_annotation = $PtoANN->{$acc};
		# pull the ortholog information
		my ($bbhfam, $mclfam, $synfam, $hmmfam);
		if (defined ($PtoB->{$acc})) { $bbhfam = $PtoB->{$acc}; }
		if (defined ($PtoM->{$acc})) { $mclfam = $PtoM->{$acc}; }
		if (defined ($PtoH->{$acc})) { $hmmfam = $PtoH->{$acc}; }
		if (defined ($PtoS->{$acc})) { $synfam = $PtoS->{$acc}; }
		
		# GenomeCount - number of genomes with an ortholog
		# must be at least 1
		my $gcount = 0;
		my $ocount = 0;
		my %genes; # holds the genes for the ortholog family
		
		print LOG "looking for orthologs to $acc (score = $ORTHEV{$acc})\n";
		# in this version, we will score potential orthologs based on their family membership
		if ($ORTHEV{$acc}) {
			my %POT_ORTH;
			# sum up evidence for each member of each ev-type
			foreach my $d (@dbs) {
				# start with the least restrictive: mcl
				foreach my $pacc (@{$MtoP->{$mclfam}->{$d}}) {
					$POT_ORTH{$d}->{$pacc} += 1 unless ($DONE{$pacc}); # a low weight for the non-specific assignment
				}
				
				# now the bbh fam
				foreach my $pacc (@{$BtoP->{$bbhfam}->{$d}}) {
					$POT_ORTH{$d}->{$pacc} += 2 unless ($DONE{$pacc});
				}
				
				# now the hmm fam
				foreach my $pacc (@{$HtoP->{$hmmfam}->{$d}}) {
					$POT_ORTH{$d}->{$pacc} += 4 unless ($DONE{$pacc});
				}
				
				# now the syn fam
				foreach my $pacc (@{$StoP->{$synfam}->{$d}}) {
					$POT_ORTH{$d}->{$pacc} += 8 unless ($DONE{$pacc});
				}
			}
			# now go through the potentials and get the best ones. If there's a tie, print them both
			foreach my $d (@dbs) {
				print LOG "\tall potentials for $d: " . join(",", keys(%{$POT_ORTH{$d}})) . "\n";
				my @best = &get_best($POT_ORTH{$d}, $ORTHEV{$acc});
				if (@best) {
					my @ids;
					foreach my $acc(@best) {
						print LOG "\tfound potential $d|$acc ($POT_ORTH{$d}->{$acc})\n";
						if ($PtoANN->{$acc}->{'locus_tag'}) {
							push @ids, $PtoANN->{$acc}->{'locus_tag'};
						} else { push @ids, $acc }
						if (! $orth_annotation ||
							$orth_annotation->{'product'} =~ /^\#/ ||
							$orth_annotation->{'product'} =~ /hypothetical protein/ ||
							$orth_annotation->{'product'} =~ /^\s*bin\=[^;]+;\s*$/) {
							$orth_annotation = $PtoANN->{$acc};
						}
					}
					$genes{$d} = join(" ", @ids);
					$DONE{$_} = 1 foreach (@best);
					$gcount++;
					$ocount += scalar(@best);
				} else {
					$genes{$d} = "";
				}
	    	}
		}
	
		if ($gcount > 0 && ! grep($acc, $genes{$db})) {
			warn "Hmm. The fam for $acc didn't include $acc. $hmmfam, mcl $mclfam, bbh $bbhfam, syn $synfam.\n";
		} elsif ($gcount == 0) {
			$genes{$db} = $acc;
			$DONE{$acc} = 1;
			$gcount++;
			$ocount++;
		}
	
		#$product =~ s/\"//g;
		print $OUT join("\t",($mclfam,
					$bbhfam,
					$synfam,
					$hmmfam,
					$gcount,
					$ocount,
					$orth_annotation->{'product'},
					$orth_annotation->{'gene'},
					$orth_annotation->{'role'},
					$orth_annotation->{'ec'},
					$orth_annotation->{'tc'},
				));
		for my $d (@dbs) { 
			print $OUT "\t$genes{$d}";
		} 
		print $OUT "\n";
    }
}

exit();

sub get_annotation {
    my $path = shift;
    my $db = shift;
    my $locus = shift;
    my $nameglob = shift;
    
    my ($product, $gene, $ec, $tc, $locus_tag);
    my $file = "$path/$db.$nameglob";
    if (! -e $file)
    { warn "Couldn't find $file in $path to pull annotation from for $db."; return() }

    my $line;
    if ($nameglob eq "gff") {
	$line = `grep "$locus;" $file | grep -w CDS`;
    } else {
	$line = `grep -w "$locus" $file`;
    }
    chomp $line;
    my @lines = split(/\n/, $line);
    if (@lines > 1) { warn "!!!!!!!!Found multiple annotation lines in $file for $locus:\n" . join("\n", @lines) . "\n"; }
    
    # if it's a fasta file...
    if ($line =~ />\s*\S+\s+(.*)/) {
	my $desc = $1;
	if ($desc =~ /product\=([^\;]+)/) {
	    $product = $1;
	    pos($desc) = 0;
	}
	if ($desc =~ /(\[EC\:(\d+(\.(\d+|\-)){3}[a-z]?)\])/ ||
	    $desc =~ /\s*(\[?EC\_number\=(\d+(\.(\d+|\-)){3}[a-z]?)\]?)/) {
	    $ec = $2;
	    my $r = $1;
	    $desc =~ s/\Q$r\E//;
	    pos($desc) = 0;
	}
	if ($desc =~ /\[?(gene\=\"?([^\"\;]+)\"?\]?)/) {
	    $gene = $2;
	    my $r = $1;
	    $desc =~ s/$r//;
	    pos($desc) = 0;
	}
	$desc =~ s/[\[\{]organism=.*[\}\]]//;
	$product = $desc if (! $product);
    }
    # assume it's a gff file
    else {
	my @l = split(/\t/, $line);
	my @tags = split(/\;/, $l[8]);
	foreach my $t (@tags) {
	    my ($k, $v) = split(/\s*=\s*/, $t);
	    if ($k eq "product") { $product = $v }
	    if ($k eq "EC_num") { $ec = $v }
	    if ($k eq "gene_sym") { $gene = $v }
	    if ($k eq "locus_tag") { $locus_tag = $v }
	}
    }
    if (defined $ec && $ec =~ /[A-Z]/) { $tc = $ec; $ec = ""; }
    return ($product, $gene, $ec, $tc, $locus_tag);
}

sub read_BBH {
    my $infile = shift;
    # open file
    open my $in, $infile or die "Can't open '$infile': $!\n";

    # first line is header
    my $line = <$in>;

    # read in rest of input file of ortholog families
    # P[rotein]toO[rthofam], O[rthofam]toP[rotein]
    our (%PtoB, %BtoP);
    while ($line = <$in>) {
	chomp $line;
	my $parafam;
	my ($orthofam, $db, $locus, $seed, $conf, $all_db, $tree_conflict) = split/\t/, $line;
	my $acc = id_from_accession($locus);
	if (defined ($BtoP{$orthofam}->{$db})) {
	    warn "Orthofam:$orthofam - $db (" . join(" ", @{$BtoP{$orthofam}->{$db}}) . ") add $acc.\n";
	}
	push @{$BtoP{$orthofam}->{$db}}, $acc;
	$PtoB{$acc} = $orthofam;
    }
    return(\%PtoB, \%BtoP);
}

sub read_SYN {
    my $infile = shift;
    open my $in, $infile or die "Can't open $infile for read: $!\n";
    # first line defines columns
    chomp(my $header = <$in>);
    my @db = split/\t/, $header;
    #shift @db; # get rid of first element which is just 'Id'
    
    my %PtoS;
    my %StoP;
    while (my $line = <$in>) {
	chomp $line;
	my @dg = split/\t/, $line;
	my $famid = $dg[0];
	for (my $i=1;$i<@dg;$i++) {
	    if (! $dg[$i]) { next }
	    my $acc = id_from_accession($dg[$i]);
	    push @{$StoP{$famid}->{$db[$i]}}, $acc;
	    $PtoS{$acc} = $famid;
	}
    }
    return (\%PtoS, \%StoP);
}

sub read_MCL {
    my $infile = shift;

    # read in clusters for superfamily names
    open my $in, $infile or die "Can't open $infile for read: $!\n";
    my %PtoC;
    my %CtoP;
    my $i = 0;
    while (my $line = <$in>) {
	$i++;
	chomp $line;
	my @f = split /\t/, $line;
	foreach my $f(@f) {
	    my $acc = id_from_accession($f);
	    push @{$CtoP{$i}->{$PtoDB->{$acc}}}, $acc if (defined $PtoDB->{$acc});
	    $PtoC{$acc} = $i; # link cluster id to gene
	}
    }
    return \%PtoC, \%CtoP;
}

sub read_HMM {
    my $infile = shift;
    open my $in, $infile or die "Can't open HMM file $infile: $!\n";
    my %PtoH;
    my %HtoP;
    while (my $line = <$in>) {
	next if ($line =~ /^#/);
	chomp $line;
	my @f = split(/\s+/, $line, 19);
	my $hmm_acc = $f[3];
	my $acc = id_from_accession($f[0]);
	$PtoH{$acc} = $hmm_acc;
	push @{$HtoP{$hmm_acc}->{$PtoDB->{$acc}}}, $acc if $PtoDB->{$acc};
    }
    return (\%PtoH, \%HtoP);	    
}

sub get_best {
    # We want this to return the best match to the input score
    # The problem is we want synteny (8) only to tie-break
    my $ref = shift;
    my $targetScore = shift;

	# get the best score
	my $bestscore = (sort {$b <=> $a} values %$ref)[0];

	# MCL alone is not enough
	if ($bestscore == 1) { return }

    my @best = ();
	# any candidate that matches the best score is returned
    foreach my $this (sort {$ref->{$b} <=> $ref->{$a}} keys %$ref) {
		if ($ref->{$this} == $bestscore) {
			push @best, $this;
		}
	}
    return (@best);
}
