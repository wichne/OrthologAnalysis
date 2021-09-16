#!/usr/bin/perl
#This is a script that makes output from the 
use strict;
use Getopt::Long;
$| = 1;

our (@dbs, $prefix, %PtoDB, %DONE);
my $bbhfile = "multiparanoid.out";
my $cluster_file = "geneFamily";
my $synfile = "synteny_ortho.out";
my $equiv_file = "all.equivalog.tblout";
my $nameglob = "gff";

&GetOptions(
    "namebase=s" => \$prefix,
    "bbh=s" => \$bbhfile,
    "syn=s" => \$synfile,
    "mcl=s" => \$cluster_file,
    "hmm=s" => \$equiv_file,
    "glob=s" => \$nameglob
    );

for my $f (glob("*.$nameglob")) { $f =~ s/\.[^\.]+$//; push @dbs, $f }
if (! @dbs) { die "Need to have a $nameglob file for each genome\n" }
$prefix or die "Please provide an orthlog family prefix with -n\n";
$bbhfile or die "Please provide multiparanoid.out file with -bbh";
#$synfile or die "Please provide synortho.out with -syn";
$cluster_file or die "Please provide out.mcl.all.I40 with -cluster";

my ($PtoB, $BtoP) = &read_BBH($bbhfile);
my ($PtoS, $StoP) = &read_SYN($synfile) if $synfile;
my ($PtoM, $MtoP) = &read_MCL($cluster_file);
my ($PtoH, $HtoP) = &read_HMM($equiv_file) if $equiv_file;

open(my $OUT, ">$prefix.ortho_out") or die "Can't open $prefix.ortho_out for write: $!\n";
select $OUT;

print "MCLFamID\tBBHFamID\tSynFamID\tHMMFam\tGenomeCount\tGeneCount\tProduct\tGene\tRole\tEC\tTC";
foreach my $d(@dbs) { print "\t$d" }
print "\n";

# go through each dataset
foreach my $db (@dbs) {
    print STDERR "Doing $db...\n";
    # to make sure we include every protein, go back to the original fasta file
    # and get every protein accession
    my $filepath;
    if (-e "input/$db.$nameglob") {
	$filepath = "input/$db.$nameglob";
    } else {
	my $in;
	until (-r $filepath) {
	    print STDERR "Enter the path to the protein fasta file used to build the orthofams for $db: ";
	    $filepath = <STDIN>;
	}
    }
    # to optimize analysis, we will first go through all the genes in a genome and count up how much
    # orthology evidence each gene has. Then we will go through all the genes in order of decreasing
    # evidence, marking 'DONE' genes that we have assigned to families
    my %ORTHEV;
    open(my $in, $filepath) or die "Can't open '$filepath' for read: $!\n";
    while (my $line = <$in>) {
	next if $line =~ /^#/;
	chomp $line;
	my $acc;
	if ($line =~ /^>\s*(\S+)\s*/) {
	    my $locus = $1;
	    # hack, hack, hack to make synteny work
	    $acc = &acc_from_locus($locus);
	} else {
	    # something in here for gff files
	    my @f = split(/\t/, $line);
	    if (@f == 9) {
		next unless $f[2] eq "CDS";
		my @n = split(";", $f[8]);
		foreach my $note (@n) {
		    if ($note =~ /locus_tag=(.+)/) { $acc = $1 }
		    elsif (! $acc && $note =~ /ID=(.+)/) { $acc = $1; $acc =~ s/cds-//; }
		    elsif (! $acc && $note =~ /Name=(.+)/) { $acc = $1 }
		    elsif ($note =~ /protein_id=(.+)/) { $acc = $1; last; }
		}
	    }
	}
	if (! $acc) { die "Sigh. No accession from $line\n"; }
	next if ($DONE{$acc});

	# pull the ortholog information
	my ($bbhfam, $mclfam, $synfam, $hmmfam, $orthflag);
	if (defined ($PtoB->{$acc})) { $bbhfam = $PtoB->{$acc}; $orthflag += 1; }
	if (defined ($PtoM->{$acc})) { $mclfam = $PtoM->{$acc}; $orthflag += 1; }
	if (defined ($PtoH->{$acc})) { $hmmfam = $PtoH->{$acc}; $orthflag += 1; }
	if (defined ($PtoS->{$acc})) { $synfam = $PtoS->{$acc}; $orthflag += 1; }
	$ORTHEV{$acc} = $orthflag;
    }

    for my $acc (sort {$ORTHEV{$b} <=> $ORTHEV{$a}} keys %ORTHEV) {
	# pull the annotation
	my ($product, $gene, $role, $ec, $tc) = &get_annotation($db, $acc, $nameglob);
	# pull the ortholog information
	my ($bbhfam, $mclfam, $synfam, $hmmfam, $orthflag);
	if (defined ($PtoB->{$acc})) { $bbhfam = $PtoB->{$acc}; $orthflag += 1; }
	if (defined ($PtoM->{$acc})) { $mclfam = $PtoM->{$acc}; $orthflag += 1; }
	if (defined ($PtoH->{$acc})) { $hmmfam = $PtoH->{$acc}; $orthflag += 1; }
	if (defined ($PtoS->{$acc})) { $synfam = $PtoS->{$acc}; $orthflag += 1; }
	
	# GenomeCount - number of genomes with an ortholog
	# must be at least 1
	my $gcount = 0;
	my $ocount = 0;
	my %genes; # holds the genes for the ortholog family
	
	# in this version, we will score potential orthologs based on their family membership
	if ($orthflag) {
	    my %POT_ORTH;
	    foreach my $d (@dbs) {
		# start with the least restrictive: mcl
		foreach my $pacc (@{$MtoP->{$mclfam}->{$d}}) {
		    $POT_ORTH{$d}->{$pacc} += 0.5 unless ($DONE{$pacc}); # a low weight for the non-specific assignment
		}
		
		# now the bbh fam
		foreach my $pacc (@{$BtoP->{$bbhfam}->{$d}}) {
		    $POT_ORTH{$d}->{$pacc} += 1 unless ($DONE{$pacc});
		}
		
		# now the hmm fam
		foreach my $pacc (@{$HtoP->{$hmmfam}->{$d}}) {
		    $POT_ORTH{$d}->{$pacc} += 1 unless ($DONE{$pacc});
		}
		
		# now the syn fam
		foreach my $pacc (@{$StoP->{$synfam}->{$d}}) {
		    $POT_ORTH{$d}->{$pacc} += 1 unless ($DONE{$pacc});
		}
	    }
	    
	    # now go through the potentials and get the best ones. If there's a tie, print them both
	    foreach my $d (@dbs) {
		my @best = &get_best($POT_ORTH{$d});
		if (@best) {
		    $genes{$d} = join(" ", @best);
		    $DONE{$_} = 1 foreach (@best);
		    $gcount++;
		    $ocount += scalar(@best);
		    if (! $product ||
			$product =~ /^\#/ ||
			$product =~ /hypothetical protein/ ||
			$product =~ /^\s*bin\=[^;]+;\s*$/) {
			($product, $gene, $role, $ec, $tc) = &get_annotation($d, $best[0], $nameglob);
		    }
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
	print join("\t",($mclfam,
			 $bbhfam,
			 $synfam,
			 $hmmfam,
			 $gcount,
			 $ocount,
			 $product,
			 $gene,
			 $role,
			 $ec,
			 $tc,
		   ));
	for my $d (@dbs) { print "\t$genes{$d}"; } 
	print "\n";
    }
}

exit();

sub get_annotation {
    my $db = shift;
    my $locus = shift;
    my $nameglob = shift;
    
    my ($product, $gene, $ec, $tc);
    my $file = "input/$db.$nameglob";
    if (! -e $file)
	#    if (-e "$db.gff") { $file = "$db.gff" }
	#    elsif (-e "$db.gtf") { $file = "$db.gtf" }
	#    elsif (-e "$db.pep") { $file = "$db.pep" }
	#    elsif (-e "$db.faa") { $file == "$db.faa" }
	#    else
    { warn "Couldn't find $file to pull annotation from for $db."; return() }
    
    my $line = `grep \"$locus" $file`;
    chomp $line;
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
	}
    }
    if (defined $ec && $ec =~ /[A-Z]/) { $tc = $ec; $ec = ""; }
    return ($product, $gene, $ec, $tc);
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
	my $acc = acc_from_locus($locus);
	if (defined ($BtoP{$orthofam}->{$db})) {
	    warn "$acc: Orthofam $orthofam already has a member(s) from $db: " . join(" ", @{$BtoP{$orthofam}->{$db}}) . ".\n";
	}
	push @{$BtoP{$orthofam}->{$db}}, $acc;
	$PtoB{$acc} = $orthofam;
	$PtoDB{$acc} = $db;
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
	    push @{$StoP{$famid}->{$db[$i]}}, $dg[$i];
	    $PtoS{$dg[$i]} = $famid;
	    $PtoDB{$dg[$i]} = $db[$i];
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
	    my $acc = acc_from_locus($f);
	    push @{$CtoP{$i}->{$PtoDB{$acc}}}, $f if (defined $PtoDB{$acc});
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
	my $acc = $f[0];
	$PtoH{$acc} = $hmm_acc;
	push @{$HtoP{$hmm_acc}->{$PtoDB{$acc}}}, $acc;
    }
    return (\%PtoH, \%HtoP);	    
}

sub acc_from_locus {
    my $locus = shift;
    my @locus = split/\|/, $locus;
    my $acc = @locus > 1 ? $locus[1] : $locus[0];
    return $acc;
}

sub get_best {
    my $ref = shift;
    my @best = ();
    my $top_score = 1;
    while (my ($this, $score) = each %$ref) {
	if ($score > $top_score) { @best = ($this); $top_score = $score; }
	elsif ($score == $top_score) { push @best, $this; }
    }
    return (@best);
}
