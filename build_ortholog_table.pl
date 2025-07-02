#!/usr/bin/env perl
#This is a script that makes output from the intermediate results files of the pipeline
use strict;
use Getopt::Long;
use List::Util qw(max);
use lib $ENV{ORTHO};
use OrthologAnalysis;

$| = 1;

our (@dbs, $prefix, %DONE);
my $bbhfile = "multiparanoid.out";
my $cluster_file = "geneFamily";
my $synfile = "synteny_ortho.out";
my $equiv_file = "all.equivalog.tblout";
my $nameglob = "gff";
my $path = ".";
my $faext = "faa";
my $OrthoFinder_file = "Orthogroups.tsv";
my $mmseqs2_file = "mmseqs2.out";
my $orthomcl_file = "groups.txt";

&GetOptions(
    "namebase=s" => \$prefix,
    "bbh=s" => \$bbhfile,
    "syn=s" => \$synfile,
    "mcl=s" => \$cluster_file,
	"orthomcl=s" => \$orthomcl_file,
    "hmm=s" => \$equiv_file,
    "orthf=s" => \$OrthoFinder_file,
    "mmseqs=s" => \$mmseqs2_file,
    "glob=s" => \$nameglob,
    "path=s" => \$path,
    "faext=s" => \$faext,
    );

if (! $prefix) { $prefix = "ORTHOLOGS"; }
# Instantiate LOG file for testing
open (my $logh, ">>$prefix.$$.log");

for my $f (glob("$path/input/*.$nameglob")) { $f =~ /\/input\/(.*)\.$nameglob$/; push @dbs, $1 }
if (! @dbs) { die "Need to have a $nameglob file for each genome in $path/input/\n" }
#$prefix or die "Please provide an orthlog family prefix with --namebase\n";
#$bbhfile or die "Please provide multiparanoid.out file with --bbh";
#$synfile or die "Please provide synortho.out with -syn";
#$cluster_file or die "Please provide out.mcl.all.I40 with --mcl";

# first we need to link all the protein identifiers with their databases
# go through each dataset and build PtoDB
# The FASTA file records
my ($PtoFAA, $FAAtoP) = &parse_protein_fasta_files("input", $faext);
my ($PtoDB, $DBtoP, $PtoANN) = &parse_gff_files("input", $nameglob, $PtoFAA, $logh);

# Make sure all accessions in fasta file and gff are associated with a db,
# even if they are not present in both files.
foreach my $acc(keys %$PtoFAA) { 
	$PtoDB->{$acc} = $PtoFAA->{$acc} if (! defined $PtoDB->{$acc});
}

my %DB; # keeps track of how many genes from each db are placed into an ortholog family. should be all of them.
foreach my $db(keys %$FAAtoP) { $DB{$db}->{in} = scalar(@{$FAAtoP->{$db}}) }

# Now we read in the different data types
# each parser returns two structures: one with protid/famid, one with famid/db/protidarray
# BBH
my ($PtoB, $BtoP);
if (-e "$path/output/$bbhfile") {
    print $logh "Parsing BBH from '$path/output/$bbhfile'...\n";
    ($PtoB, $BtoP) = &read_BBH("$path/output/$bbhfile") if ($bbhfile);
} else { warn "Couldn't find BBH file $bbhfile "}

# synteny
my ($PtoS, $StoP);
if (-e "$path/output/$synfile") {
    print $logh "Parsing synteny from '$path/output/$synfile'...\n";
    ($PtoS, $StoP) = &read_SYN("$path/output/$synfile");
} else { warn "Couldn't find synteny file $synfile "}

# OrthoMCL
my ($PtoOM, $OMtoP);
if (-e "$path/output/$orthomcl_file") {
	print $logh "Parsing OrthoMCL clusters from '$path/output/$orthomcl_file'...\n";
	($PtoOM, $OMtoP) = &read_OrthoMCL("$path/output/$orthomcl_file");
} else { warn "Couldn't find OrthoMCL file $orthomcl_file"}

# MCL
my ($PtoM, $MtoP);
if (-e "$path/output/$cluster_file") {
    print $logh "Parsing MCL clusters from '$path/output/$cluster_file'...\n";
    ($PtoM, $MtoP) = &read_MCL("$path/output/$cluster_file") if ($cluster_file);
} else { warn "Couldn't find MCL file $cluster_file "}

# OrthoFinder
my($PtoOF, $OFtoP);
if (-e "$path/output/$OrthoFinder_file") {
    print $logh "Parsing OrthoFinder clusters from '$path/output/$OrthoFinder_file'...\n";
    ($PtoOF, $OFtoP) = &read_OrthoFinder("$path/output/$OrthoFinder_file");
} else { warn "Couldn't find OrthoFinder file $OrthoFinder_file "}

# mmseqs2
my($PtoMM, $MMtoP);
if (-e "$path/output/$mmseqs2_file") {
    print $logh "Parsing mmseqs clusters from '$path/output/$mmseqs2_file'...\n";
    ($PtoMM, $MMtoP) = &read_mmseqs2("$path/output/$mmseqs2_file", $PtoDB);
} else { warn "Couldn't find mmseqs file $mmseqs2_file "}

# Equivalogs
my ($PtoH, $HtoP);
if (-e "$path/output/$equiv_file") {
    print $logh "Parsing equivalogs from '$path/output/$equiv_file'...\n";
    ($PtoH, $HtoP) = &read_HMM("$path/output/$equiv_file");
} else { warn "Couldn't find TIGRfams file $equiv_file "}

# set up the output file
my $outfile = "$path/output/$prefix.ortho.tsv";
open(my $OUT, ">$outfile") or die "Can't open $outfile for write: $!\n";
$| =1;
print $OUT "MCLFamID\tOrthoMCL\tOrthoFinder\tBBHFamID\tMMSEQS2\tSynFamID\tHMMFam\tGenomeCount\tGeneCount\tProduct\tGene\tRole\tEC\tTC";
foreach my $d(@dbs) { print $OUT "\t$d" }
print $OUT "\n";

### First go through all families of all types and figure out a score
# I have decided for now that the formula is 2*d - p
# where d is the number of databases with a family member
# and p is number of proteins in the family
# this should balance between 1 member/database and the overall number of databases
my %COUNT;
foreach my $struct ($BtoP, $MtoP, $MMtoP, $OFtoP, $StoP, $HtoP, $OMtoP) {
    foreach my $famid (keys %$struct) {
		if (defined $COUNT{$famid}) { die "Oops. famid collision between types -> $famid\n";}
		my $thisdbn = scalar(keys(%{$struct->{$famid}}));
		my $thisprotn;
		for my $db (keys (%{$struct->{$famid}})) {
			my $protn = scalar(@{$struct->{$famid}->{$db}});
			if ($protn == 0) { $thisdbn-- } # because there are no members for this db
			else { $thisprotn += $protn }
			#			$COUNT{$famid} += scalar(@{$struct->{$famid}->{$db}})
		}
		$COUNT{$famid} = (2 * $thisdbn) - $thisprotn;
    }
}

# Now go by descending order of score
for my $famid(sort { $COUNT{$b} <=> $COUNT{$a} || $a cmp $b } keys %COUNT ) {
    # For now, mark the proteins done and print out the line
    my ($type, $struct) = which_type($famid);
    my ($gcount, $ocount, $thisfam, %ann, %genes);
    for my $db (keys (%{$struct})) {
		# keep track of how many members for this genome
		my $go_n = 0;
		foreach my $p (@{$struct->{$db}}) {
			if ($DONE{$p}) { next }
			$go_n++;
			push @{$genes{$db}}, $p;
			# determine common family assignments
			$thisfam = fam_intersect($thisfam, $p);
			if (! %ann ||
			$ann{'product'} =~ /^\#/ ||
			$ann{'product'} =~ /hypothetical protein/ ||
			$ann{'product'} =~ /^\s*bin\=[^;]+;\s*$/) {
				%ann = %{$PtoANN->{$p}} if ($p && defined $PtoANN->{$p});
			}
			$DB{$db}->{'out'}++;
			$DONE{$p} = 1;
		}
		# if members > 0, increment db and add members
		if ($go_n > 0) {
			$gcount++;
			$ocount += $go_n;
		}
    }
    # screen for only fam_ids with full membership
    foreach my $t (keys %$thisfam) {
		foreach my $f (keys %{$thisfam->{$t}}) {
			if ($thisfam->{$t}->{$f} != $ocount) { delete $thisfam->{$t}->{$f} }
		}
    }
    
    # print out
    print_row($gcount, $ocount, $thisfam, \%ann, \%genes) if ($ocount > 0);
}

foreach my $db(sort keys %DB) {
    print $logh "$db had ". $DB{$db}->{in} . " proteins and " . $DB{$db}->{out} . " are in ortho fams\n";
}
exit();

########################################################################
### Subroutines
########################################################################

sub which_type {
    my $famid = shift;
    my %opts = ('mclfam' => $MtoP, 'orthomcl' => $OMtoP, 'mm' => $MMtoP, 'og' => $OFtoP,
		'synfam' => $StoP, 'bbhfam' => $BtoP, 'hmmfam' => $HtoP);
    while (my ($type, $struct) = each %opts) {
		if (defined $struct->{$famid}) {
		    return $type, $struct->{$famid};
		}
    }
}

sub fam_intersect {
    my $fams = shift;
    my $p = shift;

    if (defined $PtoB->{$p}) { $fams->{bbhfam}->{$PtoB->{$p}}++ }
    if (defined $PtoS->{$p}) { $fams->{synfam}->{$PtoS->{$p}}++ }
    if (defined $PtoOF->{$p}) { $fams->{og}->{$PtoOF->{$p}}++ }
    if (defined $PtoMM->{$p}) { $fams->{mm}->{$PtoMM->{$p}}++ }
    if (defined $PtoH->{$p}) { $fams->{hmmfam}->{$PtoH->{$p}}++ }
    if (defined $PtoM->{$p}) { $fams->{mclfam}->{$PtoM->{$p}}++ }
    if (defined $PtoOM->{$p}) { $fams->{orthomcl}->{$PtoOM->{$p}}++ }

    return $fams;
}

sub print_row {
    my $gcount = shift;
    my $ocount = shift;
    my $fam = shift;
    my $ann = shift;
    my $genes = shift;
    print $OUT join("\t",
		    (join(" ", keys(%{$fam->{mclfam}})),
		     join(" ", keys(%{$fam->{orthomcl}})),
		     join(" ", keys(%{$fam->{og}})),
		     join(" ", keys(%{$fam->{bbhfam}})),
		     join(" ", keys(%{$fam->{mm}})),
		     join(" ", keys(%{$fam->{synfam}})),
		     join(" ", keys(%{$fam->{hmmfam}})),
		     $gcount,
		     $ocount,
		     $ann->{'product'},
		     $ann->{'gene'},
		     $ann->{'role'},
		     $ann->{'ec'},
		     $ann->{'tc'},
		    ));
    for my $d (@dbs) { 
		if (defined $genes->{$d}) {
			# here we want to be sure we are returning the same identifiers as were in the input fasta file.
			#my @locus;
			#for (my $i=0; $i<@{$genes->{$d}};$i++) {
			#if ($PtoANN->{$genes->{$d}->[$i]}->{'locus_tag'}) {
			#    $locus[$i] = $PtoANN->{$genes->{$d}->[$i]}->{'locus_tag'};
			#} else {
			#    $locus[$i] = $genes->{$d}->[$i];
			#}
			#}
			print $OUT "\t" . join(" ", @{$genes->{$d}}); #@locus);
		} else { print $OUT "\t";}
    }
    print $OUT "\n";
    
}
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
	my ($bbhidx, $db, $locus, $seed, $conf, $all_db, $tree_conflict) = split/\t/, $line;
	my $fam_id = sprintf("BBH%05d", $bbhidx);
	if (! $locus) { die "No locus from $line\n"; }
	my $acc = id_from_accession($locus);
	if (defined $PtoB{$acc}) { 
	    print $logh "ERROR: BBH: $fam_id: protein $acc has already been assigned to $PtoB{$acc}\n";
	    next; 
	}
	if (defined ($BtoP{$fam_id}->{$db})) {
	    print $logh "ERROR: BBH: $fam_id: $db already contributes (" . join(" ", @{$BtoP{$fam_id}->{$db}}) . "); adding $acc.\n";
	}
	push @{$BtoP{$fam_id}->{$db}}, $acc;
	$PtoB{$acc} = $fam_id;
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
	# first column is id so start from 1
	for (my $i=1;$i<@dg;$i++) {
	    if (! $dg[$i]) { next }
	    my $acc = id_from_accession($dg[$i]);
	    if (defined $PtoS{$acc}) { 
		print $logh "SYN: $famid: protein $acc has already been assigned to $PtoS{$acc}\n";
		next; 
	    }
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
	my $fam_id = sprintf("MCL%05d", $i);
	chomp $line;
	my @f = split /\t/, $line;
	foreach my $f(@f) {
		if (! $f) { next }
	    my $acc = id_from_accession($f);
	    if (defined $PtoC{$acc}) { 
		print $logh "MCL: $fam_id: protein $acc has already been assigned to $PtoC{$acc}\n";
		next; 
	    }
	    push @{$CtoP{$fam_id}->{$PtoDB->{$acc}}}, $acc if (defined $PtoDB->{$acc});
	    $PtoC{$acc} = $fam_id; # link cluster id to gene
	}
    }
    return \%PtoC, \%CtoP;
}

sub read_OrthoFinder {
    my $infile = shift;
    # apparently, OF outputs \r's
    $/="\r\n";
    # read in the clusters
    open my $in, $infile or die "Can't open $infile for read: $!\n";
    my %PtoOF;
    my %OFtoP;
    my $i = 0;
    #first line is header
    my $header = <$in>;
    # process the rest of the lines
    while (my $line = <$in>) {
	$i++;
	chomp $line;
	my @f = split(/\t/, $line);
	my $OG = shift(@f);
	foreach my $f(@f) {
	    chomp $f;
	    my @a = split(/,\s+/, $f);
	    foreach my $id(@a) {
		if (! $id) { next }
		my $acc = id_from_accession($id);
		if (! $acc || $acc eq "") { warn "Why no id from $id?"}
		if (!defined $PtoDB->{$acc}) { die "WTF? Why is $acc not in PtoDB???\n";}
		if (! defined $PtoOF{$acc}) {
		    push @{$OFtoP{$OG}->{$PtoDB->{$acc}}}, $acc if (defined $PtoDB->{$acc});
		    $PtoOF{$acc} = $OG;

		} elsif ($PtoOF{$acc} ne $OG) {
		    print $logh "OrthoFinder: $OG: protein '$acc' has already been assigned to $PtoOF{$acc}\n";
		    next;
		}
	    }
	}
    }
    $/="\n";
    return \%PtoOF, \%OFtoP;
}

sub read_mmseqs2 {
    my $infile = shift;
    my $PtoDB = shift;
    open(my $in, $infile) or die "Can't open mmseqs2 file $infile: $!\n";
    my %PtoMM;
    my %MMtoP;
    my $index = 0;
    my $rep_seq = "";
    my $thisFam;
    #print STDERR "MMSEQS parse: ";
    while (my $line = <$in>) {
	#print STDERR ".";
	chomp $line;
	my @f = split(/\t/, $line);
	my $refacc = &id_from_accession($f[0]);
	if (! defined $PtoDB->{$refacc}) {
	    warn "$refacc isn't in a proteome faa\n";
	    next;
	}
	my $thisacc = &id_from_accession($f[1]);
	if (! defined $PtoDB->{$thisacc}) {
	    warn "$thisacc isn't in a proteome faa\n";
	    next;
	}
	# because mmseqs uses the first acc in the family as the family name,
	# we need to keep track of what that name is in $rep_seq
	if ($refacc ne $rep_seq) {
	    # we found a new family, so re-instantiate rep_seq and make a family name
	    $rep_seq = $refacc;
	    $thisFam = sprintf "MMS%05d", ++$index;
	}
	if (defined $PtoMM{$thisacc}) { # see if this member is already a member of another family
	    print $logh "MMSEQS2: $thisFam ($rep_seq): protein $thisacc has already been assigned to $PtoMM{$thisacc}\n";
	    next;
	}
	$PtoMM{$thisacc} = $thisFam; # assign the member to a family
	push @{$MMtoP{$thisFam}->{$PtoDB->{$thisacc}}}, $thisacc;
    }
    print "\n";
    return \%PtoMM, \%MMtoP;
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
	if (defined $PtoH{$acc}) { 
	    print $logh "HMM: $hmm_acc: protein $acc has already been assigned to $PtoH{$acc}\n";
	    next;
	}
	$PtoH{$acc} = $hmm_acc;
	push @{$HtoP{$hmm_acc}->{$PtoDB->{$acc}}}, $acc if $PtoDB->{$acc};
    }
    return (\%PtoH, \%HtoP);	    
}

sub get_best {
    # We want this to return the best match to the input score
    # The problem is we want synteny (8) only to tie-break
    # 2023-05-03 Have changed this so that we only want members that MATCH the target score
    # ie have the SAME family memberships
    my $ref = shift;
    my $targetScore = shift;

    # get the best score
    #my $bestscore = (sort {$b <=> $a} values %$ref)[0];

    # MCL alone is not enough
    #if ($bestscore == 1) { return }

    my @best = ();
    # any candidate that matches the best score is returned
    foreach my $this (sort {$ref->{$b} <=> $ref->{$a}} keys %$ref) {
	#		if ($ref->{$this} == $bestscore) {
	if ($ref->{$this} == $targetScore) {
	    push @best, $this;
	}
    }
    return (@best);
}

sub read_OrthoMCL {
	my $infile = shift;
	open my $in, $infile or die "Can't open $infile for read: $!\n";
    my %PtoOM;
    my %OMtoP;

    while (my $line = <$in>) {
		chomp $line;
		my @f = split(/\s+/, $line);
		my $fam_id = shift @f;
		$fam_id =~ s/://;

		foreach my $f(@f) {
			if (! $f) { next }
			my ($db, $acc) = split(/\|/, $f);
	    	if (defined $PtoOM{$acc}) { 
				print $logh "OrthoMCL: $fam_id: protein $acc has already been assigned to $PtoOM{$acc}\n";
				next; 
	    	}
	    	push @{$OMtoP{$fam_id}->{$PtoDB->{$acc}}}, $acc if (defined $PtoDB->{$acc});
	    	$PtoOM{$acc} = $fam_id; # link cluster id to gene
		}
    }
    return \%PtoOM, \%OMtoP;
}