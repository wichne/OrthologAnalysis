#!/usr/bin/perl

use strict;
use Getopt::Long;
use lib $ENV{ORTHO};
use OrthologAnalysis;

my $nameglob = "gff";
my $path = ".";
my $faext = "faa";
my %okay;
my $err = 0;

open(LOG, ">OrthologAnalysis.${$}.errorlog.txt");
select LOG;

&GetOptions(
    "glob=s" => \$nameglob,
    "path=s" => \$path,
    "faext=s" => \$faext,
    );

my (@anndb, @faadb);
for my $f (glob("$path/input/*.$nameglob")) { $f =~ /$path\/input\/(.*)\.$nameglob$/; push @anndb, $1 }
for my $f (glob("$path/input/*.$faext")) { $f =~ /$path\/input\/(.*)\.$faext$/; push @faadb, $1 }
if (@anndb != @faadb) { die "Need to have a $nameglob (" . scalar(@anndb) . ")"
			    . " and a $faext (" . scalar(@faadb) . ") file for each genome\n" }
# The FASTA file records
my ($PtoFAA, $FAAtoP) = &parse_protein_fasta_files("$path/input", $faext);
# The GFF file records
my ($PtoGFF, $GFFtoP, $PtoANN) = &parse_gff_files("$path/input", $nameglob, $PtoFAA);



foreach my $db (@faadb) {
    print LOG "checking $db\n";
    my ($faac, $gffc);
    
    # Do we have matching file names?
    if (! defined $GFFtoP->{$db} ) { die "Something wrong with input file names that is causing a disagreement in database name ('$db')\n";}

    # Do we have the same number of records?
    if (@{$GFFtoP->{$db}} != @{$FAAtoP->{$db}}) {
        $err += 1;
        warn "$db does not have same number of proteins in $nameglob and $faext files: "
             . join (",", (scalar(@{$GFFtoP->{$db}}), scalar(@{$FAAtoP->{$db}}))) . "\n";
    }

    # Do we have the same accessions in each file?
    # this is the faa to gff check
    foreach my $prot(@{$FAAtoP->{$db}}) {
        if (! exists $PtoGFF->{$prot} || $PtoGFF->{$prot} ne $db) {
	    print "$prot is in $faext file but not $nameglob file\n";
            $faac += 1;
        } else {
            $okay{$prot};
        }
    }
    if ($faac) {
        warn "\t$faac proteins in $db $faext file are not $nameglob file\n";
    }
    # this is the gff to faa check
    foreach my $prot(@{$GFFtoP->{$db}}) {
        if ($okay{$prot}) { next }
        if (! exists $PtoFAA->{$prot} || $PtoFAA->{$prot} ne $db) {
	    print "$prot is in $nameglob file but not $faext file\n";
            $gffc += 1;
        }
    }
    if ($gffc) {
        warn "$gffc proteins in $db $nameglob file are not in $faext file\n";
    }
}
if ($err == 0 ) { print "Everything look okay to go!!\n";}
exit($err)
