package OrthologAnalysis;

use strict;

use warnings;
use Exporter qw(import);
 
our @EXPORT = qw(id_from_accession parse_protein_fasta_files parse_gff_files);

sub id_from_accession {
    my $accstr = shift;
    my $acc;

    $accstr =~ s/^\s*//;
    $accstr =~ s/\s.*$//; # get rid of anything after a space
    my @acc = split(/\|/, $accstr); # split across the complex accession separator "|"

    # if the accession is complex, the first field is a descriptor, so use the second field
    if (@acc > 1) {
	if ($acc[1] ne "") { $acc = $acc[1] }
	else { $acc = $acc[0] }  # if the second field is blank, use the first field
    } else { $acc = $accstr }
    return $acc;
}

sub parse_protein_fasta_files {
    my $path = shift;
    my $filext = shift;

    print STDERR "Parsing $filext files in $path for db and accessions...\n";
    my $PtoDB = {}; # holds a mapping of the protein ids to the db name
    my $DBtoP = {};
    foreach my $fa_file (glob("$path/*.$filext")) {
	if ($fa_file =~ /\ball[\.\_]/) {
	    print STDERR "\t$fa_file looks like the agglomeration file, so skip it...\n";
	    next;
	}

	# parse the db name
	my $db;
	if ($fa_file =~ /$path\/(.*)\.$filext$/) {
	    $db = $1;
	}
	if (! $db) { die("Couldn't parse db name from '$fa_file'\n") }
	else { print STDERR "\t$fa_file -> $db\n"; }

	open(my $fa, $fa_file) or die "Can't open '$fa_file' for read: $!\n";
	while (my $h = <$fa>) {
	    if ($h =~ /^>(\S+)/) {
		my $acc = id_from_accession($1);
		if (!$acc) { die "What is up with $1 from '$h'?\n";}
		if (defined $PtoDB->{$acc}) {
		    die "Non-unique identifier '$acc' exists in both db '$db' and '"
			. $PtoDB->{$acc} . "'\n";
		} else {
		    $PtoDB->{$acc} = $db;
		    push @{$DBtoP->{$db}}, $acc;
		}
	    }
	}
    }
    print STDERR "Done.\n";
    return($PtoDB, $DBtoP);
}

sub parse_gff_files {
    my $path = shift;

    print STDERR "Parsing gff files in $path for db and accessions and annotations\n";
    my $PtoANN = {};
    my $PtoDB = {};
    my $DBtoP = {};

    foreach my $gff_file (glob("$path/*.gff")) {
	# parse the db name
	my $db;
	if ($gff_file =~ /$path\/(.*)\.gff$/) {
	    $db = $1;
	} 
	if (! $db) { die("Couldn't parse db name from '$gff_file'\n") }
	else { print STDERR "\t$gff_file -> $db\n" }

	open(my $G, $gff_file) or die "Can't open '$gff_file' for read: $!\n";
	my ($acc, $product, $ec, $gene, $locus_tag);
	while (my $line = <$G>) {
	    chomp $line;
	    my @l = split(/\t/, $line);
	    if (@l != 9) { die "gff file $gff_file has bad format?\n'$line'\n";}
	    next unless $l[2] eq "CDS";
	    my @tags = split(/\;/, $l[8]);
	    foreach my $t (@tags) {
		my ($k, $v) = split(/\s*=\s*/, $t);
		if ($k eq "ID") { $acc = id_from_accession($v) }
		if ($k eq "product") { $product = $v }
		if ($k eq "EC_num") { $ec = $v }
		if ($k eq "gene_sym") { $gene = $v }
		if ($k eq "locus_tag") { $locus_tag = $v }
	    }
	    if (! $acc) { die("Why no acc ('ID=') from '$line'?\n"); }
	    $PtoANN->{$acc} = {'product' => $product,
				   'EC_num' => $ec,
				   'gene_sym' => $gene,
				   'locus_tag' => $locus_tag };
	    if (defined $PtoDB->{$acc}) {
		die "Non-unique identifier '$acc' exists in both db '$db' and '"
		    . $PtoDB->{$acc} . "'\n";
	    } else {
		$PtoDB->{$acc} = $db;
		push @{$DBtoP->{$db}}, $acc;
	    }
	}
    }
    return ($PtoDB, $DBtoP, $PtoANN);
}

1;

