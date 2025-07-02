package OrthologAnalysis;

use strict;

use warnings;
use Exporter qw(import);
 
our @EXPORT = qw(id_from_accession parse_protein_fasta_files parse_gff_files);

sub id_from_accession {
    my $accstr = shift;
    my $acc;
	if (! $accstr) { die }
    # clean up spaces
    $accstr =~ s/^\s*//;
    $accstr =~ s/\s.*$//; # get rid of anything after a space

        # Had to add the colon because accs like fig|WP23456.1 were showing up as fig:WP23456.1 in the btab
    # Not sure if this is a btab thing or a blast thing.
    my @acc = split(/[\|:]/, $accstr); # split across the complex accession separator "|"

    # if the accession is complex, the first field is a descriptor, so use the second field
    if (@acc > 1) {
	if ($acc[1] ne "") { $acc = $acc[1] }
	else { $acc = $acc[0] }  # if the second field is blank, use the first field
    } else { $acc = $accstr }

    # get rid of 'cds-' prefix that RefSeq adds
    $acc =~ s/^cds-//;

    return $acc;
}

sub parse_protein_fasta_files {
    my $path = shift;
    my $filext = shift;

    print "Parsing $filext files in $path for db and accessions...\n";
    my $PtoDB = {}; # holds a mapping of the protein ids to the db name
    my $DBtoP = {};
    my $faa_count = 0;
    foreach my $fa_file (glob("$path/*.$filext")) {
	if ($fa_file =~ /\ball[\.\_]/) {
	    print LOG "\t$fa_file looks like the agglomeration file, so skip it...\n";
	    next;
	}

	# parse the db name
	my $db;
	if ($fa_file =~ /$path\/(.*)\.$filext$/) {
	    $db = $1;
	}
	if (! $db) { die("Couldn't parse db name from '$fa_file'\n") }
	#else { print STDERR "\t$fa_file -> $db\n"; }
	$faa_count++;
	open(my $fa, $fa_file) or die "Can't open '$fa_file' for read: $!\n";
	while (my $h = <$fa>) {
	    if ($h =~ /^>(\S+)/) {
			my $acc = id_from_accession($1);
			if (!$acc) { die "What is up with $1 from '$h'?\n";}
			if (defined $PtoDB->{$acc}) {
				warn "Non-unique identifier '$acc' exists in both db '$db' and '"
				. $PtoDB->{$acc} . "'\n";
			} else {
				$PtoDB->{$acc} = $db;
				push @{$DBtoP->{$db}}, $acc;
			}
	    }
	}
    }
    print "Processed $faa_count protein files.\n";
    return($PtoDB, $DBtoP);
}

sub parse_gff_files {
	# GFF file can be standard genome gff where the first column identifier is the contig id
	# and the location on the contig is specified in columns 4 and 5,
	# # # # OR # # #
	# GFF could be a UniProt-special protein file GFF, where the identifier in the first column
	# is the protein id.
    my $path = shift;
    my $nameglob = shift;
    my $TargetFastaIds = shift;
    my $logh = shift;
    if (! $logh) { $logh = \*STDERR; }
    
    print $logh "Parsing gff files in $path for db and accessions and annotations\n";
    my $PtoANN = {};
    my $PtoDB = {};
    my $DBtoP = {};
    my $gff_count = 0;
    foreach my $gff_file (glob("$path/*.$nameglob")) {
		# parse the db name
		my $db;
		if ($gff_file =~ /$path\/(.*)\.$nameglob$/) {
			$db = $1;
		} 
		if (! $db) { die("Couldn't parse db name from '$gff_file'\n") }
		#else { print STDERR "\t$gff_file -> $db\n" }

		$gff_count++;
		open(my $G, $gff_file) or die "Can't open '$gff_file' for read: $!\n";
		while (my $line = <$G>) {
			my ($acc, $id, $product, $ec, $gene, $locus_tag, $pseudo, $Name, $protein_id);
			next if ($line =~ /^#/);
			last if ($line =~ /^>/); # have hit the fasta section
			chomp $line;
			next if (!$line);
			my @l = split(/\t/, $line);

			# do we have the right number of fields?
			if (@l != 9) { 
				warn "gff file $gff_file has bad format (not enough fields)?\n'$line'\n";
				if ($l[2] eq "CRISPR") {
					warn "\tThis is a known problem with CRISPR elements from IMG. Skipping.";
					next;
				} else { die }
			}
			
			# We only care about protein records
			next unless $l[2] eq "CDS" or $l[2] eq "Chain";

			my @tags = split(/\;/, $l[8]);
			if ($l[2] eq "Chain") { $id = id_from_accession($l[0])}
			foreach my $t (@tags) {
				my ($k, $v) = split(/\s*=\s*/, $t);
				if ($l[2] eq "CDS" && $k eq "ID") { $id = id_from_accession($v) }
				elsif ($k eq "product") { $product = $v }
				elsif ($k eq "EC_num") { $ec = $v }
				elsif ($k eq "gene_sym") { $gene = $v }
				elsif ($k eq "locus_tag" || $k eq "PNNL") { $locus_tag = $v }
				elsif ($k eq "Name") { $Name = $v }
				elsif ($k eq "protein_id") { $protein_id = $v }
				elsif ($l[2] eq "Chain" && $k eq "Note") { $product = $v }
				elsif ($l[2] eq "pseudogene") { $pseudo = 1 }
			}
			$acc = $id;

			# if we have a target list of accessions that must be matched,
			# find the field that contains the correct accession
			if ($TargetFastaIds) {
				foreach my $a (($id, $locus_tag, $Name, $protein_id)) {
					if ($a && defined $TargetFastaIds->{$a}) { $acc = $a; last; }
				}
			} #elsif ($locus_tag) { $acc = $locus_tag }

			# Pseudogenes are likely not represented in the protein fasta file
			if ($pseudo) { warn "PSEUDO: $line\n"; next; }

			# This is a significant problem...
			if (! $acc) { warn ("FATAL ERROR: No acc (ID=) from '$line'? !!!!!!!!\n"); next; }
					
			$PtoANN->{$acc} = {'product' => $product,
					'EC_num' => $ec,
					'gene_sym' => $gene,
					'locus_tag' => $locus_tag };
			if (defined $PtoDB->{$acc}) {
				if ($PtoDB->{$acc} eq $db) {
					print $logh "WARNING: '$acc' appears multiple times in '$gff_file'. This is probably not a problem, but you have been warned.\n$line\n";
				} else {
					print $logh "FATAL: Non-unique identifier '$acc' exists in both db '$db' and '"
					. $PtoDB->{$acc} . "'\n";
				}
			} else {
				$PtoDB->{$acc} = $db;
				push @{$DBtoP->{$db}}, $acc;
			}
		}
    }
    print $logh "Processed $gff_count gff files.\n";
    return ($PtoDB, $DBtoP, $PtoANN);
}

1;

