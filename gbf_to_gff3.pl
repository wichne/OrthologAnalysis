#!/usr/bin/perl

use Bio::SeqIO;
use Bio::Seq;
#use Getopt::Std;
use Data::Dumper;
use strict;

my $args = {};
#&getopts('i:F:', $args);

#my $filename = $args->{'i'} or die "Need to provide filename with -i\n";
#my $format = $args->{'F'} ? $args->{'F'} : undef;
while (my $filename = <STDIN>) {
    my $outfile = $filename;
    $outfile =~ s/\.\w+$/\.gff/;
    open my $out, ">$outfile";
    select $out;
    my $in = Bio::SeqIO->new(-file => $filename);
    while (my $seqo = $in->next_seq) {
	my @features = $seqo->get_SeqFeatures();
	# the first feature in a genbank file is the sequence itself
	my $seqObj = shift @features;
	my $seq_id = $seqObj->seq_id;
	if (! defined $seq_id) { warn "No seq_id found in object"; next;}
	
	my ($start, $strand, $locus_tag);
	foreach my $featObj (@features) {
	    if ($featObj->primary_tag eq "CDS") {
		if ($featObj->has_tag('locus_tag')) {
		    $locus_tag = [$featObj->get_tag_values('locus_tag')]->[0];
		    if (!$locus_tag) { die "Can't find locus tag for " . $featObj->{'display_id'} . "\n"; }
		    my ($start, $end) = ($featObj->start, $featObj->end);
		    $strand = $featObj->strand < 0 ? "-" : "+";
		    
		    print "$seq_id\tGBF\tCDS\t$start\t$end\t.\t$strand\t0\t";
		    my @f = ("ID=$locus_tag");
		    my @tags = $featObj->get_all_tags();
		    foreach my $tag (@tags) {
			my @values = $featObj->get_tag_values($tag);
			push @f, "$tag=" . join(',', @values);
		    }
		    print join(';', @f);
		    print "\n";
		    next;
		}
	    }
	}
    }
}

exit();

