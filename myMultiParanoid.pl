#!/usr/bin/perl
#
# Copyright (C) Andrey Alexeyenko and Erik Sonnhammer, 2006
# For a license to use this program, contact Erik.Sonnhammer@cbc.su.se
###############################################################################
# Andrey Alexeyenko, Stockholm Bioinformatics Center, Albanova, Stockholm 2006
# mailto:avalex99@mail.ru
###############################################################################
# The program is a postprocessor of ortholog tables produced by InParanoid. It combines pairwise InParanoid outputs (say, A-B, B-C, A-C) and makes multi-species clusters of orthologs ABC. 

use strict vars;
use DBI;
use Getopt::Long;
use Cwd;
require 5.002;

our($inputdir, $makeunique,  $indirect, $Pair, 
    @fields, @cutoffs, $output, @fields_out, $out, $algo, $debug, $SQLprefix, 
    %place, $nplaces, $source_by_cluster, $source_by_gene, $source_by_taken, $seeds, $use_bootstrap, $maxINPscore, $toMainOrthologOfSelf);
my($npairs, $i, $j, $ii, $ee, $gg, $nn, $oo, $s2, $s1, $nspec, $as, $newmain, @distances, @species, $sp1, $sp2, $clusno, $specletters, $a_pair, $a_cluster, $pms, $input, $output, $inputfile, $new, $thepair);
###############################################################################

#FIRST OF ALL, set the variables $inputdir and $output below. MultiParanoid may not work without it!!!

#Parameters to submit (may be abbreviated to the first letter):
# OBLIGATORY:
# -species: a list of several species connected with "+"; e.g. mouse+cat+dog; names of the genomes exactly how they are called both in file names and inside of files. NOTE: due to phylogenetic tree conflicts, order of the genomes slightly changes the result. The input files (the InParanoid output) should thus be called:
# sqltable.dog-cat, sqltable.dog-mouse etc. 
#To use other file names, one should change the variable $SQLprefix below.

# OPTIONAL:
# -cutoff: deprecated. But one can use it to tighten the clusters (to filter out weaker in paralogs by restricting with the distance to the seed orthologs); Example: " -cutoff 0.5 " or, for doing 3 different cluster versions at the same time: 0+0.2+0.5. Default: 0 (no cutoff). 
# -algo: how to define between-paralogs distances while clustering. 
#              Options: onlymains (default and only reasonable)
#                       toallsorted
# -debug: 0 or 1; default: 0
# -input: deprecated, now defined by -species parameter, and the input files shall have conventional names
# -output: deprecated, special message will tell the name of the output file 
# -unique: 0 or 1; set to 1, if the duplicates (genes clustered twice) should be removed; default: 0

#The MultiParanoid output header:
### clusterID       species gene    is_seed_ortholog        confidence_score        species_in_cluster      tree_conflict###

#is (hopefully) self-explanatory. Though:

# "is_seed_ortholog" means the protein was a seed ortholog in at least one InParanoid cluster.
# "confidence_score" is an average InParanoid score across the input clusters 
# "tree_conflict" indicates that, from the point of view of different species, the number of inparalogs varied in at least one other species ("diff.numbers") or the numbers were the same, but the IDs differed ("diff.names"). 

#NOTE: s.c. 'distances' between cluster members are, in fact, similarities. Variables are called 'distances' for historical reasons. 
# 'main orthologs' (defined by InParanoid) are better to be called 'seed orthologs'.  

# Settings:    

$algo 		= 'onlymains';
$debug 		= 0;
$use_bootstrap = 0;
$maxINPscore = 1;
$toMainOrthologOfSelf = 0;
$indirect = 1; #this option allows using the currently default InParanoid output, where seed orthologs are recognized by the feature "IPscore = 1.00". Otherwise ($indirect = 0) a special format should be used (needs a change in the InParanoid script and is now deprecated).

open LOG, ">MultiParanoid.$$.log";

&GetOptions("inpath=s" => \$inputdir,
#	    "species=s" => \@species,
	    "debug:i" => \$debug,
	    "unique!" => \$makeunique);
$output = "$inputdir/multiparanoid.out";
rename $output, "$output.$$.old" if (-e $output);

# Format of inparanoid input
if ($indirect) {
    $place{'clusterID'} 	= 0; #these should correspond to the current InParanoid output format, but they don't
    $place{'BLASTscore'} 	= 1;
    $place{'species'} 		= 2;
    $place{'IPscore'} 		= 3;
    $place{'gene'} 		= 4;
    $place{'bootstrap'}         = 5; # bootstrap value
    $place{'main'} 		= 6; # Not in input file - instantiated in loadData
    $place{'pair'} 		= 7; # Not in input file - instantiated in loadData
}
else {
    ###The older, ava_sql, version:
    $place{'clusterID'} 	= 0;
    $place{'BLASTscore'} 	= 1;
    $place{'pair'} 			= 2;
    $place{'species'} 		= 3;
    $place{'gene'} 			= 4;
    $place{'IPscore'} 		= 5;
    $place{'main'} 		= 6;
}

my $dir = getcwd();
chdir($inputdir);
my %foo;
my @files = glob "*.ortho";
foreach my $s (@files) {
    if ($s =~ /(.*)\.(.*)\.ortho/) {
	$foo{$1}++;
	$foo{$2}++;
    } else { warn "Malformed ortho filename $s - extra '.'s" }
}
@species = sort keys(%foo);
chdir($dir);

print LOG "Species list: ".join(', ', @species)."\n";
if (@species < 3) { die "No need to run program on less than 3 species\n";}

# '-'-joined sorted list of species
$specletters = fstLetters(@species);

# Can add this to command options later - see parseParameters below
@cutoffs = (0) if !exists($cutoffs[0]);

# the number of fields in the input (why?)
$nplaces = scalar(keys(%place));

@fields_out = (clusterID, species, gene, is_seed_ortholog, confidence_score, species_in_cluster, tree_conflict);
push @fields_out, 'cutoff' if scalar(@cutoffs) > 1;

$nspec = scalar(@species);

open OUT, ">$output";
# Find the input file for each species pair and load data
# Thus input files must have format species1-species2 where species is really just the filename
print "Genome pairs found and used: \n"  if $debug > 0;
# since we don't want to assume which direction the pairs were done,
# (ie which 'species' name comes first in the input file)
# just go through all combinations
for ($s1 = 0; $s1 < $nspec; $s1++) {
    for ($s2 = $s1+1; $s2 < $nspec; $s2++) {
#	next if ($species[$s1] eq $species[$s2]);
	$Pair = $species[$s1].'.'.$species[$s2];
#	$Pair = $species[$s1].'___'.$species[$s2];
	$inputfile = "$inputdir/$Pair.ortho";
	if (-e $inputfile) {
	    loadData($inputfile);
	    print LOG "$inputfile loaded.\n";
	    $npairs++;
	} else {
	    $Pair = $species[$s2].'.'.$species[$s1];
	    $inputfile = "$inputdir/$Pair.ortho";
	    if (-e $inputfile) {
		loadData($inputfile);
		print LOG "$inputfile loaded.\n";
		$npairs++;
	    } else {
		warn "Can't find input file for $species[$s1] / $species[$s2]\n";
	    }
	}
    }
}
die "Not enough (or too many) species pair files ($npairs) for ".join(", ", @species)."\n" if $npairs !=  ($nspec * ($nspec - 1)) / 2;

$clusno = 1;
for $sp1(@species) {
    for $sp2(@species) {
	$thepair = "$sp1.$sp2";

	# $source_by_cluster is built by loadData.
	next if (!$source_by_cluster->{$thepair} or ($sp1 eq $sp2));
	print LOG "Analyzing $sp1 vs. $sp2: ";
	print LOG scalar(keys(%{$source_by_cluster->{$thepair}}))." clusters for this pair...\n";
	for $a_cluster(keys(%{$source_by_cluster->{$thepair}})) {
	    undef $new; undef $seeds;
	    $seeds = findOrthos($a_cluster, $thepair);
	    
#find recordset of orthologs of the both species (for cluster a) from the pairwise subset i-j
	    next if !defined($seeds->[0]);
	  LINE1:for $gg(@{$seeds}) {
	      next if (!$gg->[$place{'main'}] or $gg->[$nplaces]);
	      $gg->[$nplaces] = 1; # start a new field at the end of the data
	      for $a_pair(keys(%{$source_by_cluster})) {
		  # Find A-C and B-C bbh clusters
		  if (($a_pair =~ $gg->[$place{'species'}]) and ($gg->[$place{'pair'}] ne $a_pair)) {
		      # look for bbh that includes the current seed gene from the A/B-C pair
		      # and add to $new arrayref
		      $new = findNewOrtho($gg->[$place{'gene'}], $a_pair);
		      undef $newmain;
		      if (defined($new)) {
			  # add the new genes to the $seeds datastructure
			  for $nn(@{$new}) {
			      push(@{$seeds}, $nn) if isNew($seeds, $nn);
			      $newmain = 1 if $nn->[$place{'main'}];
			  }
			  # Keep looking until there's no more new genes to add
			  # and this is where the single-linkage occurs...
			  redo LINE1 if $newmain;
		      }
		  }}}
	    
	    makeClusters($seeds, \@cutoffs, \@species, $clusno++, treeConflict($seeds), $specletters);
	}}}
if ($clusno > 20) {
    print OUT join("\t", @fields_out)."\n";
    for $oo(@{$out->{'byRecord'}}) {
	print OUT join("\t", @{$oo})."\n" if defined($oo);
    }
    print "Making ortholog clusters of ".join(', ', @species)." done. The result is in file $output\n";
}
else {
    die "No output has been produced..."
}

close LOG;

###########################################################################

sub parseParameters {
    my($parameters) = @_;
    my($key, $value, $pars, %pars);

#print "$parameters\n";
    $_ = $parameters;
    while (m/\-(\w+)\s+([a-zA-Z0-9._=+]+)/g) {
	next if ((lc($2) eq 'f') or (lc($2) eq 'false') or (lc($2) eq 'no') or (lc($2) eq '0'));
	$pars{substr($1, 0, 1)} = $2;
    }
    if (scalar(keys(%pars)) < 1 or !$pars{'s'} or m/\-h/) {
	print "\nNOT ENOUGH PARAMETERS!\nOBLIGATORY:
-species: a list of several species connected with '+'; e.g. mouse+cat+dog; names of the genomes exactly how they are called both in file names and inside of files. NOTE: due to phylogenetic tree conflicts, order of the genomes slightly changes the result. The input files (the InParanoid output) should thus be called:
sqltable.dog-cat, sqltable.dog-mouse etc. 
To use other file names, one should change the variable \$SQLprefix inside the scrtipt.
\nOPTIONAL:
\n-cutoff: deprecated. But one can use it to tighten the clusters (to filter out weaker in paralogs by restricting with the distance to the seed orthologs); Example: ' -cutoff 0.5 ' or, for doing 3 different cluster versions at the same time: 0+0.2+0.5. Default: 0 (no cutoff). 
\n-algo: how to define between-paralogs distances while clustering. 
              Options: onlymains (default and only reasonable)
                       toallsorted
\n-debug: 0 or 1; default: 0
\n-input: deprecated, now defined by -species parameter and the variables \$inputdir and \$SQLprefix, and the input files shall have conventional names
\n-output: deprecated, special message will tell the name of the output file 
\n-unique: 0 or 1; set to '0' if the duplicates (genes clustered twice) should BE REMOVED, to '1' otherwise; default: 0\n\n"; exit;
    }
    while (($key, $value) = each(%pars)) {
	print "Parameter $key = $value;\n";
    }

    return(\%pars);
}

sub treeConflict { #tells if there is a discrepancy, for species B, between 2-species trees A-B and B-C 
    my($seeds) = @_;
    my($i, $ii, $jj, $kk, $sets);

    for $ii(@{$seeds}) {
	push @{$sets->{$ii->[$place{'species'}]}->{$ii->[$place{'pair'}]}}, $ii->[$place{'gene'}];
    }

    for $ii(keys(%{$sets})) {
	for $jj(keys(%{$sets->{$ii}})) {
	    for $kk(keys(%{$sets->{$ii}})) {
		next if $jj eq $kk;
		if (scalar(@{$sets->{$ii}->{$jj}}) != scalar(@{$sets->{$ii}->{$kk}})) {return('diff. numbers');}
		sort {$a <=> $b} @{$jj};
		sort {$a <=> $b} @{$kk};
		for ($i = 0; $i < scalar(@{$sets->{$ii}->{$jj}}); $i++) {
		    if ($sets->{$ii}->{$jj}->[$i] ne $sets->{$ii}->{$kk}->[$i]) {
			return('diff. names');
		    }}}}}
    return(undef);
}

sub makeClusters {
    my($seeds, $cutoffs, $species, $clusno, $conflict, $param) = @_;
    my($clo, $le, $gg, @included, $distances, $specset, $cco);

    array2StdOut($seeds, "Primers") if $debug;
    if ($algo eq 'onlymains') {
	$distances = pairWiseDist2Main($seeds);
    }
    elsif ($algo eq 'toallsorted') {
	$distances = pairWiseDist2All($seeds);
    }
    else {die "Algorithm is not set...";}
    array2StdOut(setUnique($seeds), "Unique") if $debug;
    for $cco(@{$cutoffs}) {
	undef @included;

	@included = firstMains($seeds);
	@included = @{setUnique(\@included)} if (@included);
	array2StdOut(\@included, 'Mains') if $debug;
	goto LINE3 if $cco == 1; 

	if ($algo eq 'onlymains') 		{
	    push @included, @{getOthers(
				  setUnique($seeds),
				  $distances,
				  \@included,
				  $cco)};
	}
	elsif ($algo eq 'toallsorted') 	{
	    @included 	=  	@{getClosests(
				      setUnique($seeds),
				      $distances,
				      \@included,
				      $cco)};}
      LINE3:array2StdOut(\@included, "Cluster with cut-off $cco") if $debug;
      LINE4: $specset = join('', fstLetters(sort(listSpecies(\@included, $species))));
	$conflict = (($conflict) ? $conflict : 'No');
	for $gg(@included) {
	    next if $makeunique and alreadyClustered($gg, $cco); 
	    $ii = scalar(@{$out->{'byRecord'}});
	    $out->{'byName'}->{$gg->[$place{'gene'}]}->{$cco} = $ii;
	    
	    for $ee($clusno, 
		    $gg->[$place{'species'}], 
		    $gg->[$place{'gene'}], 
		    $gg->[$place{'main'}], 
		    $gg->[$place{'IPscore'}], 
		    $specset, 
		    $conflict) {push @{$out->{'byRecord'}->[$ii]}, $ee;}
	    push @{$out->{'byRecord'}->[$ii]}, $cco if scalar(@cutoffs) > 1;
	    
	}
#$out .=  join("\t", (100 * $cco, $clusno, $gg->[$place{'species'}], $gg->[$place{'gene'}], $gg->[$place{'main'}], $gg->[$place{'IPscore'}], $specset, $conflict, $param, $algo))."\n" if !$makeunique or !alreadyClustered($gg);
    }
    return;
}


sub fstLetters {
    # This appears to first return a sorted list of species
    # The second part returns the list joined by an _. Why?
    # The sub and variable names suggest it's just taking the first letter, 
    # but its not. Must be a legacy thing.
    my(@list) = @_;
    my($letters, $ii, $firstLetter, $lastLetter);
    return(join('-', sort {$a cmp $b} @list));

    # this never executes, so I'm commenting it.
    # $firstLetter = 0; 
    # for $ii(@list) {
    # 	$lastLetter = length($ii) - $firstLetter;
    # 	$letters .= "_".substr($ii, $firstLetter, $lastLetter);
    # }
    # return $letters;
}

sub listSpecies {
    my($ar, $species) = @_;
    my($already, %already, $ii, @flags);
# checks for species present in the current cluster

    for $ii(@{$ar}) {
	push(@flags, $ii->[$place{'species'}]) if !$already{$ii->[$place{'species'}]};
	$already{$ii->[$place{'species'}]} = 1;
    }
    return(@flags);
}

sub loadData {
    my($inputfile) = @_;
    my($line, @line, $i);

    open(IN, $inputfile);
    #print "Loading $inputfile...\n";
    while (<IN>) {
	@line = split /\s+/;
#$place{'BLASTscore'}, $place{'species'}, $place{'IPscore'}, $place{'main'}
	for $i(0..scalar(keys(%place))) {
	    # making $source_by_cluster hashref with key being Pair (= species1-species2), since we're always doing indirect=1
	    if ($indirect) {
		# this will always happen. 
		# Here we're setting the 'main' field in the hashrefs (both by_gene and by_cluster) to 1 or 0 
		# depending on whether the IPscore = maxINPscore (1), which as far as I can tell, it will
		# since all the IPscores in my file are 1.
		# Does this mean we don't need 'main' in the input file? Apparently.
		if ($i == $place{'main'}) {
		    # $Pair is the species pair (s1.s2)
		    # since we already know $indirect is true (from the above 'if'), we simplify
#		    $source_by_cluster->{$indirect ? $Pair : $line[$place{'pair'}]}->{$line[$place{'clusterID'}]}->{$line[$place{'gene'}]}->[$i] = ($line[$place{'IPscore'}] == $maxINPscore) ? 1 : 0;
#		    $source_by_gene->{$indirect ? $Pair : $line[$place{'pair'}]}->{$line[$place{'gene'}]}->{$line[$place{'clusterID'}]}->[$i] = ($line[$place{'IPscore'}] == $maxINPscore) ? 1 : 0;
		    $source_by_cluster->{$Pair}->{$line[$place{'clusterID'}]}->{$line[$place{'gene'}]}->[$i] = ($line[$place{'IPscore'}] == $maxINPscore) ? 1 : 0;
		    $source_by_gene->{$Pair}->{$line[$place{'gene'}]}->{$line[$place{'clusterID'}]}->[$i] = ($line[$place{'IPscore'}] == $maxINPscore) ? 1 : 0;
		    next;
		}
		# as above, except this populates 'pair' with the Pair value ('species1-species2')
		if ($i == $place{'pair'}) {
		    $source_by_cluster->{$Pair}->{$line[$place{'clusterID'}]}->{$line[$place{'gene'}]}->[$i] = $Pair;
		    $source_by_gene->{$Pair}->{$line[$place{'gene'}]}->{$line[$place{'clusterID'}]}->[$i] = $Pair;
		    next;
		}
	    }

	    # load input line data into equivalent array space in by_cluster and by_gene hashes.
	    $source_by_cluster->{$indirect ? $Pair : $line[$place{'pair'}]}->{$line[$place{'clusterID'}]}->{$line[$place{'gene'}]}->[$i] = $line[$i];
	    $source_by_gene->{$indirect ? $Pair :$line[$place{'pair'}]}->{$line[$place{'gene'}]}->{$line[$place{'clusterID'}]}->[$i] = $line[$i];
	}
    }
    return;
}

sub findOrthos { #initially picks up orthologs in the first pair of genomes
    my($a, $thepair) = @_;
    my($set, $sm, $i, $newel);

    return(undef) if $source_by_taken->{$thepair}->{$a};
    for (keys(%{$source_by_cluster->{$thepair}->{$a}})) {
	$newel = $#{$seeds} + 1;
#	for $i(0..($nplaces - 1)) { # this seems to drop the 'pairs' field from $seeds
	for $i(0..$nplaces) {
	    $seeds->[$newel]->[$i] = $source_by_cluster->{$thepair}->{$a}->{$_}->[$i];
	}
	$source_by_taken->{$thepair}->{$a} = 1;
    }
    return $seeds;
}

sub findNewOrtho { #adds orthologs from the remaining genome pairs
    my($g, $thepair) = @_;
    my($new, $ii, $j, $jj, $sm, $newel);

    for (keys(%{$source_by_gene->{$thepair}->{$g}})) {
	next if $source_by_taken->{$thepair}->{$_};
	$newel = $#{$new} + 1;
	for $jj(@{$source_by_gene->{$thepair}->{$g}->{$_}}) {
	    $new->[$newel]->[$j] = $jj;
	    $j++;
	}
    }
    if (scalar(@{$new}) == 1) {
	$new = findOrthos($new->[0]->[$place{'clusterID'}], $new->[0]->[$place{'pair'}]);
	if (scalar(@{$new}) > 1) {
	    return($new);}
	else {print  "Unpaired  seed ortholog...$new->[0]->[$place{'gene'}]\n";}
    }
    elsif (scalar(@{$new}) > 1) {
	print  "Multiple  seed orthologs...\nPair '$thepair'\tGene '$g'";
    }
    return(undef);
}

sub isNew {
    my($seeds, $new) = @_;
    my($ii);

    for $ii(@{$seeds}) {
	if (($ii->[$place{'gene'}] eq $new->[$place{'gene'}]) and ($ii->[$place{'pair'}] eq $new->[$place{'pair'}])) {
	    return(undef);
	}
    }
    return(1);
}

sub pairWiseDist2All { #between-gene distances for algortihm 'ToAllSorted'; deprecated
    my($seeds) = @_;
    my($i, $j, $ii, $jj, $new, $old, $same, $le, $distances, $ke, $va);

    $le = scalar(@{$seeds});

    for ($i = 0; $i < $le; $i++) {
	$ii = $seeds->[$i]->[$place{'gene'}];
#next if $distances->{$ii}->{$jj};
	for ($j = 0; $j < $i; $j++) {
	    $jj = $seeds->[$j]->[$place{'gene'}];
	    undef $same; undef $old; undef $new;
	    $same = 1 if ($seeds->[$i]->[$place{'species'}] eq $seeds->[$j]->[$place{'species'}]);
	    if ($same) {
		$old = $distances->{$ii}->{$jj}->{$seeds->[$i]->[$place{'pair'}]};}
	    else {$old = $distances->{$ii}->{$jj};}

	    next if ($seeds->[$i]->[$place{'pair'}] ne $seeds->[$j]->[$place{'pair'}]);
#THEY DESCRIBE THE SAME PAIR OF GENOMES
	    if ($same) { 
#THEY ARE FROM THE SAME SPECIES 
#next if $distances->    {$ii}->{$jj}->{$seeds->[$i]->[$place{'species'}]};
		if (!$seeds->[$i]->[$place{'main'}] and !$seeds->[$j]->[$place{'main'}]) { #BOTH ARE IN-PARALOGS
		    $distances->{$ii}->{$jj}->{$seeds->[$i]->[$place{'pair'}]} = $distances->{$jj}->{$ii}->{$seeds->[$i]->[$place{'pair'}]} = $seeds->[$i]->[$place{'IPscore'}] * $seeds->[$j]->[$place{'IPscore'}];
		}
		elsif ($seeds->[$i]->[$place{'main'}]) { #THIS ONE IS IN-PARALOG
		    $distances->{$ii}->{$jj}->{$seeds->[$i]->[$place{'pair'}]} = $distances->{$jj}->{$ii}->{$seeds->[$i]->[$place{'pair'}]} = $seeds->[$j]->[$place{'IPscore'}];
		}
		elsif ($seeds->[$j]->[$place{'main'}]) { #THIS ONE IS IN-PARALOG
		    $distances->{$ii}->{$jj}->{$seeds->[$i]->[$place{'pair'}]} = $distances->{$jj}->{$ii}->{$seeds->[$i]->[$place{'pair'}]} = $seeds->[$i]->[$place{'IPscore'}];
		}
		else { #BOTH ARE SEED ORTHOLOGS
		    $distances->{$ii}->{$jj}->{$seeds->[$i]->[$place{'pair'}]} = $distances->{$jj}->{$ii}->{$seeds->[$i]->[$place{'pair'}]} = $maxINPscore;
		}
		$new = $distances->{$ii}->{$jj}->{$seeds->[$i]->[$place{'pair'}]};
	    }	
	    else { 
#THEY ARE FROM THE DIFFERENT SPECIES
		if ($seeds->[$i]->[$place{'main'}] and $seeds->[$j]->[$place{'main'}]) #BOTH ARE SEED ORTHOLOGS
{
    $distances->{$ii}->{$jj} 	= inpDist($seeds, $i, $j);
    $distances->{$jj}->{$ii} 	= inpDist($seeds, $j, $i);
}
elsif ($seeds->[$i]->[$place{'main'}] or $seeds->[$j]->[$place{'main'}]) #STILL, ONE OF THEM IS IN-PARALOG (NO MATTER WHICH)
{
    $distances->{$ii}->{$jj}	= 
	$distances->{$jj}->{$ii} 	= 
	inpDist ($seeds, $i, $j) * inpDist ($seeds, $j, $i)	;
}
else 
#BOTH ARE IN-PARALOGS 
{
    $distances->{$ii}->{$jj}	= 
	mainDist($seeds, $i, $j) 	* 
	inpDist ($seeds, $i, $j) 	* 
	inpDist ($seeds, $j, $i)	;
    $distances->{$jj}->{$ii} 	= 
	mainDist($seeds, $j, $i) 	* 
	inpDist ($seeds, $i, $j) 	* 
	inpDist ($seeds, $j, $i)	;
}
	    }
$new = $distances->{$ii}->{$jj};
if (!$same) {print "Distance $ii-$jj is $new\n" if $debug; print "Distance $jj-$ii is $new\n" if $debug;}
else {
    if ($debug) {
	print "Distance $ii-$jj is\n";
	while (($ke, $va) = each(%{$distances->{$ii}->{$jj}})) {
	    print "\tfor $ke: $va\n";
	}
}
}
print "Distance $ii-$jj is not equal to old one $old\n" if ($old and ($old != $new));
	}
    }
#for ($k = 0; $k < $le; $k++) {$distances->[$k]->[$k] = 1;}
return($distances);
}

sub pairWiseDist2Main ($) {#between-gene distances for algortihm 'OnlyMains'
    my($seeds) = @_;
    my($le, $i, $j, $ii, $jj, $new, $old, $distances);

    $le = scalar(@{$seeds});
#$le2 = scalar(@{$species});
    for ($i = 0; $i < $le; $i++) {
	$ii = $seeds->[$i]->[$place{'gene'}];
	for ($j = 0; $j < $i; $j++) {
	    $jj = $seeds->[$j]->[$place{'gene'}];
	    undef $old; undef $new;
	    $old = $distances->{$ii}->{$jj};
	    next if ($seeds->[$i]->[$place{'pair'}] ne $seeds->[$j]->[$place{'pair'}]);
	    next if ($seeds->[$i]->[$place{'main'}] and $seeds->[$j]->[$place{'main'}]); #SEED ORTHOLOGS WILL BE CLUSTERED ANYWAY 
	    if	($seeds->[$i]->[$place{'species'}] eq $seeds->[$j]->[$place{'species'}]) {
#if (!$use_bootstrap) {$distances->{$ii}->{$jj} = $distances->{$jj}->{$ii} = $maxINPscore ;next;}
		next if !$toMainOrthologOfSelf;
		$distances->{$ii}->{$jj} = $distances->{$jj}->{$ii} = 
		    $seeds->[$j]->[$place{'IPscore'}] if $seeds->[$i]->[$place{'main'}];
		$distances->{$ii}->{$jj} = $distances->{$jj}->{$ii} = 
		    $seeds->[$i]->[$place{'IPscore'}] if $seeds->[$j]->[$place{'main'}];
	    }			#THEY ARE FROM THE SAME SPECIES AND THE DISTANCE IS NOT USED
	    else { 	        #THEY ARE FROM DIFFERENT SPECIES
		if ($seeds->[$i]->[$place{'main'}] or $seeds->[$j]->[$place{'main'}]) {
#ONE OF THEM IS SEED ORTHOLOG (NO MATTER WHICH)
		    $distances->{$ii}->{$jj}	= 
			$distances->{$jj}->{$ii} 	= 
			inpDist ($seeds, $i, $j) 	* 
			inpDist ($seeds, $j, $i)	;
		}
		else {
#BOTH ARE IN-PARALOGS
		    $distances->{$ii}->{$jj}	= 
			mainDist($seeds, $i, $j) 	* 
			inpDist ($seeds, $i, $j) 	* 
			inpDist ($seeds, $j, $i)	;
		    $distances->{$jj}->{$ii} 	=
			mainDist($seeds, $j, $i) 	* 
			inpDist ($seeds, $i, $j) 	* 
			inpDist ($seeds, $j, $i)	;
		}
	    }
	    $new = $distances->{$ii}->{$jj};
	    print "Distance $ii-$jj is $new\n" if $debug;
	    print "Distance $jj-$ii is $new\n" if $debug;
	    print "Distance $ii-$jj is not equal to old one $old\n" if ($old and ($old != $new));
	}
    }
    return($distances);
}

sub mainDist ($$$) { #finds distances between two main (seed) orthologs
    my($seeds, $from, $to) = @_;
    my($ii, $pair1, $pair2, $sum, $nu, $nc);

    return($maxINPscore) if !$use_bootstrap; #if we think bootstrap values are not relevant in calculating the distance, then every Main-Main distance should bemaxINPscore = 1 (and this is the default)
    $pair1 = $seeds->[$from]->[$place{'species'}].'-'.$seeds->[$to]->[$place{'species'}];
    $pair2 = $seeds->[$to]->[$place{'species'}].'-'.$seeds->[$from]->[$place{'species'}];
    $sum = $nc = $nu = 0;

    for $ii(@{$seeds}) {
	if (
	    ($ii->[$place{'main'}]) and 
	    (($ii->[$place{'pair'}] eq $pair1) or ($ii->[$place{'pair'}] eq $pair2)) and 
	    ($ii->[$place{'species'}]) eq $seeds->[$from]->[$place{'species'}]) {
	    $nc++;
	    $sum += $ii->[$place{'IPscore'}];
	}
    }

    if (!$nc) {
	print "Distance $seeds->[$from]->[$place{'gene'}]-$seeds->[$to]->[$place{'gene'}] (pair $pair1 or $pair2) not found\n";
	return(undef);
    }
    else 		
{
    return($sum / $nc);
}
}

sub inpDist ($$$) { #finds a distance from main (seed) ortholog to in-paralog
    my($seeds, $from, $to) = @_;
    my($le, $ii, $pair1, $pair2, $sum, $nu, $nc);

    return($maxINPscore) if ($seeds->[$from]->[$place{'main'}] and !$use_bootstrap);

    $le = scalar(@{$seeds});
    $pair1 = $seeds->[$from]->[$place{'species'}].'.'.$seeds->[$to]->[$place{'species'}];
    $pair2 = $seeds->[$to]->[$place{'species'}].'.'.$seeds->[$from]->[$place{'species'}];
    $sum = $nc = $nu = 0;
    for $ii(@{$seeds}) {
	if (
	    ($ii->[$place{'gene'}] eq $seeds->[$from]->[$place{'gene'}]) and 
	    (($ii->[$place{'pair'}] eq $pair1) or 
	     ($ii->[$place{'pair'}] eq $pair2))) {
	    return($ii->[$place{'IPscore'}]);
	}
    }
    print "Distance $seeds->[$from]->[$place{'gene'}]-$seeds->[$to]->[$place{'gene'}] not found\n";
    return undef;
}

sub sortBy ($$) {
    my($num, $ar) = @_;
    my($i, $max, $maxnum, @sorted);

#array2StdOut($ar, 'Before');

LINE001: while (scalar(@{$ar}) > 0) {
    if (scalar(@{$ar}) == 1 and !$ar->[$i]->[$num]) {
	last;
    }
    $max = 0;
    undef $maxnum;
    for ($i = 0; $i < scalar(@{$ar}); $i++) {
	last if !$ar->[$i]->[$num];
	if ($ar->[$i]->[$num] > $max) {
	    $max = $ar->[$i]->[$num];
	    $maxnum = $i;
	}
    }
    if (defined($maxnum)) {
	push @sorted, splice(@{$ar}, $maxnum, 1);
	next LINE001;
    }
}
    return(\@sorted);
}

sub setUnique ($) { #removes duplicated entries from the (pre-)cluster
    my($seeds) = @_;
    my($i, $j, $unique, $sds);

    $sds->[$place{'clusterID'}] = $seeds->[$place{'clusterID'}];
    for ($i = 1; $i < scalar(@{$seeds}); $i++) {
	$unique = 1;
	for ($j = 0; $j < scalar(@{$sds}); $j++) {
	    if (($seeds->[$i]->[$place{'species'}] eq $sds->[$j]->[$place{'species'}]) and ($seeds->[$i]->[$place{'gene'}] eq $sds->[$j]->[$place{'gene'}]))
{$unique = 0; 
 last;}
	}
if ($unique) {push(@{$sds}, $seeds->[$i]);}
    }
return($sds);
}

sub isAlready ($$) {
    my($as, $included) = @_;
    my($ai);

    for $ai(@{$included}) {
	return(1) if ($ai->[4] eq $as);
    }
    return(undef);
}

sub alreadyClustered { #helps to avoid clustering same gene twice. Of two alternatives, the one is chosen where the gene is a main (seed) ortholog or was included first. 
    my($g, $cco) = @_;
    my($ai);

    return(undef) if !defined($out->{'byName'}->{$g->[$place{'gene'}]}->{$cco});

    if ($out->{'byRecord'}->[$out->{'byName'}->{$g->[$place{'gene'}]}->{$cco}]->[3]) {
	return 1;
    }
    else {
	undef $out->{'byRecord'}->[$out->{'byName'}->{$g->[$place{'gene'}]}->{$cco}];
	return(undef);
    }}

sub firstMains { #seed orthologs are included unconditionally
    my($seeds) = @_;
    my($i, $le, @included);

    $le = scalar(@{$seeds});
    for ($i = 0; $i < $le; $i++) {
# and ($seeds->[$i]->[$place{'IPscore'}] == 1)
	if ($seeds->[$i]->[$place{'main'}]) {
	    push @included, $seeds->[$i];
	}
    }
    return(@included);
}

sub getOthers ($$$$) { #retains only genes with distance (well, it is in fact a SIMILARITY) above the cut-off $cco. If it is not set, or equals to 0, the procedure just returns IPscore as a reference - to be saved in the output file.
    my($seeds, $distances, $included, $cco) = @_;
    my($i, $le, $les2, $nc, $ii, $jj, $dist, $ad, $nai, @others);

    $le = scalar(@{$seeds});
    for $ii(@{$seeds}) {
	next if isAlready($ii->[$place{'gene'}], $included);
	undef $dist; undef $nc;
	for $jj(@{$included}) {
	    next if (($jj->[$place{'species'}] eq $ii->[$place{'species'}]) and !$toMainOrthologOfSelf);
	    if (!$distances->{$ii->[$place{'gene'}]}->{$jj->[$place{'gene'}]}) {
		print "The distance between $ii->[$place{'gene'}] and $jj->[$place{'gene'}] does not exist...\n" if $debug; 
		next;
	    } 
	    if ($jj->[$place{'main'}]) {
		print "( $ii->[$place{'gene'}] - $jj->[$place{'gene'}] ) + " if $debug;
		$dist += $distances->{$ii->[$place{'gene'}]}->{$jj->[$place{'gene'}]};
		$nc++;
	    }
	}
	return(undef) if !$nc;
	print " / $nc = ".$dist/$nc."\n" if $debug;
	if ($dist / $nc >= $cco) {
	    push @others, $ii;
	    $others[$#others]->[$place{'IPscore'}] = $dist / $nc;
	}
    }
    return(\@others);
}

sub getClosests ($$$$) { #used only in the deprecated 'ToAllSorted' algorithm
    my($seeds, $distances, $included, $cco) = @_;
    my($ni, $ii, $jj, $kk, $dist, $gotnew, $toadd, $nt, $max);

    do {
	undef $gotnew;
	for $ii(@{$seeds}) {
	    next if isAlready($ii->[$place{'gene'}], $included);
	    undef $ni; undef $dist; undef $max;
	    for $jj(@{$included}) {
		undef $nt; undef $toadd;
		print "$ii->[$place{'gene'}] vs $jj->[$place{'gene'}] \n" if $debug;
		if ($ii->[$place{'species'}] ne $jj->[$place{'species'}]) {
		    $toadd = $distances->{$ii->[$place{'gene'}]}->{$jj->[$place{'gene'}]};
		}
		else {
		    for $kk(keys(%{$distances->{$ii->[$place{'gene'}]}->{$jj->[$place{'gene'}]}})) {
			$toadd += $distances->{$ii->[$place{'gene'}]}->{$jj->[$place{'gene'}]}->{$kk}; $nt++;
		    }
		    $toadd = $toadd / $nt if $nt;
		}
		$dist += $toadd; $ni++ if $toadd;
		print "A distance for $ii->[$place{'gene'}] and $jj->[$place{'gene'}] does not exist...\n" if (!$toadd and $debug);
	    }
	    if ($ni and ($dist / $ni > $max->['IPscore'])) {
		$max = $ii; $max->['IPscore'] = $dist / $ni;
	    }
	    if ($max->['IPscore'] >= $cco) {
		print "Included $max->[4] as $dist / $ni > $cco\n" if $debug;
		push(@{$included}, $max);
		$gotnew = 1;
	    }
	}
    } while ($gotnew);
    return($included);
}

sub array2StdOut { #for debug purpose only
    my($ar, $name) = @_;
    my($ii, $jj);

    print "$name: \n";
    return(0) if !scalar(@{$ar});
    for $ii(@{$ar}) {
	for $jj(@{$ii}) {
	    print "$jj\t" if defined($jj);
	}
	print "\n";
    }
    return undef;
}
