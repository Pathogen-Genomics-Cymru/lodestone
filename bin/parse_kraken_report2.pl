#!/usr/bin/env perl

use strict;
use warnings;
use JSON qw(decode_json);
use JSON::Create qw(create_json);

# REQUIREMENTS
if ( (!(defined($ARGV[0]))) or (!(defined($ARGV[1]))) or (!(defined($ARGV[2]))) or (!(defined($ARGV[3]))) )
	{ print "\n\nThis script will parse a Kraken output file to report all family/genus/species classifications in the sample, plus species complex classifications if the dominant family is Mycobacteriaceae.\n";
	  print "We require a min. coverage of x% of the total reads AND a min. number of reads PER CLASSIFICATION. Set these to 0 to report everything.\n";
	  print "USAGE:\tperl parse_kraken_report.pl [path to Kraken report] [path to output file; must end .json] [min. coverage, as %] [min. coverage, as no. of reads]\n";
	  print "E.G.:\tperl parse_kraken_report.pl report.txt out.json 1 10000\n\n\n";
	  exit 1;
	}
my $in_file = $ARGV[0]; my $out_file = $ARGV[1]; my $pct_threshold = $ARGV[2]; my $num_threshold = $ARGV[3];
if (!(-e($in_file)))
	{ die "ERROR: cannot find $in_file\n"; }
if ( (-e($in_file)) and (!(-s($in_file))) )
	{ die "ERROR: $in_file is empty\n"; }
if ($out_file !~ /^(.*?)\.json$/)
	{ die "ERROR: output file $out_file must end with suffix .json\n"; }
if ($pct_threshold !~ /^\d+\.?\d*$/)
	{ die "ERROR: $pct_threshold is not a positive number\n"; }
if ($num_threshold !~ /^\d+$/)
	{ die "ERROR: $num_threshold is not a positive integer\n"; }
if ($pct_threshold > 100)
	{ die "ERROR: $pct_threshold is a % and cannot be > 100\n"; }

# READ KRAKEN REPORT
# the Kraken report is assumed to be the standard format: 6 tab-delimited columns, with one line per taxon. This is described at https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats. We will confirm this as we parse.
my @S = (); my @G = (); my @G1 = (); my @F = ();
my $non_human_species_detected = 0;
open(IN,'<',$in_file) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  next if ($#line != 5); # in a correctly formatted Kraken report, each row will have 6 columns; we will skip those that don't
	  my $pc_frags 		   = $line[0]; # defined as "percentage of fragments covered by the clade rooted at this taxon". NOTE: for this purpose, 'fragment' is synonymous with 'read'.
	  my $num_frags_rooted = $line[1]; # defined as "number of fragments covered by the clade rooted at this taxon", i.e. all fragments at this taxonomic level AND LOWER
	  my $num_frags_direct = $line[2]; # defined as "number of fragments assigned directly to this taxon"
	  my $rank_code 	   = $line[3]; # defined as "rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies". Takes the form of one letter, optionally followed by one number.
	  my $ncbi_taxon_id    = $line[4]; # defined as "NCBI taxonomic ID number"
	  my $name 			   = $line[5]; # defined as "scientific name"
	  $pc_frags =~ s/^\s+//; $name =~ s/^\s+//;
	  next if (($pc_frags 		  < $pct_threshold) and ($name ne 'Homo sapiens')); # skip classifications not supported by a min % of fragments
	  next if (($num_frags_rooted < $num_threshold) and ($name ne 'Homo sapiens')); # skip classifications not supported by a min no. of fragments
	  if (($pc_frags =~ /^\d+\.?\d*$/) and ($num_frags_rooted =~ /^\d+$/) and ($num_frags_direct =~ /^\d+$/) and ($rank_code =~ /^\w{1}\d*$/) and ($ncbi_taxon_id =~ /^\d+$/) and ($name =~ /\w+/))
		{ if 	($rank_code eq 'S') { push(@S,[$num_frags_rooted,$pc_frags,$name,$ncbi_taxon_id]); $non_human_species_detected++ unless ($name eq 'Homo sapiens'); }
		  elsif ($rank_code eq 'G') { push(@G,[$num_frags_rooted,$pc_frags,$name,$ncbi_taxon_id]); }
		  elsif ($rank_code eq 'F') { push(@F,[$num_frags_rooted,$pc_frags,$name,$ncbi_taxon_id]); }
		  if (($name =~ /^Mycobact.*?$/) and ($rank_code eq 'G1')) # Kraken does not resolve classifications among the Mycobacteriaceae as well as Mykrobe. At best, it can detect species complexes. We shall retain these classifications to look at later, as they may indicate whether this is a mixed-mycobacterial sample.
			{ push(@G1,[$num_frags_rooted,$pc_frags,$name,$ncbi_taxon_id]); }
		}
	  else
		{ die "ERROR: malformatted Kraken report, at line $. (\"$line\")\n"; }
	}
close(IN) or die $!;

# CREATE OUTPUT FILE
open(OUT,'>',$out_file) or die $!;
my %out = ();
$out{'Thresholds'}{'reads'} = $num_threshold;
$out{'Thresholds'}{'percentage'} = $pct_threshold;
if (($#S == -1) or ($non_human_species_detected == 0))
	{ push(@{$out{'Warnings'}},"warning: no species classifications meet thresholds of > $num_threshold reads and > $pct_threshold % of total reads (human excepted)"); }
if ($#G == -1) { push(@{$out{'Warnings'}},"warning: no genus classifications meet thresholds of > $num_threshold reads and > $pct_threshold % of total reads");   }
if ($#F == -1) { push(@{$out{'Warnings'}},"warning: no family classifications meet thresholds of > $num_threshold reads and > $pct_threshold % of total reads");  }
my $top_family = ''; my $no_of_reads_assigned_to_top_family = 0;
my $top_genus = ''; my $top_species = '';
my $contaminant_species_found = 0;
for(my $x=0;$x<=2;$x++)
	{ my @arr = (); my $clade = '';
	  if    ($x == 0) { @arr = @F; $clade = 'Family';  }
	  elsif ($x == 1) { @arr = @G; $clade = 'Genus';   }
	  elsif ($x == 2) { @arr = @S; $clade = 'Species'; }
	  next if ($#arr == -1);
	  my @sorted_arr = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [$_, $_->[0]] } @arr;
	  if 	($x == 0) { $top_family  = $sorted_arr[0][2]; $no_of_reads_assigned_to_top_family = $sorted_arr[0][0]; }
	  elsif ($x == 1) { $top_genus   = $sorted_arr[0][2]; }
	  elsif ($x == 2) { $top_species = $sorted_arr[0][2]; }
	  for(my $y=0;$y<@sorted_arr;$y++)
		{ my %hash = ('reads' => $sorted_arr[$y][0], 'percentage' => $sorted_arr[$y][1], 'name' => $sorted_arr[$y][2], 'taxon_id' => $sorted_arr[$y][3]);
		  push(@{$out{$clade}},\%hash);
		  if (($x == 2) and ($y > 0) and ($sorted_arr[$y][2] ne 'Homo sapiens'))
			{ $contaminant_species_found++; } # raise a warning if a non-human species is detected that is NOT the top hit, as this indicates the sample is mixed or contaminated
		}
	}
if ($contaminant_species_found > 0)
	{ push(@{$out{'Warnings'}},'warning: sample is mixed or contaminated (contains reads from multiple non-human species)'); }
if ($top_family =~ /^Mycobact.*?$/)
	{ if (($top_genus !~ /^Mycobact.*?$/) or ($top_species !~ /^Mycobact.*?$/))
		{ push(@{$out{'Warnings'}},'warning: top family classification is mycobacterial, but this is not consistent with top genus and species classifications'); }
	}

# IF THE TOP FAMILY IS MYCOBACTERIACEAE (WHICH CAN ONLY BE THE CASE IF MINIMUM COVERAGE THRESHOLDS ARE MET), WE WILL ALSO REPORT THE KRAKEN 'G1' CLASSIFICATIONS. THESE MAY INDICATE WHETHER THIS IS A MIXED MYCOBACTERIAL SAMPLE.
if ($top_family eq 'Mycobacteriaceae')
	{ if ($no_of_reads_assigned_to_top_family < 100000)
		{ push(@{$out{'Warnings'}},"error: there are < 100k reads classified as Mycobacteriaceae");
		  $out{'Mykrobe'} = 'false';
		}
	  else
		{ $out{'Mykrobe'} = 'true'; } # as the sample is predominantly Mycobacteriaceae, we recommend the user invoke Mykrobe for higher-resolution classification. Later in the workflow, we will be using this value in a text comparison. It MUST be lower-case here otherwise it will be mistaken for a Boolean (TRUE/FALSE) instead.
	  if ($#G1 == -1)
		{ push(@{$out{'Warnings'}},"warning: top family is Mycobacteriaceae but no G1 (species complex) classifications meet thresholds of > $num_threshold reads and > $pct_threshold % of total reads (this is not necessarily a concern as not all mycobacteria have this taxonomy)"); }
	  else
		{ my @sorted_G1 = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [$_, $_->[0]] } @G1;
		  for(my $x=0;$x<@sorted_G1;$x++)
			{ my %hash = ('reads' => $sorted_G1[$x][0], 'percentage' => $sorted_G1[$x][1], 'name' => $sorted_G1[$x][2], 'taxon_id' => $sorted_G1[$x][3]);
			  push(@{$out{'Species complex'}},\%hash);
			}
		  my $num_G1 = @sorted_G1;
		  if ($num_G1 > 1)
			{ push(@{$out{'Warnings'}},'warning: sample contains multiple mycobacterial species complexes'); }
		}
	}
else
	{ push(@{$out{'Warnings'}},"error: top family is not Mycobacteriaceae");
	  $out{'Mykrobe'} = 'false';
	}
if (!(exists($out{'Warnings'})))
	{ push(@{$out{'Warnings'}},'');
	}
my $json = create_json(\%out);
print OUT JSON->new->ascii->pretty->encode(decode_json($json));
close(OUT) or die $!;
