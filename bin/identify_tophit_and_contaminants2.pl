#!/usr/bin/env perl

use strict;
use warnings;
use JSON qw(decode_json);
use JSON::Create qw(create_json);

# REQUIREMENTS
if ( (!(defined($ARGV[0]))) or (!(defined($ARGV[1]))) or (!(defined($ARGV[2]))) or (!(defined($ARGV[3]))) or (!(defined($ARGV[4]))) )
	{ print "\n\nThis script will parse the Kraken and Mykrobe output JSONs to identify (a) the dominant mycobacterial species in the sample (based on highest % coverage, determined by Mykrobe), and (b) all other species, which are considered contaminants\n";
	  print "The output will be one JSON and one txt file, both produced in rundir. The JSON will state the Mykrobe-determined dominant species in the sample, as well as list each contaminant species, showing the taxon IDs, fasta URL, and both the 'assembly level' and 'genome representation' for each\n";
	  print "The text file just contains the URLs for each contaminant fasta\n";
	  print "For each species considered a contaminant, the latest RefSeq genomes are obtained as follows:\n";
	  print "If available, obtain all NCBI 'reference genomes' of that species, provided they are defined as 'complete'\n";
	  print "Else: all 'complete genomes' of that species, regardless of whether they are the reference\n";
	  print "Else: any genome of that species, and we warn that it may not be complete (which reduces confidence in contaminant removal)\n";
	  print "A 'reference genome' is a manually-selected community standard for that species. Note that some prokaryotes can have more than one reference genome\n";
	  print "USAGE:\tperl identify_tophit_and_contaminants.pl [path to Mykrobe JSON] [path to Kraken JSON] [path to RefSeq assembly summary file] [species] [unmix myco]\n";
	  print "[species] refers to what you believe this sample to be. You will be warned if this differs from the Kraken/Mykrobe predictions\n";
	  print "By defining [species] you will automatically select this to be the genome against which reads will be aligned using Clockwork\n";
	  print "[unmix myco] is either 'yes' or 'no', given in response to the question: do you want to disambiguate mixed-mycobacterial samples by read alignment?\n";
	  print "If 'no', any contaminating mycobacteria will be recorded but NOT acted upon\n";
	  print "E.G.:\tperl identify_tophit_and_contaminants.pl mykrobe_report.json mykrobe_report.json assembly_summary_refseq.txt 1 tuberculosis yes\n\n\n";
	  exit 1;
	}
my $mykrobe_json = $ARGV[0]; my $kraken_json = $ARGV[1]; my $assembly_file = $ARGV[2]; my $supposed_species = $ARGV[3]; my $unmix_myco = $ARGV[4];
if (!(-e($mykrobe_json)))  { die "ERROR: cannot find $mykrobe_json\n";  }
if (!(-e($kraken_json)))   { die "ERROR: cannot find $kraken_json\n";   }
if (!(-e($assembly_file))) { die "ERROR: cannot find $assembly_file\n"; } # from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
if ($supposed_species ne 'null')
	{ if (($supposed_species ne 'abscessus') and ($supposed_species ne 'africanum') and ($supposed_species ne 'avium') and ($supposed_species ne 'bovis') and ($supposed_species ne 'chelonae') and ($supposed_species ne 'chimaera') and ($supposed_species ne 'fortuitum') and ($supposed_species ne 'intracellulare') and ($supposed_species ne 'kansasii') and ($supposed_species ne 'tuberculosis'))
		{ die "ERROR: if you provide a species ID, it must be one of either: abscessus|africanum|avium|bovis|chelonae|chimaera|fortuitum|intracellulare|kansasii|tuberculosis\n"; exit 1; }
	}
if ( ($unmix_myco ne 'yes') and ($unmix_myco ne 'no') ) { die "ERROR: 'unmix myco' should be either 'yes' or 'no'\n"; exit 1; }
if ( (-e($mykrobe_json))  and (!(-s($mykrobe_json)))  ) { die "ERROR: $mykrobe_json is empty\n";  exit 1; }
if ( (-e($kraken_json))   and (!(-s($kraken_json)))   ) { die "ERROR: $kraken_json is empty\n";   exit 1; }
if ( (-e($assembly_file)) and (!(-s($assembly_file))) ) { die "ERROR: $assembly_file is empty\n"; exit 1; }
my $sample_id = '';
if ($mykrobe_json =~ /^.+\/(.*?)\_mykrobe\_report\.json$/) { $sample_id = $1; } elsif ($mykrobe_json =~ /^(.*?)\_mykrobe\_report\.json$/) { $sample_id = $1; }
if ($sample_id eq '') { die "ERROR: could not identify sample ID from the filename of $mykrobe_json\n"; exit 1; }

# OUTPUT
my %out = ();
my $out_file1 = "$sample_id"."_urllist.txt";
my $out_file2 = "$sample_id"."_species_in_sample.json";
open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!;

# PARSE THE 'ASSEMBLY SUMMARY' FILE TO STORE URLs FOR EACH GENOME ASSOCIATED WITH A GIVEN TAXON ID. WE REQUIRE THAT THIS IS THE LATEST VERSION OF THE GENOME AVAILABLE IN REFSEQ. WE WILL LATER CHECK THAT (A) IT IS DEFINED AS A "COMPLETE GENOME", AND (B) IT HAS FULL GENOME REPRESENTATION.
## RefSeq defines a 'complete genome' (see ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt) as follows: "all chromosomes are gapless and have no runs of 10 or more ambiguous bases (Ns), there are no unplaced or unlocalized scaffolds, and all the expected chromosomes are present (i.e. the assembly is not noted as having partial genome representation). Plasmids and organelles may or may not be included in the assembly but if present then the sequences are gapless"
## RefSeq defines 'genome representation' (also in the above readme) as "whether the goal for the assembly was to represent the whole genome or only part of it." This takes one of two values: 'full' or 'partial'
my %urls = (); my %tax_ids = ();
open(IN,'<',$assembly_file) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  next if ($line =~ /^#/);
	  my @line = split(/\t/,$line);
	  my $refseq_category = $line[4]; my $species_taxid = $line[6]; my $species = $line[7]; my $infraspecific_name = $line[8]; my $version_status = $line[10]; my $assembly_level = $line[11]; my $genome_rep = $line[13]; my $ftp_path = $line[19];
	  next if ($version_status ne 'latest');
	  next if ($species_taxid !~ /\d+/);
	  next if ($ftp_path eq 'na');
	  $tax_ids{$species} = $species_taxid;
	  push(@{$urls{$species_taxid}},[$species,$infraspecific_name,$ftp_path,$assembly_level,$genome_rep,$refseq_category]);
	}
close(IN) or die $!;

# READ JSON FILES
open(MYK_JSON,'<',$mykrobe_json) or die $!;
local $/;
my $mykrobe = decode_json(<MYK_JSON>);
close(MYK_JSON);

open(KRK_JSON,'<',$kraken_json) or die $!;
local $/;
my $kraken = decode_json(<KRK_JSON>);
close(KRK_JSON);

# WHAT IS THE TOP HIT MYCOBACTERIAL SPECIES IN THE SAMPLE, ACCORDING TO MYKROBE?
my @species = ();
while((my $species,my $irrel)=each(%{$mykrobe->{$sample_id}{'phylogenetics'}{'species'}}))
	{ my $pc_coverage  = $mykrobe->{$sample_id}{'phylogenetics'}{'species'}{$species}{percent_coverage};
	  my $median_depth = $mykrobe->{$sample_id}{'phylogenetics'}{'species'}{$species}{median_depth};
	  $species =~ s/\_/ /g;
	  push(@species,[$pc_coverage,$median_depth,$species]);
	}
my @sorted_species = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [$_, $_->[0]] } @species;
my $top_species = $sorted_species[0][2];
$out{'top_species'}{'name'} 		= $top_species;
$out{'top_species'}{'pc_coverage'}  = $sorted_species[0][0];
$out{'top_species'}{'median_depth'} = $sorted_species[0][1];
if ($sorted_species[0][0] < 75) { push(@{$out{'Warnings'}},"warning: Mykrobe's top species hit is supported by < 75% coverage"); 		}
if ($sorted_species[0][1] < 10) { push(@{$out{'Warnings'}},"warning: Mykrobe's top species hit has a median coverage depth < 10-fold"); }

# OTHER THAN THE TOP HIT, WHAT NON-HUMAN SPECIES ARE ALSO PRESENT IN THE SAMPLE, ACCORDING TO KRAKEN?
my %other_species = ();
foreach my $key (@{$kraken->{'Species'}})
	{ my $species = ${$key}{'name'};
	  my $taxid   = ${$key}{'taxon_id'};
	  next if ($species =~ /^Mycobact.*?$/); # ignore any Kraken hits to mycobacterial species - they may be spurious. We will use only the mycobacterial classifications made by Mykrobe
	  next if ($taxid == 9606); # ignore human because we have a dedicated human read removal process elsewhere in the workflow
	  if ($species ne $top_species)
		{ $other_species{$species} = $taxid; }
	}
	
# OTHER THAN THE TOP HIT, WHAT NON-HUMAN SPECIES ARE ALSO PRESENT IN THE SAMPLE, ACCORDING TO MYKROBE?
while((my $species,my $irrel)=each(%{$mykrobe->{$sample_id}{'phylogenetics'}{'species'}}))
	{ $species =~ s/\_/ /g;
	  ## Mykrobe does not assign a taxon ID to each species, so we will need to look this up. The taxon ID is the basis on which species' genomes are downloaded - we cannot proceed without it.
	  if (!(exists($tax_ids{$species})))
		{ push(@{$out{'Warnings'}},"warning: unable to find a taxon ID for $species, which means we will not be able to locate its genome, and thereby remove it as a contaminant"); }
	  next if (!(exists($tax_ids{$species})));
	  my $taxid = $tax_ids{$species};
	  next if ($taxid == 9606); # ignore human because we have a dedicated human read removal process elsewhere in the workflow
	  if ($species ne $top_species)
		{ $other_species{$species} = $taxid; }
	}

# IDENTIFY GENOMES FOR EACH MEMBER OF THE NON-REDUNDANT LIST OF NON-TOP-HIT & NON-HUMAN SPECIES PRESENT IN THE SAMPLE
@species = ();
while((my $species,my $taxid)=each(%other_species))
	{ push(@species,$species); }
@sorted_species = sort {$a cmp $b} @species;
my %ignored_mixed_myco = ();
foreach my $species (@sorted_species)
	{ my $taxid = $other_species{$species};
	  if (!(exists($urls{$taxid})))
		{ push(@{$out{'Warnings'}},"warning: unable to find the latest RefSeq genome for $taxid, and thereby remove it as a contaminant"); }
	  next if (!(exists($urls{$taxid})));
	  
	  # for each contaminating species, we will identify the set of genomes available in RefSeq
	  # if one or more 'reference genomes' for this species are available and they are 'complete' (thereby implicitly having 'full' representation), we will use these
	  # else, we will use all complete genome available
	  # else, we will use all genomes that are available (which may simply be contigs) and warn the user that while we can proceed with contaminant removal using them, we will have reduced confidence in detecting reads from this species
	  my @arr = @{$urls{$taxid}};
	  my %reference_genomes = ();
	  for(my $x=0;$x<@arr;$x++)
		{ my $refseq_species = $arr[$x][0]; my $infraspecific_name = $arr[$x][1]; my $ftp_path = $arr[$x][2]; my $assembly_level = $arr[$x][3]; my $genome_rep = $arr[$x][4]; my $refseq_category = $arr[$x][5];
		  if (($assembly_level eq 'Complete Genome') and ($genome_rep eq 'Full') and ($refseq_category eq 'reference genome'))
			{ $reference_genomes{$ftp_path}++; }
		}
	  my %complete_genomes = ();
	  for(my $x=0;$x<@arr;$x++)
		{ my $refseq_species = $arr[$x][0]; my $infraspecific_name = $arr[$x][1]; my $ftp_path = $arr[$x][2]; my $assembly_level = $arr[$x][3]; my $genome_rep = $arr[$x][4]; my $refseq_category = $arr[$x][5];
		  next if ( (scalar keys %reference_genomes != 0) and (!(exists($reference_genomes{$ftp_path}))) ); # checkpoint: if we have found reference genomes (scalar keys %refernece_genomes != 0) then we will ignore any genomes that are not (i.e. not in the %reference_genomes hash)
		  if (($assembly_level eq 'Complete Genome') and ($genome_rep eq 'Full'))
			{ $complete_genomes{$ftp_path}++; }
		}
	  for(my $x=0;$x<@arr;$x++)
		{ my $refseq_species = $arr[$x][0]; my $infraspecific_name = $arr[$x][1]; my $ftp_path = $arr[$x][2]; my $assembly_level = $arr[$x][3]; my $genome_rep = $arr[$x][4]; my $refseq_category = $arr[$x][5];
		  my $filename = ''; if ($ftp_path =~ /^.+\/(.*?)$/) { $filename = $1; }
		  if ($filename eq '')
			{ push(@{$out{'Warnings'}},"warning: unable to parse FTP path to contaminant species, and thereby download it: $ftp_path"); }
		  next if ($filename eq '');
		  my $full_path = "$ftp_path/$filename"."_genomic.fna.gz";
#		  if ($species ne $refseq_species)
#			{ push(@{$out{'Warnings'}},"warning: RefSeq species name ($refseq_species) differs from that in the Kraken/Mykrobe reports ($species), although has the same taxon ID"); }
#		  next if ($species ne $refseq_species); # a confirmatory check that the species name in the assembly_summary file is identical to the species name from the Kraken and Mykrobe reports. This will NOT always be the case if we take the species_taxid, not the taxid, from $assembly_file (i.e. $line[6], not $line[5]);
		  next if ( (scalar keys %complete_genomes != 0) and (!(exists($complete_genomes{$ftp_path}))) ); # checkpoint: if we have found complete genomes (scalar keys %complete_genomes != 0) then we will ignore any genomes that are incomplete (i.e. not in the %complete_genomes hash)
		  my %hash = ('infraspecific_name' => $infraspecific_name, 'url' => $full_path, 'assembly_level' => $assembly_level, 'genome_representation' => $genome_rep, 'refseq_category' => $refseq_category, 'taxid' => $taxid);
		  push(@{$out{'contaminants'}{$species}},\%hash) unless (($unmix_myco eq 'no') and ($top_species =~ /^Mycobacterium/) and ($species =~ /^Mycobacterium/));
		  print OUT1 "$full_path\n"						 unless (($unmix_myco eq 'no') and ($top_species =~ /^Mycobacterium/) and ($species =~ /^Mycobacterium/));
		  if (($unmix_myco eq 'no') and ($top_species =~ /^Mycobacterium/) and ($species =~ /^Mycobacterium/))
			{ $ignored_mixed_myco{$species}++; }
		}
	  if (scalar keys %complete_genomes == 0)
		{ push(@{$out{'Warnings'}},"warning: no complete genome was found for $species. We can proceed with contaminant removal, but with reduced confidence in detecting reads from this species"); }
	}
close(OUT1) or die $!;

# ARE WE DELIBERATELY IGNORING MIXED-MYCOBACTERIAL CONTENT, AND CONSIDERING THESE *NOT* TO BE CONTAMINANTS?
while((my $name,my $irrel)=each(%ignored_mixed_myco))
	{ push(@{$out{'Warnings'}},"warning: sample is mixed-mycobacterial, its principal species ($top_species) mixed with reads from $name. Paths to these genomes are not included in this file because you have chosen NOT to take an alignment-based approach to resolving this issue"); }

# DO WE DETECT CONTAMINANTS?
if (exists($out{'contaminants'}))
	{ $out{'ContaminantsToRemove'} = 'yes'; }
else
	{ $out{'ContaminantsToRemove'} = 'no'; }

# IS THE TOP SPECIES HIT ONE OF THE 10 ACCEPTABLE POSSIBILITIES?
if ($top_species =~ /^Mycobacterium (abscessus|africanum|avium|bovis|chelonae|chimaera|fortuitum|intracellulare|kansasii|tuberculosis).*?$/)
	{ my $identified_species = $1;
	  if ($supposed_species eq 'null')
		{ $out{'AcceptableSpecies'} = 'yes'; }
	  elsif (($supposed_species ne 'null') and ($supposed_species eq $identified_species))
		{ $out{'AcceptableSpecies'} = 'yes'; }
	  elsif (($supposed_species ne 'null') and ($supposed_species ne $identified_species))
		{ push(@{$out{'Warnings'}},"warning: the top species hit is $identified_species, contrary to the expectation: $supposed_species");
		  $out{'AcceptableSpecies'} = 'no';
		}
	}
else
	{ $out{'AcceptableSpecies'} = 'no'; }

my $json = create_json(\%out);
print OUT2 JSON->new->ascii->pretty->encode(decode_json($json));
close(OUT2) or die $!;
