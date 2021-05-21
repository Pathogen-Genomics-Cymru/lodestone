#!/usr/bin/env perl

use strict;
use warnings;
use JSON qw(decode_json);
use JSON::Create qw(create_json);

# REQUIREMENTS
if ( (!(defined($ARGV[0]))) or (!(defined($ARGV[1]))) or (!(defined($ARGV[2]))) or (!(defined($ARGV[3]))) or (!(defined($ARGV[4]))) or (!(defined($ARGV[5]))) or (!(defined($ARGV[6]))) )
	{ print "\n\nThis script will parse the Kraken and Mykrobe output JSONs to identify (a) the dominant mycobacterial species in the sample (based on highest % coverage, determined by Mykrobe), and (b) all other species, which are considered contaminants\n";
	  print "The output will be one JSON and one txt file, both produced in rundir. The JSON will state the Mykrobe-determined dominant species in the sample, as well as list each contaminant species, showing the taxon IDs, fasta URL, and both the 'assembly level' and 'genome representation' for each\n";
	  print "The text file just contains the URLs for each contaminant fasta\n";
	  print "For each species considered a contaminant, the latest RefSeq genomes are obtained as follows:\n";
	  print "If available, obtain all NCBI 'reference genomes' of that species, provided they are defined as 'complete'\n";
	  print "Else: all 'complete genomes' of that species, regardless of whether they are the reference\n";
	  print "Else: any genome of that species, and we warn that it may not be complete (which reduces confidence in contaminant removal)\n";
	  print "A 'reference genome' is a manually-selected community standard for that species. Note that some prokaryotes can have more than one reference genome\n";
	  print "USAGE:\tperl identify_tophit_and_contaminants.pl [path to Mykrobe JSON] [path to Kraken JSON] [path to RefSeq assembly summary file] [species] [unmix myco] [directory containing mycobacterial reference genomes]\n";
	  print "[species] refers to what you believe this sample to be. You will be warned if this differs from the Kraken/Mykrobe predictions\n";
	  print "By defining [species] you will automatically select this to be the genome against which reads will be aligned using Clockwork\n";
	  print "[unmix myco] is either 'yes' or 'no', given in response to the question: do you want to disambiguate mixed-mycobacterial samples by read alignment?\n";
	  print "If 'no', any contaminating mycobacteria will be recorded but NOT acted upon\n";
	  print "E.G.:\tperl identify_tophit_and_contaminants.pl mykrobe_report.json mykrobe_report.json assembly_summary_refseq.txt 1 tuberculosis yes myco_dir\n\n\n";
	  exit 1;
	}
my $mykrobe_json = $ARGV[0]; my $kraken_json = $ARGV[1]; my $assembly_file = $ARGV[2]; my $supposed_species = $ARGV[3]; my $unmix_myco = $ARGV[4]; my $myco_dir = $ARGV[5]; my $prev_species_json = $ARGV[6];
if (!(-e($mykrobe_json)))  { die "ERROR: cannot find $mykrobe_json\n";  }
if (!(-e($kraken_json)))   { die "ERROR: cannot find $kraken_json\n";   }
if (!(-e($assembly_file))) { die "ERROR: cannot find $assembly_file\n"; } # from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
if (!(-d($myco_dir)))	   { die "ERROR: cannot find $myco_dir\n";		}
if ($prev_species_json ne 'null')
	{ if (!(-e($prev_species_json)))
		{ die "ERROR: cannot find $prev_species_json\n"; exit 1; }
	}
my @species = (qw/abscessus africanum avium bovis chelonae chimaera fortuitum intracellulare kansasii tuberculosis/);
foreach my $species (@species)
	{ if (!(-e("$myco_dir/$species.fasta")))
		{ die "ERROR: cannot find $myco_dir/$species.fasta"; }
	  if (!(-e("$myco_dir/$species.mmi")))
		{ die "ERROR: cannot find $myco_dir/$species.mmi"; }
	}
if ($supposed_species ne 'null')
	{ if (($supposed_species ne 'abscessus') and ($supposed_species ne 'africanum') and ($supposed_species ne 'avium') and ($supposed_species ne 'bovis') and ($supposed_species ne 'chelonae') and ($supposed_species ne 'chimaera') and ($supposed_species ne 'fortuitum') and ($supposed_species ne 'intracellulare') and ($supposed_species ne 'kansasii') and ($supposed_species ne 'tuberculosis'))
		{ die "ERROR: if you provide a species ID, it must be one of either: abscessus|africanum|avium|bovis|chelonae|chimaera|fortuitum|intracellulare|kansasii|tuberculosis\n"; exit 1; }	  
	}
if ( ($unmix_myco ne 'yes') and ($unmix_myco ne 'no') ) { die "ERROR: 'unmix myco' should be either 'yes' or 'no'\n"; exit 1; }
if ( (-e($mykrobe_json))  and (!(-s($mykrobe_json)))  ) { die "ERROR: $mykrobe_json is empty\n";  exit 1; }
if ( (-e($kraken_json))   and (!(-s($kraken_json)))   ) { die "ERROR: $kraken_json is empty\n";   exit 1; }
if ( (-e($assembly_file)) and (!(-s($assembly_file))) ) { die "ERROR: $assembly_file is empty\n"; exit 1; }
my $sample_id_MYK = ''; my $sample_id_KRK = '';
if ($mykrobe_json =~ /^.+\/(.*?)\_mykrobe\_report\.json$/) { $sample_id_MYK = $1; } elsif ($mykrobe_json =~ /^(.*?)\_mykrobe\_report\.json$/) { $sample_id_MYK = $1; }
if ($kraken_json  =~ /^.+\/(.*?)\_kraken\_report\.json$/)  { $sample_id_KRK = $1; } elsif ($kraken_json  =~ /^(.*?)\_kraken\_report\.json$/)  { $sample_id_KRK = $1; }
if ($sample_id_MYK ne $sample_id_KRK) { die "ERROR: the sample IDs of $mykrobe_json and $kraken_json are mismatched\n"; exit 1; }
my $sample_id = '';
if ($sample_id_MYK eq $sample_id_KRK) { $sample_id = $sample_id_MYK; }
if ($prev_species_json ne 'null')
	{ my $sample_id_PRE = '';
	  if ($prev_species_json =~ /^.+\/(.*?)\_species\_in\_sample\_previous\.json$/) { $sample_id_PRE = $1; } elsif ($prev_species_json =~ /^(.*?)\_species\_in\_sample\_previous\.json$/) { $sample_id_PRE = $1; }
	  if ($sample_id ne $sample_id_PRE)
		{ die "ERROR: sample ID of the previous species JSON ($prev_species_json) does not match the sample ID we have from the Kraken and Mykrobe reports ($sample_id)\n"; exit 1; }
	}
if ($sample_id eq '') { die "ERROR: could not identify sample ID from the filename of either $mykrobe_json or $kraken_json\n"; exit 1; }

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

my $prev_species = '';
if (-e($prev_species_json))
	{ open(PREV_SPE_JSON,'<',$prev_species_json) or die $!;
	  local $/;
	  $prev_species = decode_json(<PREV_SPE_JSON>);
	  close(PREV_SPE_JSON);
	}

# WHAT IS THE TOP HIT MYCOBACTERIAL SPECIES IN THE SAMPLE, ACCORDING TO MYKROBE, AND ON THE BASIS OF % COVERAGE?
@species = (); my $mykrobe_finds_nothing = 0;
while((my $species,my $irrel)=each(%{$mykrobe->{$sample_id}{'phylogenetics'}{'species'}}))
	{ my $pc_coverage  = $mykrobe->{$sample_id}{'phylogenetics'}{'species'}{$species}{percent_coverage};
	  my $median_depth = $mykrobe->{$sample_id}{'phylogenetics'}{'species'}{$species}{median_depth};
	  $species =~ s/\_/ /g;
	  push(@species,[$pc_coverage,$median_depth,$species]);
	  if ($species eq 'Unknown')
		{ $mykrobe_finds_nothing++; }
	}
my @sorted_species = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [$_, $_->[0]] } @species;
my $pc_coverage_of_top_species = $sorted_species[0][0]; my $depth_of_top_species = $sorted_species[0][1]; my $top_species = $sorted_species[0][2];
$out{'top_hit'}{'name'} = $top_species;
my $num_mykrobe_species = @sorted_species;

# WE ARE ASSUMING THE TOP MYKROBE HIT, ON THE BASIS OF % COVERAGE, IS *ALSO* THE TOP HIT ON THE BASIS OF MEDIAN DEPTH. LET'S CONFIRM THIS, AND WARN IF THIS IS NOT THE CASE.
# IT MAY BE POSSIBLE THAT A SAMPLE IS A MIXTURE OF SPECIES X (99% COVERAGE AT 10-FOLD DEPTH) AND SPECIES Y (98% COVERAGE AT 11-FOLD DEPTH). IN THIS CASE, ON WHAT BASIS DO WE CHOOSE A TOP HIT, GIVEN WE HAVE TO CHOOSE *ONE*?
for(my $x=1;$x<@sorted_species;$x++)
	{ my $pc_coverage_of_contam_species = $sorted_species[$x][0]; my $depth_of_contam_species = $sorted_species[$x][1]; my $contam_species = $sorted_species[$x][2];
	  if ($depth_of_contam_species > $depth_of_top_species)
		{ push(@{$out{'warnings'}},"warning: the top species hit ($top_species) has the highest % coverage of all Mykrobe species classifications ($pc_coverage_of_top_species) and a median depth of $depth_of_top_species, but a contaminating species ($contam_species) - although with lower coverage ($pc_coverage_of_contam_species%) - has higher depth ($depth_of_contam_species)");
		}
	}

# OTHER THAN THE TOP HIT, WHAT NON-HUMAN SPECIES ARE ALSO PRESENT IN THE SAMPLE, ACCORDING TO KRAKEN?
my $no_of_human_reads = 0;
my %other_species = ();
foreach my $key (@{$kraken->{'Species'}})
	{ my $species = ${$key}{'name'};
	  my $taxid   = ${$key}{'taxon_id'};
	  my $reads   = ${$key}{'reads'};
	  if ($taxid == 9606)
		{ $no_of_human_reads += $reads; }
	  $species =~ s/Mycobacteriodes/Mycobacterium/; # Kraken sometimes uses "Mycobacteriodes" (e.g. Mycobacteriodes abscessus) whereas Mykrobe uses "Mycobacterium" for the same. We need to standardise this to prevent downstream errors
	  next if ($species =~ /^Mycobact.*?$/); # ignore any Kraken hits to mycobacterial species - they may be spurious. We will use only the mycobacterial classifications made by Mykrobe
	  next if ($taxid == 9606); # ignore human because we have a dedicated human read removal process elsewhere in the workflow
	  if ($species ne $top_species)
		{ $other_species{$species} = $taxid; }
	}
	
# OTHER THAN THE TOP HIT, WHAT NON-HUMAN SPECIES ARE ALSO PRESENT IN THE SAMPLE, ACCORDING TO MYKROBE?
while((my $species,my $irrel)=each(%{$mykrobe->{$sample_id}{'phylogenetics'}{'species'}}))
	{ $species =~ s/\_/ /g;
	  ## Mykrobe does not assign a taxon ID to each species, so we will need to look this up. The taxon ID is the basis on which species' genomes are downloaded - we cannot proceed without it.
	  if ((!(exists($tax_ids{$species}))) and ($species ne $top_species)) # warn if we cannot find a taxon ID - but for contaminant genomes only. Note that we're not going to look up a taxon ID for the top species, because we aren't planning to remove it.
		{ push(@{$out{'warnings'}},"warning: unable to find a taxon ID for '$species', which means we will not be able to locate its genome, and thereby remove it as a contaminant. Check the Kraken report to see how this species has been reported"); }
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
my %ignored_mixed_myco = (); my %contaminant_genera = ();
foreach my $species (@sorted_species)
	{ my $taxid = $other_species{$species};
	  next if ($taxid == 9606);
	  if (!(exists($urls{$taxid})))
		{ push(@{$out{'warnings'}},"warning: unable to find the latest RefSeq genome for taxon ID $taxid, and thereby remove it as a contaminant (the Kraken report assigns this taxon ID to species '$species')"); }
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
			{ push(@{$out{'warnings'}},"warning: unable to parse FTP path to the FASTA of a contaminant species genome ('$species'), this being necessary to download it: $ftp_path"); }
		  next if ($filename eq '');
		  my $full_path = "$ftp_path/$filename"."_genomic.fna.gz";
		  next if ( (scalar keys %complete_genomes != 0) and (!(exists($complete_genomes{$ftp_path}))) ); # checkpoint: if we have found complete genomes (scalar keys %complete_genomes != 0) then we will ignore any genomes that are incomplete (i.e. not in the %complete_genomes hash)
		  my %hash = ('infraspecific_name' => $infraspecific_name, 'url' => $full_path, 'assembly_level' => $assembly_level, 'genome_representation' => $genome_rep, 'refseq_category' => $refseq_category, 'taxid' => $taxid);
		  my $contaminant_genus = ''; my $contaminant_species = '';
		  if ($species =~ /^(.*?) (.+?)$/) { $contaminant_genus = $1; $contaminant_species = $2; }
		  $contaminant_genera{$contaminant_genus}{$contaminant_species}++ unless (($unmix_myco eq 'no') and ($top_species =~ /^Mycobacterium/) and ($species =~ /^Mycobacterium/));
		  push(@{$out{'contaminants'}{$species}},\%hash) unless (($unmix_myco eq 'no') and ($top_species =~ /^Mycobacterium/) and ($species =~ /^Mycobacterium/));
		  print OUT1 "$full_path\n"						 unless (($unmix_myco eq 'no') and ($top_species =~ /^Mycobacterium/) and ($species =~ /^Mycobacterium/));
		  if (($unmix_myco eq 'no') and ($top_species =~ /^Mycobacterium/) and ($species =~ /^Mycobacterium/))
			{ $ignored_mixed_myco{$species}++; }
		}
	  if (scalar keys %complete_genomes == 0)
		{ push(@{$out{'warnings'}},"warning: no complete genome was found for the contaminant species '$species'. Alignment-based read removal will not necessarily detect all reads from this species"); }
	}
close(OUT1) or die $!;

# ARE WE DELIBERATELY IGNORING MIXED-MYCOBACTERIAL CONTENT, CONSIDERING CLASSIFICATIONS TO OTHER MYCOBACTERIAL SPECIES TO *NOT* BE CONTAMINANTS?
while((my $name,my $irrel)=each(%ignored_mixed_myco))
	{ push(@{$out{'warnings'}},"warning: sample contains a mixture of mycobacteria, its principal species ($top_species) mixed with reads from $name. As you have chosen to ignore this (--unmix_myco no), $name is NOT being counted as a contaminant"); }

# HOW MANY CONTAMINATING GENERA ARE THERE?
my $no_of_contaminant_genera = scalar keys %contaminant_genera;
	
# DO WE DETECT ONE OR MORE CONTAMINANTS?
if (exists($out{'contaminants'}))
	{ $out{'summary_questions'}{'are_there_contaminants'} = 'yes'; }
else
	{ $out{'summary_questions'}{'are_there_contaminants'} = 'no'; }
	
# *DID* WE DETECT ONE OR MORE CONTAMINANTS, WHEN (IF) WE RAN THIS SCRIPT BEFORE? NOTE THAT THIS SCRIPT CAN BE CALLED UP TO TWICE IN THE WORKFLOW.
if ($prev_species ne '')
	{ my $contam_found = 0;
	  while((my $contam_species,my $irrel)=each(%{$prev_species->{'contaminants'}}))
		{ $contam_found++;
		  push(@{$out{'summary_questions'}{'contaminants_previously_removed'}},$contam_species);
		}
	  if ($contam_found > 0)
		{ $out{'summary_questions'}{'were_contaminants_removed'} = 'yes'; }
	  elsif ($contam_found == 0)
		{ $out{'summary_questions'}{'were_contaminants_removed'} = 'no'; }
	}
else
	{ $out{'summary_questions'}{'were_contaminants_removed'} = 'no'; }
	
# IS THE TOP SPECIES HIT ONE OF THE 10 ACCEPTABLE POSSIBILITIES? IF SO, PROVIDE A LINK TO THE REFERENCE GENOME AND TO THE MYKROBE PHYLOGENETIC AND RESISTANCE PREDICTIONS.
if ($top_species =~ /^Mycobacterium (abscessus|africanum|avium|bovis|chelonae|chimaera|fortuitum|intracellulare|kansasii|tuberculosis).*?$/)
	{ my $identified_species = $1;
	  if ($supposed_species eq 'null')
		{ $out{'summary_questions'}{'is_the_top_species_appropriate'} = 'yes'; }
	  elsif (($supposed_species ne 'null') and ($supposed_species eq $identified_species))
		{ $out{'summary_questions'}{'is_the_top_species_appropriate'} = 'yes'; }
	  elsif (($supposed_species ne 'null') and ($supposed_species ne $identified_species))
		{ push(@{$out{'warnings'}},"warning: the top species hit is $identified_species, contrary to the expectation: $supposed_species");
		  $out{'summary_questions'}{'is_the_top_species_appropriate'} = 'no';
		}
	  if ($out{'summary_questions'}{'is_the_top_species_appropriate'} eq 'yes')
		{ if ($out{'summary_questions'}{'are_there_contaminants'} eq 'yes')
			{ if ($no_of_contaminant_genera == 1)
				{ my $contaminating_genus = '';
				  while((my $genus,my $irrel)=each(%contaminant_genera))
					{ $contaminating_genus = $genus; }
				  if ($contaminating_genus eq 'Mycobacterium')
					{ my @contam_species = ();
					  while((my $contam_species,my $irrel)=each(%{$contaminant_genera{'Mycobacterium'}}))
						{ push(@contam_species,$contam_species); }
					  my @sorted_contam_species = sort {$a cmp $b} @contam_species;
					  my $contam_species = join(", ",@sorted_contam_species);
					  my $num_contam_species = @sorted_contam_species;
					  push(@{$out{'warnings'}},"warning: sample contains a mixture of mycobacteria, its principal species ($identified_species) mixed with reads from $num_contam_species other Mycobacterium: $contam_species");
					}
				  else
					{ push(@{$out{'warnings'}},"warning: the top species hit ($identified_species) is supported, but the sample shows signs of contamination with species from one other genus ($contaminating_genus)"); }
				}
			  elsif ($no_of_contaminant_genera > 1)
				{ my @contam_genus = (); my $contam_with_myco = 0;
				  while((my $genus,my $irrel)=each(%contaminant_genera))
					{ if ($genus eq 'Mycobacterium') { $contam_with_myco++; }
					  next if ($genus eq 'Mycobacterium');
					  push(@contam_genus,$genus);
					}
				  my @sorted_contam_genus = sort {$a cmp $b} @contam_genus;
				  my $non_myco_contam_genus = join(", ",@sorted_contam_genus);
				  my $num_non_myco_contam_genus = @sorted_contam_genus;
				  if ($contam_with_myco == 0)
					{ if ($num_non_myco_contam_genus == 1)
						{ push(@{$out{'warnings'}},"warning: the top species hit ($identified_species) is supported, but the sample shows signs of contamination with species from $num_non_myco_contam_genus other (non-mycobacterial) genus: $non_myco_contam_genus");  }
					  elsif ($num_non_myco_contam_genus > 1)
						{ push(@{$out{'warnings'}},"warning: the top species hit ($identified_species) is supported, but the sample shows signs of contamination with species from $num_non_myco_contam_genus other (non-mycobacterial) genera: $non_myco_contam_genus"); }
					}
				  elsif ($contam_with_myco > 0)
					{ if ($num_non_myco_contam_genus == 1)
						{ push(@{$out{'warnings'}},"warning: the top species hit ($identified_species) is supported, but the sample contains multiple mycobacteria and shows signs of contamination with species from $num_non_myco_contam_genus other (non-mycobacterial) genus: $non_myco_contam_genus");  }
					  elsif ($num_non_myco_contam_genus > 1)
						{ push(@{$out{'warnings'}},"warning: the top species hit ($identified_species) is supported, but the sample contains multiple mycobacteria and shows signs of contamination with species from $num_non_myco_contam_genus other (non-mycobacterial) genera: $non_myco_contam_genus"); }
					}
				}
			  else
				{ die "ERROR: sample is considered contaminated but the number of contaminant genera is $no_of_contaminant_genera\n"; exit 1; }
			}
		  my $ref_fa  = "$myco_dir/$identified_species.fasta";
		  my $ref_dir = "$myco_dir/$identified_species";
		  my $ref_mmi = "$myco_dir/$identified_species.mmi";
		  if (!(-e($ref_fa)))  { die "ERROR: cannot find $ref_fa\n";  exit 1; }
		  if (!(-d($ref_dir))) { die "ERROR: cannot find $ref_dir\n"; exit 1; }
		  if (!(-e($ref_mmi))) { die "ERROR: cannot find $ref_mmi\n"; exit 1; }
		  $out{'top_hit'}{'file_paths'}{'ref_fa'} = $ref_fa;
		  $out{'top_hit'}{'file_paths'}{'clockwork_ref_dir'} = $ref_dir;
		  my %susceptibility = %{$mykrobe->{$sample_id}{'susceptibility'}};
		  my %phylogenetics  = %{$mykrobe->{$sample_id}{'phylogenetics'}};
		  $out{'top_hit'}{'susceptibility'} = \%susceptibility;
		  $out{'top_hit'}{'phylogenetics'}  = \%phylogenetics;
		}
	}
else
	{ $out{'summary_questions'}{'is_the_top_species_appropriate'} = 'no'; }
	
# IN THIS WORKFLOW, MYKROBE WOULD ONLY BE CALLED IF KRAKEN CLASSIFIED > 100k READS AS MYCOBACTERIACEAE, SO FOR MYKROBE TO MAKE *NO CLASSIFICATION* SOMETHING SUSPECT IS GOING ON.
# WHAT IS LIKELY TO HAVE HAPPENED IS THAT THE ALIGNMENT-BASED DECONTAMINATION PROCESS HAS TRIED TO DISAMBIGUATE A MIXTURE OF VERY SIMILAR MYCOBACTERIA AND INADVERTENTLY REMOVED TOO MANY READS. THERE WILL BE NOTHING SUBSTANTIVE LEFT FOR MYKROBE TO CLASSIFY.
if (($num_mykrobe_species == $mykrobe_finds_nothing) and ($num_mykrobe_species == 1))
	{ if ($out{'summary_questions'}{'were_contaminants_removed'} eq 'yes')
		{ push(@{$out{'warnings'}},"warning: regardless of what Kraken reports, Mykrobe did not make a species-level mycobacterial classification. If this is a mixed-mycobacterial sample, then an alignment-based contaminant-removal process may not be appropriate. Suggestion: re-run with --unmix_myco 'no'");
		}
	  elsif ($out{'summary_questions'}{'were_contaminants_removed'} eq 'no')
		{ push(@{$out{'warnings'}},"warning: regardless of what Kraken reports, Mykrobe did not make a species-level mycobacterial classification");
		}
	}
	
# IF THE TOP HIT IS APPROPRIATE AND THERE ARE NO CONTAMINANTS, WE CAN CONTINUE TO RUN CLOCKWORK
if (($out{'summary_questions'}{'is_the_top_species_appropriate'} eq 'yes') and ($out{'summary_questions'}{'are_there_contaminants'} eq 'no'))
	{ $out{'summary_questions'}{'continue_to_clockwork'} = 'yes'; }
else
	{ $out{'summary_questions'}{'continue_to_clockwork'} = 'no'; }

# IF THE SAMPLE IS *NOT CURRENTLY CONSIDERED* CONTAMINATED BUT PREVIOUSLY WAS, WE SAY SO
if (($out{'summary_questions'}{'are_there_contaminants'} eq 'no') and ($out{'summary_questions'}{'were_contaminants_removed'} eq 'yes'))
	{ push(@{$out{'warnings'}},"warning: the sample showed signs of contamination, although the contaminating reads were considered successfully removed"); }
	
# IF NO WARNINGS HAVE BEEN RAISED, WE SAY SO
if (!(exists($out{'warnings'})))
	{ push(@{$out{'warnings'}},"no warnings raised"); }

my $json = create_json(\%out);
print OUT2 JSON->new->ascii->pretty->encode(decode_json($json));
close(OUT2) or die $!;
