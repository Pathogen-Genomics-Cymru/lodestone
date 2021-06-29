#!/usr/bin/perl

use strict;
use warnings;
use JSON qw(decode_json);
use JSON::Create qw(create_json);

# REQUIREMENTS
if ( (!(defined($ARGV[0]))) or (!(defined($ARGV[1]))) ) { print "ERROR: expecting two JSONs as input\n"; exit 1; }
my $stats_json = $ARGV[0]; my $report_json = $ARGV[1];
if ( (-e($stats_json))  and (!(-s($stats_json)))  ) { die "ERROR: $stats_json is empty\n";  exit 1; }
if ( (-e($report_json)) and (!(-s($report_json))) ) { die "ERROR: $report_json is empty\n"; exit 1; }
my $sample_id_STA = ''; my $sample_id_REP = '';
if ($stats_json  =~ /^.+\/(.*?)\_alignmentStats\.json$/)      { $sample_id_STA = $1; } elsif ($stats_json  =~ /^(.*?)\_alignmentStats\.json$/)     { $sample_id_STA = $1; }
if ($report_json =~ /^.+\/(.*?)\_species\_in\_sample\.json$/) { $sample_id_REP = $1; } elsif ($report_json =~ /^(.*?)_species\_in\_sample\.json$/) { $sample_id_REP = $1; }
if ($sample_id_STA ne $sample_id_REP) { die "ERROR: the sample IDs of $stats_json and $report_json are mismatched\n"; exit 1; }
my $sample_id = '';
if ($sample_id_STA eq $sample_id_REP) { $sample_id = $sample_id_STA; }
if ($sample_id eq '') { die "ERROR: could not identify sample ID from the filename of either $stats_json or $report_json\n"; exit 1; }

# OUTPUT
my $out_json = "$sample_id"."_report.json";
open(OUT,'>',$out_json) or die $!;

# READ JSON FILES
open(STA_JSON,'<',$stats_json) or die $!;
local $/;
my $stats = decode_json(<STA_JSON>);
close(STA_JSON);

open(REP_JSON,'<',$report_json) or die $!;
local $/;
my $report = decode_json(<REP_JSON>);
close(REP_JSON);

# WHAT % OF THE REFERENCE GENOME IS COVERED AT 10-FOLD DEPTH?
my $pc_10fold = $stats->{'bam_stats'}{'coverage_breadth'}{'10'};

# HOW MANY READS WERE MAPPED TO THE REFERENCE GENOME?
my $num_mapped = $stats->{'bam_stats'}{'reads_mapped'};

# WHAT IS THE AVERAGE READ MAPPING QUALITY?
my $avg_mapq = $stats->{'bam_stats'}{'average_quality'};

# ARE THERE ANY ERRORS?
my @errors = (); my $num_errors = 0;
if ($num_mapped < 100000) { $num_errors++; push(@errors,"error: < 100k reads could be mapped to the reference genome");			    	}
if ($pc_10fold  < 50) 	  { $num_errors++; push(@errors,"error: < 50% of the reference genome is covered at 10-fold depth");			}
if ($avg_mapq   < 10)     { $num_errors++; push(@errors,"error: alignments to the reference genome have average mapping quality < 10"); }

# IF THERE ARE NO EXISTING WARNINGS BUT WE'VE NOW RULED THERE ARE ERRORS, WE REVISE THE WARNING MESSAGE TO POINT THIS OUT
my @warnings = @{$report->{'warnings'}};
my $num_warnings = @warnings; my $no_warning_message = 0;
foreach my $warning (@warnings)
	{ if ($warning eq 'no warnings raised')
		{ $no_warning_message++; }
	}
if (($no_warning_message == 1) and ($num_warnings == 1) and ($num_errors > 0))
	{ my @warnings = ();
	  if ($num_errors == 1)
		{ @warnings = ("there was $num_errors error but no warnings");   }
	  elsif ($num_errors > 1)
		{ @warnings = ("there were $num_errors errors but no warnings"); }
	  @{$report->{'warnings'}} = @warnings;
	}

# CREATE OUTPUT JSON. IF THERE ARE ERRORS AT THIS POINT, WE REPORT THEM AND UPDATE THE 'CONTINUE TO CLOCKWORK?' FLAG
my %out = %{$report};
my %bam_stats  = %{$stats->{'bam_stats'}};
$out{'bam_stats'} = \%bam_stats;
if ($num_errors > 0)
	{ $out{'errors'} = \@errors;
	  $out{'summary_questions'}{'continue_to_clockwork'} = 'no';
	}

my $json = create_json(\%out);
print OUT JSON->new->ascii->pretty->encode(decode_json($json));
close(OUT) or die $!;