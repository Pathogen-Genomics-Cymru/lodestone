#!/usr/bin/env/perl

use strict;
use warnings 'all';
use Test2::V0;
use Test::Script;

script_compiles('../../bin/identify_tophit_and_contaminants2.pl');
script_runs(['../../bin/identify_tophit_and_contaminants2.pl', 'test_mykrobe_report.json', 'test_kraken_report.json', '../../resources/assembly_summary_refseq.txt' , 'null', 'no', '../../resources', 'null']);

done_testing;
