#!/usr/bin/env/perl

use strict;
use warnings 'all';
use Test2::V0;
use Test::Script;

script_compiles('../../bin/parse_kraken_report2.pl');
script_runs(['../../bin/parse_kraken_report2.pl', 'test_kraken_report.txt', 'test_kraken_report.json', '0.5' ,'5000']);
 
done_testing;
