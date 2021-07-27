#!/usr/bin/env/perl

use strict;
use warnings 'all';
use Test2::V0;
use Test::Script;

script_compiles('../../bin/create_final_json.pl');
script_runs(['../../bin/create_final_json.pl', 'test_alignmentStats.json', 'test_species_in_sample.json']);
 
done_testing;
