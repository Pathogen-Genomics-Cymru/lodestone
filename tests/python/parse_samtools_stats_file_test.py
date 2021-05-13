import os
import pytest
import sys


# The scripts for this repo are in the bin/ dir in the root of the repo.
# We want to import the functions we're testing from the script
# parse_samtools_stats.py. Add to PYTHONPATH (is a bit hacky, but this repo
# isn't set up to use python modules)
this_file_dir = os.path.dirname(os.path.abspath(__file__))
bin_dir = os.path.join(this_file_dir, os.pardir, os.pardir, "bin")
sys.path.insert(1, bin_dir)

import parse_samtools_stats


def test_total_genome_size_from_bam():
    infile = os.path.join(
        this_file_dir, "parse_samtools_stats_file_test.genome_size_from_bam.bam"
    )
    assert parse_samtools_stats.total_genome_size_from_bam(infile) == 700


def test_coverage_hist_to_coverage_breadth():
    coverage_hist = {
        1: 2,
        2: 20,
        3: 12,
        5: 6,
        10: 1,
        11: 3,
        16: 3,
        200: 1,
    }

    genome_size = 50
    got = parse_samtools_stats.coverage_hist_to_coverage_breadth(
        coverage_hist, genome_size
    )
    expect = {
        1: 96.0,
        2: 92.0,
        3: 52.0,
        4: 28.0,
        5: 28.0,
        10: 16.0,
        15: 8.0,
        25: 2.0,
        50: 2.0,
        100: 2.0,
    }
    assert got == expect

    del coverage_hist[200]
    coverage_hist[60] = 1
    expect[100] = 0.0
    got = parse_samtools_stats.coverage_hist_to_coverage_breadth(
        coverage_hist, genome_size
    )
    assert got == expect


def test_parse_samtools_stats_file():
    infile = os.path.join(this_file_dir, "parse_samtools_stats_file_test.stats.txt")
    genome_size = 2000
    got = parse_samtools_stats.parse_samtools_stats_file(infile, genome_size)
    expect = {
        "average_quality": 35.3,
        "bases_mapped_cigar": 4050,
        "bases_trimmed": 0,
        "error_rate": 0.007407407,
        "insert_size_average": 198.6,
        "insert_size_standard_deviation": 4.1,
        "inward_oriented_pairs": 27,
        "outward_oriented_pairs": 0,
        "pairs_with_other_orientation": 0,
        "raw_total_sequences": 54,
        "reads_duplicated": 2,
        "reads_mapped": 54,
        "coverage_breadth": {
            1: 83.65,
            2: 62.15,
            3: 30.3,
            4: 17.8,
            5: 5.5,
            10: 0.0,
            15: 0.0,
            25: 0.0,
            50: 0.0,
            100: 0.0,
        },
        "genome_size": 2000,
    }
    assert got == expect
