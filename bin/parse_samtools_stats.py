#!/usr/bin/env python3

import argparse
import json
import pysam


def total_genome_size_from_bam(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    try:
        total_length = sum(bam.lengths)
    except:
        raise Exception(
            f"Error getting total of reference seq lengths from BAM file {bam_file}"
        )

    bam.close()
    return total_length


def coverage_hist_to_coverage_breadth(coverage_hist, genome_size, max_cov=100):
    cumulative_cov = {max_cov: 0}
    total_bases = 0
    for cov, depth in sorted(coverage_hist.items(), reverse=True):
        if cov >= max_cov:
            cumulative_cov[max_cov] += depth
        else:
            cumulative_cov[cov] = depth

        total_bases += depth

    depths = [1, 2, 3, 4, 5, 10, 15, 25, 50, 100]
    cov_breadth = {
        d: sum([v for k, v in cumulative_cov.items() if k >= d]) for d in depths
    }
    return {k: round(100 * v / genome_size, 2) for k, v in cov_breadth.items()}


def parse_samtools_stats_file(infile, genome_size):
    stats = {}
    wanted_keys = {
        "raw total sequences:": False,
        "reads mapped:": False,
        "reads duplicated:": False,
        "bases mapped (cigar):": False,
        "bases trimmed:": False,
        "error rate:": True,
        "average quality:": True,
        "insert size average:": True,
        "insert size standard deviation:": True,
        "inward oriented pairs:": False,
        "outward oriented pairs:": False,
        "pairs with other orientation:": False,
    }

    coverage_hist = {}

    with open(infile) as f:
        for line in f:
            if line.startswith("SN"):
                fields = line.rstrip().split("\t")
                if fields[1] in wanted_keys:
                    key = (
                        fields[1]
                        .replace(" ", "_")
                        .rstrip(":")
                        .replace("(", "")
                        .replace(")", "")
                    )
                    value = (
                        float(fields[2]) if wanted_keys[fields[1]] else int(fields[2])
                    )

                    stats[key] = value
            elif line.startswith("COV"):
                try:
                    _, _, pos, count = line.rstrip().split("\t")
                    pos = int(pos)
                    count = int(count)
                except:
                    raise Exception(
                        f"Error parsing this COV line in file {infile}: {line}"
                    )
                coverage_hist[pos] = count

    stats["coverage_breadth"] = coverage_hist_to_coverage_breadth(
        coverage_hist, genome_size
    )
    stats["genome_size"] = genome_size
    return stats


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Gather QC stats from BAM file and output from samtools stats, prints JSON to stdout",
        usage="%(prog)s <bam> <samtools_stats>",
    )
    parser.add_argument("bam", help="sorted indexed BAM file")
    parser.add_argument("samtools_stats", help="File made by samtools stats bam")
    options = parser.parse_args()
    genome_size = total_genome_size_from_bam(options.bam)
    stats = parse_samtools_stats_file(options.samtools_stats, genome_size)
    print(json.dumps(stats, sort_keys=True, indent=2))
