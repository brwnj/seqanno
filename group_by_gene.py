#!/usr/bin/env python
# encoding: utf-8
"""
total_ta: number of ta sites in that gene
dire_ta: number of ta sites > 5% and < 80%
cutoff_ta: number of ta sites per gene > cutoff 
"""
import sys
from toolshed import reader
from collections import Counter

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def counts(ta, cutoff):
    cgene = Counter()
    cdire = Counter()
    ccutoff = Counter()
    for b in reader(args.ta, header="chrom start stop name gene pos count".split()):
        if b['gene'] == ".": continue
        # all ta sites
        cgene.update([b['gene']])
        
        # all that are deemed dire
        if 80 > int(b['pos'].rstrip("%")) > 5:
            cdire.update([b['gene']])
        
        # all that passed cutoff
        if int(b['count']) >= cutoff:
            ccutoff.update([b['gene']])
    return cgene, cdire, ccutoff


def main(args):
    totalcounts, direcounts, cutoffcounts = counts(args.ta, args.cutoff)
    header = ("gene", "total_ta", "dire_ta", "cutoff_ta")
    print "\t".join(header)
    for gene in sorted(totalcounts.keys()):
        fields = (gene, totalcounts.get(gene), direcounts.get(gene), cutoffcounts.get(gene))
        print "\t".join(map(str, fields))


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('ta', help='annotated ta bed file')
    p.add_argument('--cutoff', '-c', default=5, type=int, help="inclusive read count cutoff value (default: 5)")
    args = p.parse_args()
    main(args)