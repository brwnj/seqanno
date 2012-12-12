#!/usr/bin/env python
# encoding: utf-8
"""
add uniprot annotation for each entry.
"""
import sys
from toolshed import reader, header

def todict(uniprot):
    d = {}
    h = header(uniprot)
    for u in reader(uniprot):
        for gn in u['Gene names'].split():
            d[gn] = u
    return d, h

def main(args):
    db, uheader = todict(args.UNIPROT)
    taheader = "chrom start stop name gname loc count".split()
    for ta in reader(args.TA_ANNO, header=taheader):
        # there is only one gene name per entry in the ta file
        e = db.get(ta['gname'])
        fields = []
        for h in taheader:
            fields.append(ta[h])
        if e:
            for h in uheader:
                fields.append(e[h])
        print "\t".join(map(str, fields))

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('TA_ANNO', help='result of annotate_ta.py')
    p.add_argument('UNIPROT', help='tab delimited uniprot db')
    args = p.parse_args()
    main(args)