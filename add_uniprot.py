#!/usr/bin/env python
# encoding: utf-8
"""
add uniprot annotation to gene list.
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
    db, uniprot_header = todict(args.UNIPROT)
    header = ["genename"]
    for gene in reader(args.GENELIST, header=header):
        # there is only one gene name per entry in the ta file
        e = db.get(gene['genename'])
        fields = []
        for h in header:
            fields.append(gene[h])
        if e:
            for h in uniprot_header:
                fields.append(e[h])
        print "\t".join(map(str, fields))

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('GENELIST', help='list of gene names')
    p.add_argument('UNIPROT', help='tab delimited uniprot db')
    args = p.parse_args()
    main(args)