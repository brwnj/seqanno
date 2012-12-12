#!/usr/bin/env python
# encoding: utf-8
"""
Given fasta, gff, and bam, parses for TA, annotates feature, and find coverage.
"""
import re
import sys
import tempfile
from toolshed import nopen
from subprocess import Popen, PIPE

def fasta_parser(fasta):
    name, seq = None, []
    for line in fasta:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line.lstrip(">"), []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def getlocation(gff, chrom, start, stop):
    ta = open(tempfile.mktemp(suffix=".bed"), 'w')
    ta.write('%s\t%d\t%d\n' % (chrom, start, stop))
    ta.close()

    p = Popen(["intersectBed", "-wb", "-a", ta.name, "-b", gff], stdout=PIPE, shell=False)
    location = 0
    genename = "."
    for e in p.stdout:
        e = e.strip("\r\n").split("\t")
        genestart = int(e[6])
        genestop = int(e[7])
        geneattrs = e[11]
        genename = re.findall(r'Name=([\w\.]+)', geneattrs)[0]
        location = 100 * ((start - genestart) / (genestop - genestart - 1.))
    p.wait()
    return genename, "%.0f%%" % location

def main(args):
    counter = 0
    for chrom, seq in fasta_parser(nopen(args.fasta)):
        # for each TA in this sequence
        match = [m.start() for m in re.finditer('TA', seq)]
        for m in match:
            if counter > 1 and counter % 1000 == 0:
                sys.stderr.write(">> processed %d sites\n" % counter)
            name = "TA_%d" % counter
            start = m
            stop = m + 2
            coverage = 0
            genename, location = getlocation(args.gff, chrom, start, stop)
            # count the number of reads starting at the TA site
            
            samfile = Popen(["samtools", "view", args.bam, "%s:%d-%d" % (chrom, start, stop)], stdout=PIPE, shell=False)
            for sam in samfile.stdout:
                readstart = int(sam.split("\t")[3])
                if readstart == start + 1:
                    coverage += 1
            samfile.wait()
            counter += 1
            fields = (chrom, start, stop, name, genename, location, coverage)
            print "\t".join(map(str, fields))

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('fasta')
    p.add_argument('gff')
    p.add_argument('bam')
    args = p.parse_args()
    main(args)