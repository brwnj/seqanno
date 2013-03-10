#!/usr/bin/env python
# encoding: utf-8
import re
import sys
import tempfile
import itertools
from collections import Counter
from toolshed import nopen, reader, header

def bam2bedgraph(args):
    """Convert bam to bedgraph. Args: bedgraph, bam, strand"""
    cmd = "|bedtools genomecov -bg -5 -ibam %s" % (args.bam)
    if args.strand:
        cmd = "|bedtools genomecov -bg -5 -strand %s -ibam %s" % (args.strand, args.bam)
    result_header = "chrom start stop counts".split()
    for b in reader(cmd, header=result_header):
        print "\t".join(b[r] for r in result_header)

def read_fasta(fa):
    """parses fasta filehandle and returns name and sequence."""
    for header, group in itertools.groupby(fa, lambda line: line[0] == '>'):
        if header:
            line = group.next()
            name = line[1:].strip()
        else:
            seq = ''.join(line.strip() for line in group)
            yield name, seq

def get_locs(start, genes, gene_dict):
    """For each sequence position, returns position in gene as percent."""
    locs = []
    for gene in genes.split(","):
        if gene == ".":
            locs.append("na")
            continue
        g = gene_dict.get(gene)
        loc = 100 * ((start - g['start']) / (g['stop'] - g['start'] - 1.))
        if g['strand'] == "-":
            loc = 100 - loc
        locs.append(loc)
    return ",".join(map(str, locs))

def search(args):
    """Given fasta, gff, and bam, parses for sequence, annotates feature, and
    reports coverage.
    """
    match_seq = args.seq.upper()

    # write a temp bed of sequence match sites
    site_temp = open(tempfile.mktemp(suffix=".bed"), 'wb')
    with nopen(args.fasta) as fasta:
        for chrom, seq in read_fasta(fasta):            
            if args.verbose: sys.stderr.write(">> processing %s...\n" % chrom)
            # for each sequence match
            for i, m in enumerate([s.start() for s in re.finditer(match_seq, seq)]):
                start = m
                stop = start + 2
                name = "%s_%s_%d" % (chrom, match_seq, i)
                fields = [chrom, start, stop, name]
                site_temp.write("\t".join(map(str, fields)) + "\n")
    site_temp.close()

    # convert gff to bed with gene name as bed name field
    gff_temp = open(tempfile.mktemp(suffix=".bed"), 'wb')
    result_header = "chrom source feature start stop score strand frame attributes comments".split()
    # for filtering unique and storing start and stop for each gene
    genes = {}
    if args.verbose: sys.stderr.write(">> selecting %s from gff records...\n" % args.feature)
    for g in reader(args.gff, header=result_header):
        try:
            if not g['feature'] == args.feature: continue
            # regex gene name out
            gene_name = re.findall(r'Name=([\w\.]+)', g['attributes'])[0]
            # skip already seen
            if genes.has_key(gene_name): continue
            genes[gene_name] = {'start':int(g['start']), 'stop':int(g['stop']), 'strand':g['strand']}
            fields = [g['chrom'], g['start'], g['stop'], gene_name]
            gff_temp.write("\t".join(map(str, fields)) + "\n")
        except KeyError:
            if not g['chrom'].startswith("#"):
                sys.stderr.write("ERROR parsing gff!\n")
                sys.exit(1)
    gff_temp.close()

    # sort the gene bed, map and collapse genes onto site_temp, then add counts
    if args.verbose: sys.stderr.write(">> finding relative gene location per sequence match...\n")
    result_header = "chrom start stop name gene_name counts".split()
    cmd = "|sortBed -i %s | mapBed -a %s -b - -c 4 -o collapse | mapBed -a - -b %s -c 4 -o sum"\
            % (gff_temp.name, site_temp.name, args.bedgraph)
    for b in reader(cmd, header=result_header):
        # sequence position(s) relative to gene(s) it overlaps
        locs = get_locs(int(b['start']), b['gene_name'], genes)
        fields = [b['chrom'], b['start'], b['stop'], b['name'], b['gene_name'], b['counts'], locs]
        print "\t".join(map(str, fields))

def genestat(args):
    """Reads bed-like file from annotation and returns gene level stats."""
    lower_b, upper_b = args.bounds.split(",")
    total_c = Counter()
    dire_c = Counter()
    cutoff_c = Counter()
    bedlike_header = "chrom start stop name gene_name counts loc".split()
    for b in reader(args.bed, header=bedlike_header):
        if b['counts'] == ".": b['counts'] = 0
        for gene, loc in itertools.izip(b['gene_name'].split(","), b['loc'].split(",")):
            if gene == ".": continue
            # all sequence sites
            total_c.update([gene])
            # all that are deemed disruptive
            if float(upper_b) > float(loc) > float(lower_b):
                dire_c.update([gene])
            # passed cutoff
            if int(b['counts']) >= args.cutoff:
                cutoff_c.update([gene])
    for gene in sorted(total_c.keys()):
        fields = [gene, total_c.get(gene), dire_c.get(gene), cutoff_c.get(gene)]
        print "\t".join(map(str, fields))

def get_set(filein, column):
    fset = []
    for line in reader(filein, header=False):
        fset.append(line[column - 1])
    return set(fset)

def get_dict(filein, column):
    fdict = {}
    for line in reader(filein, header=False):
        fdict[line[column - 1]] = line
    return fdict

def compare(args):
    """Compare any column of two files. Args: a, b, column"""
    seta = get_set(args.a, args.column)
    setb = get_set(args.b, args.column)
    comparisons = {'aonly':seta - setb,
                   'bonly':setb - seta,
                   'intersection': seta & setb}
    if args.mode == 'bonly':
        lookup = get_dict(args.b, args.column)
    else:
        lookup = get_dict(args.a, args.column)
    for set_items in comparisons[args.mode]:
        print "\t".join(lookup[set_items])

def uniprot(args):
    """Add Uniprot annotation to gene list."""
    uniprot_db = {}
    uniprot_header = header(args.uniprotdb)
    for entry in reader(args.uniprotdb):
        for gene in entry['Gene names'].split():
            uniprot_db[gene] = entry
    for entry in reader(args.genes, header=False):
        uniprot_fields = []
        for gene in entry[int(args.column) - 1].split(","):
            uniprot = uniprot_db.get(gene)
            if uniprot:
                for h in uniprot_header:
                    uniprot_fields.append(uniprot[h])
        print "\t".join(entry) + "\t".join(map(str, uniprot_fields))