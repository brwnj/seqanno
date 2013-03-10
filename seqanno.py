#!/usr/bin/env python
# encoding: utf-8
"""
Retrieves read counts originating at a given genomic sequence and allows
further characterization.
"""

import seqanno_functions as sqf

__version__ = '0.1'

def main(args):
    args.func(args)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    sp = p.add_subparsers(title="subcommands")
    
    # generate a bedgraph
    bg_p = sp.add_parser('bedgraph', help="Convert BAM to bedgraph.")
    bg_p.add_argument('bam', help="Alignments file in bam format.")
    bg_p.add_argument('--strand', choices=['+', '-'],
            help="Only count reads from specified strand.")
    bg_p.set_defaults(func=sqf.bam2bedgraph)
    
    # annotated a specified genomic sequence
    search_p = sp.add_parser('search',
            help="Find sequence sites, count, and annotate.")
    search_p.add_argument('bedgraph',
            help="Counts per genomic location. Obtained from `seqanno bedgraph`.")
    search_p.add_argument('fasta', help="Genomic reference sequence in fasta format.")
    search_p.add_argument('gff', help="Gene annotations in gff format. Gene \
            name will be taken from attribute 'Name'.")
    search_p.add_argument('seq', help="The sequence to characterize, e.g. TA.")
    search_p.add_argument('-f', '--feature', default="gene",
            help="Feature type to use from gff [ %(default)s ].")
    search_p.add_argument('--verbose', action='store_true', help="Maximum verbosity.")
    search_p.set_defaults(func=sqf.search)
    
    # gene level count statistics (groupby)
    genestat_p = sp.add_parser('genestat',
            description="Returns a gene's total read count, reads falling in \
            disruptive region, and reads passing cutoff.",
            help="Gene level count statistics.")
    genestat_p.add_argument('bed',
            help='Annotated bed-like file from `seqanno search`.')
    genestat_p.add_argument('cutoff', type=int,
            help="Disregard sites with fewer reads.")
    genestat_p.add_argument('bounds', help="Lower and upper bound of sequence \
            positions most likely to cause gene disruption. Expects values to \
            be comma separated, e.g. 5,80.")
    genestat_p.set_defaults(func=sqf.genestat)
    
    # compare columns of two headerless text files
    compare_p = sp.add_parser('compare', description="Compare columns within \
            any two headerless, tab-delimited text files. The result will only \
            include unique entries, so if something is listed twice, like the \
            same gene being listed more than once, one will be ignored.",
            help="Compare columns between two files.")
    compare_p.add_argument('a', help="Tab delimited text file. No header.")
    compare_p.add_argument('b', help="Tab delimited text file. No header.")
    compare_p.add_argument('mode', choices="intersection aonly bonly".split(),
            help='Desired output of the comparison.')
    compare_p.add_argument('column', type=int, help='The column to compare \
            between the files.')
    compare_p.set_defaults(func=sqf.compare)
    
    # annotate gene list using uniprot database
    uniprot_p = sp.add_parser('uniprot', description="Given Uniprot database \
            as text, appends Uniprot data for each gene.",
            help="Add Uniprot data to gene list.")
    uniprot_p.add_argument('genes', help='Headerless text file with genes.')
    uniprot_p.add_argument('uniprotdb', help="Uniprot database as text. \
            'Gene names' is one of the expected column headers.")
    uniprot_p.add_argument('column', type=int, help="Column containing gene \
            names. Their can be more than one name separated with commas in \
            the specified column.")
    uniprot_p.set_defaults(func=sqf.uniprot)

    main(p.parse_args())