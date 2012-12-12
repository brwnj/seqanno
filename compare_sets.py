#!/usr/bin/env python
# encoding: utf-8
"""
Compare any column of two files.
 
python compare_sets.py -m intersection -col gene ta_anno1.txt ta_anno2.txt
"""
from toolshed import reader, header
 
def get_set(filein, column, sep="\t"):
    fset = []
    if column.isdigit():
        column = int(column) - 1
        for line in reader(filein, header=False, sep=sep):
            fset.append(line[column])
    else:
        for line in reader(filein, header=True, sep=sep):
            fset.append(line[column])
    return set(fset)
 
def get_dict(filein, column, sep="\t"):
    fdict = {}
    if column.isdigit():
        column = int(column) - 1
        for line in reader(filein, header=False, sep=sep):
            fdict[line[column]] = line
    else:
        for line in reader(filein, header=True, sep=sep):
            fdict[line[column]] = line
    return fdict
    
def main(args):
    args = get_args()
    seta = get_set(args.filea, args.column, args.separator)
    setb = get_set(args.fileb, args.column, args.separator)
    comparisons = {'aonly':seta - setb,
                   'bonly':setb - seta,
                   'intersection': seta & setb,
                   'union':seta | setb}
    if args.mode == 'bonly':
        lookup = get_dict(args.fileb, args.column, args.separator)
        head = header(args.fileb)
    else:
        lookup = get_dict(args.filea, args.column, args.separator)
        head = header(args.fileb)
    
    if args.column.isdigit():
        for setitem in comparisons[args.mode]:
            print "\t".join([i for i in lookup[setitem]])
    else:
        for setitem in comparisons[args.mode]:
            print "\t".join([lookup[setitem][i] for i in head])
 
if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, 
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('filea')
    p.add_argument('fileb')
    p.add_argument('-m', '--mode', required=True,
                    choices="intersection union aonly bonly".split(),
                    default='intersection', help='desired output [ intersection ]')
    p.add_argument('-col', '--column', required=True, 
                    help='column name if header, number if no header')
    p.add_argument('-sep', '--separator', default="\t", 
                    help='field delimiter [ tab ]')
    args = p.parse_args()
    main(args)