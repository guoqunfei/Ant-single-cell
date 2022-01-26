#!/usr/bin/env python3
#-*- coding:utf-8 -*-


import sys
import gzip


def main():

    pep_f, cor_f = sys.argv[1:3]

    cor_dt = {}
    with open(cor_f) as f_in:
        for line in f_in:
            line = line.split()
            cor_dt[line[1]] = line[0]

    with gzip.open(pep_f, 'rt') as f_in:
        for line in f_in:
            if line.startswith('>'):
                print(line.rstrip('\n') + ' gene:' + cor_dt[line.split()[0].lstrip('>')] + ' done')
            else:
                print(line.rstrip('\n'))


if __name__ == "__main__":
    main()