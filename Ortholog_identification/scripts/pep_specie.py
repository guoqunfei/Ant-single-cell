#!/usr/bin/env python3
#-*- coding:utf-8 -*-


import sys


with open(sys.argv[1],'r') as f_in:
    for line in f_in:
        l = line.rstrip('\n')

        if l.startswith('>'):
            print('>Amel_' + l.lstrip('>'))
            continue
        print(l)
