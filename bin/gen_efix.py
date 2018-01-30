#!/usr/bin/env python3
# usage: generate preprocess (external hybridation) configuration file
import sys

help_msg = """usage: ./%s <filename> <var num> <starting value> <range>
\texample: ./gen_efix.py extfix.txt 4 5 4""" % (sys.argv[0])

if len(sys.argv) < 5:
    sys.exit(help_msg)

fname = sys.argv[1]
var_num = int(sys.argv[2])
start = int(sys.argv[3])
rge = int(sys.argv[4])

with open(fname, 'w') as fp:
    fp.write('%d %d %d\n' % (var_num, rge, start))
