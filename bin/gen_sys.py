#!/usr/bin/env python3

import argparse, sys
from random import randint

parser = argparse.ArgumentParser(description='Generate (solvable) multivariate system over GF(2).',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-n', dest='n', type=int, required=True,
          help='number of variables')
parser.add_argument('-m', dest='m', type=int, required=True,
          help='number of equations')
args = parser.parse_args()

n = args.n
m = args.m

print("""Galois Field : GF(2)
Number of variables (n) : {var_num}
Number of polynomials (m) : {eq_num}
Seed : 0
Order : graded reverse lex order

*********************""".format(var_num=n, eq_num=m))

sol = [randint(0,1) for i in range(n)]

sys.stderr.write(str(sol) + "\n")

for i in range(m):
  res = 0

  for j in range(n):
    for k in range(j+1):
      r = randint(0,1)
      print(r, end=' ')

      if r == 1:
        res ^= sol[j] & sol[k]

  for j in range(n):
    r = randint(0,1)
    print(r, end=' ')

    if r == 1:
      res ^= sol[j]

  print(str(res) + ";")

