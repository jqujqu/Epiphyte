#!/usr/bin/python2.7

import numpy as np
import argparse
import sys

def add_noise(l, mu, sigma, binary):
    l = l.split()
    s = abs(np.random.normal(mu, sigma, len(l)-2))
    for i in range(2,len(l)) :
        level = float(l[i])
        noise = min(s[i-2], 1.0)
        val =  level + noise if l[i] == '0' else level - noise
        l[i] = int(val>0.5) if binary else val
    return('\t'.join(str(e) for e in l))

def main():
    parser = argparse.ArgumentParser(description='Add noise to binary states')
    parser.add_argument('--sigma', dest='sigma', help='sigma of white noise')
    feature_parser = parser.add_mutually_exclusive_group(required=False)
    feature_parser.add_argument('--binary', dest='binary', help='return binary state', action='store_true')
    feature_parser.add_argument('--no-binary', dest='binary', help='return probability', action='store_false')
    parser.set_defaults(binary=True)
    parser.add_argument('--states', dest='states', help='input state file')
    parser.add_argument('--output', dest='output', help='output file')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    outfile = args.output
    out = open(outfile, 'w')
    sigma = float(args.sigma)
    binary = args.binary
    mu = 0.0
    with open(args.states) as f:
        first_line = f.readline()
        out.write(first_line)
        for line in f:
            s = add_noise(line, 0, sigma,binary)+ "\n"
            out.write(s)


if __name__ == "__main__":
    main()
