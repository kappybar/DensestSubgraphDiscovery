from scipy.stats import powerlaw
import argparse
import random

"""
input file format
N M
v_1 u_1
...
v_m u_m

output file format
N M
v_1 u_1 w_1
...
v_m u_m w_m

w : powerlaw distribution
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="input filename")
    args = parser.parse_args()

    input_filename = args.file
    output_filename = input_filename.replace('.txt', '_powerlaw_weighted.txt')
    output_file = open(output_filename, mode='w')
    rv = None
    i = 0

    with open(input_filename, mode='r') as input_file:
        first = True
        for line in input_file:
            if first:
                n, m = line.split()
                rv = powerlaw.rvs(a=0.1, loc=1, scale=1000, size=int(m))
                output_file.write(f'{n} {m}\n')
                first = False
            else:
                u, v = line.split()
                w = round(rv[i])
                i += 1
                output_file.write(f'{u} {v} {w}\n')

    output_file.close() 



if __name__ == "__main__":
    main()
