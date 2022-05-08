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

w : uniform(1, 1000)
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="input filename")
    args = parser.parse_args()

    input_filename = args.file
    output_filename = input_filename.replace('.txt', '_uniform_weighted.txt')
    output_file = open(output_filename, mode='w')

    with open(input_filename, mode='r') as input_file:
        first = True
        for line in input_file:
            if first:
                n, m = line.split()
                output_file.write(f'{n} {m}\n')
                first = False
            else:
                u, v = line.split()
                w = random.randint(1, 1000)
                output_file.write(f'{u} {v} {w}\n')

    output_file.close() 



if __name__ == "__main__":
    main()
