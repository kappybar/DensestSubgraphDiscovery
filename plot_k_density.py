import argparse
import matplotlib.pyplot as plt

"""
input file format
k_1 density_1
...
k_i density_i

plot (k, density) 
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="input filename")
    args = parser.parse_args()

    input_filename = args.file
    output_filename = input_filename.replace('.txt', '.png')

    ks = []
    densitys = []

    with open(input_filename, mode='r') as input_file:
        for line in input_file:
            k, density = line.split()
            k = int(float(k))
            density = float(density)
            ks.append(k)
            densitys.append(density)

    plt.plot(ks, densitys)
    index = input_filename.find('output') + len('output') + 1
    title = input_filename[index:-4]

    plt.title(title)
    plt.xlabel('k')
    plt.ylabel('triangle density')
    plt.savefig(output_filename)




if __name__ == "__main__":
    main()


