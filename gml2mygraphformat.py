import networkx as nx
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="input filename")
    args = parser.parse_args()

    G = nx.read_gml(args.file, label='id')
    output_filename = args.file.replace('gml', 'txt')

    with open(output_filename, mode='w') as f:
        minimum_id = min(G.nodes)
        N = len(G.nodes)
        M = len(G.edges)
        f.write(f'{N} {M}\n')
        for e in G.edges:
            f.write(f'{e[0] - minimum_id} {e[1] - minimum_id}\n')


if __name__ == "__main__":
    main()