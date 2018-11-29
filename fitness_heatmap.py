#!/usr/bin/env python3

import argparse
import itertools
import pandas as pd
import pickle


def heatmap_from_dataframe(dataframe, filename='heatmap.png'):
    import seaborn
    heatmap = seaborn.heatmap(dataframe, center=0)
    figure = heatmap.get_figure()
    figure.savefig(filename, dpi=300)


def fitness_dict_from_text_file(fitness_text_file):
    """parses text file of one tuple and one float per line into a dictionary

    Args:
        fitness_text_file: txt file with each line formatted as follows
            (integer, 'Capital single character amino acid'),fitness_value

    Returns:
        A dict with keys of tuples (position, amino acid char) and values of variant fitness floats

    Raises:

    """
    fitness_dict = {}
    with open(fitness_text_file, 'r') as f:
        for line in f:
            position, amino_acid, fitness = line.rstrip().split(',')
            position = int(position.split('(')[-1])
            amino_acid = amino_acid.split("'")[1]
            fitness = float(fitness)
            fitness_dict[(position, amino_acid)] = fitness
    return fitness_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script which returns some preliminary stats on mapping randomized 
        barcodes to expected library sequences""")
    parser.add_argument('-r', '--position_range',
                        help='range of positions to be plotted in heatmap. first and last position separated by dash')
    parser.add_argument('-n', '--name', help='name of heatmap image')

    required = parser.add_argument_group('required')
    required.add_argument("-p", "--pickle_file", required=True,
                          help="input pickle dictionary with keys as tuple of postion, amino acid and values "
                               "as fitness floats")
    wt_seq = 'MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVV' \
             'TGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'
    args = parser.parse_args()
    if args.pickle_file.endswith('.pkl'):
        with open(args.pickle_file, 'rb') as f:
            variant_fitness_dict = pickle.load(f)
    elif args.pickle_file.endswith('.txt'):
        variant_fitness_dict = fitness_dict_from_text_file(args.pickle_file)
    else:
        raise IOError('Please give fitness dictionary as a pickle file')

    if args.position_range:
        first_position, last_position = map(int, args.position_range.split('-'))
    else:
        first_position = 1
        last_position = 140
    position_range = range(first_position, last_position + 1)
    fitness_df = pd.DataFrame(0, index=position_range, columns=list('AVILMFYWSTNQHKRDECGP'))

    # Using Robert's fitness dictionary, which does not normalize to wt fitness
    wt_fitness = variant_fitness_dict[(0, 'WT')]
    for position, amino_acid in itertools.product(position_range, list('AVILMFYWSTNQHKRDECGP')):
        if (position, amino_acid) in variant_fitness_dict.keys():
            fitness = variant_fitness_dict[(position, amino_acid)] - wt_fitness
        elif amino_acid != wt_seq[position - 1]:
            fitness = float('NaN')
        else:
            fitness = 0
        fitness_df.loc[position, amino_acid] = fitness

    heatmap_from_dataframe(fitness_df.T, args.name)
