#!/usr/bin/env python3

import argparse
import numpy as np
import pickle

from fitness_heatmap import fitness_dict_from_text_file

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to average fitness values by position""")
    parser.add_argument('-r', '--position_range',
                        help='range of positions to to be averaged. first and last position separated by dash')

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

    # Using Robert's fitness dictionary, which does not normalize to wt fitness
    wt_fitness = variant_fitness_dict[(0, 'WT')]
    for position in position_range:
        position_fitness_list = []
        for amino_acid in 'AVILMFYWSTNQHKRDEGCP':
            if (position, amino_acid) in variant_fitness_dict.keys():
                fitness = variant_fitness_dict[(position, amino_acid)] - wt_fitness
            elif amino_acid != wt_seq[position - 1]:
                continue
            else:
                fitness = 0
            position_fitness_list.append(fitness)
        avg_fitness = np.average(position_fitness_list)
        print('{0}\t{1}'.format(position, avg_fitness))
