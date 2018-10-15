#!/usr/bin/env python3

import argparse
import numpy as np
import pickle


def compare_duplicates(fitness_dict1, fitness_dict2, max_percent_difference, name_suffix):
    variants_missing_from_rep = set(fitness_dict2.keys()) - set(fitness_dict1.keys())
    percent_diffs = []
    output_dict = {}
    print('Variants with percent difference greater than max')
    print('variant\tfitness1\tfitness2\tpercent difference')
    for variant1, fitness1 in fitness_dict1.items():
        if variant1 in fitness_dict2.keys():
            fitness2 = fitness_dict2[variant1]
            percent_diff = 100 * (abs(fitness1 - fitness2) / ((fitness1 + fitness2) / 2))
            percent_diffs.append(percent_diff)
            if percent_diff: #< max_percent_difference:
                output_dict[variant1] = np.average([fitness1, fitness2])
            else:
                print('{0}\t{1}\t{2}\t{3}'.format(variant1,
                                                  round(fitness1, 2),
                                                  round(fitness2, 2),
                                                  round(percent_diff), 2)
                      )
        else:
            variants_missing_from_rep.add(variant1)
    print()
    print('Number of variants not present in both replicates: {0}'.format(len(variants_missing_from_rep)))
    print('Variants not present in both replicates:')
    print(variants_missing_from_rep)
    print()
    print('Average percent difference: {0} +/- {1}'.format(
        round(np.average(percent_diffs), 3),
        round(np.std(percent_diffs), 3)
    ))
    if name_suffix:
        output_file = 'avg_duplicates_{0}.pkl'.format(name_suffix)
    else:
        output_file = 'avg_duplicates.pkl'
    with open(output_file, 'wb') as f:
        pickle.dump(output_dict, f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to compare fitness dictionaries between two replicates. 
    currently only supports dictionaries with a single floating point number as the key 
    and only supports two replicates""")
    required = parser.add_argument_group('required')
    required.add_argument('-f', '--fitness_pickles', nargs='*', required=True,
                          help='fitness dictionary in pickle format with keys as variant descriptors '
                               'and values as floats')
    parser.add_argument('-p', '--max_percent_difference', default=50,
                        help='output fitness values will be the average of fitness values with percent '
                             'difference less than this cutoff')
    parser.add_argument('-n', '--name_suffix')
    args = parser.parse_args()

    if len(args.fitness_pickles) != 2:
        raise IOError('two pickle files containing fitness values are reqquired')
    with open(args.fitness_pickles[0], 'rb') as f:
        fitness_dict1 = pickle.load(f)
    with open(args.fitness_pickles[1], 'rb') as f:
        fitness_dict2 = pickle.load(f)
    compare_duplicates(fitness_dict1, fitness_dict2, args.max_percent_difference, args.name_suffix)

