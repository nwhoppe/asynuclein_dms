#!/usr/bin/env python3

import argparse
import collections
import itertools
import numpy as np
import pickle
from scipy.stats import linregress


def variant_counter_from_fastqs(fastq_files, barcode_variant_dict):
    """get counter by looping over fastq and comparing to barcode dict keys"""

    variant_counter = collections.defaultdict(lambda: [0.5, 0.5, 0.5])

    for i, fastq_file in enumerate(fastq_files):
        with open(fastq_file, 'r') as f:
            for identifier, seq, spacer, quality in itertools.zip_longest(*[f] * 4, fillvalue=None):
                barcode = seq[:20]
                if barcode in barcode_variant_dict:
                    variant_counter[barcode_variant_dict[barcode]][i] += 1
    return variant_counter


def calculate_variant_fitness(variant_timepoint_counter):
    """uses method described in Doug Folwer's Enrich2 paper - fitness is the slope of linear regression line"""
    variant_fitness_dict = {}
    r2_values = []
    print('variant\tr squared\tvalues')
    for variant, count_list in variant_timepoint_counter.items():
        if variant == (0, 'WT'):
            continue
        ratio_array = np.array(count_list) / np.array(variant_timepoint_counter[(0, 'WT')])
        m_vt = np.log(ratio_array)
        slope, intercept, r_value, p_value, std_err = linregress([0, 1, 2], m_vt)

        variant_fitness_dict[variant] = slope
        r2_values.append(r_value**2)

        if r_value**2 < 0.8:
            print('{0}\t{1}\t{2}'.format(
                variant,
                round(r_value ** 2, 2),
                ' '.join(map(str, (round(y, 2) for y in m_vt)))))

    return variant_fitness_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""script to generate a pickle file containing counts of library 
    variants""")
    required = parser.add_argument_group('required')
    required.add_argument('-f', '--fastq_files', nargs='*', required=True,
                          help='illumina fastq file containing sequenced barcodes. when calculating fitness values '
                               'it is assumed that the order of fastq files corresponds to the order of time points')
    required.add_argument('-b', '--barcode_pickle', required=True,
                          help='pickle containing a dictionary with barcode as keys and library variants as values. '
                               'reads in fastq file will be searched for exact matches to dictionary keys')
    # parser.add_argument('-t', '--fitness', action='store_true',
    #                     help='calculate and output relative fitness for library variants instead of counts')
    parser.add_argument('-n', '--name_suffix')
    args = parser.parse_args()
    with open(args.barcode_pickle, 'rb') as f:
        barcode_variant_dict = pickle.load(f)

    variant_timepoint_counter = variant_counter_from_fastqs(args.fastq_files, barcode_variant_dict)
    # TODO: have script to compare replicates, test this script
    variant_fitness_dict = calculate_variant_fitness(variant_timepoint_counter)
    if args.name_suffix:
        output_file = 'variant_fitness_{0}.pkl'.format(args.name_suffix)
    else:
        output_file = 'variant_fitness.pkl'
    with open(output_file, 'wb') as f:
        pickle.dump(variant_fitness_dict, f)
