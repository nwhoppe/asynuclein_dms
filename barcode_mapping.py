#!/usr/bin/env python3

import argparse
import collections
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import numpy as np
import pickle

# TODO: histogram of number of barcodes per sequence
# TODO: number of variants present
# TODO: avg number of barcodes present plus minus std and max min - have comparison to uniform dist


def library_coverage(library_barcode_counter, expected_library_set):
    library_barcode_set = set(library_barcode_counter.keys())
    missing_library_seqs = expected_library_set.difference(library_barcode_set)
    additional_seqs = library_barcode_set.difference(expected_library_set)
    print('Expected size of library: {0}'.format(len(expected_library_set)))
    # account for wild type
    print('Number of barcoded library seqs: {0}'.format(len(library_barcode_counter.keys()) - 1))
    print('Number of wild type barcodes: {0}'.format(library_barcode_counter[(0, 'WT')]))
    print('Number of missing library seqs: {0}'.format(len(missing_library_seqs)))
    print('Missing library seqs:')
    for index, amino_acid in sorted(missing_library_seqs):
        library_barcode_counter[(index, amino_acid)] = 0
        print('{0}\t{1}'.format(index, amino_acid))

    del library_barcode_counter[(0, 'WT')]

    return library_barcode_counter


def counter_histogram(library_barcode_counter, total_element_count, xlabel):
    median = np.median(list(library_barcode_counter.values()))
    stdev = np.std(list(library_barcode_counter.values()))
    expected_per_variant = total_element_count / 2800
    bins = np.linspace(0, max(library_barcode_counter.values()), 100)
    plt.hist(list(library_barcode_counter.values()), bins=bins, density=True)
    plt.axvline(expected_per_variant, color='k', linestyle='dashed')
    plt.xlabel(xlabel)
    plt.ylabel('Count')
    plt.text(50, 0.14, 'Median: {0}\nExpected: {1}'.format(
        round(median, 2), round(expected_per_variant, 2)))
    plt.savefig('hist_barcode_per_variant.png', dpi=300)
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script which returns some preliminary stats on mapping randomized 
    barcodes to expected library sequences""")
    required = parser.add_argument_group('required')
    required.add_argument("-p", "--pickle_file", required=True,
                          help="input pickle file containing dictionary of 26 nucleotide barcodes ")
    required.add_argument("-w", "--wt_fasta", required=True,
                          help="fasta sequence of wt protein")
    args = parser.parse_args()

    with open(args.wt_fasta, 'r') as f:
        wt_seq = f.readlines()[-1].rstrip()

    with open(args.pickle_file, 'rb') as f:
        # keys: 26 randomized nucleotides; values: tuple (index, amino acid)
        barcode_mutation_dict = pickle.load(f)

    library_barcode_counter = collections.Counter()
    incorrect_length_counter = collections.Counter()
    ambiguous_barcode_counter = collections.Counter()
    for barcode, index_aa_tup in barcode_mutation_dict.items():
        if len(barcode) == 26 and 'N' not in barcode:
            library_barcode_counter[index_aa_tup] += 1
        elif len(barcode) != 26:
            incorrect_length_counter[index_aa_tup] += 1
        elif 'N' in barcode:
            ambiguous_barcode_counter[index_aa_tup] += 1
        else:
            print('Barcode does not match any of the tests')
            print(barcode)
            print(index_aa_tup)
            raise Exception('Add case to handle this example')

    # TODO: put this in a different file that would be more generally useful
    aa_name_map = {'A': 'ALA', 'P': 'PRO', 'V': 'VAL', 'L': 'LEU', 'I': 'ILE', 'M': 'MET',
                   'F': 'PHE', 'Y': 'TYR', 'W': 'TRP', 'S': 'SER', 'T': 'THR', 'C': 'CYS',
                   'K': 'LYS', 'R': 'ARG', 'H': 'HIS', 'D': 'ASP', 'E': 'GLU', 'N': 'ASN',
                   'Q': 'GLN', 'G': 'GLY'}

    possible_mutations = list(aa_name_map.keys())
    possible_mutations.extend(['*'])
    expected_library_set = set(itertools.product(range(1, 141), possible_mutations))
    for index, amino_acid in enumerate(wt_seq):
        expected_library_set.remove((index + 1, amino_acid))
    assert len(expected_library_set) == 140 * 20, 'Expected library size is not correct'

    library_barcode_counter = library_coverage(library_barcode_counter, expected_library_set)

    counter_histogram(
        library_barcode_counter,
        len(barcode_mutation_dict.keys()),
        xlabel='Number of Barcodes per Variant',
    )