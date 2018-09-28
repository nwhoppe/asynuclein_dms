#!/usr/bin/env python3

import argparse
import collections
import json
import itertools


def parse_bowtie_output(bowtie_output_file, max_mismatch):
    """Extracts read headers and fasta identifiers from bowtie output file
    makes a lot of assumptions about file format based on bowtie version 1.2.2 and output format ____

    Args:
        bowtie_output_file: a string for the path to output from Bowtie software
        max_mismatch: maximum number of mismatched positions in alignment.
            If mismatches exceed this number, the alignment is discarded

    Returns:
        A dict with keys of illumina sequencing headers and values of variant mutations from fasta header

    Raises:
        AssertionError: Read headers do not appear on sequential lines
        AssertionError: Paired reads do not align to the same variant sequence
    """
    header_variant_dict = {}
    with open(bowtie_output_file, 'r') as f:
        for line1, line2 in itertools.zip_longest(f, f, fillvalue=None):
            split_line1 = line1.rstrip().split()
            split_line2 = line2.rstrip().split()
            header1 = split_line1[0]
            header2 = split_line2[0]
            assert header1 == header2, 'Read headers do not appear on sequential lines'
            twist_tag1 = split_line1[2]
            twist_tag2 = split_line1[2]
            assert twist_tag1 == twist_tag2, 'Paired reads do not align to the same variant sequence'

            # botwie counts mismatches descending from 0
            try:
                mismatches = -1 * (int(split_line1[-10].split(':')[-1]) + int(split_line2[-10].split(':')[-1]))
            except ValueError:
                print("Read {0} missing alignment score".format(header1))
                continue

            if mismatches <= max_mismatch:
                variant = (int(twist_tag1.split('REGION01_GROUP')[-1].split(':')[0]), twist_tag1.split('-->')[-1])
                header_variant_dict[header1] = variant
    return header_variant_dict


def parse_index_fastq(index_fastq):
    """TODO: decorator
    """
    identifier_sequence_dict = {}
    with open(index_fastq, 'r') as f:
        for identifier, sequence, spacer, quality_str in itertools.zip_longest(*[f] * 4, fillvalue=None):
            identifier_sequence_dict[identifier.split()[0].split('@')[-1]] = sequence.rstrip()
    return identifier_sequence_dict


def bowtie_barcode_library_dict(bowtie_output, index_fastq, max_mismatch):
    header_barcode_dict = parse_index_fastq(index_fastq)
    header_variant_dict = parse_bowtie_output(bowtie_output, max_mismatch)

    barcode_variant_counter = collections.defaultdict(collections.Counter)
    for header, variant in header_variant_dict.items():
        barcode = header_barcode_dict[header]
        barcode_variant_counter[barcode][variant] += 1
    with open('barcode_variant_counter.json', 'w') as f:
        json.dump(barcode_variant_counter, f, sort_keys=True, indent=4)

    # TODO: check barcodes only have one variant - based on count - different script


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script with functions to take bowtie output and generate a 
    dictionary of barcodes and library variants""")
    required = parser.add_argument_group('required')
    required.add_argument('-b', '--bowtie_output', required=True,
                          help='bowtie output file, which must include the illumina read header, the twist library'
                               'fasta identifier. Best alignment only')
    required.add_argument('-i', '--index_fastq', required=True,
                          help='indexing fastq file from illumina sequencer')
    parser.add_argument('-m', '--max_mismatch', type=int, default=0,
                        help='max number of mismatches allowed between read and expected fasta sequence')
    args = parser.parse_args()
    bowtie_barcode_library_dict(args.bowtie_output, args.index_fastq, args.max_mismatch)
