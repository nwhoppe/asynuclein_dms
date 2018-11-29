#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
import sys

"""
hydrophobic = {'A': 'ALA', 'P': 'PRO', 'V': 'VAL', 'L': 'LEU', 'I': 'ILE', 'M': 'MET',
               'F': 'PHE', 'Y': 'TYR', 'W': 'TRP', 'S': 'SER', 'T': 'THR', 'C': 'CYS',
               'K': 'LYS', 'R': 'ARG', 'H': 'HIS', 'D': 'ASP', 'E': 'GLU', 'N': 'ASN',
               'Q': 'GLN', 'G': 'GLY'}
               
wt_seq = 'MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVV' \
             'TGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'
             
"""


if __name__ == '__main__':
    wt_seq = 'MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVV' \
             'TGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'

    oxidation_prone = {'C': 'CYS', 'M': 'MET', 'F': 'PHE', 'Y': 'TYR', 'W': 'TRP'}

    hydrophobic_non_aromatic = {'A': 'ALA', 'V': 'VAL', 'L': 'LEU', 'I': 'ILE'}

    sulfur = {'C': 'CYS', 'M': 'MET'}

    fitness_wt_csv = sys.argv[1]
    fitness_mel_csv = sys.argv[2]

    oxidation_wt_fit = []
    control_wt_fit = []
    oxidation_mel_fit = []
    control_mel_fit = []

    # df = pd.DataFrame(columns=['Mutation to C, M, F, Y, W', 'Mutation to A, V, I, L'])

    for fitness_csv in [fitness_wt_csv, fitness_mel_csv]:
        fitness_df = pd.read_csv(fitness_csv, header=0)
        fitness_df.dropna(inplace=True)

        for index, row in fitness_df.iterrows():
            variant_aa = row['variant_aa']
            wt_aa = wt_seq[int(row['variant_num']) - 1]
            fitness = float(row['fitness'])

            if wt_aa in hydrophobic_non_aromatic and variant_aa in oxidation_prone:
                if fitness_csv == fitness_wt_csv:
                    oxidation_wt_fit.append(fitness)
                else:
                    oxidation_mel_fit.append(fitness)
            elif variant_aa in hydrophobic_non_aromatic and wt_aa in hydrophobic_non_aromatic:
                if fitness_csv == fitness_wt_csv:
                    control_wt_fit.append(fitness)
                else:
                    control_mel_fit.append(fitness)
    df = pd.DataFrame.from_dict({
        'Untreated NA': control_wt_fit,
        'Untreated Ox': oxidation_wt_fit,
        'Treated NA': control_mel_fit,
        'Treated Ox': oxidation_mel_fit},
        orient='index'
    ).T
    ax = sns.violinplot(
                        data=df,
                        palette="muted",
                        )
    figure = ax.get_figure()
    figure.savefig('test.png', dpi=400)
