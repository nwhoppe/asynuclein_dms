#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
import sys
import matplotlib.pyplot as plt


if __name__ == '__main__':
    wt_seq = 'MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVV' \
             'TGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'

    fitness_wt_csv = sys.argv[1]
    fitness_mel_csv = sys.argv[2]

    fitness_wt = pd.read_csv(fitness_wt_csv, header=0)
    fitness_wt.set_index('full_variant', inplace=True)
    fitness_mel = pd.read_csv(fitness_mel_csv, header=0)
    fitness_mel.set_index('full_variant', inplace=True)
    fitness_df = pd.concat([fitness_wt['fitness'], fitness_mel['fitness']], axis='columns', keys=['Untreated', 'Melatonin'])
    fitness_df.dropna(inplace=True)
    ax = sns.scatterplot(x='Untreated', y='Melatonin', data=fitness_df)
    xy_line = [
        min([min(fitness_df['Untreated']), min(fitness_df['Melatonin'])]),
        max([max(fitness_wt['Untreated']), max(fitness_mel['Melatonin'])])
    ]
    sns.lineplot(x=xy_line, y=xy_line)
    figure = ax.get_figure()
    figure.savefig('xy_scatter_unadjusted.png', dpi=400)
    # sns.savefig('hist.png', dpi=400)
