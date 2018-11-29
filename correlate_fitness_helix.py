#!/usr/bin/env python3

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt


def correlate_helical_propensity(helical_propensity_df, x, y):
    """
    generate X Y plot of query df column vs helical propensity. data is averaged by amino acid

    Args:
        helical_propensity_df: pandas dataframe with helical penalty for each amino acid
        query_df: pandas dataframe containing query_column
        query_column: name of column in query_df to be plotted

    Returns:


    Raises:

    """
    sns.regplot(x=x, y=y, data=helical_propensity_df)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to generate figures analyzing alpha synuclein fitness 
    values and alpha helix propensity""")
    required = parser.add_argument_group('required')
    required.add_argument('-f', '--fitness_csv', required=True,
                          help="csv file containing at least 3 columns: fitness, variant_aa, variant_num")
    parser.add_argument('-p', '--helical_propensity_csv', default='helical_propensity.csv',
                        help='csv file containing at least 2 columns: One letter, Helical Penalty (kJ/mol)')
    args = parser.parse_args()

    fitness_df = pd.read_csv(args.fitness_csv, header=0)
    fitness_df.dropna(inplace=True)
    fitness_df = fitness_df[fitness_df['variant_aa'] != '*']

    helical_propensity_df = pd.read_csv(args.helical_propensity_csv, header=0)
    helical_propensity_df.set_index('One letter', inplace=True)

    # N terminal helix last position
    n_term_helix_end = 32

    # KxKEGV repeat last position
    repeat_end = 63

    # df_allaa_allpos = fitness_df
    # df_allaa_allpos['Group'] = 'All pos All AAs'
    #
    # df_gp_allpos = fitness_df[
    #         (fitness_df['variant_aa'] == 'P') |
    #         (fitness_df['variant_aa'] == 'G')
    #     ]
    # df_gp_allpos['Group'] = 'All pos G P'
    #
    # df_de_allpos = fitness_df[
    #     (fitness_df['variant_aa'] == 'D') |
    #     (fitness_df['variant_aa'] == 'E')
    # ]
    # df_de_allpos['Group'] = 'All pos D E'

    df_allaa_1_63 = fitness_df[fitness_df['variant_num'] < repeat_end + 1]
    df_other_1_63 = df_allaa_1_63[
        (df_allaa_1_63['variant_aa'] != 'P') |
        (df_allaa_1_63['variant_aa'] != 'G') |
        (df_allaa_1_63['variant_aa'] != 'D') |
        (df_allaa_1_63['variant_aa'] != 'E')
    ]
    df_other_1_63['Group'] = 'Repeat Other'
    df_other_1_63['AAs'] = 'Other'

    df_gp_1_63 = df_allaa_1_63[
        (df_allaa_1_63['variant_aa'] == 'P') |
        (df_allaa_1_63['variant_aa'] == 'G')
    ]
    df_gp_1_63['Group'] = 'Repeat G P'
    df_gp_1_63['AAs'] = 'G P'

    df_de_1_63 = df_allaa_1_63[
        (df_allaa_1_63['variant_aa'] == 'D') |
        (df_allaa_1_63['variant_aa'] == 'E')
        ]
    df_de_1_63['Group'] = 'Repeat D E'
    df_de_1_63['AAs'] = 'D E'

    df_allaa_140 = fitness_df[fitness_df['variant_num'] > repeat_end]
    df_other_140 = df_allaa_140[
        (df_allaa_140['variant_aa'] != 'P') |
        (df_allaa_140['variant_aa'] != 'G') |
        (df_allaa_140['variant_aa'] != 'D') |
        (df_allaa_140['variant_aa'] != 'E')
    ]
    df_other_140['Group'] = 'Past Other'
    df_other_140['AAs'] = 'Other'

    df_gp_140 = df_allaa_140[
        (df_allaa_140['variant_aa'] == 'P') |
        (df_allaa_140['variant_aa'] == 'G')
        ]
    df_gp_140['Group'] = 'Past G P'
    df_gp_140['AAs'] = 'G P'

    df_de_140 = df_allaa_140[
        (df_allaa_140['variant_aa'] == 'D') |
        (df_allaa_140['variant_aa'] == 'E')
        ]
    df_de_140['Group'] = 'Past D E'
    df_de_140['AAs'] = 'D E'

    composite_df = pd.concat(
        [df_other_1_63, df_gp_1_63, df_de_1_63, df_other_140, df_gp_140, df_de_140]
    )

    # composite_df = pd.concat(
    #     [df_allaa_allpos, df_gp_allpos, df_de_allpos, df_allaa_1_63, df_gp_1_63, df_de_1_63])

    # violin_ax = sns.violinplot(x='Group', y='fitness', data=composite_df, scale='width')
    # violin_fig = violin_ax.get_figure()
    # violin_fig.savefig('helix_groups_violin.png', dpi=500)

    box_ax = sns.boxplot(x='Group',
                         y='fitness',
                         data=composite_df,
                         showfliers=False,
                         hue='AAs'
                         )
    box_fig = box_ax.get_figure()
    box_fig.savefig('helix_groups_untreated.png', dpi=500)




    # plotting_df = pd.concat([repeat_avg_df, helical_propensity_df], axis='columns')
    # correlate_helical_propensity(plotting_df, x='Helical Penalty (kcal/mol)', y='fitness')
    # plt.savefig('test.png')



