import re
import numpy as np
from scipy.stats import hypergeom
from statsmodels.stats import multitest
from operator import itemgetter


def _count_miRNAs(mir_ref_file):
    with open(mir_ref_file) as f_in:
        mir_ref = f_in.read()

    num_miRNAs = len(re.findall(r'^>', mir_ref, flags=re.M))
    return num_miRNAs


def _calc_hypergeom_pvalue(s, rvs_dict):
    k = s['#miRNAs_share']
    M = s['#miRNAs']
    n = s['#miRNAs_target_gene']
    N = s['#miRNAs_circRNAs']

    rv = rvs_dict[(M, n, N)]
    x = np.arange(k, min(n, N))
    pvalue = rv.pmf(x).sum()

    return pvalue


def do_the_calculation_for_hypergeom_pvalue(circ_target_df):
    all_M_n_N = circ_target_df[[
        '#miRNAs',
        '#miRNAs_target_gene',
        '#miRNAs_circRNAs'
    ]].drop_duplicates().values

    all_M_n_N_tuple = [tuple(row) for row in all_M_n_N]

    rvs = np.apply_along_axis(lambda row: hypergeom(*row), 1, all_M_n_N)
    rvs_dict = dict(zip(all_M_n_N_tuple, rvs))

    p_values = circ_target_df.apply(_calc_hypergeom_pvalue, rvs_dict=rvs_dict, axis=1)

    circ_target_df_with_pv = circ_target_df.assign(p_value=p_values)

    return circ_target_df_with_pv


def do_the_hypergeometric_test(circ_mi_target_df, mir_ref_file, mir_target_db):
    num_miRNAs = _count_miRNAs(mir_ref_file)
    num_miRNAs__circRNA = circ_mi_target_df[['circ_id', 'mirna']].drop_duplicates().groupby('circ_id').agg('count')
    num_miRNAs__circRNA_target_gene = circ_mi_target_df.groupby(['circ_id', 'target_gene']).agg('count')
    num_miRNAs__target_gene = mir_target_db[['target_gene', 'mirna']].drop_duplicates().groupby('target_gene').agg('count')

    circ_target_df = num_miRNAs__circRNA_target_gene.reset_index().rename(
        {
            'mirna': '#miRNAs_share'
        },
        axis=1
    ).merge(
        num_miRNAs__circRNA,
        on='circ_id',
        how='left'
    ).rename(
        {
            'mirna': '#miRNAs_circRNAs'
        },
        axis=1
    ).assign(
        **{
            '#miRNAs': num_miRNAs
        }
    ).merge(
        num_miRNAs__target_gene,
        on='target_gene',
        how='left'
    ).rename(
        {
            'mirna': '#miRNAs_target_gene'
        },
        axis=1
    )

    # do the test
    circ_target_df_with_pv = do_the_calculation_for_hypergeom_pvalue(circ_target_df)

    # adjusting p-values
    bonferroni_corrected_p_values = circ_target_df_with_pv.groupby('circ_id')['p_value'].transform(
        lambda pvs: multitest.multipletests(pvs, method='bonferroni')[1]
    )
    bh_corrected_p_values = circ_target_df_with_pv.groupby('circ_id')['p_value'].transform(
        lambda pvs: multitest.multipletests(pvs, method='fdr_bh')[1]
    )

    circ_target_df_with_pv = circ_target_df_with_pv.assign(
        bh_corrected_p_value=bh_corrected_p_values,
        bonferroni_corrected_p_values=bonferroni_corrected_p_values,
    )

    return circ_target_df_with_pv
