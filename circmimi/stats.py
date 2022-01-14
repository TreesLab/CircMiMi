import re
import numpy as np
import concurrent.futures as cf
from scipy.stats import hypergeom
from statsmodels.stats import multitest
from operator import itemgetter


def _count_miRNAs(mir_ref_file):
    with open(mir_ref_file) as f_in:
        mir_ref = f_in.read()

    num_miRNAs = len(re.findall(r'^>', mir_ref, flags=re.M))
    return num_miRNAs


def _calc_hypergeom_pvalue(s):
    k = s['#miRNAs_share']
    M = s['#miRNAs']
    n = s['#miRNAs_target_gene']
    N = s['#miRNAs_circRNAs']

    rv = hypergeom(M, n, N)
    x = np.arange(k, min(n, N))
    pvalue = rv.pmf(x).sum()

    return pvalue


def _calc_hypergeom_pvalue_in_parallel(iterrows, num_proc=1):
    with cf.ProcessPoolExecutor(num_proc) as executor:
        p_values = executor.map(_calc_hypergeom_pvalue, iterrows)
    return p_values


def do_the_calculation_for_hypergeom_pvalue(circ_target_df, num_proc=1):
    if num_proc == 1:
        circ_target_df_with_pv = circ_target_df.assign(
            p_value=lambda sdf: sdf.apply(_calc_hypergeom_pvalue, axis=1)
        )
    else:
        rows = list(map(itemgetter(1), circ_target_df.iterrows()))
        p_values = list(_calc_hypergeom_pvalue_in_parallel(rows, num_proc=num_proc))
        circ_target_df_with_pv = circ_target_df.assign(p_value=p_values)

    return circ_target_df_with_pv


def do_the_hypergeometric_test(circ_mi_target_df, mir_ref_file, mir_target_db, num_proc=1):
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
    # circ_target_df_with_pv = circ_target_df.assign(
    #     p_value=lambda sdf: sdf.apply(lambda s: _calc_hypergeom_pvalue(s), axis=1)
    # )
    circ_target_df_with_pv = do_the_calculation_for_hypergeom_pvalue(circ_target_df, num_proc=num_proc)

    # adjusting p-values
    bonferroni_corrected_p_values = circ_target_df_with_pv.groupby('circ_id')['p_value'].transform(
        lambda pv: multitest.multipletests(pv, method='bonferroni')[1]
    )
    bh_corrected_p_values = circ_target_df_with_pv.groupby('circ_id')['p_value'].transform(
        lambda pv: multitest.multipletests(pv, method='fdr_bh')[1]
    )

    circ_target_df_with_pv = circ_target_df_with_pv.assign(
        bonferroni_corrected_p_values=bonferroni_corrected_p_values,
        bh_corrected_p_value=bh_corrected_p_values
    )

    return circ_target_df_with_pv
