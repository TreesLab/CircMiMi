import re
import numpy as np
from scipy.stats import hypergeom


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
    circ_target_df_with_pv = circ_target_df.assign(
        p_value=lambda sdf: sdf.apply(lambda s: _calc_hypergeom_pvalue(s), axis=1)
    )

    # adjusting p-values


    return circ_target_df_with_pv