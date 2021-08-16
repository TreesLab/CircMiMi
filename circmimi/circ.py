import pandas as pd
import numpy as np
from collections import Counter, defaultdict
from circmimi.annotation import Annotator
from circmimi.ambiguous import AmbiguousChecker


class CircEvents:
    INPUT_COLUMNS = ('chr', 'pos1', 'pos2', 'strand', 'circ_id')

    def __init__(self, filename):
        self._filename = filename
        self.original_df = self._read_file(self._filename)
        self.df = self._get_donor_acceptor_df(self.original_df)

        self._summary_columns = {
            'summary': [],
            'filters': []
        }

    @classmethod
    def _read_file(cls, filename):
        df = pd.read_csv(
            filename,
            sep='\t',
            names=cls.INPUT_COLUMNS,
            dtype={
                'chr': 'category',
                'pos1': 'int',
                'pos2': 'int',
                'strand': 'category'
            }
        ).rename_axis('ev_id')

        return df

    @classmethod
    def _get_donor_acceptor(cls, s):
        chr_, pos1, pos2, strand, circ_id = s[list(cls.INPUT_COLUMNS)]
        pos1, pos2 = sorted([pos1, pos2])

        if strand == '+':
            donor_site = pos2
            acceptor_site = pos1
        elif strand == '-':
            donor_site = pos1
            acceptor_site = pos2

        res = [chr_, donor_site, acceptor_site, strand, circ_id]
        res = pd.Series(res, index=cls.INPUT_COLUMNS)

        return res

    @classmethod
    def _get_donor_acceptor_df(cls, original_df):
        if original_df.empty:
            df = pd.DataFrame([], columns=cls.INPUT_COLUMNS)
        else:
            df = original_df.apply(cls._get_donor_acceptor, axis=1)

        return df

    def expand_to_all_events(self, ev_df, fillna_value):
        expanded_df = self.original_df.reset_index().rename(
            {
                'index': 'ev_id'
            },
            axis=1
        ).loc[:, ['ev_id']].merge(
            ev_df,
            on='ev_id',
            how='left'
        ).fillna(fillna_value)

        return expanded_df

    def submit_to_summary(self, summary_column, type_):
        self._summary_columns[type_].append(summary_column)

    @staticmethod
    def _get_host_genes(anno_df):

        def get_gene_symbol(anno_df):
            return anno_df['transcript'].apply(
                lambda t: t.gene.gene_symbol
            )

        host_gene_df = anno_df.assign(
            host_gene=get_gene_symbol
        )[[
            'ev_id',
            'host_gene'
        ]].drop_duplicates(
        ).groupby(
            'ev_id'
        ).agg(lambda genes: ','.join(sorted(genes))).reset_index()

        return host_gene_df

    @staticmethod
    def _get_circ_ids(host_genes):
        genes = host_genes['host_gene']
        gene_count = Counter(genes)
        idx_dict = defaultdict(int)

        circ_ids = []
        for gene in genes:
            if gene_count[gene] == 1:
                circ_ids.append(f'circ{gene}')
            else:
                idx_dict[gene] += 1
                circ_ids.append(f'circ{gene}_{idx_dict[gene]}')

        circ_ids_df = host_genes.assign(
            circ_id=pd.Series(circ_ids)
        )[['ev_id', 'circ_id']]

        return circ_ids_df

    def _update_circ_ids(self, circ_ids):
        all_ev_with_circ_id = self.expand_to_all_events(circ_ids, np.nan)
        self.original_df = self.original_df.assign(circ_id=all_ev_with_circ_id)
        self.df = self.df.assign(circ_id=all_ev_with_circ_id)

    def check_annotation(self, anno_db_file):
        self._annotator = Annotator(anno_db_file)
        self.anno_df, anno_status = self.df.pipe(self._annotator.annotate)

        self.host_genes = self._get_host_genes(self.anno_df)
        self.circ_ids = self._get_circ_ids(self.host_genes)
        self._update_circ_ids(self.circ_ids)

        self.submit_to_summary(anno_status, type_='filters')

    def check_ambiguous(self,
                        anno_db_file,
                        ref_file,
                        other_ref_file,
                        work_dir='.',
                        num_proc=1):

        self.checker = AmbiguousChecker(
            anno_db_file,
            ref_file,
            other_ref_file,
            work_dir=work_dir,
            num_proc=num_proc
        )
        checking_result = self.checker.check(self.df)

        self.submit_to_summary(checking_result, type_='filters')

    @staticmethod
    def _merge_columns(df, column_dfs):
        merged_df = df

        for col_df in column_dfs:
            merged_df = merged_df.merge(col_df, on='ev_id', how='left')

        return merged_df

    def get_summary(self):
        filters_df = pd.DataFrame(self.df.index).pipe(
            self._merge_columns,
            self._summary_columns['filters']
        )

        pass_column = filters_df.fillna('2').set_index(
            'ev_id'
        ).astype(
            int
        ).agg(
            'sum',
            axis=1
        ).eq(
            0
        ).replace(
            {
                True: 'yes',
                False: 'no'
            }
        ).reset_index().rename({0: 'pass'}, axis=1)

        columns_to_be_merged = [self.host_gene, self.circ_ids]
        columns_to_be_merged += self._summary_columns['summary']
        columns_to_be_merged += [pass_column, filters_df.fillna('NA')]

        summary_df = self.original_df.reset_index().pipe(
            self._merge_columns,
            columns_to_be_merged
        ).set_index('ev_id')

        return summary_df

    @property
    def clear_df(self):
        pass_df = self.get_summary()[lambda df: df['pass'] == 'yes'].reset_index()[['ev_id']]
        clear_df = self.original_df.merge(pass_df, on='ev_id')
        return clear_df

    @property
    def clear_anno_df(self):
        return self.anno_df.merge(
            self.clear_df.reset_index(),
            on='ev_id',
            how='inner'
        )

    @property
    def clear_df_with_gene(self):
        return self.clear_df.merge(self.host_gene, on='ev_id', how='left').set_index('ev_id')
