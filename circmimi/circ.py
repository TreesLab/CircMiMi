import pandas as pd
from circmimi.annotation import Annotator
from circmimi.ambiguous import AmbiguousChecker


class CircEvents:
    def __init__(self, filename):
        self._filename = filename
        self.original_df = self._read_file(self._filename)
        self.df = self._get_donor_acceptor_df(self.original_df)

        self._summary_columns = {
            'summary': [],
            'filters': []
        }

    @staticmethod
    def _read_file(filename):
        df = pd.read_csv(
            filename,
            sep='\t',
            names=['chr', 'pos1', 'pos2', 'strand'],
            dtype={
                'chr': 'category',
                'pos1': 'int',
                'pos2': 'int',
                'strand': 'category'
            }
        ).rename_axis('ev_id')

        return df

    @staticmethod
    def _get_donor_acceptor(s):
        chr_, pos1, pos2, strand = s[['chr', 'pos1', 'pos2', 'strand']]
        pos1, pos2 = sorted([pos1, pos2])

        if strand == '+':
            donor_site = pos2
            acceptor_site = pos1
        elif strand == '-':
            donor_site = pos1
            acceptor_site = pos2

        res = [chr_, donor_site, acceptor_site, strand]
        res = pd.Series(res, index=['chr', 'donor', 'acceptor', 'strand'])

        return res

    @classmethod
    def _get_donor_acceptor_df(cls, original_df):
        if original_df.empty:
            df = pd.DataFrame([], columns=['chr', 'donor', 'acceptor', 'strand'])
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

    def check_annotation(self, anno_db_file):
        self._annotator = Annotator(anno_db_file)
        self.anno_df, anno_status = self.df.pipe(self._annotator.annotate)

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

        summary_df = self.original_df.reset_index().pipe(
            self._merge_columns,
            [pass_column] + self._summary_columns['summary'] + [filters_df.fillna('NA')]
        ).set_index('ev_id')

        return summary_df

    @property
    def clear_df(self):
        return self.get_summary()[lambda df: df['pass'] == 'yes'].loc[:, :'strand']

    @property
    def clear_anno_df(self):
        return self.anno_df.merge(
            self.clear_df.reset_index(),
            on='ev_id',
            how='inner'
        )
