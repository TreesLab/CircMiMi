import pandas as pd
from itertools import cycle
from bisect import bisect
from circmimi.bed import IntersectBED, Bed


class RBPBindingSites:
    _titles = Bed.BED_TITLE \
        + tuple(f'{title}_rbp' for title in Bed.BED_TITLE[:6]) \
        + ('overlap',)

    def __init__(self, rbp_file, bedtools_bin="bedtools"):
        self.rbp_file = rbp_file
        self.bedtools_bin = bedtools_bin
        self.intersect_bed = IntersectBED(
            bedtools_bin=self.bedtools_bin,
            options='-wo'
        )

    @staticmethod
    def _recover_blocks(bed_line):
        bed_start = bed_line['start'] + 1
        blockStarts = list(map(int, bed_line['blockStarts'].split(',')))
        blockSizes = list(map(int, bed_line['blockSizes'].split(',')))

        blocks = [Bed._region(bed_start + start, bed_start + start + size - 1)
                  for start, size in zip(blockStarts, blockSizes)]

        return blocks

    @staticmethod
    def _shrink_site_to_region(r, s, offset=0):
        bisect_result = bisect(r, s)

        if bisect_result == 0:
            s_in_r = r.start + offset - 1
        elif bisect_result == 1:
            s_in_r = s
        elif bisect_result == 2:
            s_in_r = r.end + offset

        return s_in_r

    @classmethod
    def _calculate_two_region_overlap(cls, r1, r2):
        start_in_r1 = cls._shrink_site_to_region(r1, r2.start, offset=1)
        end_in_r1 = cls._shrink_site_to_region(r1, r2.end)
        overlap = end_in_r1 - start_in_r1 + 1
        return overlap

    @classmethod
    def _calculate_overlap(cls, bed_line):
        blocks = cls._recover_blocks(bed_line)
        rbp_block = Bed._region(bed_line['start_rbp'] + 1, bed_line['end_rbp'])

        overlap = 0
        for b in blocks:
            overlap += cls._calculate_two_region_overlap(b, rbp_block)

        return overlap

    @classmethod
    def _get_real_overlap(cls, bed_line):
        blockCount = bed_line['blockCount']

        if blockCount == 1:
            overlap = bed_line['overlap']
        elif blockCount > 1:
            overlap = cls._calculate_overlap(bed_line)

        return overlap

    @staticmethod
    def _get_split_rbp_name(idx):
        return lambda df: df['name_rbp'].apply(lambda name: name.split('_')[idx])

    def overlap(self, bed_df):
        intersect_result_df = pd.DataFrame(
            self.intersect_bed.intersect(bed_df, self.rbp_file),
            columns=self._titles,
            dtype='object'
        ).astype(
            {
                'start': int,
                'end': int,
                'blockCount': int,
                'start_rbp': int,
                'end_rbp': int,
                'overlap': int
            }
        ).astype(
            'object'
        ).assign(
            sample_id=self._get_split_rbp_name(0),
            RBP=self._get_split_rbp_name(1)
        ).assign(
            real_overlap=lambda df: df.apply(self._get_real_overlap, axis=1)
        ).assign(
            rbp_region_len=lambda df: df['end_rbp'] - df['start_rbp'],
            total_blocks_len=lambda df: df['blockSizes'].apply(
                lambda sizes: sum(map(int, sizes.split(',')))
            )
        )

        return intersect_result_df

    @staticmethod
    def _get_joined_overlap(group_df):
        joined_overlap = sum(group_df['real_overlap'])
        group_elts = tuple(group_df.index)
        joined_overlap_df = group_df.assign(
            joined_overlap=joined_overlap,
            group_elements=','.join(map(str, group_elts))
        )

        return joined_overlap_df

    @classmethod
    def append_joined_overlap(cls, df):
        groupby_df = df.groupby(
            [
                'name',
                'sample_id'
            ]
        ).apply(
            cls._get_joined_overlap
        ).reset_index(
            [
                'name',
                'sample_id'
            ],
            drop=True
        ).sort_index()

        return groupby_df


class RBPBindingSitesFilters:
    def AGO_overlap_filter(df):
        coverage_df = df[[
            'name',
            'sample_id',
            'RBP',
            'joined_overlap',
            'total_blocks_len'
        ]].drop_duplicates(
        ).assign(
            coverage=lambda df: (df['joined_overlap'] / df['total_blocks_len']).apply(lambda n: round(n, 2))
        )[
            lambda df: df['coverage'] >= 0.8
        ]

        return coverage_df

    def RBP_overlap_filter(df):
        coverage_df = df[[
            'name',
            'chr_rbp',
            'start_rbp',
            'end_rbp',
            'strand_rbp',
            'sample_id',
            'RBP',
            'real_overlap',
            'rbp_region_len'
        ]].assign(
            coverage=lambda df: (df['real_overlap']/df['rbp_region_len']).apply(lambda n: round(n, 2))
        )[
            lambda df: df['coverage'] >= 0.8
        ]

        return coverage_df


class PosMap:
    def __init__(self, genomic_regions):
        self._regions = genomic_regions

        self._init_lens(genomic_regions)

    def _init_lens(self, regions):
        self._lens = list(map(lambda r: r[2] - r[1] + 1, regions))

    @staticmethod
    def get_real_pos(rel_pos, genomic_region):
        chr_, start, end, strand = genomic_region

        if strand == '+':
            pos = start + rel_pos - 1
        elif strand == '-':
            pos = end - rel_pos + 1

        return (chr_, pos, strand)

    def regions_iter(self):
        return cycle(zip(self._lens, self._regions))

    @staticmethod
    def _to_region(pos1, pos2):
        c1, p1, s1 = pos1
        c2, p2, s2 = pos2

        assert c1 == c2
        assert s1 == s2

        p1, p2 = sorted([p1, p2])

        return (c1, p1, p2, s1)

    def get_real_blocks(self, rel_start, rel_end):
        real_blocks = []
        regions = self.regions_iter()

        l, r = next(regions)

        while rel_start > l:
            rel_start -= l
            rel_end -= l
            l, r = next(regions)

        start_pos = self.get_real_pos(rel_start, r)

        if rel_end <= l:
            end_pos = self.get_real_pos(rel_end, r)
            real_blocks.append(self._to_region(start_pos, end_pos))
        else:
            while rel_end > l:
                end_pos = self.get_real_pos(l - 1, r)
                real_blocks.append(self._to_region(start_pos, end_pos))

                rel_end -= l
                l, r = next(regions)
                start_pos = self.get_real_pos(0, r)

            end_pos = self.get_real_pos(rel_end, r)
            real_blocks.append(self._to_region(start_pos, end_pos))

        return tuple(real_blocks)
