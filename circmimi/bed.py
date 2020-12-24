import pandas as pd
import io
import subprocess as sp
import tempfile as tp
from collections import namedtuple
from operator import itemgetter


class Bed:
    BED_TITLE = (
        'chr',
        'start',
        'end',
        'name',
        'score',
        'strand',
        'thickStart',
        'thickEnd',
        'itemRGB',
        'blockCount',
        'blockSizes',
        'blockStarts'
    )
    _region = namedtuple('Region', ('start', 'end'))

    def __init__(self, regions, name='.', union=False):
        self.regions = regions
        self.name = name
        self.union = union

        if not self.union:
            self.regions = self.reverse_regions_if_minus_strand(self.regions)
        else:
            self.regions = self.get_union_regions(self.regions)

        self._parse_regions(self.regions)

    @staticmethod
    def reverse_regions_if_minus_strand(regions):
        strand = regions[0][3]
        if strand == '-':
            regions = list(reversed(regions))

        return regions

    def _parse_regions(self, regions):
        self._parse_region(regions[0])
        for region in regions[1:]:
            self._add_block(region)

    def _parse_region(self, region):
        self.chrom = region[0]
        self.start = int(region[1]) - 1
        self.end = int(region[2])
        self.strand = region[3]

        self.block_count = 1
        self.block_sizes = [self.end - self.start]
        self.block_starts = [0]

    def _add_block(self, region):
        other = Bed([region])

        assert self.strand == other.strand

        self.block_count += other.block_count
        self.block_sizes += other.block_sizes

        all_starts = [
            start + self.start for start in self.block_starts
        ] + [other.start]

        if other.start < self.start:
            self.start = other.start

        self.block_starts = [start - self.start for start in all_starts]

        if other.end > self.end:
            self.end = other.end

    def get_data(self, all_fields=False):
        fields = [
            self.chrom,
            self.start,
            self.end,
            self.name,
            '.',
            self.strand
        ]

        if all_fields:
            fields += [
                self.start,
                self.end,
                0,
                self.block_count,
                self._list_to_str(self.block_sizes, sep=','),
                self._list_to_str(self.block_starts, sep=',')
            ]

        return fields

    def to_string(self, all_fields=False):
        data = self.get_data(all_fields=all_fields)
        bed_txt = self._list_to_str(data)
        return bed_txt

    @staticmethod
    def _list_to_str(list_, sep='\t', end=''):
        with io.StringIO() as tmp:
            print(*list_, sep=sep, end=end, file=tmp)
            return tmp.getvalue()

    @classmethod
    def get_union_regions(cls, regions):
        assert len(set(map(itemgetter(0), regions))) == 1, \
            "Not all regions in the same chromosome!"
        assert len(set(map(itemgetter(3), regions))) == 1, \
            "Not all regions at the same strand!"

        chr_ = regions[0][0]
        strand = regions[0][3]

        regions = sorted(
            map(
                lambda r: cls._region(r[1], r[2]),
                regions
            ),
            key=lambda r: r.start
        )

        union_regions = []
        r1 = regions[0]
        for r2 in regions[1:]:
            if r1.end < r2.start:
                union_regions.append(r1)
                r1 = r2
            else:
                if r1.end < r2.end:
                    r1 = cls._region(r1.start, r2.end)
        else:
            union_regions.append(r1)

        union_regions = tuple((chr_, start, end, strand)
                              for start, end in union_regions)
        return union_regions


class BedUtils:
    @classmethod
    def to_regions_df(cls, exons_df):
        if exons_df.empty:
            regions_df = pd.DataFrame([], columns=['regions_id', 'regions'])
        else:
            regions_df = exons_df[['exons_id', 'exons']].assign(
                regions_id=lambda df: df.exons_id,
                regions=lambda df: df.exons.apply(cls._exons_to_regions)
            )[['regions_id', 'regions']]

        return regions_df

    @staticmethod
    def _exons_to_regions(exons):
        regions = list(map(lambda exon: exon.region, exons))
        return regions

    @classmethod
    def to_bed_df(cls, regions_df, union=False):
        if regions_df.empty:
            bed_df = pd.DataFrame([], columns=Bed.BED_TITLE)
        else:
            bed_df = regions_df.apply(cls._to_bed, union=union, axis=1)
            bed_df.columns = Bed.BED_TITLE

        return bed_df

    @classmethod
    def _to_bed(cls, s, union=False):
        bed = Bed(s.regions, name=s.regions_id, union=union)
        return pd.Series(bed.get_data(all_fields=True))


class IntersectBED:
    _dummy_tmpfile = namedtuple('DummyTmpFile', ['name'])

    def __init__(self, bedtools_bin='bedtools', options=''):
        self._bedtools_bin = bedtools_bin

        if isinstance(options, str):
            self._options = options.split()
        else:
            self._options = options

    @classmethod
    def _save_temp(cls, bed):
        if isinstance(bed, str):
            tmp_file = cls._dummy_tmpfile(bed)
        else:
            tmp_file = tp.NamedTemporaryFile(dir='.', prefix='intersectBED.')
            bed.to_csv(tmp_file.name, sep='\t', index=False)

        return tmp_file

    def _generate_cmd(self, bed_file_1, bed_file_2):
        cmd = [self._bedtools_bin, 'intersect'] + self._options
        cmd += ['-a', bed_file_1]
        cmd += ['-b', bed_file_2]

        return cmd

    def intersect(self, bed_1, bed_2):
        bed_tmp_file_1 = self._save_temp(bed_1)
        bed_tmp_file_2 = self._save_temp(bed_2)

        cmd = self._generate_cmd(bed_tmp_file_1.name, bed_tmp_file_2.name)

        with sp.Popen(cmd, stdout=sp.PIPE, encoding='utf-8') as p:
            for line in p.stdout:
                yield line.rstrip('\n').split('\t')
