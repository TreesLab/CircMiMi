import pandas as pd
import io


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

    def __init__(self, regions, name='.'):
        self.regions = regions
        self.name = name

        self._parse_regions(regions)

    def _parse_regions(self, regions):
        self._parse_region(regions[0])
        for region in regions[1:]:
            self._add_block(region)

        if self.strand == '-':
            self._reverse_region()

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

    def _reverse_region(self):
        self.block_sizes = list(reversed(self.block_sizes))
        self.block_starts = list(reversed(self.block_starts))

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
    def to_bed_df(cls, regions_df):
        if regions_df.empty:
            bed_df = pd.DataFrame([], columns=Bed.BED_TITLE)
        else:
            bed_df = regions_df.apply(cls._to_bed, axis=1)
            bed_df.columns = Bed.BED_TITLE

        return bed_df

    @classmethod
    def _to_bed(cls, s):
        bed = Bed(s.regions, name=s.regions_id)
        return pd.Series(bed.get_data(all_fields=True))
