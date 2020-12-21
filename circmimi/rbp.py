from itertools import cycle


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
            pos = start + rel_pos
        elif strand == '-':
            pos = end - rel_pos

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

        while rel_start >= l:
            rel_start -= l
            rel_end -= l
            l, r = next(regions)

        start_pos = self.get_real_pos(rel_start, r)

        if rel_end < l:
            end_pos = self.get_real_pos(rel_end, r)
            real_blocks.append(self._to_region(start_pos, end_pos))
        else:
            while rel_end >= l:
                end_pos = self.get_real_pos(l - 1, r)
                real_blocks.append(self._to_region(start_pos, end_pos))

                rel_end -= l
                l, r = next(regions)
                start_pos = self.get_real_pos(0, r)

            end_pos = self.get_real_pos(rel_end, r)
            real_blocks.append(self._to_region(start_pos, end_pos))

        return tuple(real_blocks)
