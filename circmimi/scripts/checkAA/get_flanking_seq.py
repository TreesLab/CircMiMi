#! /usr/bin/env python

"""
A python script to get the sequences consist of the flanking +-100bp
 sequences on circRNA from the circRNA junction.

Note:
  The `bedtools` is required.
  We use `bedtools getfasta` to get the needed sequences from the genome.

Usage:
  get_flanking_seq.py [Genome] [CircRNAs] [Output]

  Genome:
    The reference genome file.

  CircRNAs:
    CircRNAs data in the 4-columns format with (chr, pos1, pos2, strand).

  Output:
    The output file in FASTA format.
"""

import sys
import io
import re
import tempfile as tp
import subprocess as sp
import logging


logging.basicConfig(
    format="{asctime} - {message}",
    level=logging.INFO,
    style='{'
)


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


def get_fasta(bed_file, ref_file, use_blocks=False, bedtools_bin='bedtools'):
    cmd = [bedtools_bin, 'getfasta']
    cmd += ['-fi', ref_file]
    cmd += ['-bed', bed_file]
    cmd += ['-fo', '/dev/stdout']
    cmd += ['-name', '-s']

    if use_blocks:
        cmd += ['-split']

    res = sp.run(cmd, stdout=sp.PIPE, encoding='utf-8')

    return res.stdout


class CircRNAs:
    def __init__(self, chr_, pos1, pos2, strand, circ_id=None):
        self.chr_ = chr_
        self.pos1 = int(pos1)
        self.pos2 = int(pos2)
        self.strand = strand
        self.circ_id = circ_id

        self.get_donor_acceptor()

    def get_donor_acceptor(self):
        pos1, pos2 = sorted([self.pos1, self.pos2])

        if self.strand == '+':
            self.donor_site = pos2
            self.acceptor_site = pos1
        elif self.strand == '-':
            self.donor_site = pos1
            self.acceptor_site = pos2

    def get_flanking_region(self):
        donor_flanking = self._get_flanking_region_of_donor()
        acceptor_flanking = self._get_flanking_region_of_acceptor()
        return [donor_flanking, acceptor_flanking]

    def _get_flanking_region_of_donor(self):
        if self.strand == '+':
            donor_start = self.donor_site - 99
            donor_end = self.donor_site

        elif self.strand == "-":
            donor_start = self.donor_site
            donor_end = self.donor_site + 99

        donor_flanking = [
            self.chr_,
            donor_start,
            donor_end,
            self.strand
        ]

        return donor_flanking

    def _get_flanking_region_of_acceptor(self):
        if self.strand == '+':
            acceptor_start = self.acceptor_site
            acceptor_end = self.acceptor_site + 99

        elif self.strand == "-":
            acceptor_start = self.acceptor_site - 99
            acceptor_end = self.acceptor_site

        acceptor_flanking = [
            self.chr_,
            acceptor_start,
            acceptor_end,
            self.strand
        ]

        return acceptor_flanking

    @property
    def id(self):
        if self.circ_id:
            return self.circ_id
        else:
            return f"{self.chr_}:{self.donor_site}|{self.acceptor_site}({self.strand})"


def get_flanking_seq(genome_file, circRNAs_file, out_file):
    logging.info("Start to get the flanking sequences.")

    tmp_bed = tp.NamedTemporaryFile(
        prefix='flanking_region.',
        suffix='.bed',
        dir='.'
    )

    with open(circRNAs_file) as circRNA_data, \
            open(tmp_bed.name, 'w') as bed_out:
        for line in circRNA_data:
            circRNA = CircRNAs(*line.rstrip('\n').split('\t'))
            flanking_region = circRNA.get_flanking_region()
            region_bed = Bed(flanking_region, name=circRNA.id)
            print(region_bed.to_string(all_fields=True), file=bed_out)

    fasta_data = get_fasta(tmp_bed.name, genome_file, use_blocks=True)

    with open(out_file, 'w') as out:
        out.write(
            re.sub(r'(?<=\([+-]\))\([+-]\)$', '', fasta_data, flags=re.M)
        )

    tmp_bed.close()

    logging.info("Process completed!")


def print_usage():
    print(__doc__, file=sys.stderr)


def cli():
    if len(sys.argv) < 4:
        print_usage()
        exit(1)

    get_flanking_seq(*sys.argv[1:])


if __name__ == "__main__":
    cli()
