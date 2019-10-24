import tempfile as tp
import subprocess as sp
import re
import pandas as pd


class Seq:
    @classmethod
    def get_seq(cls, bed_df, ref_file):
        with tp.NamedTemporaryFile(dir='.') as tmp_bed_file:

            with open(tmp_bed_file.name, 'w') as bed_out:
                bed_str = bed_df.to_csv(sep='\t', header=False, index=False)
                bed_out.write(bed_str)

            fasta_txt = get_fasta(tmp_bed_file.name, ref_file, use_blocks=True)

            fasta_df = pd.DataFrame(
                parse_fasta(fasta_txt),
                columns=['name', 'seq']
            )

            return fasta_df

    @staticmethod
    def to_fasta(fasta_df):
        fa_txt = ""
        for name, seq in fasta_df.values:
            fa_txt += ">{}\n{}\n".format(name, seq)
        return fa_txt

    @staticmethod
    def extend_seq_for_circ_js(s):
        seq = s['seq']
        extended_seq = seq + seq[:30]
        return extended_seq

    @classmethod
    def get_extended_seq(cls, bed_df, ref_file):
        fasta_df = cls.get_seq(bed_df, ref_file)
        if not fasta_df.empty:
            fasta_df['seq'] = fasta_df.apply(
                cls.extend_seq_for_circ_js,
                axis=1
            )

        return fasta_df


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


def parse_fasta(fasta_txt):
    for m in re.finditer(r">(.+?)(?:\([+-]\))?\n([^>]+)", fasta_txt):
        fa_id = m.group(1)
        fa_seq = ''.join(m.group(2).rstrip('\n').split('\n'))
        yield [fa_id, fa_seq]
