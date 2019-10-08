'''
Use case:

from .resource import (EnsemblGenome,
                               EnsemblAnnotation,
                               GencodeGenome,
                               GencodeAnnotation)

reference_genome = EnsemblGenome('homo_sapiens', '97')
reference_genome.download('./refs')
genome_filename = reference_genome.filename

'''


import subprocess as sp
import ftplib
import os
import re


class ResourceBase:
    def __init__(self, url):
        self.url = url
        self.filename = os.path.basename(url)

    def download(self, dir_='.'):
        res = sp.run(['wget', self.url, '-P', dir_])
        return res


class FTPResource(ResourceBase):
    ftp_site = ""
    release_dir = ""
    release_pattern = r""
    file_pattern = r""

    def __init__(self, species, version):
        self.species = species
        self.version = version
        url = self.get_url()
        super().__init__(url)

    def get_url(self):
        with ftplib.FTP(self.ftp_site) as ftp:
            ftp.login()

            if self.version == "current":
                self.version = self.get_current_release(ftp)

            dir_path = self.get_dir_path()

            ftp.cwd(dir_path)
            dir_list = ftp.nlst()

            filename = self.get_the_filename(dir_list)

            url = "ftp://{}{}".format(
                self.ftp_site,
                os.path.join(dir_path, filename)
            )

            return url

    def get_current_release(self, ftp, max_key=float):
        ftp.cwd(self.release_dir.format(species=self.species))
        dir_list = ftp.nlst()
        current_release = max(map(self.get_release_number, dir_list), key=max_key)
        return current_release

    @classmethod
    def get_release_number(cls, name):
        m = re.search(cls.release_pattern, name)
        if m:
            return m.group(1)
        else:
            return '-1'

    def get_dir_path(self):
        raise NotImplementedError

    @classmethod
    def get_the_filename(cls, dir_list):
        for name in dir_list:
            m = re.search(cls.file_pattern, name)
            if m:
                return m.group(0)
        else:
            return ""


class EnsemblResource(FTPResource):
    ftp_site = "ftp.ensembl.org"
    release_dir = "/pub"
    release_pattern = r"^release-([0-9]+)$"


class EnsemblGenome(EnsemblResource):
    file_pattern = r".+\.dna\.primary_assembly\.fa\.gz$"

    def get_dir_path(self):
        dir_path = os.path.join(
            '/',
            'pub',
            'release-{}'.format(self.version),
            'fasta',
            self.species,
            'dna'
        )

        return dir_path


class EnsemblAnnotation(EnsemblResource):
    file_pattern = r".+\.chr\.gtf\.gz$"

    def get_dir_path(self):
        dir_path = os.path.join(
            '/',
            'pub',
            'release-{}'.format(self.version),
            'gtf',
            self.species
        )

        return dir_path


class GencodeResource(FTPResource):
    ftp_site = "ftp.ebi.ac.uk"
    release_dir = "/pub/databases/gencode/Gencode_{species}/"
    release_pattern = r"^release_(M?[0-9]+)$"

    def get_current_release(self, ftp, max_key=float):
        if self.species == "mouse":
            max_key = self.get_digit_part

        return super().get_current_release(ftp, max_key)

    @staticmethod
    def get_digit_part(version):
        m = re.search(r'[^0-9]*([0-9]+)', version)
        if m:
            digit = float(m.group(1))
        else:
            digit = -1

        return digit

    def get_dir_path(self):
        dir_path = os.path.join(
            '/',
            'pub/databases/gencode',
            'Gencode_{}'.format(self.species),
            'release_{}'.format(self.version)
        )

        return dir_path


class GencodeGenome(GencodeResource):
    file_pattern = r".+\.primary_assembly\.genome\.fa\.gz$"


class GencodeAnnotation(GencodeResource):
    file_pattern = r"^gencode.v[^.]+.annotation.gtf.gz$"


class MiRBaseResource(FTPResource):
    ftp_site = "mirbase.org"
    release_dir = "/pub/mirbase/"
    release_pattern = r"([0-9]+(?:\.[0-9]+)?)"

    def get_dir_path(self):
        dir_path = os.path.join(
            '/',
            'pub/mirbase',
            '{}'.format(self.version)
        )

        return dir_path


class MiRBaseMiRNA(MiRBaseResource):
    file_pattern = r"mature\.fa\.gz"


class MiRTarBaseResource(ResourceBase):
    url_templ = "http://mirtarbase.mbc.nctu.edu.tw/cache/download/{version}/miRTarBase_MTI.xlsx"

    def __init__(self, species, version):
        self.species = species
        self.version = version
        url = self.get_url()
        super().__init__(url)

    def get_url(self):
        url = self.url_templ.format(species=self.species, version=self.version)
        return url
