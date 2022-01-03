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
import urllib.request
import urllib.error
import os
import re


class Resource:
    wget_options = []

    def __init__(self, url):
        self.url = url
        self.filename = os.path.basename(url)

    def download(self, dir_='.'):
        if self.filename == '':
            return ''

        dest_path = os.path.join(dir_, self.filename)

        if os.path.exists(dest_path):
            self.filename = dest_path

        elif os.path.exists(dest_path.rstrip('.gz')):
            self.filename = dest_path.rstrip('.gz')

        else:
            res = sp.run(['wget', self.url, '-P', dir_] + self.wget_options)

            if res.returncode == 0:
                self.filename = dest_path

        return self.filename

    def unzip(self):
        if self.filename.endswith('.gz'):
            res = sp.run(['gunzip', self.filename])

            if res.returncode == 0:
                self.filename = re.sub(r'\.gz$', '', self.filename)

        return self.filename

    def rename(self, dest_name):
        if os.path.exists(dest_name):
            raise FileExistsError

        os.renames(self.filename, dest_name)
        self.filename = dest_name

        return self.filename


class FTPResource(Resource):
    ftp_site = ""
    release_dir = ""
    release_pattern = r""
    file_pattern = r""
    another_pats = ()

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

            try:
                ftp.cwd(dir_path)
            except ftplib.error_perm:
                raise URLError(dir_path) from None

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
            m1 = re.search(cls.file_pattern, name)
            if m1:
                return m1.group(0)
        else:
            for pat in cls.another_pats:
                for name in dir_list:
                    m2 = re.search(pat, name)
                    if m2:
                        return m2.group(0)
            else:
                return ""


class EnsemblResource(FTPResource):
    ftp_site = "ftp.ensembl.org"
    release_dir = "/pub"
    release_pattern = r"^release-([0-9]+)$"


class EnsemblGenome(EnsemblResource):
    file_pattern = r".+\.dna\.primary_assembly\.fa\.gz$"
    another_pats = (r".+\.dna\.toplevel\.fa\.gz$",)

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


class EnsemblCDna(EnsemblResource):
    file_pattern = r".+\.cdna\.all\.fa\.gz$"

    def get_dir_path(self):
        dir_path = os.path.join(
            '/',
            'pub',
            'release-{}'.format(self.version),
            'fasta',
            self.species,
            'cdna'
        )

        return dir_path


class EnsemblNCRna(EnsemblResource):
    file_pattern = r".+\.ncrna\.fa\.gz$"

    def get_dir_path(self):
        dir_path = os.path.join(
            '/',
            'pub',
            'release-{}'.format(self.version),
            'fasta',
            self.species,
            'ncrna'
        )

        return dir_path


class EnsemblAnnotation(EnsemblResource):
    file_pattern = r".+\.chr\.gtf\.gz$"
    another_pats = (r".+\.[0-9]+\.gtf\.gz$",)

    def get_dir_path(self):
        dir_path = os.path.join(
            '/',
            'pub',
            'release-{}'.format(self.version),
            'gtf',
            self.species
        )

        return dir_path


class EnsemblSisterResource(FTPResource):
    ftp_site = "ftp.ensemblgenomes.org"
    release_dir = "/pub"
    release_pattern = r"^release-([0-9]+)$"

    def __init__(self, field, species, version):
        self.field = field
        super().__init__(species, version)


class EnsemblSisterGenome(EnsemblSisterResource):
    file_pattern = r".+\.dna\.primary_assembly\.fa\.gz$"
    another_pats = (r".+\.dna\.toplevel\.fa\.gz$",)

    def get_dir_path(self):
        dir_path = os.path.join(
            '/',
            'pub',
            'release-{}'.format(self.version),
            self.field,
            'fasta',
            self.species,
            'dna'
        )

        return dir_path


class EnsemblSisterCDna(EnsemblSisterResource):
    file_pattern = r".+\.cdna\.all\.fa\.gz$"

    def get_dir_path(self):
        dir_path = os.path.join(
            '/',
            'pub',
            'release-{}'.format(self.version),
            self.field,
            'fasta',
            self.species,
            'cdna'
        )

        return dir_path


class EnsemblSisterNCRna(EnsemblSisterResource):
    file_pattern = r".+\.ncrna\.fa\.gz$"

    def get_dir_path(self):
        dir_path = os.path.join(
            '/',
            'pub',
            'release-{}'.format(self.version),
            self.field,
            'fasta',
            self.species,
            'ncrna'
        )

        return dir_path


class EnsemblSisterAnnotation(EnsemblSisterResource):
    file_pattern = r".+\.chr\.gtf\.gz$"
    another_pats = (r".+\.[0-9]+\.gtf\.gz$",)

    def get_dir_path(self):
        dir_path = os.path.join(
            '/',
            'pub',
            'release-{}'.format(self.version),
            self.field,
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


class GencodeProteinCodingTranscripts(GencodeResource):
    file_pattern = r".+\.pc_transcripts\.fa\.gz$"


class GencodeLongNonCodingTranscripts(GencodeResource):
    file_pattern = r".+\.lncRNA_transcripts\.fa\.gz$"


class GencodeAnnotation(GencodeResource):
    file_pattern = r"^gencode.v[^.]+.annotation.gtf.gz$"


class MiRBaseResource(Resource):
    url_templ = "https://www.mirbase.org/ftp/{version}/{file_pattern}"
    wget_options = ['--no-check-certificate']

    def __init__(self, species, version):
        self.species = species
        self.version = version
        url = self.get_url()
        super().__init__(url)

    def get_url(self):
        url = self.url_templ.format(
            version=self.version,
            file_pattern=self.file_pattern
        )
        return url

    def rename_with_version(self):
        dirname, filename = os.path.split(self.filename)
        name_with_ver = self._add_version_to_name(filename, self.version)
        dest_name = os.path.join(dirname, name_with_ver)
        self.rename(dest_name)

        return self.filename

    @staticmethod
    def _add_version_to_name(filename, version):
        sep_name = filename.split('.')
        sep_name.insert(1, "v{}".format(version))
        return '.'.join(sep_name)


class MiRBaseMiRNA(MiRBaseResource):
    file_pattern = "mature.fa.gz"


class MiRBaseDiff(MiRBaseResource):
    file_pattern = "miRNA.diff.gz"


class MiRBaseDat(MiRBaseResource):
    file_pattern = "miRNA.dat.gz"


class MiRTarBaseResource(Resource):
    url_templ = "https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2019/cache/download/{version}/miRTarBase_MTI.xlsx"
    url_templ_2 = "ftp://treeslab1.genomics.sinica.edu.tw/CircMiMi/miRTarBase_backup/{version}/miRTarBase_MTI.xlsx"

    def __init__(self, species, version):
        self.species = species
        self.version = version
        url = self.get_url()
        super().__init__(url)

    def get_url(self):
        url = self.url_templ.format(species=self.species, version=self.version)

        try:
            urllib.request.urlopen(url, timeout=10)

        except urllib.error.URLError:
            url = self.url_templ_2.format(species=self.species, version=self.version)

        return url


class MiRDBResource(Resource):
    url_templ = "http://www.mirdb.org/download/miRDB_v{version}_prediction_result.txt.gz"
    url_templ_2 = "http://www.mirdb.org/download/MirTarget2_v{version}_prediction_result.txt.gz"

    def __init__(self, species, version):
        self.species = species
        self.version = version
        url = self.get_url()
        super().__init__(url)

    def get_url(self):
        if self.version in ['1.0', '2.0', '3.0', '4.0']:
            url_templ = self.url_templ_2
        else:
            url_templ = self.url_templ

        url = url_templ.format(
            species=self.species,
            version=self.version
        )

        return url


class Gene2AccessionResource(Resource):
    url = "ftp://ftp.ncbi.nih.gov/gene/DATA/gene2accession.gz"

    def __init__(self):
        super().__init__(self.url)


class MiRDBData(Resource):
    url_templ = "ftp://treeslab1.genomics.sinica.edu.tw/CircMiMi/miRDB_data/v{version}/miRDB_v{version}_{species}.tsv.gz"
    available_species = ('hsa', 'mmu', 'rno', 'cfa', 'gga')

    def __init__(self, species, version):

        if species in self.available_species:
            self.species = species
        else:
            self.species = "others"

        self.version = version
        url = self.get_url()
        super().__init__(url)

    def get_url(self):
        url = self.url_templ.format(species=self.species, version=self.version)
        return url


class EncoriRBPData(Resource):
    url_templ = "ftp://treeslab1.genomics.sinica.edu.tw/CircMiMi/ENCORI_RBP/ENCORI_RBP_binding_sites.{}.{}.bed.gz"
    url_templ_2 = "ftp://treeslab1.genomics.sinica.edu.tw/CircMiMi/ENCORI_RBP/ENCORI_RBP_binding_sites.{}.no_chr.{}.bed.gz"

    def __init__(self, species, source, version, only_AGO=False):
        self.species = species
        self.source = source
        self.version = version
        self.only_AGO = only_AGO

        url = self.get_url()
        super().__init__(url)

    def get_url(self):
        if self.species == 'hsa':
            assembly = 'hg38'
        elif self.species == 'mmu':
            if self.version == 'current':
                assembly = 'mm39'
            else:
                if self.source == 'gencode':
                    if GencodeResource.get_digit_part(self.version) >= 26:
                        assembly = 'mm39'
                    else:
                        assembly = 'mm10'
                elif self.source == 'ensembl':
                    if int(self.version) >= 103:
                        assembly = 'mm39'
                    else:
                        assembly = 'mm10'
        else:
            return ''

        if self.source == 'gencode':
            url_templ = self.url_templ
        else:
            url_templ = self.url_templ_2

        if self.only_AGO:
            url = url_templ.format(assembly, 'AGO')
        else:
            url = url_templ.format(assembly, 'all')

        return url


class EncoriRBPTargetData(Resource):
    url_templ = "ftp://treeslab1.genomics.sinica.edu.tw/CircMiMi/ENCORI_RBP/ENCORI_RBP_Target.grouped.tsv.gz"
    url_templ_2 = "ftp://treeslab1.genomics.sinica.edu.tw/CircMiMi/ENCORI_RBP/ENCORI_mouse_RBP_Target.grouped.tsv.gz"

    def __init__(self, species):
        self.species = species

        url = self.get_url()
        super().__init__(url)

    def get_url(self):
        if self.species == 'hsa':
            url = self.url_templ.format()
        elif self.species == 'mmu':
            url = self.url_templ_2.format()
        else:
            return ''

        return url


class EncoriMiRNATargetData(Resource):
    url_templ = "ftp://treeslab1.genomics.sinica.edu.tw/CircMiMi/ENCORI_miRNA/mir_target_ref.ENCORI.miRBase_v22.tsv.gz"
    url_templ_2 = "ftp://treeslab1.genomics.sinica.edu.tw/CircMiMi/ENCORI_miRNA/mir_target_ref.ENCORI_mouse.miRBase_v22.tsv.gz"

    def __init__(self, species):
        self.species = species

        url = self.get_url()
        super().__init__(url)

    def get_url(self):
        if self.species == 'hsa':
            url = self.url_templ.format()
        elif self.species == 'mmu':
            url = self.url_templ_2.format()
        else:
            return ''

        return url


class Error(Exception):
    pass


class URLError(Error):
    pass


class SourceNotSupportError(Error):
    pass
