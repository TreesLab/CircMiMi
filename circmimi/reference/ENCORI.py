"""
- https://starbase.sysu.edu.cn/api/ref/starBase_referenceData.zip
    - refData/hg19_clip_ref.txt
    - refData/mm10_clip_ref.txt
- https://starbase.sysu.edu.cn/api/bindingSite/?assembly= *** hg19 *** &datasetID= *** SBDH27 ***
    - 
"""


import requests
import io
import csv
import os
from operator import itemgetter
from tqdm import tqdm
from circmimi.reference.utils import cwd, read_zip_file

import logging

logger = logging.getLogger(__name__)

class EncoriClipData:
    ASSEMBLY = {
        'hsa': 'hg19',
        'mmu': 'mm10',
        'cel': 'ce10',
        'dme': 'dm6',
        'sce': 'sacCer3',
    }
    ref_data = "https://starbase.sysu.edu.cn/api/ref/starBase_referenceData.zip"
    clip_file_templ = "refData/{}_clip_ref.txt"
    url_templ = "https://starbase.sysu.edu.cn/api/bindingSite/?assembly={}&datasetID={}"


    def __init__(self, species):
        logger.debug("Initializing")

        self.species = species
        self.assembly = self.ASSEMBLY.get(species, 'unknown')
        logger.debug(f"species: {self.species}, assembly: {self.assembly}")

        logger.debug('creating headers')
        self._headers = self._create_headers()

        logger.debug('getting all datasets info')
        self.all_datasets_info = self.get_all_datasets_info(self.assembly, self._headers)

        logger.info(
            f"species: {self.species}, assembly: {self.assembly}, "
            f"num_datasets: {len(self.all_datasets_info)}"
        )

    @staticmethod
    def _create_headers():
        headers = requests.utils.default_headers()
        headers.update(
            {
                'User-Agent': 'My User Agent 1.0',
            }
        )
        return headers

    @classmethod
    def get_all_datasets_info(cls, assembly, headers):
        logger.debug("making requests for datasets infos")
        r = requests.get(cls.ref_data, headers=headers)
        logger.debug("reading the zip file")
        info_data = read_zip_file(
            io.BytesIO(r.content),
            cls.clip_file_templ.format(assembly)
        )
        reader = csv.DictReader(info_data, delimiter='\t')

        get_info = itemgetter('datasetID', 'GeneSymbol')
        all_datasets_info = [get_info(data) for data in reader]
        logger.debug("returning the datasets info")
        return all_datasets_info

    @classmethod
    def get_data(cls, assembly, dataset_id, headers):
        logger.debug(f"fetching data of ({assembly}, {dataset_id})")
        try:
            r = requests.get(cls.url_templ.format(assembly, dataset_id), headers=headers)
        except requests.exceptions.ConnectionError as e:
            return b""

        logger.debug("checking response")

        if r.ok:
            logger.debug(f"returning content with len = {len(r.content)}")
            return r.content
        else:
            return b""


    def download(self, dir_='.'):
        self._raw_data_dir = dir_

        logger.info("start to download data files")
        os.makedirs(dir_, exist_ok=True)

        with cwd(dir_):
            for dataset_id, gene_symbol in tqdm(self.all_datasets_info):
                out_filename = f"{dataset_id}_{gene_symbol}.txt"

                if os.path.exists(out_filename):
                    continue

                data = self.get_data(self.assembly, dataset_id, self._headers)

                with open(out_filename, 'wb') as out:
                    out.write(data)

    @staticmethod
    def _rename_id(old_id, RBP):
        dataset_id, idx = old_id.split('-')
        new_id = '_'.join([dataset_id, RBP, idx])
        return new_id

    def merge_files(self, out_file):
        self.merged_file = out_file

        files_to_be_merged = [
            os.path.join(self._raw_data_dir, f"{dataset_id}_{gene_symbol}.txt")
            for dataset_id, gene_symbol in self.all_datasets_info
        ]

        with open(out_file, 'w') as out:
            for file_ in files_to_be_merged:
                gene_symbol = os.path.basename(file_).rstrip('.txt').split('_')[1]

                with open(file_) as f_in:
                    for line in f_in:
                        if not line.startswith('#'):
                            res = line.rstrip('\n').split('\t')
                            res[3] = self._rename_id(res[3], gene_symbol)
                            print(*res, sep='\t', file=out)
