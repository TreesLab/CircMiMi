import os
import configparser


DEFAULT_REF_CONFIG = 'refs.cfg'

CONFIG_TEMPL = """
    [info]
    species =
    source =
    version =

    [refs]
    anno_db =
    ref_file =
    mir_ref =
    mir_target =
    other_transcripts =
    AGO_data =
    RBP_data =
    RBP_target =
"""


class RefConfig:
    def __init__(self):
        self.config = configparser.ConfigParser()
        self.config.optionxform = str
        self.config.read_string(CONFIG_TEMPL)

    def read(self, cfg_file):
        self.config.read(cfg_file)

    def write(self, ref_dir, cfg_file_name=DEFAULT_REF_CONFIG):
        cfg_file = os.path.join(ref_dir, cfg_file_name)
        with open(cfg_file, 'w') as config_file:
            self.config.write(config_file)

    def __getitem__(self, key):
        return self.config[key]


def get_refs(ref_dir):
    cfg_file = os.path.join(ref_dir, DEFAULT_REF_CONFIG)

    configObj = RefConfig()
    configObj.read(cfg_file)
    config = configObj.config

    anno_db = prepend_dirname_to_file(ref_dir, config['refs']['anno_db'])
    ref_file = prepend_dirname_to_file(ref_dir, config['refs']['ref_file'])
    mir_ref = prepend_dirname_to_file(ref_dir, config['refs']['mir_ref'])
    mir_target = prepend_dirname_to_file(ref_dir, config['refs']['mir_target'])
    other_transcripts = prepend_dirname_to_file(ref_dir, config['refs']['other_transcripts'])
    AGO_data = prepend_dirname_to_file(ref_dir, config['refs']['AGO_data'])
    RBP_data = None #prepend_dirname_to_file(ref_dir, config['refs']['RBP_data'])
    RBP_target = None #prepend_dirname_to_file(ref_dir, config['refs']['RBP_target'])

    return (anno_db,
            ref_file,
            mir_ref,
            mir_target,
            other_transcripts,
            AGO_data,
            RBP_data,
            RBP_target)


def prepend_dirname_to_file(ref_dir, filename):
    if filename == '':
        return None
    else:
        return os.path.join(ref_dir, filename)
