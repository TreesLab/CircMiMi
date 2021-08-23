#! /usr/bin/env python
import click
import os
import logging


logger = logging.getLogger(__name__)


@click.group()
@click.version_option()
@click.option('--debug', 'debug_mode', is_flag=True)
def cli(debug_mode):
    """
    A toolset for investigating the interactions between circRNA, miRNA, and mRNA.
    """

    root_logger = logging.getLogger()

    ch = logging.StreamHandler()

    if debug_mode:
        logger.warning('Under debug mode')
        root_logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s (%(name)s) [%(levelname)s] %(message)s')
    else:
        root_logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(message)s')

    ch.setFormatter(formatter)
    root_logger.addHandler(ch)


@cli.command('interactions')
@click.option('-r', '--ref', 'ref_dir', type=click.Path(), metavar="REF_DIR", required=True)
@click.option('-i', '--circ', 'circ_file', metavar="CIRC_FILE", required=True)
@click.option('-o', '--out-prefix', 'out_prefix', default='./', metavar="OUT_PREFIX")
@click.option('-p', '--num_proc', default=1, type=click.INT, metavar="NUM_PROC",
    help="Number of processes")
@click.option('--checkAA', 'checkAA', is_flag=True,
    help="Check if the circRNA has ambiguous alignments.", hidden=True)
@click.option('--miranda-sc', 'sc', metavar='S', type=click.FLOAT)
@click.option('--miranda-en', 'en', metavar='-E', type=click.FLOAT)
@click.option('--miranda-scale', 'scale', metavar='Z', type=click.FLOAT)
@click.option('--miranda-strict', 'strict', is_flag=True)
@click.option('--miranda-go', 'go', metavar='-X', type=click.FLOAT)
@click.option('--miranda-ge', 'ge', metavar='-Y', type=click.FLOAT)
def predict_interactions(circ_file, ref_dir, out_prefix, num_proc, checkAA, **miranda_options):
    """
    Predict the interactions.

    Predict the interactions between circRNA-miRNA-mRNA.
    """

    logger.info('Preparing ...')

    from circmimi.utils import add_prefix

    output_dir = os.path.dirname(out_prefix)
    if output_dir == '':
        output_dir = '.'

    if output_dir != '.':
        os.makedirs(output_dir, exist_ok=True)

    from circmimi.reference.config import get_refs
    anno_db, ref_file, mir_ref, mir_target, other_transcripts, AGO_data, RBP_data, RBP_target = get_refs(ref_dir)

    if checkAA:
        other_ref_file = other_transcripts
    else:
        other_ref_file = None

    miranda_options_list = []
    for k, v in miranda_options.items():
        if v is not None:
            cmd = ['-' + k]
            if k != 'strict':
                cmd.append(str(v))

            miranda_options_list.extend(cmd)

    from circmimi.circmimi import Circmimi
    circmimi_result = Circmimi(
        anno_db,
        ref_file,
        mir_ref,
        mir_target,
        AGO_data,
        RBP_data,
        RBP_target,
        other_ref_file,
        work_dir=output_dir,
        num_proc=num_proc,
        miranda_options=miranda_options_list
    )

    logger.info('Starting the main pipeline.')
    circmimi_result.run(circ_file)
    logger.info('Pipeline completed.')

    logger.info('Saving results ...')
    res_file = add_prefix('all_interactions.miRNA.tsv', out_prefix)
    circmimi_result.save_result(res_file)
    logger.info('miRNA part ... done')

    # RBP_res_file = add_prefix('all_interactions.RBP.tsv', out_prefix)
    # circmimi_result.save_RBP_result(RBP_res_file)
    # logger.info('RBP part ... done')

    summary_file = add_prefix('summary_list.tsv', out_prefix)
    circmimi_result.save_circRNAs_summary(summary_file)
    logger.info('summary file ... done')
    logger.info('All results are saved.')

    logger.info('Process completed.')


@cli.command('genref')
@click.option('--species', 'species', metavar="SPECIES_KEY", required=True, help="eg. hsa, mmu, rno, ...")
@click.option('--source', 'source', metavar="SOURCE", required=True,
    help="""
        Available sources are "gencode", "ensembl", "ensembl_plants", and "ensembl_metazoa"
    """)
@click.option('--gencode', 'source', flag_value='gencode', type=click.STRING, help="--source gencode", hidden=True)
@click.option('--ensembl', 'source', flag_value='ensembl', type=click.STRING, help="--source ensembl", hidden=True)
@click.option('--version', 'version', default='current', metavar="VERSION",
    help="""
        The release version. If not specified, it will be automatically set to the latest version of the SOURCE.
    """)
@click.option('--init', 'init', is_flag=True, help="Create an init template ref_dir.", hidden=True)
@click.argument('ref_dir')
def generate_references(species, source, version, ref_dir, init):
    """
    Generate the references.                                          

    This command is used to generate the needed references for CircMiMi.
    It would download the genome ane annotation files automatically, with user-specified version.

    The generated reference files will be saved in the REF_DIR.
    """

    os.makedirs(ref_dir, exist_ok=True)

    from circmimi.reference.config import RefConfig
    config = RefConfig()

    if not init:
        from circmimi.reference import genref
        info, ref_files = genref.generate(species, source, version, ref_dir)

        config['info'].update(info)
        config['refs'].update(ref_files)

    config.write(ref_dir)


@cli.command('gendb', hidden=True)
@click.argument('gtf_path')
@click.argument('db_path', metavar='OUT_PATH')
def generate_annotation_database(gtf_path, db_path):
    from circmimi.reference import gendb

    gendb.generate(gtf_path, db_path)


@cli.command('genmirdb', hidden=True)
@click.option('--species', 'species', metavar="SPECIES_KEY", required=True)
@click.option('--version', 'version', default='current', metavar="VERSION", required=True)
@click.option('-r', '--ref', 'ref_dir', type=click.Path(), metavar="REF_DIR", required=True)
@click.option('-o', '--out_file', 'out_file', metavar="OUT_FILE", required=True)
@click.option('-a', '--show-accession', 'show_accession', is_flag=True)
def generate_miRNA_database(species, version, ref_dir, out_file, show_accession):
    os.makedirs(ref_dir, exist_ok=True)

    from circmimi.reference import genmirdb

    genmirdb.generate(species, version, ref_dir, out_file, show_accession)


@cli.group(hidden=True)
def check():
    pass


@check.command('annotation')
@click.argument('anno_db')
@click.argument('circ_file')
@click.argument('out_file')
def check_annotation(circ_file, anno_db, out_file):
    from circmimi.circ import CircEvents

    circ_events = CircEvents(circ_file)
    circ_events.check_annotation(anno_db)

    circ_events.get_filters_results(True).to_csv(out_file, sep='\t', index=False)


@check.command('ambiguous')
@click.argument('ref_file')
@click.argument('other_ref_file')
@click.argument('circ_file')
@click.argument('out_file')
@click.option('-p', '--num_proc', default=1, type=click.INT,
              metavar="NUM_PROC", help="Number of processes")
def check_ambiguous(ref_file, other_ref_file, circ_file, out_file, num_proc):
    import subprocess as sp
    import tempfile as tp

    tmp_fa = tp.NamedTemporaryFile(
        dir='.',
        prefix='flanking_seq.',
        suffix='.fa'
    )

    sp.run(['get_flanking_seq.py', ref_file, circ_file, tmp_fa.name])
    sp.run(['checkAA_reads.py', '-rg', ref_file, '-ro', other_ref_file,
            tmp_fa.name, out_file, '-p', str(num_proc)])


@check.command('RCS')
@click.argument('ref_file')
@click.argument('circ_file')
@click.argument('out_file')
@click.option('-d', '--dist', default=10000, type=click.INT)
@click.option('-m', '--min-matches', default=80, type=click.FLOAT)
@click.option('-l', '--min-aln-len', default=50, type=click.INT)
@click.option('-b', '--min-bitscore', default=100, type=click.FLOAT)
@click.option('-p', '--num_proc', default=1, type=click.INT,
              metavar="NUM_PROC", help="Number of processes")
def check_RCS(ref_file,
              circ_file,
              out_file,
              dist,
              min_matches,
              min_aln_len,
              min_bitscore,
              num_proc):

    import subprocess as sp

    cmd1 = ['get_RCS.py', ref_file, circ_file, '--dist', dist, '-p', num_proc]
    cmd2 = ['RCS_filter.py', '-', '-m', min_matches, '-l', min_aln_len, '-b', min_bitscore]
    cmd3 = ['get_RCS_summary.py', '--circ_file', circ_file, '-']

    cmd1 = [str(c) for c in cmd1]
    cmd2 = [str(c) for c in cmd2]
    cmd3 = [str(c) for c in cmd3]

    with open(out_file, 'w') as out:
        p1 = sp.Popen(cmd1, stdout=sp.PIPE, encoding='UTF-8')
        p2 = sp.Popen(cmd2, stdin=p1.stdout, stdout=sp.PIPE, encoding='UTF-8')
        p3 = sp.Popen(cmd3, stdin=p2.stdout, stdout=out, encoding='UTF-8')

        p1.wait()
        p2.wait()
        p3.wait()


@cli.command('checking')
@click.option('-i', '--circ', 'circ_file', metavar="CIRC_FILE", required=True)
@click.option('-r', '--ref', 'ref_dir', type=click.Path(), metavar="REF_DIR",
              required=True)
@click.option('-o', '--out-prefix', 'out_prefix', default='./',
              metavar="OUT_PREFIX")
@click.option('-p', '--num_proc', default=1, type=click.INT,
              metavar="NUM_PROC", help="Number of processes")
@click.option('-d', '--dist', default=10000, type=click.INT)
@click.pass_context
def default_checking(ctx, circ_file, ref_dir, out_prefix, num_proc, dist):
    """
    Check the circRNAs.

    Check the circRNAs if there are
    (i) ambiguous alignments for the fragment around the circRNA junction;
    (ii) RCS pair in the flanking sequences.
    """

    import tempfile as tp
    from circmimi.reference.config import get_refs
    from circmimi.utils import add_prefix

    output_dir = os.path.dirname(out_prefix)
    if output_dir == '':
        output_dir = '.'

    if output_dir != '.':
        os.makedirs(output_dir, exist_ok=True)

    anno_db, ref_file, _, _, other_transcripts, _, _, _ = get_refs(ref_dir)

    annotation_result_file = tp.NamedTemporaryFile(
        dir=output_dir,
        prefix='check.annotation.',
        suffix='.tsv'
    )

    ambiguous_result_file = tp.NamedTemporaryFile(
        dir=output_dir,
        prefix='check.AA.',
        suffix='.tsv'
    )
    RCS_result_file = tp.NamedTemporaryFile(
        dir=output_dir,
        prefix='check.RCS.',
        suffix='.tsv'
    )
    out_file = add_prefix('checking.results.tsv', out_prefix)

    ctx.invoke(
        check_annotation,
        circ_file=circ_file,
        anno_db=anno_db,
        out_file=annotation_result_file.name
    )

    ctx.invoke(
        check_ambiguous,
        ref_file=ref_file,
        other_ref_file=other_transcripts,
        circ_file=circ_file,
        out_file=ambiguous_result_file.name,
        num_proc=num_proc
    )

    ctx.invoke(
        check_RCS,
        ref_file=ref_file,
        circ_file=circ_file,
        out_file=RCS_result_file.name,
        num_proc=num_proc,
        dist=dist
    )

    import pandas as pd

    ann_df = pd.read_csv(annotation_result_file.name, sep='\t', dtype='object')
    amb_df = pd.read_csv(ambiguous_result_file.name, sep='\t', dtype='object')
    rcs_df = pd.read_csv(RCS_result_file.name, sep='\t', dtype='object')
    result_df = pd.concat(
        [
            ann_df,
            amb_df.drop('circRNA_id', axis=1),
            rcs_df.drop('circRNA_id', axis=1)
        ],
        axis=1
    )
    result_df.to_csv(out_file, sep='\t', index=False)


@cli.group(hidden=True)
def network():
    """
    Command to create the network file for Cytoscape.
    """

    pass


@network.command('create')
@click.argument('in_file')
@click.argument('out_file')
@click.option('-1', 'idx_circRNA', type=int, default=1,
              help='column key for circRNAs.')
@click.option('-2', 'idx_mediator', type=int, default=2,
              help='column key for mediators.')
@click.option('-3', 'idx_mRNA', type=int, default=3,
              help='column key for mRNAs.')
@click.option('-f', '--format', 'format_', default='xgmml',
              help="Assign the format of the OUT_FILE.", hidden=True)
def create_network(in_file, out_file, idx_circRNA, idx_mediator, idx_mRNA, format_):
    """
    Create the network file.

    This command would transform the IN_FILE to XGMML format,
    and save to OUT_FILE.
    """

    from circmimi.network.network import CyNetwork, Layout, Style

    network = CyNetwork()
    network.load_data(
        in_file,
        k1=idx_circRNA,
        k2=idx_mediator,
        k3=idx_mRNA
    )

    layout = Layout()
    style = Style()

    network.apply_layout(layout)
    network.apply_style(style)

    if format_ == "xgmml":
        network.to_xgmml(out_file)


@cli.command('visualize')
@click.argument('in_file')
@click.argument('out_file')
@click.option('-1', 'idx_circRNA', type=int, default=6,
              help='column key for circRNAs.')
@click.option('-2', 'idx_mediator', type=int, default=7,
              help='column key for mediators.')
@click.option('-3', 'idx_mRNA', type=int, default=13,
              help='column key for mRNAs.')
@click.option('-f', '--format', 'format_', default='xgmml',
              help="Assign the format of the OUT_FILE.", hidden=True)
@click.option('--header', 'header', flag_value=True, type=click.BOOL,
              default=True, hidden=True,
              help="Use this option to specify if there are titles in the IN_FILE.")
@click.option('--no-header', 'header', flag_value=False, type=click.BOOL,
              help="Use this option to specify if there are no titles in the IN_FILE.")
def visualize_interactions(in_file,
                           out_file,
                           idx_circRNA,
                           idx_mediator,
                           idx_mRNA,
                           header,
                           format_):

    """
    Visualize the interactions.

    This command is used to convert the file of interactions into XGMML format,
     which then can be loaded into Cytoscape for visualization and further analysis.
    """

    from circmimi.network.network import CyNetwork, Layout, Style

    network = CyNetwork()
    network.load_data(
        in_file,
        k1=idx_circRNA,
        k2=idx_mediator,
        k3=idx_mRNA,
        header=header
    )

    layout = Layout()
    style = Style()

    network.apply_layout(layout)
    network.apply_style(style)

    if format_ == "xgmml":
        network.to_xgmml(out_file)


@cli.group('update_mirna', hidden=True)
def update_mirna():
    pass


@update_mirna.command('genmaps')
@click.option('--from', 'from_', required=True)
@click.option('--to', 'to_', required=True)
@click.option('--species')
@click.option('-o', '--out-prefix', 'out_prefix', default='./', metavar="OUT_PREFIX")
def generate_miRNAs_mappings(from_, to_, species, out_prefix):
    from circmimi.reference.mirbase import MatureMiRNAUpdater
    from circmimi.utils import add_prefix

    output_dir, prefix_name = os.path.split(out_prefix)

    if output_dir == '':
        output_dir = '.'

    if output_dir != '.':
        os.makedirs(output_dir, exist_ok=True)

    updater = MatureMiRNAUpdater(from_, to_, species)
    updater.create(output_dir)

    out_file = 'miRNA.maps_{}_to_{}.tsv'.format(from_, to_)
    out_file = add_prefix(out_file, prefix_name)
    out_file = os.path.join(output_dir, out_file)

    updater.save(out_file)


@update_mirna.command('update')
@click.argument('in_file')
@click.argument('out_file')
@click.option('-m', '--maps', 'mapping_file', metavar="MAPPING_FILE", required=True)
@click.option('-k', '--col-key', 'column_key', type=click.INT, default=1,
              help="The column number of miRNA IDs.")
@click.option('-i', '--inplace', is_flag=True)
@click.option('-R', '--remove-deleted', is_flag=True)
def update_miRNAs(in_file, out_file, mapping_file, column_key, inplace, remove_deleted):
    column_key = column_key - 1

    from circmimi.reference.mirbase import MatureMiRNAUpdater
    updater = MatureMiRNAUpdater(None, None, None)
    updater.load_maps(mapping_file)
    updater.update_file(in_file, out_file, column_key, inplace, remove_deleted)


if __name__ == "__main__":
    cli()
