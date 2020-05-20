#! /usr/bin/env python
import click
import os


@click.group()
@click.version_option()
def cli():
    """
    A toolset for investigating the interactions between circRNA, miRNA, and mRNA.
    """

    pass


@cli.command()
@click.option('-r', '--ref', 'ref_dir', type=click.Path(), metavar="REF_DIR", required=True)
@click.option('-i', '--circ', 'circ_file', metavar="CIRC_FILE", required=True)
@click.option('-o', '--out-prefix', 'out_prefix', default='./out/', metavar="OUT_PREFIX")
@click.option('-p', '--num_proc', default=1, type=click.INT, metavar="NUM_PROC",
    help="Number of processes")
@click.option('--checkAA', 'checkAA', is_flag=True, help="Check if the circRNA has ambiguous alignments.")
@click.option('--header', 'header', flag_value=True, type=click.BOOL,
              default=True, hidden=True)
@click.option('--no-header', 'header', flag_value=False, type=click.BOOL, hidden=True)
@click.option('--miranda-sc', 'sc', metavar='S', type=click.FLOAT)
@click.option('--miranda-en', 'en', metavar='-E', type=click.FLOAT)
@click.option('--miranda-scale', 'scale', metavar='Z', type=click.FLOAT)
@click.option('--miranda-strict', 'strict', is_flag=True)
@click.option('--miranda-go', 'go', metavar='-X', type=click.FLOAT)
@click.option('--miranda-ge', 'ge', metavar='-Y', type=click.FLOAT)
def run(circ_file, ref_dir, out_prefix, num_proc, header, checkAA, **miranda_options):
    """
    This command is the main pipeline of CircMiMi.
    """

    from circmimi.utils import add_prefix

    output_dir = os.path.dirname(out_prefix)
    if output_dir == '':
        output_dir = '.'

    if output_dir != '.':
        os.makedirs(output_dir, exist_ok=True)

    from circmimi.config import get_refs
    anno_db, ref_file, mir_ref, mir_target, other_transcripts = get_refs(ref_dir)

    summary_file = add_prefix('summary_list.tsv', out_prefix)
    clear_file = add_prefix('circ.clear.tsv', out_prefix)

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
        other_ref_file,
        work_dir=output_dir,
        num_proc=num_proc,
        miranda_options=miranda_options_list
    )

    circmimi_result.run(circ_file)

    result_table = circmimi_result.get_result_table()
    res_file = add_prefix('all_interactions.tsv', out_prefix)
    result_table.to_csv(res_file, sep='\t', index=False, header=header)
    circmimi_result.save_circRNAs_summary(summary_file)
    circmimi_result.save_clear_circRNAs(clear_file)


@cli.command()
@click.option('--species', 'species', metavar="SPECIES_KEY", required=True)
@click.option('--source', 'source', metavar="SOURCE", required=True,
    help="""
        Available sources are "gencode", "ensembl", "ensembl_plants", and "ensembl_metazoa"
    """)
@click.option('--gencode', 'source', flag_value='gencode', type=click.STRING, help="--source gencode", hidden=True)
@click.option('--ensembl', 'source', flag_value='ensembl', type=click.STRING, help="--source ensembl", hidden=True)
@click.option('--version', 'version', default='current', metavar="VERSION",
    help="""
        The release version. If not assigned, it will be automatically set to the latest version of the SOURCE.
    """)
@click.option('--init', 'init', is_flag=True, help="Create an init template ref_dir.", hidden=True)
@click.argument('ref_dir')
def genref(species, source, version, ref_dir, init):
    """
    Generate the references.                                          
    
    The generated reference files will be saved in the REF_DIR.
    """

    os.makedirs(ref_dir, exist_ok=True)

    from circmimi.config import RefConfig
    config = RefConfig()

    if not init:
        from circmimi.reference import genref
        info, ref_files = genref.generate(species, source, version, ref_dir)

        config['info'].update(info)
        config['refs'].update(ref_files)

    config.write(ref_dir)


@cli.command(hidden=True)
@click.argument('gtf_path')
@click.argument('db_path', metavar='OUT_PATH')
def gendb(gtf_path, db_path):
    from circmimi.reference import gendb

    gendb.generate(gtf_path, db_path)


@cli.command(hidden=True)
@click.option('--species', 'species', metavar="SPECIES_KEY", required=True)
@click.option('--version', 'version', default='current', metavar="VERSION", required=True)
@click.option('-r', '--ref', 'ref_dir', type=click.Path(), metavar="REF_DIR", required=True)
@click.option('-o', '--out_file', 'out_file', metavar="OUT_FILE", required=True)
@click.option('-a', '--show-accession', 'show_accession', is_flag=True)
def genmirdb(species, version, ref_dir, out_file, show_accession):
    os.makedirs(ref_dir, exist_ok=True)

    from circmimi.reference import genmirdb

    genmirdb.generate(species, version, ref_dir, out_file, show_accession)


@cli.group(hidden=True)
def check():
    pass


@check.command()
@click.argument('circ_file')
@click.option('-r', '--ref', 'ref_dir', type=click.Path(), metavar="REF_DIR", required=True)
@click.option('-o', '--out-prefix', 'out_prefix', default='./out/', metavar="OUT_PREFIX")
def annotation(circ_file, ref_dir, out_prefix):
    from circmimi.config import get_refs
    from circmimi.circ import CircEvents
    from circmimi.utils import add_prefix

    output_dir = os.path.dirname(out_prefix)
    if output_dir == '':
        output_dir = '.'

    if output_dir != '.':
        os.makedirs(output_dir, exist_ok=True)

    anno_db, _, _, _, _ = get_refs(ref_dir)

    circ_events = CircEvents(circ_file)
    circ_events.check_annotation(anno_db)
    res_file = add_prefix('circ.summary_list.tsv', out_prefix)
    circ_events.get_summary().to_csv(res_file, sep='\t', index=False)


@check.command()
@click.argument('circ_file')
@click.option('-r', '--ref', 'ref_dir', type=click.Path(), metavar="REF_DIR", required=True)
@click.option('-o', '--out', 'output_dir', metavar="OUT_DIR", required=True)
@click.option('-p', '--num_proc', default=1, type=click.INT, metavar="NUM_PROC",
    help="Number of processes")
def ambiguous(circ_file, ref_dir, output_dir, num_proc):
    os.makedirs(output_dir, exist_ok=True)

    from circmimi.config import get_refs
    anno_db, ref_file, _, _, other_transcripts = get_refs(ref_dir)

    from circmimi.circ import CircEvents
    circ_events = CircEvents(circ_file)
    circ_events.check_ambiguous(
        anno_db,
        ref_file,
        other_transcripts,
        work_dir=output_dir,
        num_proc=num_proc
    )

    result_file = os.path.join(output_dir, 'circ.summary_list.tsv')
    circ_events.get_summary().to_csv(result_file, sep='\t', index=False)


@cli.group()
def network():
    """
    Command to create the network file for Cytoscape.
    """

    pass


@network.command()
@click.argument('in_file')
@click.argument('out_file')
@click.option('-f', '--format', 'format_', default='xgmml',
              help="Assign the format of the OUT_FILE.", hidden=True)
def create(in_file, out_file, format_):
    """
    Create the network file.

    This command would transform the IN_FILE to XGMML format,
    and save to OUT_FILE.
    """

    from circmimi.network.network import CyNetwork, Layout, Style

    network = CyNetwork()
    network.load_data(in_file)

    layout = Layout()
    style = Style()

    network.apply_layout(layout)
    network.apply_style(style)

    if format_ == "xgmml":
        network.to_xgmml(out_file)


@cli.group(hidden=True)
def update_mirna():
    pass


@update_mirna.command()
@click.option('--from', 'from_', required=True)
@click.option('--to', 'to_', required=True)
@click.option('--species')
@click.option('-o', '--out-prefix', 'out_prefix', default='./', metavar="OUT_PREFIX")
def genmaps(from_, to_, species, out_prefix):
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


@update_mirna.command()
@click.argument('in_file')
@click.argument('out_file')
@click.option('-m', '--maps', 'mapping_file', metavar="MAPPING_FILE", required=True)
@click.option('-k', '--col-key', 'column_key', type=click.INT, default=1,
              help="The column number of miRNA IDs.")
@click.option('-i', '--inplace', is_flag=True)
@click.option('-R', '--remove-deleted', is_flag=True)
def update(in_file, out_file, mapping_file, column_key, inplace, remove_deleted):
    column_key = column_key - 1

    from circmimi.reference.mirbase import MatureMiRNAUpdater
    updater = MatureMiRNAUpdater(None, None, None)
    updater.load_maps(mapping_file)
    updater.update_file(in_file, out_file, column_key, inplace, remove_deleted)


if __name__ == "__main__":
    cli()
