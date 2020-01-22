#! /usr/bin/env python
import click
import os


@click.group(help="""
    A toolset for investigating the interactions between circRNA, miRNA, and mRNA.

    Example.                                                    
    circmimi_tools genref --species hsa --source ensembl --version 98 ./refs
    circmimi_tools run -r ./refs -p 10 circ_events.tsv > out.tsv
    """)
def cli():
    pass


@cli.command(help="""
    Main pipeline.

    Example.                                                    
    circmimi_tools run -r ./refs -p 10 circ_events.tsv > out.tsv
    """)
@click.argument('circ_file')
@click.option('-r', '--ref', 'ref_dir', type=click.Path(), metavar="REF_DIR", required=True)
@click.option('-o', '--out', 'output_dir', metavar="OUT_DIR", required=True)
@click.option('-p', '--num_proc', default=1, type=click.INT, metavar="NUM_PROC",
    help="Number of processes")
@click.option('--checkAA', 'checkAA', is_flag=True, help="Check if the circRNA has ambiguous alignments.")
@click.option('--header', 'header', flag_value=True, type=click.BOOL,
              default=True, hidden=True)
@click.option('--no-header', 'header', flag_value=False, type=click.BOOL)
def run(circ_file, ref_dir, output_dir, num_proc, header, checkAA):
    os.makedirs(output_dir, exist_ok=True)

    from circmimi.config import get_refs
    anno_db, ref_file, mir_ref, mir_target, other_transcripts = get_refs(ref_dir)

    if checkAA:
        from circmimi.ambiguous import AmbiguousChecker
        checker = AmbiguousChecker(
            ref_file=ref_file,
            other_ref_file=other_transcripts,
            work_dir=output_dir,
            num_proc=num_proc
        )
        checker.check(circ_file)
        checker.save_result()
        circ_file = checker.clear_circ_file

    from circmimi.circmimi import Circmimi
    circmimi_result = Circmimi(work_dir=output_dir)
    circmimi_result.run(
        circ_file,
        anno_db,
        ref_file,
        mir_ref,
        mir_target,
        num_proc=num_proc
    )

    result_table = circmimi_result.get_result_table()
    res_file = os.path.join(output_dir, 'out.tsv')
    result_table.to_csv(res_file, sep='\t', index=False, header=header)


@cli.command(help="""
    Generate index and references.                                          
    The generated files would be saved in the directory REF_DIR.

    Example.                                                    
    circmimi_tools genref --species hsa --source ensembl --version 98 ./refs
    
    ---------------------                                                    
    | Available Species |                                                   
    --------------------------------------------------------------------------
    | Key | Name                    | Ensembl | Gencode | Alternative Source |
    | --- | ----------------------- | ------- | ------- | --------------------
    | ath | Arabidopsis thaliana    |    *    |         | Ensembl Plants     |
    | bmo | Bombyx mori             |    *    |         | Ensembl Metazoa    |
    | bta | Bos taurus              |    V    |         |                    |
    | cel | Caenorhabditis elegans  |    V    |         | Ensembl Metazoa    |
    | cfa | Canis familiaris        |    V    |         |                    |
    | dre | Danio rerio             |    V    |         |                    |
    | dme | Drosophila melanogaster |    V    |         |                    |
    | gga | Gallus gallus           |    V    |         |                    |
    | hsa | Homo sapiens            |    V    |    V    |                    |
    | mmu | Mus musculus            |    V    |    V    |                    |
    | osa | Oryza sativa            |    *    |         | Ensembl Plants     |
    | ola | Oryzias latipes         |    V    |         |                    |
    | oar | Ovis aries              |    V    |         |                    |
    | rno | Rattus norvegicus       |    V    |         |                    |
    | ssc | Sus scrofa              |    V    |         |                    |
    | tgu | Taeniopygia guttata     |    V    |         |                    |
    | xtr | Xenopus tropicalis      |    V    |         |                    |
    --------------------------------------------------------------------------
    * Only in the alternative source
    """)
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
@click.option('--mirdb', 'use_miRDB', is_flag=True)
@click.option('--init', 'init', is_flag=True, help="Create an init template ref_dir.", hidden=True)
@click.argument('ref_dir')
def genref(species, source, version, ref_dir, use_miRDB, init):
    os.makedirs(ref_dir, exist_ok=True)

    from circmimi.config import RefConfig
    config = RefConfig()

    if not init:
        from circmimi.reference import genref
        info, ref_files = genref.generate(species, source, version, ref_dir, use_miRDB)

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


@cli.command(hidden=True)
@click.argument('circ_file')
@click.option('-r', '--ref', 'ref_dir', type=click.Path(), metavar="REF_DIR", required=True)
@click.option('-o', '--out', 'output_dir', metavar="OUT_DIR", required=True)
@click.option('-p', '--num_proc', default=1, type=click.INT, metavar="NUM_PROC",
    help="Number of processes")
def checkaa(circ_file, ref_dir, output_dir, num_proc):
    os.makedirs(output_dir, exist_ok=True)

    from circmimi.config import get_refs
    _, ref_file, _, _, other_transcripts = get_refs(ref_dir)

    from circmimi.ambiguous import AmbiguousChecker
    checker = AmbiguousChecker(
        ref_file=ref_file,
        other_ref_file=other_transcripts,
        work_dir=output_dir,
        num_proc=num_proc
    )
    checker.check(circ_file)
    checker.save_result()


if __name__ == "__main__":
    cli()
