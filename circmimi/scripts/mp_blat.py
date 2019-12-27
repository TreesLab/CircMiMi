#! /usr/bin/env python
from __future__ import print_function

import argparse
import subprocess as sp
import os
import sys
import shutil
import random
import time
import re


def mp_blat(reference_file, fasta_file, output_file, num_of_process=1, tmp_path=".", blat_bin="blat", blat_opt=""):

    blat_opt = re.findall("([^ ]+)", blat_opt)

    if num_of_process == 1:
        # do blat directly
        sp.call([blat_bin, reference_file, fasta_file, output_file] + blat_opt)
    elif num_of_process > 1:

        # split file

        print("Start to split the input into {} part ... ".format(num_of_process), file=sys.stderr)
        split_start = time.time()

        fasta_file_size = os.stat(fasta_file).st_size
        with open(fasta_file) as fa_data:

            # get the number of all reads and their start pos in this file
            all_reads_pos = []
            while fa_data.tell() < fasta_file_size:
                next_char = fa_data.read(1)
                if next_char == ">":
                    all_reads_pos.append(fa_data.tell()-1)
            num_of_reads = len(all_reads_pos)

            fa_data.seek(0)

            if num_of_reads < num_of_process:
                splice_positions = all_reads_pos[1:]
            else:
                # eg. 70/20 = 3 ... 10
                #   ----> [3, 3, 3, ..., 3, 3, ..., 3]
                #   ----> [4, 4, 4, ..., 4, 3, ..., 3]
                #          ^^^^^^^^^^^^^^^
                #                10
                #          ^^^^^^^^^^^^^^^^^^^^^^^^^^
                #                      20
                #         
                #   ----> [4, 8, 12, ..., 40, 43, ..., 67]
                #          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                #                       19
                #
                #  Apply the index to "all_reads_pos"
                #   ----> [x[4], x[8], ..., x[67]]
                #   ----> [x[0]] + [x[4], x[8], ..., x[67]] + [file_size]
                #
                #  Get the size for each file
                #   ----> [x[4] - x[0], x[8] - x[4], ..., file_size - x[67]]
                #

                num_of_reads_for_each_blat = [num_of_reads // num_of_process]*num_of_process
                for i in range(num_of_reads % num_of_process):
                    num_of_reads_for_each_blat[i] += 1

                splice_positions = []
                idx = 0
                for n in num_of_reads_for_each_blat[:-1]:
                    idx += n
                    splice_positions.append(all_reads_pos[idx])

            edge_pos = [all_reads_pos[0]] + splice_positions + [fasta_file_size]
            read_size_for_each_blat = [edge_pos[1:][i] - edge_pos[:-1][i] for i in range(len(edge_pos)-1)]

            del all_reads_pos

            # create tmp files
            tmp_dir = "{}/mp_blat.tmp.{}".format(tmp_path.rstrip('/'), str(random.random())[2:8])
            os.mkdir(tmp_dir)

            tmp_fa_files = []
            for i, size in enumerate(read_size_for_each_blat):
                tmp_file_name = "{}/blat_tmp.part_{}.fa".format(tmp_dir, i+1)

                tmp_file = open(tmp_file_name, 'w')
                tmp_file.write(fa_data.read(size))
                tmp_file.close()

                tmp_fa_files.append(tmp_file_name)

        split_stop = time.time()
        print("Time cost of split_file = {} sec".format(split_stop - split_start), file=sys.stderr)


        # multi-processing blat

        print("Start to do the mp_blat using {} processes ... ".format(num_of_process), file=sys.stderr)
        blat_start = time.time()

        tmp_res_files = []
        for i, fa_file in enumerate(tmp_fa_files):
            tmp_file_name = "{}/blat_tmp.part_{}.psl".format(tmp_dir, i+1)
            tmp_res_files.append(tmp_file_name)

        proc = []
        for i in range(len(tmp_fa_files)):
            proc.append(sp.Popen([blat_bin, reference_file, tmp_fa_files[i], tmp_res_files[i]] + blat_opt))

        for p in proc:
            p.wait()

        blat_stop = time.time()
        print("Time cost of mp_blat = {} sec".format(blat_stop - blat_start), file=sys.stderr)

        # merge results
        print("Merging results ... ", file=sys.stderr)
        merge_start = time.time()

        with open(output_file, 'w') as output:
            all_res_files = [open(res_file) for res_file in tmp_res_files]

            # titles
            for i in range(5):
                output.write(all_res_files[0].readline())

            after_title_pos = all_res_files[0].tell()

            for res_file in all_res_files:
                res_file.seek(after_title_pos)
                output.write(res_file.read())
                res_file.close()

        merge_stop = time.time()
        print("Time cost of merge_result = {} sec".format(merge_stop - merge_start), file=sys.stderr)

        # clean tmp files and dir
        shutil.rmtree(tmp_dir)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("fasta")
    parser.add_argument("output")
    parser.add_argument("-p", "--process", type=int, default=1)
    parser.add_argument("--tmp_path", type=str, default=".")
    parser.add_argument("--blat_bin", default="blat")
    parser.add_argument("--blat_opt", type=str, default="")
    args = parser.parse_args()

    start = time.time()
    mp_blat(args.reference, args.fasta, args.output, args.process, args.tmp_path, args.blat_bin, args.blat_opt)
    stop = time.time()
    print("Total time cost = {} sec".format(stop - start), file=sys.stderr)


if __name__ == "__main__":
    cli()
