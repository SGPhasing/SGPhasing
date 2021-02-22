# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.reader read fastx file."""

import gzip
from pathlib import Path
from sys import exit

import mappy as mp

from SGPhasing.sys_output import Output


def open_fastx(input_fastx: str) -> tuple:
    """Check input and open."""
    if input_fastx.endswith(('fastq', 'fq')):
        open_input_fastx = open(input_fastx, 'r')
        input_format = 'fastq'
    elif input_fastx.endswith(('fastq.gz', 'fq.gz')):
        open_input_fastx = gzip.open(input_fastx, 'rb')
        input_format = 'fastq'
    elif input_fastx.endswith(('fasta', 'fa')):
        open_input_fastx = open(input_fastx, 'r')
        input_format = 'fasta'
    elif input_fastx.endswith(('fasta.gz', 'fa.gz')):
        open_input_fastx = gzip.open(input_fastx, 'rb')
        input_format = 'fasta'
    else:
        output = Output()
        output.error('input error: input format must be '
                     'fastq, fasta, fq, fa, or gzipped file.')
        exit()
    return open_input_fastx, input_format


def check_index(input_ref: str, threads: int = 1) -> None:
    """Build index for minimap2 if not exists.

    Args:
        input_ref (str): input reference fasta file path str.
        threads (int): threads using for mappy to index.
    """
    index_file = input_ref + '.mmi'
    index_path = Path(index_file)
    if not index_path.exists():
        output = Output()
        output.info('preparing genome index for minimap2.')
        opened_aligner = mp.Aligner(input_ref, preset='splice:hq',
                                    k=17, best_n=100,
                                    n_threads=threads, fn_idx_out=index_file)
