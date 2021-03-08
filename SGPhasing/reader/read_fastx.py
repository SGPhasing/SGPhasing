# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.reader read fasta/q file.

Functions:
  - open_fastx
  - get_seq_dict
  - check_index
  - faidx
"""

import gzip
from pathlib import Path
from sys import exit

import mappy as mp
import pysam

from SGPhasing.sys_output import Output


def open_fastx(input_fastx: str) -> tuple:
    """Check input and open.

    Args:
        input_fastx (str): input fasta/q file path string.

    Returns:
        opened_fastx (TextIOWrapper): opened input fastx handle.
        input_format (str): input file in fasta or fastq format.
    """
    if input_fastx.endswith(('fastq', 'fq')):
        opened_fastx = open(input_fastx, 'r')
        input_format = 'fastq'
    elif input_fastx.endswith(('fastq.gz', 'fq.gz')):
        opened_fastx = gzip.open(input_fastx, 'rb')
        input_format = 'fastq'
    elif input_fastx.endswith(('fasta', 'fa')):
        opened_fastx = open(input_fastx, 'r')
        input_format = 'fasta'
    elif input_fastx.endswith(('fasta.gz', 'fa.gz')):
        opened_fastx = gzip.open(input_fastx, 'rb')
        input_format = 'fasta'
    else:
        output = Output()
        output.error('Input error: input format must be '
                     'fastq, fasta, fq, fa, or gzipped file.')
        exit()
    return opened_fastx, input_format


def get_seq_dict(input_fastx: str) -> dict:
    """Get sequences dict from fasta/q file.

    Args:
        input_fastx (str): input fasta/q file path string.

    Returns:
        id_seq_dict (dict): id as key and sequence as value.
    """
    id_seq_dict = {}
    opened_fastx, input_format = open_fastx(input_fastx)
    seq_id, sequence, seq_line = '', '', True
    if input_format == 'fasta':
        identifier = '>'
    else:
        identifier = '@'
    for eachline in opened_fastx:
        if eachline[0] == identifier:
            if seq_id:
                id_seq_dict.update({seq_id: sequence})
            seq_id = eachline.strip().split()[0][1:]
            sequence = ''
            seq_line = True
        elif eachline[0] == '+':
            seq_line = False
        elif seq_line:
            sequence += eachline.strip()
    id_seq_dict.update({seq_id: sequence})
    return id_seq_dict


def check_index(reference: str, threads: int = 1) -> None:
    """Build index for minimap2 if not exists.

    Args:
        reference (str): reference fasta file path string.
        threads (int): threads using for mappy to index, default 1.
    """
    index_file = reference + '.mmi'
    index_path = Path(index_file)
    if not index_path.exists():
        output = Output()
        output.info('Preparing genome index for minimap2')
        opened_aligner = mp.Aligner(reference, preset='splice:hq',
                                    k=17, best_n=100,
                                    n_threads=threads, fn_idx_out=index_file)


def faidx(reference: str) -> None:
    """Build samtools fasta index."""
    pysam.faidx('--length', '70', reference)
