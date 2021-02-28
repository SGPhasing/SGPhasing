# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.processor minimap2."""

from subprocess import PIPE, Popen


def splice_mapper(input_ref: str,
                  input_fastx: str,
                  output_sam: str,
                  threads: int = 1) -> str:
    """Call minimap2 for splice mapping.

    Args:
        input_ref (str): input reference fasta file path str.
        input_fastx (str): input sequences fasta/q file path str.
        output_sam (str): output sam file path str.
        threads (int): threads using for minimap2.

    Returns:
        (str): minimap2 print.
    """
    with Popen(['minimap2', '-k', '17', '-uf', '-a', '-o', output_sam,
                '--MD', '-t', str(threads), '-x', 'splice:hq',
                input_ref+'.mmi', input_fastx], stdout=PIPE) as proc:
        return proc.stdout.read().decode('utf-8').strip()


def genomic_mapper(input_ref: str,
                   input_fastx: str,
                   output_sam: str,
                   preset: str,
                   threads: int = 1) -> str:
    """Call minimap2 for splice mapping.

    Args:
        input_ref (str): input reference fasta file path str.
        input_fastx (str): input sequences fasta/q file path str.
        output_sam (str): output sam file path str.
        preset (str): minimap2 preset.
        threads (int): threads using for minimap2.

    Returns:
        (str): minimap2 print.
    """
    with Popen(['minimap2', '-k', '17', '-a', '-o', output_sam,
                '--MD', '-t', str(threads), '-x', preset,
                input_ref, input_fastx], stdout=PIPE) as proc:
        return proc.stdout.read().decode('utf-8').strip()
