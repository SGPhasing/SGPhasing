# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.processor minimap2.

Functions:
  - splice_mapper
  - genomic_mapper
"""

from subprocess import PIPE, Popen, STDOUT


def splice_mapper(reference: str,
                  input_fastx: str,
                  output_sam: str,
                  threads: int = 1) -> str:
    """Call minimap2 for splice mapping.

    Args:
        reference (str): input reference fasta file path string.
        input_fastx (str): input sequences fasta/q file path string.
        output_sam (str): output sam file path string.
        threads (int): threads using for minimap2, default = 1.

    Returns:
        (str): minimap2 print.
    """
    with Popen(['minimap2', '-k', '17', '-uf', '-a', '-o', output_sam,
                '--MD', '-t', str(threads), '-x', 'splice:hq',
                reference+'.mmi', input_fastx],
               stdout=PIPE, stderr=STDOUT) as proc:
        return proc.stdout.read().decode('utf-8')


def genomic_mapper(reference: str,
                   input_fastx: str,
                   output_sam: str,
                   preset: str,
                   threads: int = 1) -> str:
    """Call minimap2 for genomic mapping.

    Args:
        reference (str): input reference fasta file path string.
        input_fastx (str): input sequences fasta/q file path string.
        output_sam (str): output sam file path string.
        preset (str): minimap2 preset.
        threads (int): threads using for minimap2, default = 1.

    Returns:
        (str): minimap2 print.
    """
    with Popen(['minimap2', '-k', '17', '-a', '-o', output_sam,
                '--MD', '-t', str(threads), '-x', preset,
                reference, input_fastx], stdout=PIPE, stderr=STDOUT) as proc:
        return proc.stdout.read().decode('utf-8')
