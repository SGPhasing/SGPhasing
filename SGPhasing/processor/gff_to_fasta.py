# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.processor convert gff3 to fasta.

Functions:
  - gff_to_fasta
"""

from subprocess import PIPE, Popen, STDOUT


def gff_to_fasta(input_gff: str,
                 reference: str,
                 output_fasta: str) -> str:
    """Convert gff3 to fasta.

    Args:
        input_gff (str): input gff3 file path string.
        reference (str): reference fasta file path string.
        output_fasta (str): output fasta file path string.

    Returns:
        (str): gffread print.
    """
    with Popen(['gffread', input_gff, '-g', reference,
                '-w', output_fasta], stdout=PIPE, stderr=STDOUT) as proc:
        return proc.stdout.read().decode('utf-8')
