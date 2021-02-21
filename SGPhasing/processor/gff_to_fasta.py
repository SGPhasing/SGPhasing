# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.processor convert gff3 to fasta."""

from subprocess import PIPE, Popen


def gff_to_fasta(input_gff: str,
                 input_ref_fasta: str,
                 output_fasta: str) -> str:
    with Popen(['gffread', input_gff, '-g', input_ref_fasta,
                '-w', output_fasta], stdout=PIPE) as proc:
        return proc.stdout.read().decode('utf-8').strip()
