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


def splice_mapper(input_ref: str, input_fastx: str, threads: int):
    with Popen(['minimap2', '-k', '17', '-uf', '-a', '-t', str(threads),
                '-x', 'splice:hq', input_ref+'.mmi',
                input_fastx], stdout=PIPE) as proc:
        return proc.stdout.read().decode('utf-8').strip()


def genomic_mapper(input_ref: str,
                   input_fastx: str,
                   threads: int,
                   preset: str):
    with Popen(['minimap2', '-k', '17', '-a', '-t', str(threads),
                '-x', preset, input_ref+'.mmi',
                input_fastx], stdout=PIPE) as proc:
        return proc.stdout.read().decode('utf-8').strip()
