# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.writer write fasta/q file.

Functions:
  - write_partial_fastx
"""

from io import TextIOWrapper


def write_partial_fastx(opened_input_fastx: TextIOWrapper,
                        opened_output_fastx: TextIOWrapper,
                        fastx_format: str,
                        limit_reads_set: set) -> None:
    """Write fastx for limit reads.

    Args:
        opened_input_fastx (TextIOWrapper): opened input fastx handle.
        opened_output_fastx (TextIOWrapper): opened output fastx handle.
        fastx_format (str): input fasta/q file format.
        limit_reads_set (set): limited reads id set.
    """
    if fastx_format == 'fastq':
        unit = 4
    else:
        unit = 2
    read_unit = []
    line_id, is_write = 0, False
    for eachline in opened_input_fastx:
        line_id += 1
        if line_id % unit == 1:
            if is_write:
                opened_output_fastx.write(''.join(read_unit))
            read_id = eachline.strip().split()[0][1:]
            is_write = True if read_id in limit_reads_set else False
            read_unit = [eachline]
        else:
            read_unit.append(eachline)
    opened_output_fastx.close()
