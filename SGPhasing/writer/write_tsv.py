# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.writer write tsv file.

Functions:
  - write_reads_bases_matrix
"""


def write_reads_bases_matrix(output_file: str,
                             positions_list: list,
                             reads_id_list: list,
                             reads_bases_matrix: list) -> None:
    """Write read base matrix in tsv file.

    Args:
        output_file (str): output file path string.
        positions_list (list): position list for each base.
        reads_id_list (list): read_id list for each read.
        reads_bases_matrix (list): bases matrix for each read.
    """
    with open(output_file, 'w') as opened_tsv:
        HEADER = [str(pos) for pos in sorted(positions_list)]
        HEADER.insert(0, '')
        opened_tsv.write('\t'.join(HEADER)+'\n')
        for read_index, read_id in enumerate(reads_id_list):
            outline = reads_bases_matrix[read_index][::]
            outline.insert(0, read_id)
            opened_tsv.write('\t'.join(outline)+'\n')
