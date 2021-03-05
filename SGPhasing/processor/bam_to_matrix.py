# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.processor convert bam to matrix.

Functions:
  - bam_to_matrix
"""

from pathlib import Path

from SGPhasing.reader import read_xam
from SGPhasing.writer import write_tsv


def bam_to_matrix(input_type: str,
                  link_floder_path: Path,
                  positions_list: list,
                  expand_lalign_bam: str) -> list:
    """Extract reads base matrix and write into tsv.

    Args:
        input_type (str): input fasta/q type string,
                          in ('reference', 'index', 'phase').
        link_floder_path (Path): linked region floder Path.
        positions_list (list): position list for each base.
        expand_lalign_bam (str): mapped to expand reference and
                                 left aligned bam path string.

    Returns:
        reads_bases_matrix (list): extracted bases matrix for each read.
    """
    reads_id_list, reads_bases_matrix = (
        read_xam.extract_read_matrix(expand_lalign_bam, positions_list))
    reads_bases_matrix_path = (
        link_floder_path / f'linked_region.{input_type}.reads_bases.tsv')
    if input_type != 'reference':
        reads_id_list, reads_bases_matrix = read_xam.remove_blank(
            reads_id_list, reads_bases_matrix)
    write_tsv.write_reads_bases_matrix(
        str(reads_bases_matrix_path), positions_list,
        reads_id_list, reads_bases_matrix)
    return reads_bases_matrix
