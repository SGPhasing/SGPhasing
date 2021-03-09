# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.writer write index file.

Functions:
  - write_index
  - seq_to_array
"""

from pathlib import Path

from h5py import File
import numpy as np

from SGPhasing.reader.read_fastx import get_seq_dict


def write_index(output_index: str,
                link_id_list: list,
                process_link_returns: list,
                tmp_floder_path: Path) -> None:
    """Write hdf5 index file.

    Args:
        output_index (str): output hdf5 index file path string.
        link_id_list (list): linked regions id in list.
        process_link_returns (list): (region_id_main_dict, positions_list)
                                     in list.
        tmp_floder_path (Path): temporary folder Path.
    """
    opened_index = File(output_index, 'w')
    for link_id, link_tuple in zip(link_id_list, process_link_returns):
        region_id_main_dict, positions_list = link_tuple
        link_group = opened_index.create_group(link_id)
        link_fasta_path = (
            tmp_floder_path / link_id / 'expanded_primary.reference.fasta')
        link_id_seq_dict = get_seq_dict(str(link_fasta_path))
        link_group.create_dataset(
            'reference', data=seq_to_array(link_id_seq_dict.get(link_id, '')))
        link_group.create_dataset(
            'positions', data=np.array(positions_list, dtype=np.uint32))
        for region_id, region_main_list in region_id_main_dict.items():
            chrom, start, end, matrix = region_main_list
            region_group = link_group.create_group(region_id)
            region_group.attrs.create('chromosome', chrom)
            region_group.attrs.create('start', start)
            region_group.attrs.create('end', end)
            region_group.create_dataset('index',
                                        data=np.array(matrix, np.uint8))
    opened_index.close()


def seq_to_array(sequence: str) -> np.ndarray:
    """Convert sequence to numpy array.

    Args:
        sequence (str): sequence string.

    Returns:
        (ndarray): numpy array.
    """
    BASE_INT = {
        'A': 0, 'a': 0,
        'C': 1, 'c': 1,
        'G': 2, 'g': 2,
        'T': 3, 't': 3,
        'M': 4, 'm': 4,
        'R': 5, 'r': 5,
        'W': 6, 'w': 6,
        'S': 7, 's': 7,
        'Y': 8, 'y': 8,
        'K': 9, 'k': 9,
        'V': 10, 'v': 10,
        'H': 11, 'h': 11,
        'D': 12, 'd': 12,
        'B': 13, 'b': 13,
        'N': 14, 'n': 14}
    return np.array([BASE_INT[base] for base in sequence], dtype=np.uint8)
