# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.processor phase each linked region.

Functions:
  - phase_each_link
"""

from gc import collect

from SGPhasing.processor.bam_to_matrix import bam_to_matrix
from SGPhasing.processor.gatk4 import create_sequence_dictionary
from SGPhasing.processor.region_to_bam import region_to_bam
from SGPhasing.reader.read_fastx import faidx
from SGPhasing.Regions import Linked_Region, Region


def phase_each_link(args_tuple: tuple):
    """Phase each linked region.

    Args:
        link_id (str): linked_region id.
        positions_list (list): position list for each base.
        region_id_main_dict (dict): region id as key and
                                    [chr, start, end, matrix] as value.
        tmp_dir (Path): temporary folder Path.
        index_xam (str): input high quality bam/cram file for indexing.
        index_fastx (str): input high quality fasta/q file for indexing.
        threads (int): threads using for minimap2, pysam and gatk4, default 1.

    Returns:

    """
    (link_id, positions_list, region_id_main_dict, tmp_dir,
     index_xam, index_fastx, threads) = args_tuple
    link_floder_path = tmp_dir / link_id
    opened_minimap2_log = (link_floder_path / 'minimap2.log').open('a')
    opened_gatk4_log = (link_floder_path / 'gatk4.log').open('a')

    expand_fasta_path = link_floder_path / 'expanded_primary.reference.fasta'
    expand_fasta_fai_path = (
        link_floder_path / 'expanded_primary.reference.fasta.fai')
    if not expand_fasta_fai_path.exists():
        faidx(str(expand_fasta_path))
    expand_fasta_dict_path = (
        link_floder_path / 'expanded_primary.reference.dict')
    if not expand_fasta_dict_path.exists():
        opened_gatk4_log.write(
            create_sequence_dictionary(str(expand_fasta_path)))

    secondary_regions_list = []
    for region_id, region_main_list in region_id_main_dict.items():
        chrom, start, end, prototype_bases_matrix = region_main_list
        if region_id.split('.')[-1] == '0':
            primary_region = Region(chrom, start, end)
        else:
            secondary_regions_list.append(Region(chrom, start, end))
    linked_region = Linked_Region(primary_region)
    linked_region.update_secondary(secondary_regions_list)
    del primary_region, secondary_regions_list
    collect()

    phase_expand_lalign_bam, reads_num = region_to_bam(
        'phase', link_id, linked_region, link_floder_path,
        index_xam, index_fastx, str(expand_fasta_path),
        opened_minimap2_log, opened_gatk4_log, threads)
    phase_reads_id_list, phase_reads_bases_matrix = bam_to_matrix(
        'phase', link_floder_path, positions_list, phase_expand_lalign_bam)
