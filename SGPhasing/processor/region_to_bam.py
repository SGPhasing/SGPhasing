# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.processor convert region to bam.

Functions:
  - region_to_bam
  - fastx_to_bam
"""

from io import TextIOWrapper
from pathlib import Path

from SGPhasing.processor.gatk4 import add_or_replace_read_groups
from SGPhasing.processor.gatk4 import left_align_indels
from SGPhasing.processor.minimap2 import genomic_mapper
from SGPhasing.reader.read_fastx import open_fastx
from SGPhasing.reader.read_xam import open_xam
from SGPhasing.Regions import Linked_Region
from SGPhasing.writer import write_fastx, write_xam


def region_to_bam(input_type: str,
                  link_id: str,
                  linked_region: Linked_Region,
                  link_floder_path: Path,
                  input_xam: str,
                  input_fastx: str,
                  reference: str,
                  opened_minimap2_log: TextIOWrapper,
                  opened_gatk4_log: TextIOWrapper,
                  threads: int = 1) -> tuple:
    """Extract reads and run fastx_to_bam.

    Args:
        input_type (str): input fasta/q type string,
                          in ('reference', 'index', 'phase').
        link_id (str): link_id string.
        linked_region (Linked_Region): linked region.
        link_floder_path (Path): linked region floder Path.
        input_xam (str): input bam/cram file path string.
        input_fastx (str): input fasta/q file path string.
        reference (str): reference fasta file path string.
        opened_minimap2_log (TextIOWrapper): opened minimap2 log file handle.
        opened_gatk4_log (TextIOWrapper): opened gatk4 log file handle.
        threads (int): threads using for minimap2 and pysam, default 1.

    Returns:
        input_expand_lalign_bam (str): input fasta/q left aligned bam
                                       file path string.
        reads_num (int): extracted reads number.
    """
    link_reads_set = set()
    opened_input_xam, input_xam_format = open_xam(input_xam)
    for read in opened_input_xam.fetch(linked_region.Primary_Region.chrom,
                                       linked_region.Primary_Region.start,
                                       linked_region.Primary_Region.end):
        link_reads_set.add(read.query_name)
    for region in linked_region.Secondary_Regions_list:
        for read in opened_input_xam.fetch(
                region.chrom, region.start, region.end):
            link_reads_set.add(read.query_name)
    opened_input_xam.close()

    opened_input_fastx, input_fastx_format = open_fastx(input_fastx)
    input_fastx_path = (
        link_floder_path /
        (f'linked_region.{input_type}.' + input_fastx_format))
    write_fastx.write_partial_fastx(
        opened_input_fastx, input_fastx_path.open('w'),
        input_fastx_format, link_reads_set)
    opened_input_fastx.close()

    input_expand_lalign_bam = fastx_to_bam(
        input_type, link_id, link_floder_path, str(input_fastx_path),
        reference, opened_minimap2_log, opened_gatk4_log, threads)
    return input_expand_lalign_bam, len(link_reads_set)


def fastx_to_bam(input_type: str,
                 link_id: str,
                 link_floder_path: Path,
                 input_fastx: str,
                 reference: str,
                 opened_minimap2_log: TextIOWrapper,
                 opened_gatk4_log: TextIOWrapper,
                 threads: int = 1) -> str:
    """Run minimap2 & gatk from fasta/q to left aligned bam.

    Args:
        input_type (str): input fasta/q type string,
                          in ('reference', 'index', 'phase').
        link_id (str): link_id string.
        link_floder_path (Path): linked region floder Path.
        input_fastx (str): input fasta/q file path string.
        reference (str): reference fasta file path string.
        opened_minimap2_log (TextIOWrapper): opened minimap2 log file handle.
        opened_gatk4_log (TextIOWrapper): opened gatk4 log file handle.
        threads (int): threads using for minimap2 and pysam, default 1.

    Returns:
        input_expand_lalign_bam (str): input fasta/q left aligned bam
                                       file path string.
    """
    if input_type in ('reference', 'index'):
        minimap2_preset = 'asm20'
    elif input_type == 'phase':
        minimap2_preset = 'map-ont'
    else:
        exit()
    input_expand_sam_path = (
        link_floder_path / f'linked_region.{input_type}.minimap2_expand.sam')
    opened_minimap2_log.write(genomic_mapper(
        reference, input_fastx, str(input_expand_sam_path),
        minimap2_preset, threads))

    input_expand_bam_path = (
        link_floder_path / f'linked_region.{input_type}.minimap2_expand.bam')
    write_xam.sam_to_bam(str(input_expand_sam_path), reference,
                         str(input_expand_bam_path), threads)

    input_expand_group_bam_path = (
        link_floder_path /
        f'linked_region.{input_type}.minimap2_expand.group.bam')
    opened_gatk4_log.write(add_or_replace_read_groups(
        str(input_expand_bam_path), str(input_expand_group_bam_path),
        input_type, 'SGPhasing', link_id, '0'))

    input_expand_lalign_bam_path = (
        link_floder_path /
        f'linked_region.{input_type}.minimap2_expand.lalign.bam')
    opened_gatk4_log.write(left_align_indels(
        str(input_expand_group_bam_path), str(input_expand_lalign_bam_path),
        reference))
    return str(input_expand_lalign_bam_path)
