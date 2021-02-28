# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.processor process each linked region.

Functions:
  - process_each_link
"""

from SGPhasing.processor.gff_to_fasta import gff_to_fasta
from SGPhasing.processor.minimap2 import genomic_mapper
from SGPhasing.reader import read_fastx, read_xam
from SGPhasing.writer import write_fastx


def process_each_link(args_tuple: tuple):
    """Process each linked region.

    Args:
        link_id (str): linked_region id.
        linked_region (Linked_Region): linked_region.
        tmp_dir (PosixPath): temporary folder PosixPath.
        reference (str): input reference fasta path string.
        main_xam (str): input high quality bam/cram file.
        main_fastx (str): input high quality fasta/q file.
        threads (int):  threads using for minimap2, default = 1.
    """
    (link_id, linked_region, tmp_dir,
     reference, main_xam, main_fastx, threads) = args_tuple
    link_floder_path = tmp_dir / link_id
    if not link_floder_path.is_dir():
        link_floder_path.mkdir()
    link_gff_path = (
        link_floder_path / 'linked_region.minimap2_reference.gff3')
    linked_region.write_gff(str(link_gff_path))
    link_fasta_path = (
        link_floder_path / 'linked_region.reference.fasta')
    gff_to_fasta(str(link_gff_path), reference, str(link_fasta_path))

    expanded_primary_region = linked_region.extend_primary()
    expanded_primary_region.update_info_id(link_id)
    opened_expand_gff = (
        link_floder_path /
        'expanded_primary.minimap2_reference.gff3').open('w')
    expanded_primary_region.write_gff(opened_expand_gff)
    opened_expand_gff.close()
    expand_fasta_path = (
        link_floder_path / 'expanded_primary.reference.fasta')
    gff_to_fasta(opened_expand_gff.name, reference, str(expand_fasta_path))

    link_reads_set = set()
    opened_main_xam, main_xam_format = read_xam.open_xam(main_xam)
    for read in opened_main_xam.fetch(linked_region.Primary_Region.chrom,
                                      linked_region.Primary_Region.start,
                                      linked_region.Primary_Region.end):
        link_reads_set.add(read.query_name)
    for region in linked_region.Secondary_Regions_list:
        for read in opened_main_xam.fetch(
                region.chrom, region.start, region.end):
            link_reads_set.add(read.query_name)
    opened_main_xam.close()

    opened_main_fastx, main_fastx_format = read_fastx.open_fastx(main_fastx)
    main_fastx_path = (
        link_floder_path / ('linked_region.main.' + main_fastx_format))
    write_fastx.write_partial_fastx(
        opened_main_fastx, main_fastx_path.open('w'),
        main_fastx_format, link_reads_set)
    opened_main_fastx.close()

    link_expand_sam_path = (
        link_floder_path / 'linked_region.reference.minimap2_expand.sam')
    genomic_mapper(str(expand_fasta_path), str(link_fasta_path),
                   str(link_expand_sam_path), 'asm20', threads)
    main_expand_sam_path = (
        link_floder_path / 'linked_region.main.minimap2_expand.sam')
    genomic_mapper(str(expand_fasta_path), str(main_fastx_path),
                   str(main_expand_sam_path), 'asm20', threads)
