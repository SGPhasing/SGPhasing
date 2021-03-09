# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.processor index each linked region.

Functions:
  - index_each_link
"""

from SGPhasing.processor.bam_to_matrix import bam_to_matrix
from SGPhasing.processor.gatk4 import create_sequence_dictionary
from SGPhasing.processor.gatk4 import haplotype_caller
from SGPhasing.processor.gff_to_fasta import gff_to_fasta
from SGPhasing.processor.region_to_bam import fastx_to_bam, region_to_bam
from SGPhasing.reader.read_fastx import faidx
from SGPhasing.reader.read_vcf import get_alt_positions
from SGPhasing.threader.thread_haplotypes import onehot_decoder
from SGPhasing.threader.thread_haplotypes import thread_haplotypes


def index_each_link(args_tuple: tuple) -> tuple:
    """Index each linked region.

    Args:
        link_id (str): linked_region id.
        linked_region (Linked_Region): linked region.
        tmp_dir (Path): temporary folder Path.
        reference (str): reference fasta file path string.
        index_xam (str): input high quality bam/cram file for indexing.
        index_fastx (str): input high quality fasta/q file for indexing.
        threads (int): threads using for minimap2, pysam and gatk4, default 1.

    Returns:
        link_reads_id_main_dict (dict): region id as key and
                                        [chr, start, end, matrix] as value.
        positions_list (list): position list for each base.
    """
    (link_id, linked_region, tmp_dir,
     reference, index_xam, index_fastx, threads) = args_tuple
    link_floder_path = tmp_dir / link_id
    if not link_floder_path.is_dir():
        link_floder_path.mkdir()
    opened_gffread_log = (link_floder_path / 'gffread.log').open('w')
    opened_minimap2_log = (link_floder_path / 'minimap2.log').open('w')
    opened_gatk4_log = (link_floder_path / 'gatk4.log').open('w')

    link_gff_path = (
        link_floder_path / 'linked_region.minimap2_reference.gff3')
    linked_region.write_gff(str(link_gff_path))
    link_fasta_path = link_floder_path / 'linked_region.reference.fasta'
    opened_gffread_log.write(gff_to_fasta(
        str(link_gff_path), reference, str(link_fasta_path)))

    expanded_primary_region = linked_region.extend_primary()
    expanded_primary_region.update_info_id(link_id)
    opened_expand_gff = (
        link_floder_path/'expanded_primary.minimap2_reference.gff3').open('w')
    expanded_primary_region.write_gff(opened_expand_gff)
    opened_expand_gff.close()
    expand_fasta_path = link_floder_path / 'expanded_primary.reference.fasta'
    opened_gffread_log.write(gff_to_fasta(
        opened_expand_gff.name, reference, str(expand_fasta_path)))
    opened_gffread_log.close()

    faidx(str(expand_fasta_path))
    opened_gatk4_log.write(create_sequence_dictionary(str(expand_fasta_path)))

    link_expand_lalign_bam = fastx_to_bam(
        'reference', link_id, link_floder_path, str(link_fasta_path),
        str(expand_fasta_path), opened_minimap2_log, opened_gatk4_log, threads)
    index_expand_lalign_bam, reads_num = region_to_bam(
        'index', link_id, linked_region, link_floder_path,
        index_xam, index_fastx, str(expand_fasta_path),
        opened_minimap2_log, opened_gatk4_log, threads)
    opened_minimap2_log.close()

    ploidy = len(linked_region.Secondary_Regions_list) + 1
    index_expand_vcf_path = (
        link_floder_path / 'linked_region.index.minimap2_expand.hapcal.vcf')
    opened_gatk4_log.write(haplotype_caller(
        index_expand_lalign_bam, str(index_expand_vcf_path),
        str(expand_fasta_path), reads_num, 20, ploidy, threads))
    opened_gatk4_log.close()

    positions_list = get_alt_positions(str(index_expand_vcf_path), ploidy)
    if positions_list:
        link_reads_id_list, link_reads_bases_matrix = bam_to_matrix(
            'reference', link_floder_path,
            positions_list, link_expand_lalign_bam)
        index_reads_id_list, index_reads_bases_matrix = bam_to_matrix(
            'index', link_floder_path,
            positions_list, index_expand_lalign_bam)

        (clusters_indexes_array, sample_cluster_indexes,
         new_prototypes_array) = thread_haplotypes(
            link_reads_bases_matrix, index_reads_bases_matrix,
            len(positions_list))

        link_reads_id_main_dict = {}
        link_reads_id_main_dict.update(
            linked_region.Primary_Region.get_main_dict())
        for region in linked_region.Secondary_Regions_list:
            link_reads_id_main_dict.update(region.get_main_dict())

        for read_id, prototype_bases_matrix in zip(
                link_reads_id_list,
                onehot_decoder(new_prototypes_array, False)):
            link_reads_id_main_dict[read_id].append(prototype_bases_matrix)
        return link_reads_id_main_dict, positions_list
    else:
        return {}, []
