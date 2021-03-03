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

from SGPhasing.indexer.thread import thread_haplotypes
from SGPhasing.processor.gatk4 import add_or_replace_read_groups
from SGPhasing.processor.gatk4 import create_sequence_dictionary
from SGPhasing.processor.gatk4 import haplotype_caller
from SGPhasing.processor.gatk4 import left_align_indels
from SGPhasing.processor.gff_to_fasta import gff_to_fasta
from SGPhasing.processor.minimap2 import genomic_mapper
from SGPhasing.reader import read_fastx, read_vcf, read_xam
from SGPhasing.writer import write_fastx, write_tsv, write_xam


def process_each_link(args_tuple: tuple) -> None:
    """Process each linked region.

    Args:
        link_id (str): linked_region id.
        linked_region (Linked_Region): linked_region.
        tmp_dir (PosixPath): temporary folder PosixPath.
        reference (str): input reference fasta file path string.
        index_xam (str): input high quality bam/cram file for indexing.
        index_fastx (str): input high quality fasta/q file for indexing.
        threads (int): threads using for minimap2, default = 1.
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
    link_fasta_path = (
        link_floder_path / 'linked_region.reference.fasta')
    opened_gffread_log.write(gff_to_fasta(
        str(link_gff_path), reference, str(link_fasta_path)))

    expanded_primary_region = linked_region.extend_primary()
    expanded_primary_region.update_info_id(link_id)
    opened_expand_gff = (
        link_floder_path /
        'expanded_primary.minimap2_reference.gff3').open('w')
    expanded_primary_region.write_gff(opened_expand_gff)
    opened_expand_gff.close()
    expand_fasta_path = (
        link_floder_path / 'expanded_primary.reference.fasta')
    opened_gffread_log.write(gff_to_fasta(
        opened_expand_gff.name, reference, str(expand_fasta_path)))
    opened_gffread_log.close()

    link_reads_set = set()
    opened_index_xam, index_xam_format = read_xam.open_xam(index_xam)
    for read in opened_index_xam.fetch(linked_region.Primary_Region.chrom,
                                       linked_region.Primary_Region.start,
                                       linked_region.Primary_Region.end):
        link_reads_set.add(read.query_name)
    for region in linked_region.Secondary_Regions_list:
        for read in opened_index_xam.fetch(
                region.chrom, region.start, region.end):
            link_reads_set.add(read.query_name)
    opened_index_xam.close()

    opened_index_fastx, index_fastx_format = read_fastx.open_fastx(index_fastx)
    index_fastx_path = (
        link_floder_path / ('linked_region.index.' + index_fastx_format))
    write_fastx.write_partial_fastx(
        opened_index_fastx, index_fastx_path.open('w'),
        index_fastx_format, link_reads_set)
    opened_index_fastx.close()

    link_expand_sam_path = (
        link_floder_path / 'linked_region.reference.minimap2_expand.sam')
    opened_minimap2_log.write(genomic_mapper(
        str(expand_fasta_path), str(link_fasta_path),
        str(link_expand_sam_path), 'asm20', threads))
    index_expand_sam_path = (
        link_floder_path / 'linked_region.index.minimap2_expand.sam')
    opened_minimap2_log.write(genomic_mapper(
        str(expand_fasta_path), str(index_fastx_path),
        str(index_expand_sam_path), 'asm20', threads))
    opened_minimap2_log.close()

    link_expand_bam_path = (
        link_floder_path / 'linked_region.reference.minimap2_expand.bam')
    write_xam.sam_to_bam(str(link_expand_sam_path), str(expand_fasta_path),
                         str(link_expand_bam_path), threads)
    index_expand_bam_path = (
        link_floder_path / 'linked_region.index.minimap2_expand.bam')
    write_xam.sam_to_bam(str(index_expand_sam_path), str(expand_fasta_path),
                         str(index_expand_bam_path), threads)

    read_fastx.faidx(str(expand_fasta_path))
    opened_gatk4_log.write(create_sequence_dictionary(str(expand_fasta_path)))
    link_expand_group_bam_path = (
        link_floder_path / 'linked_region.reference.minimap2_expand.group.bam')
    opened_gatk4_log.write(add_or_replace_read_groups(
        str(link_expand_bam_path), str(link_expand_group_bam_path),
        'reference', 'SGPhasing', link_id, '0'))
    index_expand_group_bam_path = (
        link_floder_path / 'linked_region.index.minimap2_expand.group.bam')
    opened_gatk4_log.write(add_or_replace_read_groups(
        str(index_expand_bam_path), str(index_expand_group_bam_path),
        'index', 'SGPhasing', link_id, '0'))

    link_expand_lalign_bam_path = (
        link_floder_path/'linked_region.reference.minimap2_expand.lalign.bam')
    opened_gatk4_log.write(left_align_indels(
        str(link_expand_group_bam_path), str(link_expand_lalign_bam_path),
        str(expand_fasta_path)))
    index_expand_lalign_bam_path = (
        link_floder_path / 'linked_region.index.minimap2_expand.lalign.bam')
    opened_gatk4_log.write(left_align_indels(
        str(index_expand_group_bam_path), str(index_expand_lalign_bam_path),
        str(expand_fasta_path)))

    ploidy = len(linked_region.Secondary_Regions_list) + 1
    index_expand_vcf_path = (
        link_floder_path / 'linked_region.index.minimap2_expand.hapcal.vcf')
    opened_gatk4_log.write(haplotype_caller(
        str(index_expand_lalign_bam_path), str(index_expand_vcf_path),
        str(expand_fasta_path), len(link_reads_set), 20, ploidy, threads))
    opened_gatk4_log.close()

    positions_list = read_vcf.get_alt_positions(
        str(index_expand_vcf_path), ploidy)
    if positions_list:
        link_reads_id_list, link_reads_bases_matrix = (
            read_xam.extract_read_matrix(
                str(link_expand_lalign_bam_path), positions_list))
        link_expand_reads_bases_matrix_path = (
            link_floder_path / 'linked_region.reference.reads_bases.tsv')
        write_tsv.write_reads_bases_matrix(
            str(link_expand_reads_bases_matrix_path), positions_list,
            link_reads_id_list, link_reads_bases_matrix)
        index_reads_id_list, index_reads_bases_matrix = (
            read_xam.extract_read_matrix(
                str(index_expand_lalign_bam_path), positions_list))
        index_expand_reads_bases_matrix_path = (
            link_floder_path / 'linked_region.index.reads_bases.tsv')
        write_tsv.write_reads_bases_matrix(
            str(index_expand_reads_bases_matrix_path), positions_list,
            index_reads_id_list, index_reads_bases_matrix)

        (clusters_indexes_array, sample_cluster_indexes,
         new_prototypes_array) = thread_haplotypes(
            link_reads_bases_matrix, index_reads_bases_matrix)
    else:
        return None
