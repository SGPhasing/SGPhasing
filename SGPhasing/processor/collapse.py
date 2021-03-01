# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.processor collapse isoforms by sam.

Functions:
  - collapse_isoforms_by_sam
  - get_most_supported_isoforms
"""

from pathlib import Path

from cupcake.tofu.utils import check_ids_unique
from cupcake.tofu.branch.branch_simple2 import BranchSimple


def collapse_isoforms_by_sam(input_sam: str,
                             input_fastx: str,
                             is_fq: bool,
                             tmp_dir: Path) -> tuple:
    """Collapse isoforms by sam.

    Args:
        input_sam (str): input sam file path string.
        input_fastx (str): input fasta/q file path string.
        is_fq (bool): input_fastx is in fastq format,
        tmp_dir (PosixPath): temporary folder PosixPath.

    Returns:
        gff_path (str): primary_reference.collapsed.gff file path string.
        most_iso_id_list (list): most supported isoforms id list.
    """
    opened_gff = (tmp_dir/'primary_reference.collapsed.gff').open('w')
    opened_group = (tmp_dir/'primary_reference.collapsed.group.txt').open('w')
    opened_ignore = (tmp_dir/'primary_reference.ignored_ids.txt').open('w')
    check_ids_unique(input_fastx, is_fq=is_fq)
    branch_simple = BranchSimple(
        input_fastx, cov_threshold=1, min_aln_coverage=0.85,
        min_aln_identity=0.85, is_fq=is_fq, max_5_diff=1000, max_3_diff=100)
    iter = branch_simple.iter_gmap_sam(input_sam, opened_ignore)
    for recs in iter:
        # recs is {'+': list of list of records, '-': list of list of records}
        for recs_values in recs.values():
            for recs_value in recs_values:
                if len(recs_value) > 0:
                    branch_simple.process_records(
                        recs_value, False, False,
                        opened_gff, opened_gff, opened_group)
    opened_gff.close()
    opened_group.close()
    opened_ignore.close()
    return opened_gff.name, get_most_supported_isoforms(opened_group.name)


def get_most_supported_isoforms(input_group: str) -> list:
    """Get most supported isoforms.

    Args:
        input_group (str): primary_reference.collapsed.group.txt
                           file path string.

    Returns:
        most_iso_id_list (list): most supported isoforms id list.
    """
    most_iso_id_list = []
    gene_iso_flncnum_dict = {}
    with open(input_group, 'r') as opened_group:
        for eachline in opened_group:
            full_iso_id, reads_id = eachline.strip().split()
            pb, gene_id, iso_id = full_iso_id.split('.')
            flnc_num = len(reads_id.split(','))
            gene_iso_flncnum_dict.setdefault(
                pb+'.'+gene_id, {}).update({iso_id: flnc_num})
    for gene_id, gene_dict in gene_iso_flncnum_dict.items():
        max_flnc_num, most_iso_id = 0, ''
        for iso_id, flnc_num in gene_dict.items():
            if flnc_num > max_flnc_num:
                max_flnc_num = flnc_num
                most_iso_id = iso_id
        most_iso_id_list.append(gene_id+'.'+most_iso_id)
    return most_iso_id_list
