# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.reader read gff file.

Functions:
  - read_gff
  - get_info_dict
"""

from gc import collect
from io import TextIOWrapper

from SGPhasing.Regions import Linked_Region, Region
from SGPhasing.Regions import check_two_linked_regions
from SGPhasing.Regions import merge_two_linked_regions


def read_gff(opened_gff: TextIOWrapper) -> dict:
    """Read gff file and extract region information.

    Args:
        opened_gff (TextIOWrapper): opened input gff handle.

    Returns:
        gene_id_linked_region (dict): gene_id as key and
                                      linked_region as value.
    """
    gene_id_linked_region = {}
    for eachline in opened_gff:
        if eachline[0] != '#':
            sp = eachline.strip().split()
            if sp[2] == 'mRNA':
                gene_id = get_info_dict(sp[8])['Parent']
                gene_region = Region(sp[0], int(sp[3]),
                                     int(sp[4]), sp[6], sp[8])
                if gene_id in gene_id_linked_region:
                    gene_id_linked_region[
                        gene_id].append_secondary(gene_region)
                else:
                    gene_id_linked_region.update({
                        gene_id: Linked_Region(gene_region)})
            elif sp[2] == 'exon':
                gene_id = get_info_dict(sp[8])['Parent']
                exon_region = Region(sp[0], int(sp[3]),
                                     int(sp[4]), sp[6], sp[8])
                if gene_id_linked_region[gene_id].Secondary_Regions_list:
                    gene_id_linked_region[gene_id].Secondary_Regions_list[
                        -1].child_list.append(exon_region)
                else:
                    gene_id_linked_region[
                        gene_id].Primary_Region.child_list.append(exon_region)
    for gene_id, linked_region in gene_id_linked_region.items():
        linked_region.Primary_Region.check_child()
        for secondary_region in linked_region.Secondary_Regions_list:
            secondary_region.check_child()
    merged_linked_region, repeat_id_set = {}, set()
    gene_id_list = list(gene_id_linked_region.keys())
    genes_num = len(gene_id_list)
    for index1, gene_id in enumerate(gene_id_list):
        for index2 in range(index1+1, genes_num):
            linked_region1 = gene_id_linked_region.get(gene_id_list[index1])
            linked_region2 = gene_id_linked_region.get(gene_id_list[index2])
            if check_two_linked_regions(linked_region1, linked_region2):
                if index1 in repeat_id_set:
                    merged_linked_region.update({
                        gene_id_list[index2]: merge_two_linked_regions(
                            linked_region2, linked_region1)})
                else:
                    merged_linked_region.update({
                        gene_id: merge_two_linked_regions(
                            linked_region1, linked_region2)})
                    repeat_id_set.add(index2)
    for id in repeat_id_set:
        del gene_id_linked_region[gene_id_list[id]]
    gene_id_linked_region.update(merged_linked_region)
    del merged_linked_region
    collect()
    return gene_id_linked_region


def get_info_dict(anno_info: str) -> dict:
    """Get information from gff3 9th column.

    Args:
        anno_info (str): gff3 9th column string.

    Returns:
        info_dict (dict): gff3 9th column information dictionary.
    """
    info_dict = {}
    info_sp = anno_info.split(';')
    for each_info in info_sp:
        each_info_sp = each_info.split('=')
        info_dict.update({each_info_sp[0]: each_info_sp[1]})
    return info_dict
