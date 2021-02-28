# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""xxxx.

What's here:

xxxx.
----------------------------------------------------------

Classes:
  - Region
  - Linked_Region
"""

from SGPhasing.sys_output import Output


class Region(object):
    """The Region class.

    Attributes:
      - chrom: .
      - start: .
      - end: .
      - strand: .
      - info: .
      - child_list: .
    """

    def __init__(self,
                 chrom: str,
                 start: int,
                 end: int,
                 strand: str,
                 info: str = '') -> None:
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.info = info
        self.child_list = []
        super().__init__()

    def copy(self):
        new_region = Region(self.chrom, self.start, self.end,
                            self.strand, self.info)
        new_region.update_child_list(self.child_list)
        return new_region

    def update_child_list(self, child_list: list) -> None:
        self.child_list = child_list
        if not self.check_child():
            self.child_list = []
            output = Output()
            output.warning('Region.child_list does not updated. '
                           'Please check your child_list input.')

    def check_child(self) -> bool:
        return all([each_child.chrom == self.chrom and
                    each_child.strand == self.strand and
                    self.start <= each_child.start < each_child.end <= self.end
                    for each_child in self.child_list])

    def update_info_id(self, new_id: str) -> None:
        self.info = update_info_str_id(self.info, new_id)
        for child_id, each_child in enumerate(self.child_list):
            self.child_list[child_id].info = update_info_str_id(
                each_child.info, new_id)

    def to_gff_string(self) -> str:
        return '\t'.join([
            self.chrom, 'sgphasing_tmp', 'mRNA' if self.child_list else 'exon',
            str(self.start), str(self.end), '.', self.strand, '.', self.info])

    def write_gff(self, opened_gff) -> None:
        opened_gff.write(self.to_gff_string() + '\n')
        for each_child in self.child_list:
            opened_gff.write(each_child.to_gff_string() + '\n')


class Linked_Region(object):
    """The Linked_Region class.

    Attributes:
      - Primary_Region: .
      - Secondary_Regions_list: .
    """

    def __init__(self, Primary_Region: Region) -> None:
        self.Primary_Region = Primary_Region
        self.Secondary_Regions_list = []
        super().__init__()

    def update_secondary(self, Secondary_Regions_list: list) -> None:
        self.Secondary_Regions_list = Secondary_Regions_list

    def append_secondary(self, Secondary_Region: Region) -> None:
        self.Secondary_Regions_list.append(Secondary_Region)

    def update_info_id(self, new_id) -> None:
        self.Primary_Region.update_info_id(new_id+'.0')
        for region_id in range(len(self.Secondary_Regions_list)):
            self.Secondary_Regions_list[region_id].update_info_id(
                new_id+'.'+str(region_id+1))

    def extend(self) -> list:
        new_region_list = self.Secondary_Regions_list[::]
        new_region_list.insert(0, self.Primary_Region)
        return new_region_list

    def expand_primary(self, length: int = 1000) -> Region:
        new_region = self.Primary_Region.copy()
        new_region.start -= length
        new_region.end += length
        return new_region

    def write_gff(self, output_gff: str) -> None:
        with open(output_gff, 'w') as opened_gff:
            self.Primary_Region.write_gff(opened_gff)
            for region in self.Secondary_Regions_list:
                region.write_gff(opened_gff)


def merge_two_linked_regions(Linked_Region1: Linked_Region,
                             Linked_Region2: Linked_Region,
                             threshold_coverage: float = 0.5) -> Linked_Region:
    region2_list = Linked_Region2.extend()
    region2_id_list = []
    for region2_id, region2 in enumerate(region2_list):
        for region1 in Linked_Region1.extend():
            if coverage_two_regions(region1, region2) > threshold_coverage:
                region2_id_list.append(region2_id)
    for region2_id in region2_id_list:
        Linked_Region1.append_secondary(region2_list[region2_id])
    return Linked_Region1


def check_two_linked_regions(Linked_Region1: Linked_Region,
                             Linked_Region2: Linked_Region,
                             threshold_coverage: float = 0.5) -> bool:
    return any(
        coverage > threshold_coverage for coverage in
        [coverage_two_regions(Linked_Region1.Primary_Region, secondary_region)
         for secondary_region in Linked_Region2.Secondary_Regions_list] +
        [coverage_two_regions(Linked_Region2.Primary_Region, secondary_region)
         for secondary_region in Linked_Region1.Secondary_Regions_list])


def coverage_two_regions(Region1: Region, Region2: Region) -> float:
    if Region1.chrom == Region2.chrom and Region1.strand == Region2.strand:
        gap_list = [Region1.end - Region2.start, Region2.end - Region1.start]
        if all([gap > 0 for gap in gap_list]):
            return min(gap_list)/max(gap_list)
        else:
            return 0
    else:
        return 0


def update_info_str_id(info_str: str, new_id: str) -> str:
    info_dict = {}
    info_sp = info_str.split(';')
    for each_info in info_sp:
        each_info_sp = each_info.split('=')
        info_dict.update({each_info_sp[0]: each_info_sp[1]})
    new_name = new_id + '.' + info_dict.get('Name', '').split('.')[-1]
    info_dict.update({'ID': new_name, 'Name': new_name, 'Parent': new_id})
    return ';'.join(
        [info_k+'='+info_v for info_k, info_v in info_dict.items()])
