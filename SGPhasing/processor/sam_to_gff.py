# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.processor convert sam to gff3.

Functions:
  - sam_to_gff
  - sam_record_to_gff_record
"""

from sys import stderr

from BCBio import GFF as BCBio_GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from cupcake.io.BioReaders import GMAPSAMReader


def sam_to_gff(input_sam: str,
               opened_output_gff,
               input_fasta: str = '') -> None:
    """Convert sam to gff3.

    Args:
        input_sam (str): input sam file path string.
        opened_output_gff: opened output gff3 file handle.
        input_fasta (str): input reads fasta file.
    """
    read_id_len_dict = {}
    if input_fasta:
        read_id_len_dict = {
            read.id: len(read.seq)
            for read in SeqIO.parse(open(input_fasta), 'fasta')}
    records = [sam_record_to_gff_record(record, 'sgphasing_tmp')
               for record in GMAPSAMReader(input_sam, True,
                                           query_len_dict=read_id_len_dict)]
    BCBio_GFF.write([record for record in records if record is not None],
                    opened_output_gff)


def sam_record_to_gff_record(sam_record, source: str = 'sgphasing_tmp'):
    """Convert sam record to gff3 record.

    Args:
        sam_record: GMAPSAMRecord record.
        source (str): gff3 source.

    Returns:
        gff_record: SeqRecord ready to be written as GFF3.
    """
    if sam_record.sID == '*':
        print(f'Skipping {sam_record.qID} because unmapped.', file=stderr)
        return None
    strand = 1 if sam_record.flag.strand == '+' else -1
    matches_num = sam_record.num_mat_or_sub - sam_record.num_nonmatches
    mismatches_num = sam_record.num_nonmatches
    indels_num = sam_record.num_ins + sam_record.num_del
    coverage_str = '{0:.2f}'.format(
        sam_record.qCoverage*10**2) if sam_record.qCoverage else 'NA'
    identity_str = '{0:.2f}'.format(sam_record.identity*10**2)

    gene_qualifiers = {'source': source,
                       'ID': sam_record.qID,
                       'Name': sam_record.qID}  # for gene record
    mRNA_qualifiers = {'source': source,
                       'ID': sam_record.qID + '.mRNA',
                       'Name': sam_record.qID + '.mRNA',
                       'Parent': sam_record.qID,
                       'coverage': coverage_str,
                       'identity': identity_str,
                       'matches': matches_num,
                       'mismatches': mismatches_num,
                       'indels': indels_num}

    # gene line, one per record
    top_feature = SeqFeature(
        FeatureLocation(sam_record.sStart, sam_record.sEnd),
        type='gene', strand=strand, qualifiers=gene_qualifiers)
    # mRNA line, one per record
    top_feature.sub_features = [
        SeqFeature(
            FeatureLocation(sam_record.sStart, sam_record.sEnd),
            type="mRNA", strand=strand, qualifiers=mRNA_qualifiers)]

    # exon lines, as many exons per record
    for exon_id, exon in enumerate(sam_record.segments):
        exon_id_str = '{0}.exon{1}'.format(sam_record.qID, exon_id+1)
        exon_qual = {'source': source, 'ID': exon_id_str, 'Name': exon_id_str}
        top_feature.sub_features.append(
            SeqFeature(
                FeatureLocation(exon.start, exon.end),
                type='exon', strand=strand, qualifiers=exon_qual))

    length = sum(exon.end - exon.start for exon in sam_record.segments)
    # DO NOT CARE since sequence is not written in GFF3
    fake_sequence = Seq('A'*length)
    gff_record = SeqRecord(fake_sequence, sam_record.sID)
    gff_record.features = [top_feature]
    return gff_record
