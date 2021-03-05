# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.indexer thread haplotypes using read_base_matrix.

Functions:
  - thread_haplotypes
  - onehot_encode
  - distance_cluster
  - update_prototype
  - hardmax
"""

import numpy as np


def thread_haplotypes(ref_reads_bases_matrix: list,
                      index_reads_bases_matrix: list,
                      bases_num: int) -> tuple:
    """Threads haplotypes through the clusters.

    Args:
        ref_reads_bases_matrix (list): reference bases matrix for each read.
        index_reads_bases_matrix (list): index bases matrix for each read.

    Returns:
        clusters_indexes_array (list): [indexes, array] list for each cluster.
        sample_cluster_indexes (list): cluster index for each sample.
        new_prototypes_array (list): prototype array for each cluster.
    """
    ref_reads_bases_indexes, ref_reads_bases_array = onehot_encode(
        ref_reads_bases_matrix, True)
    index_reads_bases_indexes, index_reads_bases_array = onehot_encode(
        index_reads_bases_matrix)

    clusters_indexes_array, sample_cluster_indexes = [], []
    new_prototypes_array = []
    prototypes_array = ref_reads_bases_array[::]
    if_converge = False
    while not if_converge:
        clusters_indexes_array, sample_cluster_indexes = distance_cluster(
            prototypes_array, index_reads_bases_array,
            index_reads_bases_indexes)
        new_prototypes_array = update_prototype(
            clusters_indexes_array, prototypes_array, bases_num)
        if_converge = np.all(
            np.array(new_prototypes_array) == np.array(prototypes_array))
        prototypes_array = new_prototypes_array[::]
    return clusters_indexes_array, sample_cluster_indexes, new_prototypes_array


def onehot_encode(reads_bases_matrix: list, padding: bool = False) -> tuple:
    """Convert reads bases to one-hot matrix.

    Args:
        reads_bases_matrix (list): bases matrix for each read.
        padding (pool): if padding null value bases, default False.

    Returns:
        reads_bases_indexes (list): bases index for each read.
        reads_bases_array (list): bases numpy array for each read.
    """
    ONE_HOT = {
        '-': [0., 0., 0., 0., 1.],
        'A': [1., 0., 0., 0., 0.],
        'C': [0., 1., 0., 0., 0.],
        'G': [0., 0., 1., 0., 0.],
        'T': [0., 0., 0., 1., 0.],
        'M': [0.5, 0.5, 0., 0., 0.],
        'R': [0.5, 0., 0.5, 0., 0.],
        'W': [0.5, 0., 0., 0.5, 0.],
        'S': [0., 0.5, 0.5, 0., 0.],
        'Y': [0., 0.5, 0., 0.5, 0.],
        'K': [0., 0., 0.5, 0.5, 0.],
        'V': [1/3, 1/3, 1/3, 0., 0.],
        'H': [1/3, 1/3, 0., 1/3, 0.],
        'D': [1/3, 0., 1/3, 1/3, 0.],
        'B': [0., 1/3, 1/3, 1/3, 0.],
        'N': [0.25, 0.25, 0.25, 0.25, 0.]}
    reads_bases_indexes, reads_bases_array = [], []
    for read_bases_matrix in reads_bases_matrix:
        base_indexes, base_array = [], []
        for base_id, base in enumerate(read_bases_matrix):
            if base:
                base_indexes.append(base_id)
                base_array.append(ONE_HOT[base])
            elif padding:
                base_indexes.append(base_id)
                base_array.append([0., 0., 0., 0., 0.])
        reads_bases_indexes.append(base_indexes)
        reads_bases_array.append(np.array(base_array, dtype=np.float16))
    return reads_bases_indexes, reads_bases_array


def distance_cluster(prototypes_array: list,
                     samples_array: list,
                     samples_indexes: list) -> tuple:
    """Calculate distance between every prototype and every sample.

    Args:
        prototypes_array (list): prototype array for each cluster.
        samples_array (list): read bases array for each sample.
        samples_indexes (list): read bases index for each sample.

    Returns:
        clusters_indexes_array (list): [indexes, array] list for each cluster.
        sample_cluster_indexes (list): cluster index for each sample.
    """
    clusters_indexes_array = [[[], []] for _ in range(len(prototypes_array))]
    sample_cluster_indexes = []
    for sample_indexes, sample_array in zip(samples_indexes, samples_array):
        distance_list = [
            np.linalg.norm(sample_array - prototype_array[sample_indexes])
            for prototype_array in prototypes_array]
        min_distance_index = np.argmin(distance_list)
        sample_cluster_indexes.append(min_distance_index)
        clusters_indexes_array[min_distance_index][0].append(sample_indexes)
        clusters_indexes_array[min_distance_index][1].append(sample_array)
    return clusters_indexes_array, sample_cluster_indexes


def update_prototype(clusters_indexes_array: list,
                     prototypes_array: list,
                     bases_num: int) -> list:
    """Calculate new prototype array for each cluster.

    Args:
        clusters_indexes_array (list): [indexes, array] list for each cluster.
        prototypes_array (list): prototype array for each cluster.
        bases_num (int): prototype array bases number.

    Returns:
        new_prototypes_array (list): new prototype array for each cluster.
    """
    new_prototypes_array = []
    for cluster_id, cluster_indexes_array in enumerate(clusters_indexes_array):
        bases_count = [0 for _ in range(bases_num)]
        new_prototype_array = np.zeros((bases_num, 5), dtype=np.float16)
        clusters_indexes, cluster_array = cluster_indexes_array
        for sample_indexes, sample_array in zip(
                clusters_indexes, cluster_array):
            for base_index, base_array in zip(sample_indexes, sample_array):
                bases_count[base_index] += 1
                new_prototype_array[base_index] += base_array
        for base_id, base_count, base_array in zip(
                range(bases_num), bases_count, new_prototype_array):
            if base_count:
                new_prototype_array[base_id] = base_array / base_count
            else:
                new_prototype_array[base_id] = prototypes_array[
                    cluster_id][base_id][::]
        new_prototypes_array.append(hardmax(new_prototype_array))
    return new_prototypes_array


def hardmax(input_array: np.ndarray) -> np.ndarray:
    """Hardmax activation function.

    Args:
        input_array (ndarray): numpy array in (m, n) shape.

    Returns:
        (ndarray): numpy array in (m, n) shape.
    """
    new_array = np.zeros(input_array.shape, dtype=np.float16)
    for array_index, each_array in enumerate(input_array):
        new_array[array_index][np.argmax(each_array)] = 1.
    return new_array
