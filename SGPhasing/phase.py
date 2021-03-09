# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
""".

What's here:

.
----------------------------------------------------------

Classes:
  - Phase
"""

from gc import collect
from logging import getLogger
from multiprocessing import Pool
from pathlib import Path

from SGPhasing.processor.phase_each_link import phase_each_link
from SGPhasing.reader.read_index import read_index
from SGPhasing.reader.read_xam import check_index
from SGPhasing.sys_output import Output

logger = getLogger(__name__)  # pylint: disable=invalid-name


class Phase(object):
    """The Phase process.

    Attributes:
        args: Arguments.
        output: Output info, warning and error.
    """

    def __init__(self, arguments) -> None:
        self.args = arguments
        self.output = Output()
        if not self.args.tmp:
            self.args.tmp = self.args.input + '.sgphasing.tmp'
        self.output.info(
            f'Initializing {self.__class__.__name__}: (args: {arguments}')
        logger.debug(
            f'Initializing {self.__class__.__name__}: (args: {arguments}')
        super().__init__()

    def prepare(self) -> None:
        """Check genome index."""
        self.tmp_floder_path = Path(self.args.tmp)
        if not self.tmp_floder_path.is_dir():
            self.tmp_floder_path.mkdir()
            self.output.info(f'Creating temporary folder at {self.args.tmp}')
        self.opened_log_file = (
            self.tmp_floder_path / 'sgphasing.log').open('a')
        self.args.input = check_index(self.args.input)

    def read_sgp_index(self) -> None:
        """Read SGPhasing index."""
        (self.link_id_list, self.positions_list_list,
         self.region_id_main_dict_list) = read_index(
            self.args.index, self.tmp_floder_path)

    def process_links(self) -> None:
        """Using multiply threads phase each linked region."""
        in_pool_threads = int(self.args.threads / len(self.link_id_list))
        in_pool_threads = in_pool_threads if in_pool_threads else 1
        out_pool_threads = int(self.args.threads / in_pool_threads)
        self.output.info(
            f'Calling {out_pool_threads} pools to phase each region')
        process_link_args = []
        for link_id, positions_list, region_id_main_dict in zip(
                self.link_id_list, self.positions_list_list,
                self.region_id_main_dict_list):
            process_link_args.append((
                link_id, region_id_main_dict, self.tmp_floder_path,
                self.args.input, self.args.fastx, in_pool_threads))
        del self.positions_list_list, self.region_id_main_dict_list
        collect()
        with Pool(processes=out_pool_threads) as pool:
            pool.map(phase_each_link, process_link_args)
        del process_link_args
        collect()

    def process(self) -> None:
        """Call the Phase object."""
        self.output.info('Starting Phase Process')
        logger.debug('Starting Phase Process')
        self.prepare()
        self.read_sgp_index()
        self.process_links()
        self.opened_log_file.close()
        self.output.info('Completed Phase Process')
        logger.debug('Completed Phase Process')
