# -*- coding: utf-8 -*-
"""setup script for SGPhasing.

This uses setuptools which is now the standard python mechanism for
installing packages. If you have downloaded and uncompressed the
SGPhasing source code, or fetched it from git, for the simplest
installation just type the command:

    python setup.py install

However, you would normally install the latest SGPhasing release from
the PyPI archive with:

    pip install SGPhasing

Or, you can also install it via conda's bioconda channel:

    conda install SGPhasing --channel bioconda

For more in-depth instructions, see the Install section of the
SGPhasing manual, linked to from:

    https://github.com/SGPhasing/SGPhasing/wiki

If all else fails, feel free to post on GitHub Issue:

    https://github.com/SGPhasing/SGPhasing/issues

or contact Shang Xie and ask for help:

    xieshang0608@gmail.com
"""

import ctypes
import locale
import os
import platform
import re
from subprocess import CalledProcessError, PIPE, Popen, run
import sys

from SGPhasing.sys_output import Output

INSTALL_FAILED = False


class Environment(object):
    """The current install environment.

    Attributes:
      - output: Output info, warning and error.
      - is_installer (bool): Whether to enter the installer mode.
      - missing_packages (list): Missing packages from current environment.
      - conda_missing_packages (list): Missing packages from current
        conda environment.

    """

    def __init__(self) -> None:
        """Initialize an Environment on current environment."""
        self.output = Output()
        self.is_installer = False
        self.conda_required_packages = [('bcbio-gff', 'bioconda'),
                                        ('cdna_cupcake', 'bioconda'),
                                        ('colorlover', 'plotly'),
                                        ('minimap2', 'bioconda'),
                                        ('numpy',),
                                        ('plotly', 'plotly'),
                                        ('pysam', 'bioconda'),
                                        ('samtools', 'bioconda'),
                                        ('seqkit', 'bioconda')]
        self.missing_packages = []
        self.conda_missing_packages = []

        self.process_arguments()
        self.check_permission()
        self.check_system()
        self.check_python()
        self.output_runtime_info()
        self.check_pip()
        self.upgrade_pip()
        self.get_installed_packages()
        self.get_installed_conda_packages()
        self.get_required_packages()

    @property
    def is_admin(self) -> bool:
        """Check whether user is admin."""
        try:
            retval = os.getuid() == 0
        except AttributeError:
            retval = ctypes.windll.shell32.IsUserAnAdmin() != 0
        return retval

    @property
    def os_version(self) -> tuple:
        """Get OS Verion."""
        return platform.system(), platform.release()

    @property
    def py_version(self):
        """Get Python Verion."""
        return platform.python_version(), platform.architecture()[0]

    @property
    def is_conda(self) -> bool:
        """Check whether using Conda."""
        return bool('conda' in sys.version.lower())

    @property
    def is_virtualenv(self) -> bool:
        """Check whether this is a virtual environment."""
        if not self.is_conda:
            retval = (hasattr(sys, 'real_prefix') or
                      (hasattr(sys, 'base_prefix') and
                       sys.base_prefix != sys.prefix))
        else:
            prefix = os.path.dirname(sys.prefix)
            retval = (os.path.basename(prefix) == 'envs')
        return retval

    @property
    def encoding(self):
        """Get system encoding."""
        return locale.getpreferredencoding()

    def process_arguments(self) -> None:
        """Process any cli arguments."""
        argv = [arg for arg in sys.argv]
        for arg in argv:
            if arg == 'install':
                self.is_installer = True

    def check_permission(self) -> None:
        """Check for Admin permissions."""
        if self.is_admin:
            self.output.info('Running as Root/Admin')
        else:
            self.output.warning('Running without root/admin privileges')

    def check_system(self) -> None:
        """Check the system."""
        self.output.info('The tool provides tips for installation\n'
                         'and installs required python packages')
        self.output.info(f'Setup in {self.os_version[0]} {self.os_version[1]}')
        if not self.os_version[0] in ['Windows', 'Linux', 'Darwin']:
            self.output.error(
                f'Your system {self.os_version[0]} is not supported!')
            sys.exit(1)

    def check_python(self) -> None:
        """Check python and virtual environment status."""
        self.output.info(
            f'Installed Python: {self.py_version[0]} {self.py_version[1]}')
        if not (self.py_version[0].split('.')[0] == '3' and
                self.py_version[0].split('.')[1] in ('6', '7', '8', '9')):
            self.output.error('Please run this script with Python version '
                              '3.6, 3.7, 3.8 or 3.9 and try again.')
            sys.exit(1)

    def output_runtime_info(self) -> None:
        """Output runtime info."""
        if self.is_conda:
            self.output.info('Running in Conda')
        if self.is_virtualenv:
            self.output.info('Running in a Virtual Environment')
        self.output.info(f'Encoding: {self.encoding}')

    def check_pip(self) -> None:
        """Check installed pip version."""
        try:
            import pip  # noqa pylint:disable=unused-import
        except ImportError:
            self.output.error(
                'Import pip failed. Please Install python3-pip and try again')
            sys.exit(1)

    def upgrade_pip(self) -> None:
        """Upgrade pip to latest version."""
        if not self.is_conda:
            # Don't do this with Conda, as we must use conda's pip
            self.output.info('Upgrading pip...')
            pipexe = [sys.executable, '-m', 'pip']
            pipexe.extend(['install', '--no-cache-dir', '-qq', '--upgrade'])
            if not self.is_admin and not self.is_virtualenv:
                pipexe.append('--user')
            pipexe.append('pip')
            run(pipexe)
        import pip
        self.output.info(f'Installed pip: {pip.__version__}')

    def get_installed_packages(self) -> None:
        """Get currently installed packages."""
        self.installed_packages = {}
        chk = Popen(f'"{sys.executable}" -m pip freeze',
                    shell=True, stdout=PIPE)
        installed = chk.communicate()[0].decode(self.encoding).splitlines()

        for pkg in installed:
            if '==' not in pkg:
                continue
            item = pkg.split('==')
            self.installed_packages.update({item[0]: item[1]})

    def get_installed_conda_packages(self) -> None:
        """Get currently installed conda packages."""
        if not self.is_conda:
            return
        self.installed_conda_packages = {}
        chk = os.popen('conda list').read()
        installed = [re.sub(' +', ' ', line.strip())
                     for line in chk.splitlines() if not line.startswith('#')]
        for pkg in installed:
            item = pkg.split(' ')
            self.installed_conda_packages.update({item[0]: item[1]})

    def get_required_packages(self) -> None:
        """Load requirements list."""
        self.required_packages = []
        pypath = os.path.dirname(os.path.realpath(__file__))
        requirements_file = os.path.join(pypath, 'requirements.txt')
        with open(requirements_file) as req:
            for package in req.readlines():
                package = package.strip()
                if package and (not package.startswith('#')):
                    self.required_packages.append(package)


class Install(object):
    """Install the requirements.

    Attributes:
      - output: Output info, warning and error.
      - env: Environment.

    """

    def __init__(self, environment) -> None:
        """Initialize an Install.

        Args:
          - environment: Store an environment.

        """
        self.output = Output()
        self.env = environment

        if not self.env.is_installer:
            self.ask_continue()

        self.install_missing_dep()
        self.output.info('All python3 dependencies are met.\r\n'
                         'You are good to go.\r\n\r\n'
                         'Enter:  python sgphasing.py -h to see the options')

    def ask_continue(self) -> None:
        """Ask Continue with Install."""
        inp = input(
            'Please ensure your System Dependencies are met. Continue? [y/N] ')
        if inp in ('', 'N', 'n'):
            self.output.error('Please install system dependencies to continue')
            sys.exit(1)

    def check_missing_dep(self) -> None:
        """Check for missing dependencies."""
        for pkg in self.env.required_packages:
            pkgs = pkg.split('==')
            key = pkgs[0]
            if key not in self.env.installed_packages:
                self.env.missing_packages.append(pkg)
                continue
            else:
                if len(pkgs) > 1:
                    if pkgs[1] != self.env.installed_packages.get(key):
                        self.env.missing_packages.append(pkg)
                        continue

    def check_conda_missing_dep(self) -> None:
        """Check for conda missing dependencies."""
        if not self.env.is_conda:
            return
        for pkg in self.env.conda_required_packages:
            pkgs = pkg.split('==')
            key = pkgs[0]
            if key not in self.env.installed_packages:
                self.env.conda_missing_packages.append(pkg)
                continue
            else:
                if len(pkgs) > 1:
                    if pkgs[1] != self.env.installed_conda_packages.get(key):
                        self.env.conda_missing_packages.append(pkg)
                        continue

    def install_missing_dep(self) -> None:
        """Install missing dependencies."""
        # Install conda packages first
        self.check_conda_missing_dep()
        if self.env.conda_missing_packages:
            self.install_conda_packages()
        self.check_missing_dep()
        if self.env.missing_packages:
            self.install_python_packages()

    def install_conda_packages(self) -> None:
        """Install required conda packages."""
        self.output.info(
            'Installing Required Conda Packages. This may take some time...')
        for pkg in self.env.conda_missing_packages:
            channel = None if len(pkg) != 2 else pkg[1]
            self.conda_installer(pkg[0], channel=channel, conda_only=True)

    def conda_installer(self, package: str,
                        channel: bool = None,
                        verbose: bool = False,
                        conda_only: bool = False) -> bool:
        """Install a conda package."""
        success = True
        condaexe = ['conda', 'install', '-y']
        if not verbose:
            condaexe.append('-q')
        if channel:
            condaexe.extend(['-c', channel])
        condaexe.append(package)
        self.output.info(f'Installing {package}')
        try:
            if verbose:
                run(condaexe, check=True)
            else:
                with open(os.devnull, 'w') as devnull:
                    run(condaexe, stdout=devnull, stderr=devnull, check=True)
        except CalledProcessError:
            if not conda_only:
                self.output.info(
                    f'Couldn\'t install {package} with Conda. Trying pip')
            else:
                self.output.warning(f'Couldn\'t install {package} with Conda. '
                                    'Please install this package manually')
            success = False
        return success

    def install_python_packages(self, verbose: bool = False) -> None:
        """Install required pip packages."""
        self.output.info(
            'Installing Required Python Packages. This may take some time...')
        for pkg in self.env.missing_packages:
            self.pip_installer(pkg)

    def pip_installer(self, package: str) -> None:
        """Install a pip package."""
        pipexe = [sys.executable, '-m', 'pip']
        # hide info/warning and fix cache hang
        pipexe.extend(['install', '-qq', '--no-cache-dir'])
        # install as user to solve perm restriction
        if not self.env.is_admin and not self.env.is_virtualenv:
            pipexe.append('--user')
        msg = f'Installing {package}'
        self.output.info(msg)
        pipexe.append(package)
        try:
            run(pipexe, check=True)
        except CalledProcessError:
            self.output.warning(f'Couldn\'t install {package} with pip. '
                                'Please install this package manually')


def main() -> None:
    """Create an environment and install missing packages."""
    ENV = Environment()
    Install(ENV)


if __name__ == '__main__':
    main()
