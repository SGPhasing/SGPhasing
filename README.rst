SGPhasing
=========

.. class:: no-web no-pdf

    .. image:: https://github.com/SGPhasing/SGPhasing/blob/main/images/sgphasing_logo.svg

.. class:: no-web no-pdf

    |language| |version| |update|

.. contents::

.. section-numbering::

Description
-----------

**SGPhasing** is a python3 package for Phasing sequences of Similar Genes.

Getting Started
---------------

Requirements
~~~~~~~~~~~~

SGPhasing is a python3 package. To use SGPhasing, python version 3.7 or higher is required.

Installation
~~~~~~~~~~~~

.. code-block:: shell

    git clone https://github.com/SGPhasing/SGPhasing.git
    cd SGPhasing
    conda create -n sgphasing -y python=3.8
    python setup.py install

or

.. code-block:: shell

    git clone https://github.com/SGPhasing/SGPhasing.git
    cd SGPhasing
    conda env create -f SGPhasing.yml

or

.. code-block:: shell

    git clone https://github.com/SGPhasing/SGPhasing.git
    conda create -n sgphasing -y python=3.8 bcbio-gff cdna_cupcake gatk4 h5py mappy minimap2 numpy pysam rich

Running the tests
~~~~~~~~~~~~~~~~~

Usage
-----

Support
-------

Contributing
------------

Citation
--------

License
-------

Changelog
---------

.. |language| image:: https://img.shields.io/badge/language-python-blue.svg

.. |version| image:: https://img.shields.io/badge/version-v0.0.1a-green.svg

.. |update| image:: https://img.shields.io/badge/last%20updated-08%20Mar%202021-orange.svg