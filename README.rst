SGPhasing
=========

.. image:: https://github.com/SGPhasing/SGPhasing/blob/main/Doc/images/sgphasing_logo.svg

.. image:: https://img.shields.io/badge/language-python-blue.svg
.. image:: https://img.shields.io/badge/version-v0.0.1a-green.svg
.. image:: https://img.shields.io/badge/last%20updated-24%20Feb%202021-orange.svg

Table of contents
-----------------

Description
-----------

**SGPhasing** is a python3 package for Phasing sequences of Similar Genes.

Getting Started
---------------

Requirements
~~~~~~~~~~~~

SGPhasing is a python3 package. To use SGPhasing, python version 3.7 is required.

Installation
~~~~~~~~~~~~

.. code-block:: shell

    git clone https://github.com/SGPhasing/SGPhasing.git
    conda create -n sgphasing -y --no-default-packages python=3.8 bcbio-gff cdna_cupcake gatk4 mappy minimap2 numpy pysam samtools

or

.. code-block:: shell

    git clone https://github.com/SGPhasing/SGPhasing.git
    cd SGPhasing
    conda create -n sgphasing -y --no-default-packages python=3.8
    python setup.py

or

.. code-block:: shell

    git clone https://github.com/SGPhasing/SGPhasing.git
    cd SGPhasing
    conda env create -f SGPhasing.yml

Running the tests
~~~~~~~~~~~~~~~~~

.. code-block:: shell

    cd SGPhasing
    python test.py

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