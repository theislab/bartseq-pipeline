.. image:: https://zenodo.org/badge/139153315.svg
   :target: https://doi.org/10.5281/zenodo.3251773

BART-Seq pipeline
=================

This is a pipeline for BART-Seq implemented with Snakemake

The primer design code and an older version lives in `theislab/bartSeq <https://github.com/theislab/bartSeq>`_.

Entry point
-----------
The pipeline can be run via ``snakemake [-j 4] [-s …/bartseq/Snakefile] [-d …/mydata]``,
where ``-j`` specifies the number of threads,
and the other parameters default to ``./Snakefile`` and ``.``, respectively.

Within the data directory, the following structure is expected:

- ``in/``
  
  - ``reads/``
    
    - ``<libname>_R1_001.fastq.gz``
    - ``<libname>_R2_001.fastq.gz``
    - optionally additional libraries…
  
  - ``amplicons.fa`` or ``amplicons/<libname>.fa`` for all libraries
  - ``barcodes.fa`` or ``barcodes/<libname>.fa`` for all libraries

- ``config.yml`` – A file with the defaults

  .. code:: yaml
  
     amplicon-min-length: null  # You can set an integer like 70
     allow-mismatch:      True  # You can set this to False

Through the way Snakemake works, you need to create this file.
leave it empty to use the defaults.

Output
------
The pipeline creates a ``process`` and an ``out`` directory.

The ``out`` directory contains plots and summary spreadsheets.

``process/3-tagged`` contains the tagged FASTQ reads without alignment,
but you can add a tag for the mapped amplicon by executing e.g.
``python -m bartseq browse NGS16 Lib1_S1_L001 | gzip >Lib1.fq.gz``

Data and statistics
-------------------

Read structure
~~~~~~~~~~~~~~
- trash
- 3nt protection CCA ()
- 8nt barcode (known from set)
- Linker (one for left bcs, one for right bcs)
- Primer + Rest of Amplicon

Interesting Statistics
~~~~~~~~~~~~~~~~~~~~~~
Make statistics: How many reads have a barcode, ...

from reads tagged with info:

- Barcode available?
- Trash before bc?
- Where bc?
- Concatamere? (bc[-bc-bc…]-linker-primer)
- Which nucleotides where bc should be?
- Amplicon maps to which gene?

Possible Problems
~~~~~~~~~~~~~~~~~
- No Amplicons: Only bc and linker
- Amplicon quality bad at the end
- Trash at the beginning
- Barcodes can have mismatches
