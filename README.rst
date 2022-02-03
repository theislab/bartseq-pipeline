.. image:: https://zenodo.org/badge/139153315.svg
   :target: https://doi.org/10.5281/zenodo.3251773

BART-Seq pipeline
=================

This is a pipeline for BART-Seq implemented with Snakemake

The primer design code and an older version lives in `theislab/bartSeq <https://github.com/theislab/bartSeq>`_.

Entry point
-----------
- Delete old environment in case there is one

``conda deactivate
conda env remove -n bartseq_snakemake``

- Build new environment from yml file

``conda env create -f Environment.yml``

Activate the new environment

``conda activate [newenvironmentname]``

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

Command line interface
----------------------

``python -m bartseq tag [<options>] [in_1] [out_1]``

in_1
   Read1 file to read from. Supported compression: see --in-compression
out_1
   Read1 file to write to. Supported compression: see --out-compression

--in-2 IN_2                                    Read2 file to read from. Supported compression: see --in-compression
--out-2 OUT_2                                  Read2 file to write to. Supported compression: see --out-compression
--bc-file=BC_FILE, -b BC_FILE                  Barcode file in the format ``<ID> <Sequence>`` (with header)
--stats-file=STATS_FILE, -s STATS_FILE         File to write final stats to (in JSON format)
--bc-table=BC_TABLE, -B BC_TABLE               File name for the HTML table of barcode mismatches
--total=TOTAL, -t TOTAL                        Number of fastq records in file. “0” means no progressbar
--len-primer=LEN_PRIMER, -p LEN_PRIMER         Primer length for stats
--len-linker=LEN_LINKER, -l LEN_LINKER         Linker length to cut out
--in-compression=<gz|xz|bz2>, -i <gz|xz|bz2>   Specify compression if reading from stdin or a file with unusual suffix
--out-compression=<gz|xz|bz2>, -o <gz|xz|bz2>  Specify compression if writing to stdout or a file with unusual suffix
--dry-run, -n                                  Only print what would be done and exit

``python -m bartseq count [<options>] data_dir [library]``

data_dir
   Data directory to read from. Needs to have the directories “./process/{3-tagged,4-mapped}” filled.
library
   Library name. E.g. “Lib1_S1_L001” for input files named “Lib1_S1_L001_R{12}_001.fastq.gz”. Omittable if only one library exists.

--no-mismatch  Ignore barcodes with mismatches while counting.
--both         Print the count results for both to stdout. Default: Write to “./process/5-counts” instead
--one          Print the count results for one to stdout. Default: Write to “./process/5-counts” instead

``python -m bartseq browse [<options>] data_dir [library] [out]``

data_dir
   Data directory to read from. Needs to have the directories “./process/{3-tagged,4-mapped}” filled.
library
   Library name. E.g. “Lib1_S1_L001” for input files named “Lib1_S1_L001_R{12}_001.fastq.gz”. Omittable if only one library exists.
out
   FASTA file to write to. Supported compression: see --out-compression

--out-compression <gz|xz|bz2>, -o <gz|xz|bz2>  Specify compression if writing to stdout or a file with unusual suffix

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

Notes
~~~~~~~~~~~~~~~~~
For MacOS, please use Snakefile_for_MacOS.

Contributors
----------------------
This repository is forked from theislab/bartseq-pipeline.

Achim Kramer Lab provided the Environment.yml and Snakefile_for_MacOS.

- Bert Maier, Environment.yml 

- Merve Busra Duman, Snakefile_for_MacOS
