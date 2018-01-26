#!/bin/sh

pypy3 -m bartseq \
    -b 'data/ngs15/barcodes/barcodes.txt' \
    -t 5466115 \
    -s 'data/ngs15/stats.json' \
    --in-2 'data/ngs15/rawdata/NGS15_SC_Rep1_S1_R2_001.fastq.gz' \
    --out-2 'data/ngs15/tagged/NGS15_SC_Rep1_R2.fastq.gz' \
    'data/ngs15/rawdata/NGS15_SC_Rep1_S1_R1_001.fastq.gz' \
    'data/ngs15/tagged/NGS15_SC_Rep1_R1.fastq.gz'
