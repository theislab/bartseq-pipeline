import re
from collections import Counter
from pathlib import Path
from typing import Optional, Tuple

from tqdm import tqdm

from ..io import transparent_open


bc_re = re.compile(r'barcode=(\w+)')


def get_barcodes(fastq_file):
	for header in fastq_file:
		next(fastq_file)
		next(fastq_file)
		next(fastq_file)
		yield bc_re.search(header).group(1)


def count(
	data_dir: Path,
	library: str,
	*,
	total: Optional[int] = None,
	amp_min: Optional[int] = None,
) -> Tuple[Counter, Counter]:
	reads = [f'{data_dir}/process/3-tagged/{library}_R{read}.fastq.gz' for read in [1, 2]]
	mappings = [f'{data_dir}/process/4-mapped/{library}_R{read}.tsv' for read in [1, 2]]
	
	counts = Counter()
	counts_one = Counter()
	with \
		transparent_open(reads[0]) as r1, open(mappings[0]) as a1, \
		transparent_open(reads[1]) as r2, open(mappings[1]) as a2:
		
		bcs1 = get_barcodes(r1)
		bcs2 = get_barcodes(r2)
		amps1 = (a.strip().split('\t') for a in a1)
		amps2 = (a.strip().split('\t') for a in a2)
		for bc1, bc2, (amp1, amp1s), (amp2, amp2s) in tqdm(zip(bcs1, bcs2, amps1, amps2), total=total):
			bc1, bc2 = sorted([bc1, bc2])
			if amp_min is not None:
				if len(amp1s) < amp_min: amp1 = '*'
				if len(amp2s) < amp_min: amp2 = '*'
			if amp1 == amp2:
				if amp1 == '*':
					counts[bc1, bc2, '-unmapped'] += 1
				else:
					counts[bc1, bc2, amp1] += 1
					counts_one[bc1, bc2, amp1] += 1
			elif amp1 == '*':
				counts[bc1, bc2, '-one-mapped'] += 1
				counts_one[bc1, bc2, amp2] += 1  # amp1 == '*'
			else:
				counts[bc1, bc2, '-mismatch'] += 1
	
	return counts, counts_one
