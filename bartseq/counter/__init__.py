import json
import re
import collections
from pathlib import Path
from typing import Optional, Tuple, Counter, Generator, Iterator

from tqdm import tqdm

from ..io import transparent_open


bc_re = re.compile(r'barcode=(\w+)')
bc_mm_re = re.compile(r'barcode-mismatch=(True|False)')


def get_barcodes(fastq_file: Iterator[str], allow_mismatch: bool) -> Generator[Tuple[str, Optional[bool]], None, None]:
	for header in fastq_file:
		next(fastq_file)
		next(fastq_file)
		next(fastq_file)
		mm = None if allow_mismatch else (bc_mm_re.search(header).group(1) == 'True')
		yield bc_re.search(header).group(1), mm


def count(
	data_dir: Path,
	library: str,
	*,
	allow_mismatch: bool = True,
	total: Optional[int] = None,
	amp_min: Optional[int] = None,
) -> Tuple[Counter[Tuple[str, str, str]], Counter[Tuple[str, str, str]]]:
	reads = [f'{data_dir}/process/3-tagged/{library}_R{read}.fastq.gz' for read in [1, 2]]
	mappings = [f'{data_dir}/process/4-mapped/{library}_R{read}.tsv' for read in [1, 2]]
	
	if total is None:
		try:
			total = json.loads(Path(f'{data_dir}/process/3-tagged/{library}_stats.json').read_bytes())['n_both_regular']
		except Exception:
			pass
	
	counts_both = collections.Counter()
	counts_one = collections.Counter()
	with \
		transparent_open(reads[0]) as r1, open(mappings[0]) as a1, \
		transparent_open(reads[1]) as r2, open(mappings[1]) as a2:
		
		bcs1 = get_barcodes(r1, allow_mismatch)
		bcs2 = get_barcodes(r2, allow_mismatch)
		amps1 = (a.strip().split('\t') for a in a1)
		amps2 = (a.strip().split('\t') for a in a2)
		for (bc1, bc1mm), (bc2, bc2mm), (amp1, amp1s), (amp2, amp2s) in tqdm(zip(bcs1, bcs2, amps1, amps2), total=total):
			# If we don’t allow mismatches in barcodes, we skip this read pair
			if not allow_mismatch and bc1mm or bc2mm:
				continue
			# We don’t know which read is “the left one”, so e.g. (L3,R4) == (R4,L3)
			bc1, bc2 = sorted([bc1, bc2])
			if amp_min is not None:
				if len(amp1s) < amp_min: amp1 = '*'
				if len(amp2s) < amp_min: amp2 = '*'
			if amp1 == amp2:
				if amp1 == '*':
					counts_both[bc1, bc2, '-unmapped'] += 1
				else:
					counts_both[bc1, bc2, amp1] += 1
					counts_one[bc1, bc2, amp1] += 1
			elif amp1 == '*':
				counts_both[bc1, bc2, '-one-mapped'] += 1
				counts_one[bc1, bc2, amp2] += 1  # amp1 == '*'
			else:
				counts_both[bc1, bc2, '-mismatch'] += 1
	
	return counts_both, counts_one
