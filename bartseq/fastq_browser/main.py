import re
from itertools import cycle
from pathlib import Path
from typing import TextIO, Union

from tqdm import tqdm

from ..io import transparent_open, iter_fq


fields = [
	'barcode', 'linker', 'multi-bc', 'just-primer',
	'other-bcs', 'barcode-mismatch', 'junk',
]
RE_FIELDS = re.compile(f'(?P<field>{"|".join(fields)})=(?P<value>[^ ]*)')


def interleave(*iterables):
	unexhausted = set(range(len(iterables)))
	for i, it in cycle(enumerate(map(iter, iterables))):
		if not unexhausted:
			break
		try:
			yield i, next(it)
		except StopIteration:
			unexhausted -= {i}


def main(data_dir: Path, library: str, out: Union[Path, str, TextIO], out_compression: str):
	dir_process = data_dir / 'process'
	dir_tagged = dir_process / '3-tagged'
	dir_mapped = dir_process / '4-mapped'
	
	paths_fsq = [dir_tagged / f'{library}_R{r}.fastq.gz' for r in [1, 2]]
	paths_map = [dir_mapped / f'{library}_R{r}.tsv' for r in [1, 2]]
	path_count = dir_process / '1-index' / f'{library}.count.txt'
	
	with path_count.open() as c_f:
		n_pairs = int(c_f.read())
	
	with \
			transparent_open(paths_fsq[0]) as fsq_r1, \
			transparent_open(paths_fsq[1]) as fsq_r2, \
			transparent_open(paths_map[0]) as map_r1, \
			transparent_open(paths_map[1]) as map_r2, \
			transparent_open(out, 'wt', suffix=out_compression, ensure_parentdir=True) as f_out:
		print(
			'read',
			'header', 'read_seq', 'quality_seq',
			'amplicon', 'match_len',
			'barcode', 'linker', 'has_multiple_bcs',
			'is_just_primer', 'other_bcs', 'has_bc_mismatch', 'junk',
			sep='\t', file=f_out,
		)
		
		for read_side, ((header, read, qual), mapping) in tqdm(interleave(
			zip(iter_fq(fsq_r1), map_r1),
			zip(iter_fq(fsq_r2), map_r2),
		), total=2 * n_pairs):
			fields = list(RE_FIELDS.finditer(header))
			amp, match = mapping.strip().split('\t')
			
			print(
				read_side + 1,
				header[:fields[0].start() - 1],
				read, qual,
				amp, len(match),
				*[v if v != 'None' else '' for v in (m.group('value') for m in fields)],
				sep='\t', file=f_out,
			)
