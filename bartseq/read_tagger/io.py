import json
from pathlib import Path
from typing import Union, Optional

from ..io import read_fasta, transparent_open


def write_bc_table(path_bc_file: Union[Path, str], path_bc_table: Union[Path, str]):
	"""This is mainly independent of the rest, so do simple duplicate work to be able to create this separately"""
	from . import get_tagger
	bc_table = get_tagger(read_fasta(path_bc_file)).get_barcode_table()
	Path(path_bc_table).parent.mkdir(parents=True, exist_ok=True)
	with transparent_open(path_bc_table, 'wt') as f_bc:
		f_bc.write(bc_table)


def write_stats(stats_file: Union[Path, str], n_reads: int, n_both_regular: Optional[int], stats1: dict, stats2: Optional[dict]=None):
	with transparent_open(stats_file, 'wt') as f_s:
		stats = dict(n_reads=n_reads, read1=stats1)
		if n_both_regular is not None:
			stats['n_both_regular'] = n_both_regular
		if stats2:
			stats['read2'] = stats2
		json.dump(stats, f_s, indent='\t')
