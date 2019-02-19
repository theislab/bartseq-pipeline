import sys
from pathlib import Path
from typing import Optional, Counter, Tuple
from typing.io import TextIO

from . import count


def print_counter(counter: Counter[Tuple[str, str, str]], of: Optional[TextIO] = None):
	if of is None:
		of = sys.stdout
	print('bc_l', 'bc_r', 'amp', 'count', sep='\t', file=of)
	for fields, c in counter.items():
		print(*fields, c, sep='\t', file=of)


def main(
	data_dir: Path,
	library: str,
	*,
	allow_mismatch: bool = True,
	both: Optional[bool] = None,
	total: Optional[int] = None,
	amp_min: Optional[int] = None,
):
	counts_both, counts_one = count(data_dir, library, allow_mismatch=allow_mismatch, total=total, amp_min=amp_min)
	
	if both is None:
		for counter, counting in [(counts_both, 'both'), (counts_one, 'one')]:
			with open(f'process/5-counts/{counting}/{library}.tsv', 'w') as of:
				print_counter(counter, of)
	else:
		print_counter(counts_both if both else counts_one)
