from contextlib import contextmanager
from typing import Generator, Union, Iterable, Optional

from tqdm import tqdm

from . import defaults, ReadTagger, get_tagger
from .io import write_bc_tables, write_stats
from ..io import transparent_open, iter_fq, read_fasta
from ..logging import init_logging


def run(
	in_1: Union[str, Iterable[str]],
	out_1: Union[str, Iterable[str]],
	*,
	in_2: Union[str, Iterable[str]],
	out_2: Union[str, Iterable[str]],
	linker_file: Optional[str],
	bc_file: str,
	stats_file: str,
	bc_table: Optional[str] = None,
	total: int = defaults.total,
	len_primer: int = defaults.len_primer,
	len_linker: int = defaults.len_linker,
	in_compression: Optional[str] = None,
	out_compression: Optional[str] = None,
	dry_run=False,
	log_init=True
):
	has_two_reads = bool(in_2)
	
	if dry_run:
		if bc_table:
			print('Would write bc table to', bc_table)
		print('Would write', 'two read files:' if has_two_reads else 'one read file:')
		print('\tWould read from', in_1, f'and {in_2}' if has_two_reads else '')
		print('\tWould write to', out_1, f'and {out_2}' if has_two_reads else '')
		print('Would write stats to', stats_file)
		return
	
	if log_init:
		init_logging()
	
	bcs_all = list(read_fasta(bc_file))
	if linker_file:  # If a linker file is passed, we match the full thing
		linkers = dict(read_fasta(linker_file))
		bcs_l = [(h, bc+linkers['Left' ]) for h, bc in bcs_all if h[0] == 'L']
		bcs_r = [(h, bc+linkers['Right']) for h, bc in bcs_all if h[0] == 'R']
		bcs_all = bcs_l + bcs_r
		len_linker = 0
	
	if bc_table:
		write_bc_tables([bc_file], bc_table)
	
	# Two taggers to get two sets of statistics
	tagger1 = get_tagger(bcs_all, len_linker, len_primer)
	tagger2 = get_tagger(bcs_all, len_linker, len_primer) if has_two_reads else None
	
	with tqdm(total=total) if total != 0 else ctx_dummy() as pb, \
			transparent_open(in_1, 'rt', suffix=in_compression) as f_in_1, \
			transparent_open(out_1, 'wt', suffix=out_compression, ensure_parentdir=True) as f_out_1, \
			transparent_open(in_2, 'rt', suffix=in_compression) \
				if has_two_reads else ctx_dummy() as f_in_2, \
			transparent_open(out_2, 'wt', suffix=out_compression, ensure_parentdir=True) \
				if has_two_reads else ctx_dummy() as f_out_2:
		
		n_reads = 0
		n_both_regular = 0
		if has_two_reads:
			for r, (parts1, parts2) in enumerate(zip(iter_fq(f_in_1), iter_fq(f_in_2))):
				read1 = tagger1.tag_read(*parts1)
				read2 = tagger2.tag_read(*parts2)
				
				try:
					if read1.is_regular and read2.is_regular:
						n_both_regular += 1
						f_out_1.write(str(read1))
						f_out_2.write(str(read2))
				except BrokenPipeError:
					break
				
				update_pb(pb, tagger1, r)
				n_reads = r
		# TODO: maybe both?
		else:
			for r, (header, seq_read, seq_qual) in enumerate(iter_fq(f_in_1)):
				read = tagger1.tag_read(header, seq_read, seq_qual)
				
				try:
					f_out_1.write(str(read))
				except BrokenPipeError:
					break
				
				update_pb(pb, tagger1, r)
				n_reads = r
		
		if pb:
			pb.close()
	
	write_stats(
		stats_file, n_reads, n_both_regular if has_two_reads else None,
		tagger1.stats, tagger2.stats if has_two_reads else None,
	)


def update_pb(pb: Optional[tqdm], tgr: ReadTagger, i: int):
	if not pb:
		return
	pb.update(1)
	if i % 10000 == 0:
		pb.set_description(', '.join(f'{k}: {v}' for k, v in tgr.stats.items()), refresh=False)


@contextmanager
def ctx_dummy() -> Generator[None, None, None]:
	yield None
