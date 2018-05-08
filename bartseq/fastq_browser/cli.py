import re
from argparse import ArgumentParser, Namespace
from functools import wraps
from itertools import cycle

from pathlib import Path

from tqdm import tqdm

from ..io import openers, transparent_open, iter_fq
from ..cli_helpers import CLI, t_out_file


fields = [
	'barcode', 'linker', 'multi-bc', 'just-primer',
	'other-bcs', 'barcode-mismatch', 'junk',
]

RE_FIELDS = re.compile(f'(?P<field>{"|".join(fields)})=(?P<value>[^ ]+)')
RE_READ_FILE = re.compile(r'(?P<lib>.+)_R(?P<read>[12])_001\.fastq\.gz')


def interleave(*iterables):
	unexhausted = set(range(len(iterables)))
	for i, it in cycle(enumerate(map(iter, iterables))):
		if not unexhausted:
			break
		try:
			yield i, next(it)
		except StopIteration:
			unexhausted -= {i}


def clean_kbdinterrupt(callback):
	@wraps(callback)
	def decorated(*args, **kw):
		try:
			return callback(*args, **kw)
		except KeyboardInterrupt:
			pass
	return decorated


class FastqBrowserCLI(CLI):
	@staticmethod
	def populate_parser(parser: ArgumentParser) -> ArgumentParser:
		parser.add_argument(
			'data_dir',
			help='Data directory to read from. Needs to have the directories “./process/{3-tagged,4-mapped}” filled.')
		parser.add_argument(
			'library',
			help='Library name. E.g. “Lib1_S1_L001” for input files named “Lib1_S1_L001_R{12}_001.fastq.gz”')
		parser.add_argument(
			'out', nargs='?', default='-', type=t_out_file,
			help='FASTA file to write to. Supported compression: see --out-compression')
		parser.add_argument(
			'--out-compression', '-o', choices=openers.keys(),
			help='Specify compression if writing to stdout or a file with unusual suffix')
		return parser
	
	@staticmethod
	@clean_kbdinterrupt
	def run(parser: ArgumentParser, args: Namespace):
		dir_process = Path(args.data_dir) / 'process'
		dir_tagged = dir_process / '3-tagged'
		dir_mapped = dir_process / '4-mapped'
		
		paths_fsq = [dir_tagged / f'{args.library}_R{r}_001.fastq.gz' for r in [1, 2]]
		paths_map = [dir_mapped / f'{args.library}_R{r}_001.txt' for r in [1, 2]]
		path_count = dir_process / '1-index' / f'{args.library}_001.count.txt'
		
		with path_count.open() as c_f:
			n_pairs = int(c_f.read())
		
		with \
			transparent_open(paths_fsq[0]) as fsq_r1, \
			transparent_open(paths_fsq[1]) as fsq_r2, \
			transparent_open(paths_map[0]) as map_r1, \
			transparent_open(paths_map[1]) as map_r2, \
			transparent_open(args.out, 'wt', suffix=args.out_compression, ensure_parentdir=True) as out:
			
			print(
				'read',
				'header', 'read_seq', 'quality_seq',
				'amplicon',
				'barcode', 'linker', 'has_multiple_bcs',
				'is_just_primer', 'other_bcs', 'has_bc_mismatch', 'junk',
				sep='\t', file=out,
			)
			
			for read_side, ((header, read, qual), amp) in tqdm(interleave(
				zip(iter_fq(fsq_r1), map_r1),
				zip(iter_fq(fsq_r2), map_r2),
			), total=2*n_pairs):
				fields = list(RE_FIELDS.finditer(header))
				
				print(
					read_side + 1,
					header[:fields[0].start()-1],
					read, qual,
					amp.strip(),
					*[v if v != 'None' else '' for v in (m.group('value') for m in fields)],
					sep='\t', file=out,
				)


cli = FastqBrowserCLI()
