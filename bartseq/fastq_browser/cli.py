import re
from argparse import ArgumentParser, Namespace
from functools import wraps
from pathlib import Path

from .main import main
from ..io import openers
from ..cli_helpers import CLI, t_out_file


RE_READ_FILE = re.compile(r'(?P<lib>.+)_R(?P<read>[12])_001\.fastq\.gz')


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
			'data_dir', type=Path,
			help='Data directory to read from. Needs to have the directories “./process/{3-tagged,4-mapped}” filled.')
		parser.add_argument(
			'library', nargs='?', default=None, help=(
				'Library name. E.g. “Lib1_S1_L001” for input files named “Lib1_S1_L001_R{12}_001.fastq.gz”. '
				'Omittable if only one library exists.'))
		parser.add_argument(
			'out', nargs='?', default='-', type=t_out_file,
			help='FASTA file to write to. Supported compression: see --out-compression')
		parser.add_argument(
			'--out-compression', '-o', choices=openers.keys(),
			help='Specify compression if writing to stdout or a file with unusual suffix')
		return parser
	
	@staticmethod
	def check_args(parser: ArgumentParser, args: Namespace):
		if args.library is None:
			libraries = {
				RE_READ_FILE.fullmatch(p.name)['lib']
				for p in (args.data_dir / 'in' / 'reads').iterdir()
			}
			if len(libraries) == 1:
				args.library, = libraries
			else:
				parser.error(f'You have to specify a library from: {", ".join(libraries)}')
	
	@staticmethod
	@clean_kbdinterrupt
	def run(parser: ArgumentParser, args: Namespace):
		kwargs = vars(args)
		del kwargs['func']
		main(**kwargs)


cli = FastqBrowserCLI()
