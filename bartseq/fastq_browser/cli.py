from argparse import ArgumentParser, Namespace

from pathlib import Path

from ..io import openers, transparent_open, iter_fq
from ..cli_helpers import CLI, t_out_file


class FastqBrowserCLI(CLI):
	@staticmethod
	def populate_parser(parser: ArgumentParser) -> ArgumentParser:
		parser.add_argument(
			'fq_in',  # needs to be a path not stdin
			help='FASTQ file to read from. Supported compression: see --in-compression')
		parser.add_argument(
			'fa_out', nargs='?', default='-', type=t_out_file,
			help='FASTA file to write to. Supported compression: see --out-compression')
		parser.add_argument(
			'--out-compression', '-o', choices=openers.keys(),
			help='Specify compression if writing to stdout or a file with unusual suffix')
		return parser
	
	@staticmethod
	def run(parser: ArgumentParser, args: Namespace):
		path_fq_in = Path(args.fq_in)
		if path_fq_in.suffix in openers:
			stem = path_fq_in.with_suffix('').with_suffix('').name
		else:
			stem = path_fq_in.with_suffix('').name
		path_map_in = (path_fq_in.parent.parent / '4-mapped' / stem).with_suffix('.txt')
		
		with \
				transparent_open(path_fq_in) as fq_in, \
				transparent_open(path_map_in) as map_in, \
				transparent_open(args.fa_out, suffix=args.out_compression, ensure_parentdir=True) as fa_out:
			for (header, seq_read, seq_qual), amp in zip(iter_fq(fq_in), map_in):
				fa_out.write(f'>{header} amplicon={amp.strip()}\n{seq_read}\n')


cli = FastqBrowserCLI()
