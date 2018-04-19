import sys
from argparse import ArgumentParser, _SubParsersAction, Namespace
from typing import Sequence, Optional

from .read_tagger.cli import ReadTaggerCLI

SUBCMDS = [
	ReadTaggerCLI,
]


def run_cli(argv: Optional[Sequence[str]] = None):
	if argv is None:
		argv = sys.argv[1:]
	
	parser = ArgumentParser()
	subparsers = parser.add_subparsers()  # type: _SubParsersAction
	
	for subcmd in SUBCMDS:
		subcmd.register_as_subcmd(subparsers)
	
	args = parser.parse_args(argv)  # type: Namespace
	if 'func' in args:
		return args.func(args)
	elif len(argv) == 0:
		parser.print_help()
	else:
		assert False, 'This should never happen'
