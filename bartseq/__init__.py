from pathlib import Path

from .read_tagger import ReadTagger


def main(args):
	path_in = Path(args[1])
	
	with path_in.open() as f_in:
		tagger = ReadTagger(['TODO'])
		tagger.tag_reads(l.strip() for l in f_in)
