from collections import OrderedDict
from typing import NamedTuple, Iterable, FrozenSet, Tuple, Optional

from ahocorasick import Automaton


LEN_PRIMER = NotImplemented


_TaggedReadBase = NamedTuple('_TaggedReadBase', [
	('junk', Optional[str]),
	('barcode', Optional[str]),
	('amplicon', str),
	('other_barcodes', FrozenSet[str]),
])

class TaggedRead(_TaggedReadBase):
	@property
	def has_multiple_barcodes(self):
		return len(self.other_barcodes) > 0
	
	@property
	def is_just_primer(self):
		return self.barcode is not None and len(self.amplicon) > LEN_PRIMER


BASES = set('ATGC')


def get_mismatches(barcode):
	yield barcode
	for i in range(len(barcode)):
		for mismatch in BASES - {barcode[i]}:
			yield f'{barcode[:i]}{mismatch}{barcode[i+1:]}'


class ReadTagger:
	def __init__(self, barcodes: Iterable[str]):
		self.barcodes = barcodes
		self.automaton = Automaton()
		for barcode in barcodes:
			for pattern in get_mismatches(barcode):
				self.automaton.add_word(pattern, barcode)
		self.automaton.make_automaton()
	
	def search_barcode(self, read) -> Tuple[int, int, str]:
		for end, barcode in self.automaton.iter(read):
			start = end - len(barcode) + 1
			yield start, end + 1, barcode
	
	def tag_read(self, read: str) -> TaggedRead:
		# as ordered set
		matches = OrderedDict((match, None) for match in self.search_barcode(read))
		
		match_iter = iter(matches)
		bc_start, bc_end, barcode = next(match_iter, (None, None, None))
		other_barcodes = frozenset(set(bc for _, _, bc in match_iter) - {barcode})
		
		junk = read[:bc_start] if bc_start else None
		amplicon = read[bc_end:] if bc_end else read
		# TODO: no amplicon if only primer
		return TaggedRead(junk, barcode, amplicon, other_barcodes)
	
	def tag_reads(self, reads: Iterable[str]):
		for read in reads:
			yield self.tag_read(read)
