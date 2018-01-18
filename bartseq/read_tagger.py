from collections import OrderedDict
from typing import NamedTuple, Iterable, FrozenSet, Tuple, Optional, Generator

from ahocorasick import Automaton

from .logging import log


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
	
	def is_just_primer(self, len_primer):
		return self.barcode is not None and len(self.amplicon) <= len_primer


BASES = set('ATGC')


def get_mismatches(barcode: str, *, max_mm: int=1) -> Generator[str, None, None]:
	yield barcode
	if max_mm != 1:
		raise NotImplemented
	if max_mm == 1:
		for i in range(len(barcode)):
			for mismatch in BASES - {barcode[i]}:
				yield f'{barcode[:i]}{mismatch}{barcode[i+1:]}'


def get_all_barcodes(barcodes: Iterable[str], *, max_mm: int=1):
	found = {}
	blacklist = set()
	
	for barcode in barcodes:
		for pattern in get_mismatches(barcode, max_mm=max_mm):
			previous_barcode = found.get(pattern)
			if previous_barcode:
				blacklist.add(pattern)
				log.warning(
					'Barcodes with one mismatch are ambiguous: '
					f'Modification {pattern} encountered in '
					f'barcode {previous_barcode} and {barcode}'
				)
			found[pattern] = barcode
	
	for pattern in blacklist:
		del found[pattern]
	
	return found


class ReadTagger:
	def __init__(self, barcodes: Iterable[str], *, max_mm: int=1):
		self.barcodes = barcodes
		self.automaton = Automaton()
		for pattern, barcode in get_all_barcodes(barcodes, max_mm=max_mm).items():
			self.automaton.add_word(pattern, barcode)
		self.automaton.make_automaton()
	
	def search_barcode(self, read: str) -> Tuple[int, int, str]:
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
