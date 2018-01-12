from pytest import raises

from bartseq.read_tagger import get_mismatches, ReadTagger, TaggedRead


def test_get_mismatches_1():
	mms = list(get_mismatches('a'))
	expected = ['a', 'A', 'T', 'G', 'C']
	assert sorted(expected) == sorted(mms), set(mms) ^ set(expected)


expected_2 = [
	'ab',
	'Ab', 'Tb', 'Gb', 'Cb',
	'aA', 'aT', 'aG', 'aC',
]


def test_get_mismatches_2():
	mms = list(get_mismatches('ab'))
	assert sorted(expected_2) == sorted(mms), set(mms) ^ set(expected_2)


def test_get_mismatches_3():
	mms = list(get_mismatches('abc'))
	expected = [
		'abc',
		'Abc', 'Tbc', 'Gbc', 'Cbc',
		'aAc', 'aTc', 'aGc', 'aCc',
		'abA', 'abT', 'abG', 'abC',
	]
	assert sorted(expected) == sorted(mms), set(mms) ^ set(expected)


def test_search_barcode():
	tagger = ReadTagger(['ab'])
	assert {'ab'} == set(tagger.automaton.values())
	assert set(expected_2) == set(tagger.automaton.keys())
	
	matches = list(tagger.search_barcode('_ab,Gb,Xb'))
	assert [(1, 3, 'ab'), (4, 6, 'ab')] == matches


# searches


def test_tagger_ambiguous_barcodes():
	with raises(ValueError, match='ambiguous.*in barcode ab and ac'):
		ReadTagger(['ab', 'ac'])


def test_tag_read_find1():
	tagger = ReadTagger(['ab'])
	
	assert tagger.tag_read('XXabblah') == TaggedRead('XX', 'ab', 'blah', frozenset())
	assert tagger.tag_read('XXbvblah') == TaggedRead(None, None, 'XXbvblah', frozenset())
	# two occurrences
	assert tagger.tag_read('XXabblab') == TaggedRead('XX', 'ab', 'blab', frozenset())


def test_tag_read_find_mismatch():
	tagger = ReadTagger(['ab'])
	
	assert tagger.tag_read('XXaGblah') == TaggedRead('XX', 'ab', 'blah', frozenset())
	# two occurrences
	assert tagger.tag_read('XXaGblab') == TaggedRead('XX', 'ab', 'blab', frozenset())


def test_tag_read_find2():
	tagger = ReadTagger(['ab', 'bx'])
	assert tagger.tag_read('XXabxblah') == TaggedRead('XX', 'ab', 'xblah', frozenset({'bx'}))
