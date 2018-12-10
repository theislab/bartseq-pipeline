from pytest import warns

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

qual_9 = '!' * 9


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
	tagger = ReadTagger(dict(ab='A'), 1, 1)
	assert {'ab'} == set(tagger.automaton.values())
	assert set(expected_2) == set(tagger.automaton.keys())
	
	matches = list(tagger.search_barcode('_ab,Gb,Xb'))
	assert [(1, 3, 'ab'), (4, 6, 'ab')] == matches


# searches


def test_tagger_ambiguous_barcodes():
	with warns(UserWarning, match='ambiguous.*in barcode ab and ac'):
		ReadTagger(dict(ab='A', ac='B'), 1, 1)


def test_tag_read_find1():
	tagger = ReadTagger(dict(ab='A'), 1, 1)
	
	assert tagger.tag_read('a', 'XXabLblah', qual_9) == TaggedRead('a', qual_9, 1, 'XX', 'A', 'L', 'blah', frozenset(), False)  # noqa
	assert tagger.tag_read('b', 'abLblahXX', qual_9) == TaggedRead('b', qual_9, 1, None, 'A', 'L', 'blahXX', frozenset(), False)  # noqa
	assert tagger.tag_read('c', 'XXbvblahL', qual_9) == TaggedRead('c', qual_9, 1, None, None, None, 'XXbvblahL', frozenset(), False)  # noqa
	# two occurrences
	assert tagger.tag_read('d', 'XXabLblab', qual_9) == TaggedRead('d', qual_9, 1, 'XX', 'A', 'L', 'blab', frozenset(), False)  # noqa
	
	assert tagger.stats == dict(
		n_only_primer=0,
		n_multiple_bcs=0,
		n_no_barcode=1,
		n_barcode_mismatch=0,
		n_junk=2,
		n_regular=3,
	)


def test_tag_read_find_mismatch():
	tagger = ReadTagger(dict(ab='A'), 1, 1)
	
	assert tagger.tag_read('a', 'XXaGLblah', qual_9) == TaggedRead('a', qual_9, 1, 'XX', 'A', 'L', 'blah', frozenset(), True)  # noqa
	# two occurrences
	assert tagger.tag_read('b', 'XXaGLblab', qual_9) == TaggedRead('b', qual_9, 1, 'XX', 'A', 'L', 'blab', frozenset(), True)  # noqa
	
	assert tagger.stats == dict(
		n_only_primer=0,
		n_multiple_bcs=0,
		n_no_barcode=0,
		n_barcode_mismatch=2,
		n_junk=2,
		n_regular=2,
	)


def test_tag_read_find2():
	tagger = ReadTagger(dict(ab='A', bL='B'), 1, 1)
	assert tagger.tag_read('a', 'XXabLxblah', qual_9) == TaggedRead('a', qual_9, 1, 'XX', 'A', 'L', 'xblah', frozenset({'B'}), False)  # noqa
	
	assert tagger.stats == dict(
		n_only_primer=0,
		n_multiple_bcs=1,
		n_no_barcode=0,
		n_barcode_mismatch=0,
		n_junk=1,
		n_regular=1,
	)
