from bartseq.read_tagger.dna import reverse_complement


def test_reverse_complement():
	assert 'ATATATAAATTTTTTTTCCCC' == reverse_complement('GGGGAAAAAAAATTTATATAT')