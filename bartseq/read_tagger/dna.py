COMPLEMENT = str.maketrans('ATGC', 'TACG')


def reverse_complement(seq: str):
	return seq[::-1].translate(COMPLEMENT)
