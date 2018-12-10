import numpy as np
import pandas as pd

from plotnine import (
	ggplot, aes,
	geom_tile,
	scale_fill_cmap,
	coord_fixed,
	labs, theme,
	element_text,
)


def bc_range(bc_series: pd.Series):
	prefixes = set(bc[0] for bc in bc_series)
	r = []
	for prefix in sorted(prefixes):
		bc_ns = set(int(bc[1:]) for bc in bc_series[bc_series.str.startswith(prefix)])
		r = r + [f'{prefix}{i:02}' for i in range(min(bc_ns), max(bc_ns) + 1)]
	return r


def plot_counts(counts: pd.DataFrame):
	counts.bc_l = pd.Categorical(counts.bc_l, bc_range(counts.bc_l))
	counts.bc_r = pd.Categorical(counts.bc_r, bc_range(counts.bc_r))
	log = np.log
	return (
		ggplot(counts, aes('bc_r', 'bc_l', fill='log(Count)'))
		+ geom_tile()
		# + scale_y_reverse()
		+ coord_fixed()
		+ scale_fill_cmap('inferno')
		+ theme(axis_text_x=element_text(angle=90, vjust=.5))
		+ labs(x='Right Barcode', y='Left Barcode')
	)
