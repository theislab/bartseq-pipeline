from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl

from plotnine import (
	ggplot, aes,
	geom_tile,
	scale_fill_cmap, scale_y_reverse, scale_x_discrete, scale_y_discrete,
	coord_fixed, facet_wrap,
	labs, theme,
	element_text,
)


def get_count_files(count_dir):
	for d in count_dir.iterdir():
		if not d.is_dir(): continue
		path = (count_dir / d.name / d.name).with_suffix('.tsv')
		yield d.name, pd.read_csv(path, '\t')


def bc_range(bc_series):
	prefix = next(iter(bc_series))[0]
	n_min = int(bc_series.unique().min()[1:])
	n_max = int(bc_series.unique().max()[1:])
	return [f'{prefix}{i:02}' for i in range(n_min, n_max+1)]


def get_counts(count_dir):
	count_dir = Path(count_dir)
	counts_long = pd.concat(dict(get_count_files(count_dir)), names=['Amplicon']).reset_index(0)
	counts_na = counts_long.melt(['Amplicon', 'bc_l'], var_name='bc_r', value_name='Count')
	counts = counts_na[~counts_na.Count.isna()].copy()
	counts.bc_l = pd.Categorical(counts.bc_l, bc_range(counts.bc_l))
	counts.bc_r = pd.Categorical(counts.bc_r, bc_range(counts.bc_r))
	return counts


def plot_counts(data_dir, nrow):
	count_dir = Path(data_dir) / 'out/counts'
	counts = get_counts(count_dir)
	
	gg = (ggplot(counts, aes('bc_r', 'bc_l', fill='np.log(Count)')) +
		facet_wrap('~Amplicon', nrow) + 
		geom_tile() +
		#scale_y_reverse() +
		coord_fixed() +
		scale_fill_cmap('inferno') +
		theme(
			axis_text_x=element_text(size=1.8, angle=90, vjust=.5),
			axis_text_y=element_text(size=1.8),
		) +
		labs(x='Right Barcode', y='Left Barcode'))
	
	gg.save(str(count_dir / 'all.png'), dpi=300, verbose=False)
	return gg


if __name__ == '__main__':
	plot_counts('ngs15', 5)
	plot_counts('ngs23', 6)
