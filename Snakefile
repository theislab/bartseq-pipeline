# usage example: snakemake -d data/ngs15 -j 4
import re
import json
from collections import Counter
from functools import reduce
from operator import add
from pathlib import Path

from snakemake.utils import min_version, listfiles
import pandas as pd
from tqdm import tqdm
import matplotlib
matplotlib.rcParams['backend'] = 'agg'  # make pypy work without Qt
from plotnine import facet_wrap, theme, element_text

from bartseq.io import write_bc_table, transparent_open
from bartseq.heatmaps import plot_counts


min_version('4.5.1')


with open('in/amplicons.fa') as a_f:
	amplicons = [line.lstrip('>').strip() for line in a_f.readlines() if line.startswith('>')]
	amplicons += ['unmapped', 'one-mapped']

dir_qc = 'out/qc'
amplicon_index_stem = 'process/1-index/amplicons'
amplicon_index_files = expand('{stem}.{n}.ht2', stem=[amplicon_index_stem], n=range(1, 9))
all_reads_in = [(w.readname, w.read) for _, w in listfiles('in/reads/{readname}_R{read,[12]}_001.fastq.gz')]
lib_names = [readname for readname, read in all_reads_in if read == '1']

def get_read_path(prefix, name, read, suffix='.fastq.gz'):
	return '{prefix}/{name}_R{read}_001{suffix}'.format_map(locals())

def get_read_paths(prefix, *suffixes):
	if len(suffixes) == 0:
		suffixes = ['.fastq.gz']
	return [
		get_read_path(prefix, n, r, suffix)
		for n, r in all_reads_in
		for suffix in suffixes
	]

reads_raw = get_read_paths('rawdata')

wildcard_constraints:
    which = '.*'

rule all:
	input:
		get_read_paths('out/qc', '_fastqc.html', '_fastqc.zip'),
		expand('out/counts/{amplicon}/{amplicon}{which}-log.png', amplicon=amplicons, which=['-all', '']),
		expand('out/counts/{amplicon}/bylib/{amplicon}-{lib_name}{which}-log.png', amplicon=amplicons, lib_name=lib_names, which=['-all', '']),
		'out/counts/all.png',
		'out/barcodes.htm',

rule get_qc:
	input:
		'in/reads/{name_full}.fastq.gz'
	output:
		expand('{dir_qc}/{{name_full}}_fastqc{suffix}', dir_qc=[dir_qc], suffix=['.html', '.zip'])
	shell:
		'fastqc {input:q} -o {dir_qc:q}'

rule get_read_count:
	input:
		'in/reads/{name}_R1_001.fastq.gz'
	output:
		'process/1-index/{name}_001.count.txt'
	shell:
		'''
		lines=$(zcat {input:q} | wc -l)
		echo "$lines / 4" | bc > {output:q}
		'''

rule trim_quality:
	input:
		expand('in/reads/{{name}}_R{read}_001.fastq.gz', read=[1,2]),
	output:
		expand('process/2-trimmed/{{name}}_R{read}_001.fastq.gz', read=[1,2]),
		single = 'process/2-trimmed/{name}_single_001.fastq.gz',
	shell:
		'''
		sickle pe --gzip-output --qual-type=sanger \
			--pe-file1={input[0]:q} --output-pe1={output[0]:q} \
			--pe-file2={input[1]:q} --output-pe2={output[1]:q} \
			--output-single={output.single:q}
		'''

rule bc_table:
	input:
		'in/barcodes.fa'
	output:
		'out/barcodes.htm'
	run:
		write_bc_table(input[0], output[0])

rule tag_reads:
	input:
		expand('process/2-trimmed/{{name}}_R{read}_001.fastq.gz', read=[1,2]),
		count_file = 'process/1-index/{name}_001.count.txt',
		bc_file = 'in/barcodes.fa',
	output:
		expand('process/3-tagged/{{name}}_R{read}_001.fastq.gz', read=[1,2]),
		stats_file='process/3-tagged/{name}_stats.json',
	run:
		from bartseq.main import run
		with open(input.count_file) as c_f:
			total = int(c_f.read())
		run(
			in_1=input[0], out_1=output[0],
			in_2=input[1], out_2=output[1],
			bc_file=input.bc_file,
			stats_file=output.stats_file,
			total=total,
		)

rule tag_stats:
	input:
		expand('process/3-tagged/{name}_stats.json', name=lib_names)
	run:
		for path in input:
			with open(path) as f:
				stats = json.load(f)
				print(path)
				print('', 'n_reads', stats['n_reads'], sep='\t')
				print('', 'both_regular', '{:.1%}'.format(stats['n_both_regular'] / stats['n_reads']), sep='\t')
				for read in ['read1', 'read2']:
					print('', read, sep='\t')
					for stat, count in stats[read].items():
						print('', '', stat[2:], '{:.1%}'.format(count / stats['n_reads']), sep='\t')

# Helper rule for Lukas’ pipeline file
rule amplicon_fa:
	input:
		'in/amplicons.txt'
	output:
		'in/amplicons.fa'
	shell:
		r"cat {input:q} | tail -n +2 | cut -f 1,4 | sed 's/^ */>/;s/\t/\n/' > {output:q}"

rule build_index:
	input:
		'in/amplicons.fa'
	output:
		idx = amplicon_index_files,
	threads: 4
	shell:
		'hisat2-build -p {threads} {input:q} {amplicon_index_stem:q}'

rule map_reads:
	input:
		amplicons = amplicon_index_files,
		read = 'process/3-tagged/{name_full}.fastq.gz',
	output:
		map = 'process/4-mapped/{name_full}.txt',
		summary = 'process/4-mapped/{name_full}_summary.txt'
	threads: 4
	shell:
		'''
		hisat2 \
			--threads {threads} \
			--reorder \
			-k 1 \
			-x {amplicon_index_stem:q} \
			--new-summary --summary-file {output.summary:q} \
			-q -U {input.read:q} | \
			grep -v "^@" - | \
			cut -f3 > {output.map:q}
		'''

rule count:
	input:
		reads = expand('process/3-tagged/{{name}}_R{read}_001.fastq.gz', read=[1,2]),
		mappings = expand('process/4-mapped/{{name}}_R{read}_001.txt', read=[1,2]),
		stats_file = 'process/3-tagged/{name}_stats.json'
	output:
		'process/5-counts/{name}_001.tsv'
	run:
		bc_re = re.compile(r'barcode=(\w+)')
		with open(input.stats_file) as s_f:
			total = json.load(s_f)['n_both_regular']
		def get_barcodes(fastq_file):
			for header in fastq_file:
				next(fastq_file)
				next(fastq_file)
				next(fastq_file)
				yield bc_re.search(header).group(1)
		counts = Counter()
		with \
			transparent_open(input.reads[0]) as r1, open(input.mappings[0]) as a1, \
			transparent_open(input.reads[1]) as r2, open(input.mappings[1]) as a2:
			
			bcs1 = get_barcodes(r1)
			bcs2 = get_barcodes(r2)
			amps1 = (a.strip() for a in a1)
			amps2 = (a.strip() for a in a2)
			for bc1, bc2, amp1, amp2 in tqdm(zip(bcs1, bcs2, amps1, amps2), total=total):
				bc1, bc2 = sorted([bc1, bc2])
				if amp1 == amp2:
					counts[bc1, bc2, 'unmapped' if amp1 == '*' else amp1] += 1
				else:
					counts[bc1, bc2, 'one-mapped'] += 1
			
			with open(output[0], 'w') as of:
				print('bc_l', 'bc_r', 'amp', 'count', sep='\t', file=of)
				for fields, c in counts.items():
					print(*fields, c, sep='\t', file=of)

rule amplicon_counts_lib_all:
	input:
		'process/5-counts/{name}_001.tsv'
	output:
		'out/counts/{amplicon}/bylib/{amplicon}-{name}-all.tsv'
	run:
		entries_lib = pd.read_csv(input[0], '\t')
		entries = entries_lib[entries_lib.amp == wildcards.amplicon].drop(columns=['amp'])
		table = entries.pivot('bc_l', 'bc_r', 'count')
		table.to_csv(output[0], '\t')

rule amplicon_counts_all:
	input:
		expand('out/counts/{{amplicon}}/bylib/{{amplicon}}-{name}-all.tsv', name=lib_names)
	output:
		'out/counts/{amplicon}/{amplicon}-all.tsv'
	run:
		table = reduce(add, [pd.read_csv(f, '\t', index_col='bc_l') for f in input])
		table.to_csv(output[0], '\t')

rule amplicon_counts:
	input:
		'out/counts/{amplicon}/{matrix}-all.tsv'
	output:
		'out/counts/{amplicon}/{matrix}.tsv'
	run:
		table_all = pd.read_csv(input[0], '\t', index_col='bc_l')
		if table_all.shape == (0, 0):  # here, the empty index can’t use str methods
			table = table_all
		else:
			table = table_all.loc[table_all.index.str.match('L.*'), table_all.columns.str.match('R.*')]
		table.to_csv(output[0], '\t')

rule plot_counts:
	input:
		'out/counts/{amplicon}/{matrix}{which}.tsv'
	output:
		'out/counts/{amplicon}/{matrix}{which}-log.png'
	run:
		counts_long = pd.read_csv(input[0], '\t')
		counts_na = counts_long.melt(['bc_l'], var_name='bc_r', value_name='Count')
		counts = counts_na[~counts_na.Count.isna()].copy()
		
		if counts.shape[0] <= 1:
			with open(output[0], 'wb') as empty:
				empty.write(b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01\x08\x04\x00\x00\x00\xb5\x1c\x0c\x02\x00\x00\x00\x0bIDATx\xdacd`\x00\x00\x00\x06\x00\x020\x81\xd0/\x00\x00\x00\x00IEND\xaeB`\x82')
			return
		
		gg = plot_counts(counts) + \
			theme(axis_text_x=element_text(size=4), axis_text_y=element_text(size=4))
		gg.save(output[0], dpi=300, verbose=False)

rule plot_counts_all:
	input:
		expand('out/counts/{amplicon}/{amplicon}.tsv', amplicon=amplicons)
	output:
		expand('out/counts/all.png')
	run:
		tables = {Path(i).parent.name: pd.read_csv(i, '\t') for i in input}
		table = pd.concat(tables, names=['Amplicon']).reset_index(0)
		counts_na = table.melt(['Amplicon', 'bc_l'], var_name='bc_r', value_name='Count')
		counts = counts_na[~counts_na.Count.isna()].copy()
		
		gg = plot_counts(counts) + \
			facet_wrap('~Amplicon', ncol=4) + \
			theme(axis_text_x=element_text(size=1.8), axis_text_y=element_text(size=1.8))
		gg.save(output[0], dpi=300, verbose=False)

#Needs https://bitbucket.org/snakemake/snakemake/pull-requests/264
rule dag:
	shell:
		'snakemake -s {__file__} --dag | dot -Tsvg | gwenview /dev/stdin'
