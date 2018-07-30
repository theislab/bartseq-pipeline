# usage example: snakemake -d data/ngs15 -j 4
import re
import json
from collections import Counter
from functools import reduce
from operator import add
from pathlib import Path
from typing import Dict

from snakemake.utils import min_version, listfiles, format
import pandas as pd
from tqdm import tqdm
import matplotlib
matplotlib.rcParams['backend'] = 'agg'  # make pypy work without Qt
from plotnine import facet_wrap, theme, element_text

from bartseq.io import transparent_open
from bartseq.read_tagger.io import write_bc_tables
from bartseq.read_tagger.defaults import len_linker
from bartseq.heatmaps import plot_counts

min_version('4.5.1')

EMPTY_PNG = b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01\x08\x04\x00\x00\x00\xb5\x1c\x0c\x02\x00\x00\x00\x0bIDATx\xdacd`\x00\x00\x00\x06\x00\x020\x81\xd0/\x00\x00\x00\x00IEND\xaeB`\x82'

dir_qc = 'out/qc'
amplicon_index_stem = 'process/1-index/amplicons'
amplicon_index_files = expand('{stem}/{{lib_name}}.{n}.ht2', stem=amplicon_index_stem, n=range(1, 9))
all_reads_in = [(w.readname, w.read) for _, w in listfiles('in/reads/{readname}_R{read,[12]}_001.fastq.gz')]
lib_names = sorted([readname for readname, read in all_reads_in if read == '1'])

def get_input_seq_paths(dir_):
	by_lib = Path('in', dir_).is_dir()
	for lib in lib_names:
		if by_lib:
			yield lib, Path('in', dir_, lib + '.fa')
		else:
			yield lib, Path('in', dir_+'.fa')

amplicons = {
	lib: [line.lstrip('>') for line in path.read_text().splitlines() if line.startswith('>')] + ['-unmapped', '-one-mapped', '-mismatch']
	for lib, path in get_input_seq_paths('amplicons')
}

len_barcode = max(
	len(line)
	for _, path in get_input_seq_paths('barcodes')
	for line in path.read_text().splitlines()
	if not line.startswith('>')
)

len_protection = 3
len_3prime_junk = len_linker + len_barcode + len_protection

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

re_amplicon = '({})'.format('|'.join(re.escape(a) for a in set(a for amps in amplicons.values() for a in amps)))
re_lib_name = '({})'.format('|'.join(re.escape(ln) for ln in lib_names))
re_matrix = '(both|one|{re_lib_name}/({re_lib_name}|{re_amplicon}/{re_lib_name}-{re_amplicon}))'.format(re_amplicon=re_amplicon, re_lib_name=re_lib_name)


wildcard_constraints:
	which = '(-all|)',
	counting = '(both|one)',
	amplicon = re_amplicon,
	lib_name = re_lib_name,
	matrix = re_matrix,
	read = '[12]'


rule all:
	input:
		get_read_paths('out/qc', '_fastqc.html', '_fastqc.zip'),
		[
			path
			for lib_name, amps in amplicons.items()
			for path in expand(
				'out/counts/{counting}/{lib_name}/{amplicon}/{lib_name}-{amplicon}{which}-log.png',
				counting=['both', 'one'], lib_name=lib_name, amplicon=amps, which=['-all', ''],
			)
		],
		expand('out/counts/{counting}/{lib_name}/{lib_name}{which}-log.png', counting=['both', 'one'], lib_name=lib_names, which=['-all', '']),
		#expand('out/counts/{counting}/{counting}{which}-log.png', counting=['both', 'one'], which=['-all', '']),
		expand('out/counts/{counting}/{counting}.xlsx', counting=['both', 'one']),
		'out/barcodes.htm',

rule get_qc:
	input:
		'in/reads/{name_full}.fastq.gz'
	output:
		expand('{dir_qc}/{{name_full}}_fastqc{suffix}', dir_qc=dir_qc, suffix=['.html', '.zip'])
	shell:
		'fastqc {input:q} -o {dir_qc:q}'

rule get_read_count:
	input:
		'in/reads/{lib_name}_R1_001.fastq.gz'
	output:
		'process/1-index/{lib_name}.count.txt'
	shell:
		'''
		lines=$(zcat {input:q} | wc -l)
		echo "$lines / 4" | bc > {output:q}
		'''

rule seqs_by_lib:
	input:  'in/{seqs_type}/{lib_name}.fa'
	output: 'process/1-index/{seqs_type}/{lib_name}.fa'
	shell:  'cp -T {input:q} {output:q}'

rule seqs_universal:
	input:  'in/{seqs_type}.fa'
	output: 'process/1-index/{seqs_type}/{lib_name}.fa'
	shell:  'cp -T {input:q} {output:q}'

rule trim_quality:
	input:
		expand('in/reads/{{lib_name}}_R{read}_001.fastq.gz', read=[1,2]),
	output:
		expand('process/2-trimmed/{{lib_name}}_R{read}.fastq.gz', read=[1,2]),
		single = 'process/2-trimmed/{lib_name}_single.fastq.gz',
	shell:
		'''
		sickle pe --gzip-output --qual-type=sanger \
			--pe-file1={input[0]:q} --output-pe1={output[0]:q} \
			--pe-file2={input[1]:q} --output-pe2={output[1]:q} \
			--output-single={output.single:q}
		'''

rule bc_table:
	input:
		expand('process/1-index/barcodes/{lib_name}.fa', lib_name=lib_names)
	output:
		'out/barcodes.htm'
	run:
		write_bc_tables(input, output[0])

rule tag_reads:
	input:
		expand('process/2-trimmed/{{lib_name}}_R{read}.fastq.gz', read=[1,2]),
		count_file = 'process/1-index/{lib_name}.count.txt',
		bc_file = 'process/1-index/barcodes/{lib_name}.fa',
	output:
		expand('process/3-tagged/{{lib_name}}_R{read}.fastq.gz', read=[1,2]),
		stats_file='process/3-tagged/{lib_name}_stats.json',
	run:
		from bartseq.read_tagger.main import run
		total = int(Path(input.count_file).read_text('utf-8'))
		run(
			in_1=input[0], out_1=output[0],
			in_2=input[1], out_2=output[1],
			bc_file=input.bc_file,
			stats_file=output.stats_file,
			total=total,
		)

rule tag_stats:
	input:
		expand('process/3-tagged/{lib_name}_stats.json', lib_name=lib_names)
	run:
		for path in input:
			stats = json.loads(Path(path).read_bytes())
			print(path)
			print('', 'n_reads', stats['n_reads'], sep='\t')
			print('', 'both_regular', '{:.1%}'.format(stats['n_both_regular'] / stats['n_reads']), sep='\t')
			for read in ['read1', 'read2']:
				print('', read, sep='\t')
				for stat, count in stats[read].items():
					print('', '', stat[2:], '{:.1%}'.format(count / stats['n_reads']), sep='\t')

rule build_index:
	input:
		'process/1-index/amplicons/{lib_name}.fa'
	output:
		idx = amplicon_index_files,
	threads: 4
	shell:
		'hisat2-build -p {threads} {input:q} {amplicon_index_stem:q}/{wildcards.lib_name:q}'

rule map_reads:
	input:
		amplicons = amplicon_index_files,
		read = 'process/3-tagged/{lib_name}_R{read}.fastq.gz',
	output:
		map = 'process/4-mapped/{lib_name}_R{read}.txt',
		summary = 'process/4-mapped/{lib_name}_R{read}_summary.txt'
	threads: 4
	shell:
		'''
		hisat2 \
			--threads {threads} \
			--reorder \
			-k 1 \
			-3 {len_3prime_junk} \
			-x {amplicon_index_stem:q}/{wildcards.lib_name:q} \
			--new-summary --summary-file {output.summary:q} \
			-q -U {input.read:q} | \
			grep -v "^@" - | \
			cut -f3 > {output.map:q}
		'''

rule count:
	input:
		reads = expand('process/3-tagged/{{lib_name}}_R{read}.fastq.gz', read=[1,2]),
		mappings = expand('process/4-mapped/{{lib_name}}_R{read}.txt', read=[1,2]),
		stats_file = 'process/3-tagged/{lib_name}_stats.json'
	output:
		expand('process/5-counts/{counting}/{{lib_name}}.tsv', counting=['both', 'one'])
	run:
		bc_re = re.compile(r'barcode=(\w+)')
		total = json.loads(Path(input.stats_file).read_bytes())['n_both_regular']
		def get_barcodes(fastq_file):
			for header in fastq_file:
				next(fastq_file)
				next(fastq_file)
				next(fastq_file)
				yield bc_re.search(header).group(1)
		counts = Counter()
		counts_one = Counter()
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
					if amp1 == '*':
						counts[bc1, bc2, '-unmapped'] += 1
					else:
						counts[bc1, bc2, amp1] += 1
						counts_one[bc1, bc2, amp1] += 1
				elif amp1 == '*':
					counts[bc1, bc2, '-one-mapped'] += 1
					counts_one[bc1, bc2, amp2] += 1  # amp1 == '*'
				else:
					counts[bc1, bc2, '-mismatch'] += 1
			
			for path, counter in zip(output, [counts, counts_one]):
				with open(path, 'w') as of:
					print('bc_l', 'bc_r', 'amp', 'count', sep='\t', file=of)
					for fields, c in counter.items():
						print(*fields, c, sep='\t', file=of)

rule amplicon_counts_all:
	input:
		'process/5-counts/{counting}/{lib_name}.tsv'
	output:
		'out/counts/{counting}/{lib_name}/{amplicon}/{lib_name}-{amplicon}-all.tsv'
	run:
		entries_lib = pd.read_csv(input[0], '\t')
		entries = entries_lib[entries_lib.amp == wildcards.amplicon].drop(columns=['amp'])
		table = entries.pivot('bc_l', 'bc_r', 'count')
		table.to_csv(output[0], '\t')

rule lib_counts_all:
	input:
		lambda wildcards: expand(
			'out/counts/{counting}/{lib_name}/{amplicon}/{lib_name}-{amplicon}-all.tsv',
			counting=wildcards.counting,
			lib_name=wildcards.lib_name,
			amplicon=amplicons[wildcards.lib_name],
		)
	output:
		'out/counts/{counting}/{lib_name}/{lib_name}-all.tsv'
	run:
		table = reduce(add, [pd.read_csv(f, '\t', index_col='bc_l') for f in input])
		table.to_csv(output[0], '\t')

# TODO: maybe merge with above somehow?
rule counts_all:
	input:
		expand('out/counts/{{counting}}/{lib_name}/{lib_name}-all.tsv', lib_name=lib_names)
	output:
		'out/counts/{counting}/{counting}-all.tsv'
	run:
		table = reduce(add, [pd.read_csv(f, '\t', index_col='bc_l') for f in input])
		table.to_csv(output[0], '\t')

rule amplicon_counts:
	input:
		'out/counts/{counting}/{matrix}-all.tsv'
	output:
		'out/counts/{counting}/{matrix}.tsv'
	run:
		table_all = pd.read_csv(input[0], '\t', index_col='bc_l')
		if table_all.shape == (0, 0):  # here, the empty index can’t use str methods
			table = table_all
		else:
			table = table_all.loc[table_all.index.str.match('L.*'), table_all.columns.str.match('R.*')]
		table.to_csv(output[0], '\t')

rule plot_counts:
	input:
		'out/counts/{counting}/{lib_name}/{amplicon}/{lib_name}-{amplicon}{which}.tsv'
	output:
		'out/counts/{counting}/{lib_name}/{amplicon}/{lib_name}-{amplicon}{which}-log.png'
	run:
		counts_long = pd.read_csv(input[0], '\t')
		counts_na = counts_long.melt(['bc_l'], var_name='bc_r', value_name='Count')
		counts = counts_na[~counts_na.Count.isna()].copy()
		
		if counts.shape[0] <= 2:
			Path(output[0]).write_bytes(EMPTY_PNG)
			return
		
		gg = plot_counts(counts) + \
			theme(axis_text_x=element_text(size=4), axis_text_y=element_text(size=4))
		gg.save(output[0], dpi=300, verbose=False)

def amp_tables_to_counts(tables: Dict[str, pd.DataFrame]):
	table = pd.concat(tables, names=['Amplicon']).reset_index(0)
	counts_na = table.melt(['Amplicon', 'bc_l'], var_name='bc_r', value_name='Count')
	return counts_na[~counts_na.Count.isna()].copy()

rule plot_counts_lib:
	input:
		lambda wildcards: expand(
			'out/counts/{counting}/{lib_name}/{amplicon}/{lib_name}-{amplicon}{which}.tsv',
			counting=wildcards.counting,
			lib_name=wildcards.lib_name,
			amplicon=amplicons[wildcards.lib_name],
			which=wildcards.which,
		)
	output:
		'out/counts/{counting}/{lib_name}/{lib_name}{which}-log.png'
	run:
		tables = {Path(i).parent.name: pd.read_csv(i, '\t') for i in input}
		counts = amp_tables_to_counts(tables)
		
		gg = plot_counts(counts) + \
			facet_wrap('~Amplicon', ncol=4) + \
			theme(axis_text_x=element_text(size=1.8), axis_text_y=element_text(size=1.8))
		gg.save(output[0], dpi=300, verbose=False)

re_summary = re.compile(r'''HISAT2 summary stats:
	Total reads: (?P<total>\d+)
		Aligned 0 time: (?P<zero>\d+) \(\d+\.\d+%\)
		Aligned 1 time: (?P<one>\d+) \(\d+\.\d+%\)
		Aligned >1 times: (?P<more>\d+) \(\d+\.\d+%\)
	Overall alignment rate: \d+\.\d+%''')

rule spreadsheet:
	input:
		summaries = expand('process/4-mapped/{lib_name}_R{read}_summary.txt', lib_name=lib_names, read=[1,2]),
		counts = [
			path
			for lib_name in lib_names
			for path in expand('out/counts/{{counting}}/{lib_name}/{amplicon}/{lib_name}-{amplicon}.tsv', amplicon=amplicons[lib_name], lib_name=lib_name)
		],
	output:
		'out/counts/{counting}/{counting}.xlsx'
	run:
		import openpyxl
		from openpyxl.utils.dataframe import dataframe_to_rows
		
		re_name = re.compile(format('{re_lib_name}-{re_amplicon}\.tsv'))
		tables_all = {
			tuple(re_name.fullmatch(Path(i).name).groups()):
			pd.read_csv(i, '\t') for i in input.counts
		}
		
		wb = openpyxl.Workbook()
		wb.active.title = 'Statistics'
		
		stats = ['total', 'zero', 'one', 'more']
		wb.active.cell(1, 1).value = 'library'
		wb.active.cell(1, 2).value = 'read'
		for c, stat in enumerate(stats, 3):
			wb.active.cell(1, c).value = stat
		for r, path in enumerate(input.summaries, 2):
			path = Path(path)
			stats_txt = path.read_text('utf-8')
			stats_match = re_summary.match(stats_txt)
			wb.active.cell(r, 1).value = re.sub(r'_R[12]_summary.txt$', '', path.name)
			wb.active.cell(r, 2).value = int(re.fullmatch(r'.*_R([12])_summary.txt', path.name)[1])
			for c, stat in enumerate(stats, 3):
				wb.active.cell(r, c).value = stats_match[stat]
		
		for lib in lib_names:
			tables_lib = {amp: table for (l, amp), table in tables_all.items() if l == lib}
			counts = amp_tables_to_counts(tables_lib)
			sheet = counts.pivot_table('Count', ['bc_l', 'bc_r'], 'Amplicon').reset_index().fillna(0)
			
			ws = wb.create_sheet(lib)
			for r in dataframe_to_rows(sheet, index=False, header=True):
				if not (len(r) == 1 and r[0] is None):
					ws.append(r)
		wb.save(output[0])

#Needs https://bitbucket.org/snakemake/snakemake/pull-requests/264
rule dag:
	shell:
		'snakemake -s {__file__} --dag | dot -Tsvg | gwenview /dev/stdin'
