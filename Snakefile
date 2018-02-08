# usage example: snakemake -d data/ngs15 -j 4
import re
from collections import Counter

from snakemake.utils import min_version, listfiles
import pandas as pd
import seaborn as sns
from tqdm import tqdm

from bartseq.io import write_bc_table, transparent_open


min_version('4.5.1')

read_file_names = [(w.readname, w.read) for _, w in listfiles('rawdata/{readname}_R{read,[12]}_001.fastq.gz')]
lib_names = [readname for readname, read in read_file_names if read == '1']
count_files = expand('counts/counts{which}{ext}', which=['', '-all'], ext=['.tsv', '-log.svg'])

def get_read_path(prefix, name, read, suffix='.fastq.gz'):
	return '{prefix}/{name}_R{read}_001{suffix}'.format_map(locals())

def get_read_paths(prefix, *suffixes):
	if len(suffixes) == 0:
		suffixes = ['.fastq.gz']
	return [
		get_read_path(prefix, n, r, suffix)
		for n, r in read_file_names
		for suffix in suffixes
	]

reads_raw = get_read_paths('rawdata')


rule all:
	input:
		get_read_paths('qc', '_fastqc.html', '_fastqc.zip'),
		count_files,
		'barcodes/barcodes.htm',

rule get_qc:
	input:
		'rawdata/{name_full}.fastq.gz'
	output:
		expand('qc/{{name_full}}_fastqc{suffix}', suffix=['.html', '.zip'])
	shell:
		'fastqc {input:q} -o qc'

rule get_read_count:
	input:
		'rawdata/{name}_R1_001.fastq.gz'
	output:
		'rawdata/{name}_001.count.txt'
	shell:
		'''
		lines=$(zcat {input[0]:q} | wc -l)
		echo "$lines / 4" | bc > {output[0]:q}
		'''

rule trim_quality:
	input:
		expand('rawdata/{{name}}_R{read}_001.fastq.gz', read=[1,2]),
	output:
		expand('trimmed/{{name}}_R{read}_001.fastq.gz', read=[1,2]),
		single='trimmed/{name}_single_001.fastq.gz',
	shell:
		'''
		sickle pe --gzip-output --qual-type=sanger \
			--pe-file1={input[0]:q} --output-pe1={output[0]:q} \
			--pe-file2={input[1]:q} --output-pe2={output[1]:q} \
			--output-single={output.single:q}
		'''

rule bc_table:
	input:
		'barcodes/barcodes.fa'
	output:
		'barcodes/barcodes.htm'
	run:
		write_bc_table(input[0], output[0])

rule tag_reads:
	input:
		expand('trimmed/{{name}}_R{read}_001.fastq.gz', read=[1,2]),
		count_file='rawdata/{name}_001.count.txt',
		bc_file='barcodes/barcodes.fa',
	output:
		expand('tagged/{{name}}_R{read}_001.fastq.gz', read=[1,2]),
		stats_file='tagged/{name}_stats.json',
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

rule amplicon_fa:
	input:
		'amplicons/amplicons.txt'
	output:
		'amplicons/amplicons.fa'
	shell:
		r"cat {input} | tail -n +2 | cut -f 1,4 | sed 's/^ */>/;s/\t/\n/' > {output}"

rule build_index:
	input:
		'amplicons/amplicons.fa'
	output:
		idx = expand('amplicons/amplicons.{n}.ht2', n=range(1, 9)),
	threads: 4
	shell:
		r'''
		hisat2-build -p {threads} {input} amplicons/amplicons
		'''

rule map_reads:
	input:
		amplicons = expand('amplicons/amplicons.{n}.ht2', n=range(1, 9)),
		read = 'tagged/{name_full}.fastq.gz',
	output:
		map = 'mapped/{name_full}.txt',
		summary = 'mapped/{name_full}_summary.txt'
	threads: 4
	shell:
		'''
		hisat2 \
			--threads {threads} \
			--reorder \
			-k 1 \
			-x amplicons/amplicons \
			--new-summary --summary-file {output.summary} \
			-q -U {input.read} | \
			grep -v "^@" - | \
			cut -f3 > {output.map}
		'''

rule count:
	input:
		reads = expand('tagged/{{name}}_R{read}_001.fastq.gz', read=[1,2]),
		mappings = expand('mapped/{{name}}_R{read}_001.txt', read=[1,2]),
		count_file = 'rawdata/{name}_001.count.txt',
	output:
		'counts/{name}_001.tsv'
	run:
		bc_re = re.compile(r'barcode=(\w+)')
		with open(input.count_file) as c_f:
			total = int(c_f.read())
		def get_barcodes(fastq_file):
			for header in fastq_file:
				next(fastq_file)
				next(fastq_file)
				next(fastq_file)
				yield bc_re.search(header).group(1)
		counts = Counter()
		with \
			transparent_open(input.reads[0]) as r1, open(input.mappings[0]) as m1, \
			transparent_open(input.reads[1]) as r2, open(input.mappings[1]) as m2:
			
			bcs1 = get_barcodes(r1)
			bcs2 = get_barcodes(r2)
			for bc1, bc2, amp1, amp2 in tqdm(zip(bcs1, bcs2, m1, m2), total=total):
				bc1, bc2 = sorted([bc1, bc2])
				if amp1 == amp2:
					counts[bc1, bc2] += 1
			
			with open(output[0], 'w') as of:
				print('bc_l', 'bc_r', 'count', sep='\t', file=of)
				for (bc_l, bc_r), c in counts.items():
					print(bc_l, bc_r, c, sep='\t', file=of)

rule combine_counts:
	input:
		expand('counts/{name}_001.tsv', name=lib_names)
	output:
		count_files
	run:
		entries_all = pd.concat([pd.read_csv(f, '\t') for f in input])
		entries_useful = entries_all[
			entries_all.bc_l.str.match('L.*') &
			entries_all.bc_r.str.match('R.*')]
		table_useful = entries_useful.pivot('bc_l', 'bc_r', 'count')
		table_all = entries_all.pivot('bc_l', 'bc_r', 'count')
		for o in output:
			table: pd.DataFrame = table_all if '-all' in o else table_useful
			if o.endswith('.svg'):
				plot = sns.heatmap(table.transform(pd.np.log1p))
				plot.set(xlabel='', ylabel='')
				plot.get_figure().savefig(o)
			else:
				table.to_csv(o, '\t')


#Needs https://bitbucket.org/snakemake/snakemake/pull-requests/264
rule dag:
	shell:
		'snakemake -s {__file__} --dag | dot -Tsvg | gwenview /dev/stdin'
