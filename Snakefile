# usage example: snakemake -d data/ngs15 -j 4

from snakemake.utils import min_version, listfiles


min_version('4.5.1')

read_file_names = [(w.readname, w.read) for _, w in listfiles('rawdata/{readname}_R{read,[12]}_001.fastq.gz')]

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
		get_read_paths('tagged'),
		'barcodes/barcodes.htm',

rule get_qc:
	input:
		'rawdata/{name_full}.fastq.gz'
	output:
		expand('qc/{{name_full}}_fastqc{suffix}', suffix=['.html', '.zip'])
	shell:
		'fastqc {input:q} -o qc'

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
		'barcodes/barcodes.txt'
	output:
		'barcodes/barcodes.htm'
	run:
		from bartseq.io import write_bc_table
		write_bc_table(input[0], output[0])

rule tag_reads:
	input:
		expand('trimmed/{{name}}_R{read}_001.fastq.gz', read=[1,2]),
		bc_file='barcodes/barcodes.txt',
	output:
		expand('tagged/{{name}}_R{read}_001.fastq.gz', read=[1,2]),
		stats_file='tagged/{name}_stats.json',
	run:
		from bartseq.main import run
		run(
			in_1=input[0], out_1=output[0],
			in_2=input[1], out_2=output[1],
			bc_file=input.bc_file,
			stats_file=output.stats_file,
		)

#Needs https://bitbucket.org/snakemake/snakemake/pull-requests/264
rule dag:
	shell:
		'snakemake -s {__file__} --dag | dot -Tsvg | gwenview /dev/stdin'
