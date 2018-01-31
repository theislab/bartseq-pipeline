configfile: 'Snakefile-example-config.json'

from pathlib import Path

from snakemake.utils import listfiles

dataset = Path(config["dataset"])

read_file_names = [(w.readname, w.read) for _, w in listfiles(f'{dataset}/rawdata/{{readname}}_R{{read,[12]}}_001.fastq.gz')]

def get_read_path(prefix, name, read, suffix='.fastq.gz'):
	return f'{dataset}/{prefix}/{name}_R{read}_001{suffix}'

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

rule get_qc:
	input:
		f'{dataset}/rawdata/{{name_full}}.fastq.gz',
	output:
		expand(str(dataset / 'qc/{name_full}_fastqc{suffix}'), suffix=['.html', '.zip']),

