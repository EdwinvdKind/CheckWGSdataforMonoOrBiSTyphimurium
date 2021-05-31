import yaml

SAMPLES = {}
with open("sample_sheet.yaml") as sample_sheet_file:
    SAMPLES = yaml.safe_load(sample_sheet_file) 

configfile: 'config/parameters.yaml'

rule all:
	input:
		'salmonella_multi_report.csv'

rule Seqsero2:
	input:
		r1 = lambda wildcards: SAMPLES[wildcards.sample]['R1'],
		r2 = lambda wildcards: SAMPLES[wildcards.sample]['R2']
	output:
		'{sample}_serotype/SeqSero_result.tsv'
	params:
		output_dir = "{sample}_serotype"
	shell:
		"SeqSero2_package.py -t 2 -i {input} -d {params.output_dir}"

rule checkamplicons:
	input:
		"{sample}_serotype/SeqSero_result.tsv",
		r1 = lambda wildcards: SAMPLES[wildcards.sample]['R1'],
		r2 = lambda wildcards: SAMPLES[wildcards.sample]['R2'],
		r3 = "amplicon/FFLIB_FFLIA.fasta",
		r4 = "amplicon/sense_59_antisense_83.fasta"		
	output:
		"{sample}_combinedresult.tsv"
	params:
		config["readthreshold"]["styphi"],
		config["readthreshold"]["biphasic"]
	threads:
		config["threads"]["checkamplicons"]
	script:
		"bin/checkamplicons.py"

rule salmonella_serotype_multireport:
	input:
		expand('{sample}_combinedresult.tsv', sample=SAMPLES),
	output:
		'salmonella_multi_report.csv'
	shell:
		'python bin/seqsero2_multireport.py -i {input} -o {output}'

