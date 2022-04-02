configfile: "config.yaml"

#--------------------------------------------------------------------------------
# TargetRule cONTsensus
#--------------------------------------------------------------------------------

rule cONTsensus:
	input:
		expand("OTU_Processing/{project}.fastq",project=config["project"]),
		expand("OTU_Processing/{project}.fastq.barcode_trimmed.porechop",project=config["project"]),


#--------------------------------------------------------------------------------
# Init: Initializing files and folder
#--------------------------------------------------------------------------------
rule Init:
	input:
		raw_fastq=config["raw_fastq"],
	output:
		expand("OTU_Processing/{project}.fastq",project=config["project"]),
	params:
		project=config["project"],
	shell:
		"""
		cp {input.raw_fastq} OTU_Processing/{params.project}.fastq
		"""

#--------------------------------------------------------------------------------
# Porechop: Remove Adapters from nanopore fastq file
#--------------------------------------------------------------------------------
rule Porechop:
	input:
		rules.Init.output,
	output:
		expand("OTU_Processing/{project}.fastq.barcode_trimmed.porechop",project=config["project"]),
	params:
		n_threads=config["n_threads"],
		project=config["project"],
		porechop_runner=config["porechop_runner"],	
	shell:
		"""
		echo STARTING TRIMMING OF FASTQ VIA PORECHOP...
		python {params.porechop_runner} -i {input} -o OTU_Processing/{params.project}.fastq.barcode_trimmed.porechop --threads={params.n_threads}
		"""
