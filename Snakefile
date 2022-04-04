configfile: "config.yaml"

#--------------------------------------------------------------------------------
# TargetRule cONTsensus
#--------------------------------------------------------------------------------

rule cONTsensus:
	input:
		expand("OTU_Processing/Input_Files/{project}.fastq",project=config["project"]),
		expand("OTU_Processing/Trimmed_Fastq_Porechop/{project}.fastq.barcode_trimmed.porechop",project=config["project"]),
		expand("OTU_Processing/Trimmed_Fastq_Guppy_Barcoder/{project}.fastq.barcode_trimmed.guppy_barcoder",project=config["project"]),
		expand("OTU_Processing/ProwlerTrimmer_QF_{quality_filter}/{project}.fastq.ProwlerTrimmer.QF_{quality_filter}",project=config["project"],quality_filter=config["quality_filter"]),

#--------------------------------------------------------------------------------
# Init: Initializing files and folder
#--------------------------------------------------------------------------------
rule Init:
	input:
		raw_fastq=config["raw_fastq"],
	output:
		expand("OTU_Processing/Input_Files/{project}.fastq",project=config["project"]),
	params:
		project=config["project"],
	shell:
		"""
		snakemake --dag | dot -Tsvg > dag.svg
		cp {input.raw_fastq} OTU_Processing/Input_Files/{params.project}.fastq
		"""

#--------------------------------------------------------------------------------
# Porechop: Remove Adapters from nanopore fastq file method #1
#--------------------------------------------------------------------------------
rule Porechop:
	input:
		rules.Init.output,
	output:
		expand("OTU_Processing/Trimmed_Fastq_Porechop/{project}.fastq.barcode_trimmed.porechop",project=config["project"]),
	params:
		n_threads=config["n_threads"],
		project=config["project"],
		porechop_runner=config["porechop_runner"],	
	shell:
		"""
		echo STARTING TRIMMING OF FASTQ VIA PORECHOP...
		python {params.porechop_runner} -i {input} -o OTU_Processing/Trimmed_Fastq_Porechop/{params.project}.fastq.barcode_trimmed.porechop --threads={params.n_threads}
		"""

#--------------------------------------------------------------------------------
# Guppy_Barcoder: Remove Adapters from nanopore fastq file method #2
#--------------------------------------------------------------------------------
rule Guppy_Barcoder:
	input:
		rules.Init.output,
	output:
		expand("OTU_Processing/Trimmed_Fastq_Guppy_Barcoder/{project}.fastq.barcode_trimmed.guppy_barcoder",project=config["project"]),
	params:
		n_threads=config["n_threads"],
		project=config["project"],
		guppy_barcoder=config["guppy_barcoder"],	
	shell:
		"""
		echo STARTING TRIMMING OF FASTQ VIA GUPPY_BARCODER...
		{params.guppy_barcoder} -i OTU_Processing/Input_Files/ -s OTU_Processing/Trimmed_Fastq_Guppy_Barcoder --trim_barcodes --detect_barcodes 
		cat $(find OTU_Processing/Trimmed_Fastq_Guppy_Barcoder/ -type f -name "*.fastq") > OTU_Processing/Trimmed_Fastq_Guppy_Barcoder/{params.project}.fastq.barcode_trimmed.guppy_barcoder
		"""

#--------------------------------------------------------------------------------
# ProwlerTrimmer: Remove Adapters from nanopore fastq file method #3
#--------------------------------------------------------------------------------
rule ProwlerTrimmer:
	input:
		rules.Guppy_Barcoder.output,
	output:
		expand("OTU_Processing/ProwlerTrimmer_QF_{quality_filter}/{project}.fastq.ProwlerTrimmer.QF_{quality_filter}",project=config["project"],quality_filter=config["quality_filter"]),
	params:
		project=config["project"],
		ProwlerTrimmer=config["ProwlerTrimmer"],
		quality_filter=config["quality_filter"],
	shell:
		"""
		QF_LIST=$(echo {params.quality_filter})

		for Q_filter in $QF_LIST
			do

			echo STARTING FILTERING VIA PROWLER_TRIMMER AT QUALITY FILTER OF: $Q_filter ...
			python3 {params.ProwlerTrimmer} -f {params.project}.fastq.barcode_trimmed.guppy_barcoder -i OTU_Processing/Trimmed_Fastq_Guppy_Barcoder/ -o OTU_Processing/ProwlerTrimmer_QF_$Q_filter/{params.project} -w 200 -l 600 -g F1 -m D -q $Q_filter 

			cd OTU_Processing/ProwlerTrimmer_QF_$Q_filter
			mv $(find . -type f -name "*.fastq") {params.project}.fastq.ProwlerTrimmer.QF_$Q_filter
			
			cd ../..
			done

		"""