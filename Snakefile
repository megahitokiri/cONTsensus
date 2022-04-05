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
		expand("OTU_Processing/Filtlong_QF_{quality_filter}/{project}.fastq.Filtlong.QF_{quality_filter}",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_Processing/BLAST/{project}.QF_{quality_filter}.fasta",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_Processing/BLAST/{project}.QF_{quality_filter}.db_16S_ribosomal_RNA.unfilter_top_hits.txt",project=config["project"],quality_filter=config["quality_filter"]),

		
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

#-----------------------------------------------------------------------------------
# ProwlerTrimmer: Asses Quality filtering using window approach after Guppy_Barcoder
#-----------------------------------------------------------------------------------
rule ProwlerTrimmer:
	input:
		rules.Guppy_Barcoder.output,
	output:
		expand("OTU_Processing/ProwlerTrimmer_QF_{quality_filter}/{project}.fastq.ProwlerTrimmer.QF_{quality_filter}",project=config["project"],quality_filter=config["quality_filter"]),
	params:
		project=config["project"],
		ProwlerTrimmer=config["ProwlerTrimmer"],
		quality_filter=config["quality_filter"],
		L_filt=config["L_filt"],
	shell:
		"""
		QF_LIST=$(echo {params.quality_filter})

		for Q_filter in $QF_LIST
			do

			echo STARTING FILTERING VIA PROWLER_TRIMMER AT QUALITY FILTER OF: $Q_filter ...
			python3 {params.ProwlerTrimmer} -f {params.project}.fastq.barcode_trimmed.guppy_barcoder -i OTU_Processing/Trimmed_Fastq_Guppy_Barcoder/ -o OTU_Processing/ProwlerTrimmer_QF_$Q_filter/{params.project} -w 200 -l {params.L_filt} -g F1 -m D -q $Q_filter 

			cd OTU_Processing/ProwlerTrimmer_QF_$Q_filter
			mv $(find . -type f -name "*.fastq") {params.project}.fastq.ProwlerTrimmer.QF_$Q_filter
			
			cd ../..
			done

		"""
		
#-----------------------------------------------------------------------------------
# Filtlong: Asses Quality filtering using window approach after Porechop
#-----------------------------------------------------------------------------------		
rule Filtlong:
	input:
		rules.Porechop.output,
	output:
		expand("OTU_Processing/Filtlong_QF_{quality_filter}/{project}.fastq.Filtlong.QF_{quality_filter}",project=config["project"],quality_filter=config["quality_filter"]),
	params:
		project=config["project"],
		Filtlong=config["Filtlong"],
		quality_filter=config["quality_filter"],
		L_filt=config["L_filt"],
	shell:
		"""
		QF_LIST=$(echo {params.quality_filter})

		for Q_filter in $QF_LIST
			do
			vQF=$(python -c "Result=(1-(10**(-(float($Q_filter)/10))))*100; print(Result)")
			
			echo STARTING FILTERING VIA PROWLER_TRIMMER AT QUALITY FILTER OF $Q_filter: corresponding to vQ of: $vQF ...
			{params.Filtlong} OTU_Processing/Trimmed_Fastq_Porechop/{params.project}.fastq.barcode_trimmed.porechop --min_mean_q $vQF --min_length {params.L_filt} > OTU_Processing/Filtlong_QF_$Q_filter/{params.project}.fastq.Filtlong.QF_$Q_filter
			
			done

		"""

#-----------------------------------------------------------------------------------
# BLAST: Compares the sequences against a fasta or blast database
#-----------------------------------------------------------------------------------		
rule BLAST:
	input:
		rules.Filtlong.output,
	output:
		blast_fasta="OTU_Processing/BLAST/{project}.QF_{quality_filter}.fasta",
		db_16S_ribosomal_RNA="OTU_Processing/BLAST/{project}.QF_{quality_filter}.db_16S_ribosomal_RNA.unfilter_top_hits.txt"
	params:
		n_threads=config["n_threads"],
		quality_filter=config["quality_filter"],
		seqtk=config["seqtk"],
		blastn=config["blastn"],
		db_16S_ribosomal_RNA=config["db_16S_ribosomal_RNA"],
		max_target_seqs=config["max_target_seqs"],
		evalue=config["evalue"],
	shell:
		"""
		echo CONVERTING FASTQ TO FASTA FILE...
		{params.seqtk} seq -a {input} > OTU_Processing/BLAST/{wildcards.project}.QF_{wildcards.quality_filter}.fasta

		echo INITIATING BLASTn IN DATABASE: {params.db_16S_ribosomal_RNA} ON FASTA FILE: {wildcards.project}.QF_{wildcards.quality_filter}.fasta
		{params.blastn} -db {params.db_16S_ribosomal_RNA} -query OTU_Processing/BLAST/{wildcards.project}.QF_{wildcards.quality_filter}.fasta -num_threads {params.n_threads} -outfmt '6 qseqid evalue qcovhsp salltitles pident' -max_target_seqs {params.max_target_seqs} -evalue {params.evalue} > OTU_Processing/BLAST/{wildcards.project}.QF_{wildcards.quality_filter}.db_16S_ribosomal_RNA.unfilter_top_hits.txt
		"""
		