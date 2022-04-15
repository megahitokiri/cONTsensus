configfile: "config.yaml"

#--------------------------------------------------------------------------------
# TargetRule cONTsensus
#--------------------------------------------------------------------------------

rule cONTsensus:
	input:
		expand("OTU_Processing/Input_Files/{project}.fastq",project=config["project"]),
		expand("OTU_Processing/Trimmed_Fastq_Porechop/{project}.fastq.barcode_trimmed.porechop",project=config["project"]),
		expand("OTU_Processing/Trimmed_Fastq_Guppy_Barcoder/{project}.fastq.barcode_trimmed.guppy_barcoder",project=config["project"]),
		expand("OTU_Processing/ProwlerTrimmer_QF_{quality_filter}/{project}.ProwlerTrimmer.QF_{quality_filter}.fastq",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_Processing/Filtlong_QF_{quality_filter}/{project}.fastq.Filtlong.QF_{quality_filter}",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_Processing/BLAST/{project}.QF_{quality_filter}.fasta",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_Processing/BLAST/{project}.QF_{quality_filter}.db_16S_ribosomal_RNA.BLAST.unfilter_top_hits.txt",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_Processing/16S_ppm_results/Out16S_{project}.QF_{quality_filter}.db_16S_ribosomal_RNA.BLAST_species_counts.txt",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_Processing/GUPPY_ALIGNER_QF_{quality_filter}/{project}.ProwlerTrimmer.QF_{quality_filter}.bam",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_Processing/GUPPY_ALIGNER_QF_{quality_filter}/{project}.Guppy_aligner.QF_{quality_filter}.summary",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_Processing/16S_ppm_results/Out16S_{project}.QF_{quality_filter}.db_16S_ribosomal_RNA.GUPPY_ALIGNER_species_counts.txt",project=config["project"],quality_filter=config["quality_filter"]),
		
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
		expand("OTU_Processing/ProwlerTrimmer_QF_{quality_filter}/{project}.ProwlerTrimmer.QF_{quality_filter}.fastq",project=config["project"],quality_filter=config["quality_filter"]),
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
			mv $(find . -type f -name "*.fastq") {params.project}.ProwlerTrimmer.QF_$Q_filter.fastq
			
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
# BLAST_NCBI: Compares the sequences against NCBI db_16S_ribosomal_RNA blast database
#-----------------------------------------------------------------------------------		
rule BLAST_NCBI:
	input:
		rules.Filtlong.output,
	output:
		blast_fasta="OTU_Processing/BLAST/{project}.QF_{quality_filter}.fasta",
		db_16S_ribosomal_RNA="OTU_Processing/BLAST/{project}.QF_{quality_filter}.db_16S_ribosomal_RNA.BLAST.unfilter_top_hits.txt"
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
		{params.seqtk} seq -a OTU_Processing/Filtlong_QF_{wildcards.quality_filter}/{wildcards.project}.fastq.Filtlong.QF_{wildcards.quality_filter} > OTU_Processing/BLAST/{wildcards.project}.QF_{wildcards.quality_filter}.fasta

		echo INITIATING BLASTn IN DATABASE: {params.db_16S_ribosomal_RNA} ON FASTA FILE: {wildcards.project}.QF_{wildcards.quality_filter}.fasta
		{params.blastn} -db {params.db_16S_ribosomal_RNA} -query OTU_Processing/BLAST/{wildcards.project}.QF_{wildcards.quality_filter}.fasta -num_threads {params.n_threads} -outfmt '6 qseqid evalue qcovhsp salltitles pident' -max_target_seqs {params.max_target_seqs} -evalue {params.evalue} > OTU_Processing/BLAST/{wildcards.project}.QF_{wildcards.quality_filter}.db_16S_ribosomal_RNA.BLAST.unfilter_top_hits.txt
		"""

#------------------------------------------------------------------------------------------------------------------------
# Guppy_Aligner_NCBI: Compares the sequences against NCBI db_16S_ribosomal_RNA blast database (minimap2_based)
# sed -e 's/\s\+/_/g' 16S_ribosomal_RNA.fasta  to 16S_ribosomal_RNA.fasta.names_corrected
#------------------------------------------------------------------------------------------------------------------------
rule Guppy_Aligner_NCBI:
	input:
		rules.ProwlerTrimmer.output,
	output:
		guppy_16S_ribosomal_RNA="OTU_Processing/GUPPY_ALIGNER_QF_{quality_filter}/{project}.ProwlerTrimmer.QF_{quality_filter}.bam",
		summary_Guppy_aligner="OTU_Processing/GUPPY_ALIGNER_QF_{quality_filter}/{project}.Guppy_aligner.QF_{quality_filter}.summary"
	params:
		n_threads=config["n_threads"],
		project=config["project"],
		guppy_aligner=config["guppy_aligner"],	
		db_16S_ribosomal_RNA=config["db_16S_ribosomal_RNA"],
	shell:
		"""
		echo STARTING ALIGNING OF FASTQ VIA GUPPY_BARCODER IN DATABASE: {params.db_16S_ribosomal_RNA} ON FASTA FILE: ProwlerTrimmer_QF_{wildcards.quality_filter}
		{params.guppy_aligner} -i OTU_Processing/ProwlerTrimmer_QF_{wildcards.quality_filter}/ -s OTU_Processing/GUPPY_ALIGNER_QF_{wildcards.quality_filter} --align_ref {params.db_16S_ribosomal_RNA}.fasta.names_corrected --align_type auto -t {params.n_threads} --bam_out
		cp OTU_Processing/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/alignment_summary.txt OTU_Processing/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/{wildcards.project}.Guppy_aligner.QF_{wildcards.quality_filter}.summary
		"""

#---------------------------------------------------------------------------------------
# NCBI_16S_ppm: Transform blast and guppy_aligner results into a percentage content file
#----------------------------------------------------------------------------------------
rule NCBI_16S_ppm:
	input:
		BLAST_results=rules.BLAST_NCBI.output.db_16S_ribosomal_RNA,
	output:
		"OTU_Processing/16S_ppm_results/Out16S_{project}.QF_{quality_filter}.db_16S_ribosomal_RNA.BLAST_species_counts.txt",
	params:
		python_16S_ppm=config["16S_ppm"],
	shell:
		"""
		cp {input.BLAST_results} {input.BLAST_results}.bak
		python {params.python_16S_ppm} {input.BLAST_results} OTU_Processing/BLAST/
		mv {input.BLAST_results}.bak {input.BLAST_results}
		"""
		
#---------------------------------------------------------------------------------------
# GUPPY_ALIGNER_16S_ppm: Transform blast and guppy_aligner results into a percentage content file
#----------------------------------------------------------------------------------------
rule GUPPY_ALIGNER_16S_ppm:
	input:
		GUPPY_ALIGNER_results=rules.Guppy_Aligner_NCBI.output.summary_Guppy_aligner,
	output:
		"OTU_Processing/16S_ppm_results/Out16S_{project}.QF_{quality_filter}.db_16S_ribosomal_RNA.GUPPY_ALIGNER_species_counts.txt",
	params:
		python_16S_ppm=config["16S_ppm"],
	shell:
		"""
		BASEDIR=$PWD
		Rscript --vanilla scripts/guppy_aligner_parser.R $BASEDIR/{input} $BASEDIR/OTU_Processing/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/{wildcards.project}.QF_{wildcards.quality_filter}.db_16S_ribosomal_RNA.GUPPY_ALIGNER.unfilter_top_hits.txt
		python {params.python_16S_ppm} $BASEDIR/OTU_Processing/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/{wildcards.project}.QF_{wildcards.quality_filter}.db_16S_ribosomal_RNA.GUPPY_ALIGNER.unfilter_top_hits.txt OTU_Processing/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/
		"""		
