configfile: "config.yaml"

#--------------------------------------------------------------------------------
# TargetRule cONTsensus
#--------------------------------------------------------------------------------

rule cONTsensus:
	input:
		expand("OTU_{project}/Input_Files/{project}.fastq",project=config["project"]),
		expand("OTU_{project}/Trimmed_Fastq_Porechop/{project}.fastq.barcode_trimmed.porechop",project=config["project"]),
		expand("OTU_{project}/Trimmed_Fastq_Guppy_Barcoder/{project}.fastq.barcode_trimmed.guppy_barcoder",project=config["project"]),
		expand("OTU_{project}/ProwlerTrimmer_QF_{quality_filter}/{project}.ProwlerTrimmer.QF_{quality_filter}.fastq",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_{project}/Filtlong_QF_{quality_filter}/{project}.fastq.Filtlong.QF_{quality_filter}",project=config["project"],quality_filter=config["quality_filter"]),

		expand("OTU_{project}/BLAST/{project}.QF_{quality_filter}.fasta",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_{project}/BLAST/{project}.QF_{quality_filter}.db_16S_ribosomal_RNA.BLAST.unfilter_top_hits.txt",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_{project}/BLAST/{project}.QF_{quality_filter}.db_GGC.BLAST.unfilter_top_hits.txt",project=config["project"],quality_filter=config["quality_filter"]),

		expand("OTU_{project}/16S_ppm_results/Out16S_{project}.QF_{quality_filter}.db_16S_ribosomal_RNA.BLAST_species_counts.txt",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_{project}/GUPPY_ALIGNER_QF_{quality_filter}/{project}.ProwlerTrimmer.QF_{quality_filter}.bam",project=config["project"],quality_filter=config["quality_filter"]),
		
		expand("OTU_{project}/GUPPY_ALIGNER_QF_{quality_filter}/{project}.Guppy_aligner.QF_{quality_filter}.summary",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_{project}/GUPPY_ALIGNER_QF_{quality_filter}/{project}.Guppy_aligner.QF_{quality_filter}.db_GGC.summary",project=config["project"],quality_filter=config["quality_filter"]),
		
		expand("OTU_{project}/16S_ppm_results/Out16S_{project}.QF_{quality_filter}.db_16S_ribosomal_RNA.GUPPY_ALIGNER_species_counts.txt",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_{project}/16S_ppm_results/Out16S_{project}.QF_{quality_filter}.db_GGC.BLAST_species_counts.txt",project=config["project"],quality_filter=config["quality_filter"]),		
		expand("OTU_{project}/16S_ppm_results/Out16S_{project}.QF_{quality_filter}.db_GGC.GUPPY_ALIGNER_species_counts.txt",project=config["project"],quality_filter=config["quality_filter"]),
		
		expand("OTU_{project}/OTU_Analysis/QF_{quality_filter}_All_Combined_species_counts.txt",project=config["project"],quality_filter=config["quality_filter"]),
		expand("OTU_{project}/OTU_Analysis/QF_{quality_filter}_All_Combined_species_counts.expanded",project=config["project"],quality_filter=config["quality_filter"]),
		
#--------------------------------------------------------------------------------
# Init: Initializing files and folder
#--------------------------------------------------------------------------------
rule Init:
	input:
		raw_fastq=config["raw_fastq"],
	output:
		expand("OTU_{project}/Input_Files/{project}.fastq",project=config["project"]),
	params:
		project=config["project"],
	shell:
		"""
		snakemake --dag | dot -Tsvg > dag.svg
		cp {input.raw_fastq} OTU_{params.project}/Input_Files/{params.project}.fastq
		"""

#--------------------------------------------------------------------------------
# Porechop: Remove Adapters from nanopore fastq file method #1
#--------------------------------------------------------------------------------
rule Porechop:
	input:
		rules.Init.output,
	output:
		expand("OTU_{project}/Trimmed_Fastq_Porechop/{project}.fastq.barcode_trimmed.porechop",project=config["project"]),
	params:
		n_threads=config["n_threads"],
		project=config["project"],
		porechop_runner=config["porechop_runner"],	
	shell:
		"""
		echo STARTING TRIMMING OF FASTQ VIA PORECHOP...
		python {params.porechop_runner} -i {input} -o OTU_{params.project}/Trimmed_Fastq_Porechop/{params.project}.fastq.barcode_trimmed.porechop --threads={params.n_threads}
		"""

#--------------------------------------------------------------------------------
# Guppy_Barcoder: Remove Adapters from nanopore fastq file method #2
#--------------------------------------------------------------------------------
rule Guppy_Barcoder:
	input:
		rules.Init.output,
	output:
		expand("OTU_{project}/Trimmed_Fastq_Guppy_Barcoder/{project}.fastq.barcode_trimmed.guppy_barcoder",project=config["project"]),
	params:
		n_threads=config["n_threads"],
		project=config["project"],
		guppy_barcoder=config["guppy_barcoder"],	
	shell:
		"""
		echo STARTING TRIMMING OF FASTQ VIA GUPPY_BARCODER...
		{params.guppy_barcoder} -i OTU_{params.project}/Input_Files/ -s OTU_{params.project}/Trimmed_Fastq_Guppy_Barcoder --trim_barcodes --detect_barcodes 
		cat $(find OTU_{params.project}/Trimmed_Fastq_Guppy_Barcoder/ -type f -name "*.fastq") > OTU_{params.project}/Trimmed_Fastq_Guppy_Barcoder/{params.project}.fastq.barcode_trimmed.guppy_barcoder
		"""

#-----------------------------------------------------------------------------------
# ProwlerTrimmer: Asses Quality filtering using window approach after Guppy_Barcoder
#-----------------------------------------------------------------------------------
rule ProwlerTrimmer:
	input:
		rules.Guppy_Barcoder.output,
	output:
		expand("OTU_{project}/ProwlerTrimmer_QF_{quality_filter}/{project}.ProwlerTrimmer.QF_{quality_filter}.fastq",project=config["project"],quality_filter=config["quality_filter"]),
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
			python3 {params.ProwlerTrimmer} -f {params.project}.fastq.barcode_trimmed.guppy_barcoder -i OTU_{params.project}/Trimmed_Fastq_Guppy_Barcoder/ -o OTU_{params.project}/ProwlerTrimmer_QF_$Q_filter/{params.project} -w 200 -l {params.L_filt} -g F1 -m D -q $Q_filter 

			cd OTU_{params.project}/ProwlerTrimmer_QF_$Q_filter
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
		expand("OTU_{project}/Filtlong_QF_{quality_filter}/{project}.fastq.Filtlong.QF_{quality_filter}",project=config["project"],quality_filter=config["quality_filter"]),
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
			{params.Filtlong} OTU_{params.project}/Trimmed_Fastq_Porechop/{params.project}.fastq.barcode_trimmed.porechop --min_mean_q $vQF --min_length {params.L_filt} > OTU_{params.project}/Filtlong_QF_$Q_filter/{params.project}.fastq.Filtlong.QF_$Q_filter
			
			done

		"""

#---------------------------------------------------------------------------------------
# BLAST: Compares the sequences against NCBI db_16S_ribosomal_RNA blast database and GCC
#---------------------------------------------------------------------------------------		
rule BLAST:
	input:
		rules.Filtlong.output,
	output:
		blast_fasta="OTU_{project}/BLAST/{project}.QF_{quality_filter}.fasta",
		db_16S_ribosomal_RNA="OTU_{project}/BLAST/{project}.QF_{quality_filter}.db_16S_ribosomal_RNA.BLAST.unfilter_top_hits.txt",
		db_GGC="OTU_{project}/BLAST/{project}.QF_{quality_filter}.db_GGC.BLAST.unfilter_top_hits.txt"
	params:
		project=config["project"],
		n_threads=config["n_threads"],
		quality_filter=config["quality_filter"],
		seqtk=config["seqtk"],
		blastn=config["blastn"],
		db_16S_ribosomal_RNA=config["db_16S_ribosomal_RNA"],
		db_GGC=config["db_ggc"],
		max_target_seqs=config["max_target_seqs"],
		evalue=config["evalue"],
	shell:
		"""
		echo CONVERTING FASTQ TO FASTA FILE...
		{params.seqtk} seq -a OTU_{params.project}/Filtlong_QF_{wildcards.quality_filter}/{wildcards.project}.fastq.Filtlong.QF_{wildcards.quality_filter} > OTU_{params.project}/BLAST/{wildcards.project}.QF_{wildcards.quality_filter}.fasta

		echo INITIATING BLASTn IN DATABASE: {params.db_16S_ribosomal_RNA} ON FASTA FILE: {wildcards.project}.QF_{wildcards.quality_filter}.fasta
		{params.blastn} -db {params.db_16S_ribosomal_RNA} -query OTU_{params.project}/BLAST/{wildcards.project}.QF_{wildcards.quality_filter}.fasta -num_threads {params.n_threads} -outfmt '6 qseqid evalue qcovhsp salltitles pident' -max_target_seqs {params.max_target_seqs} -evalue {params.evalue} > OTU_{params.project}/BLAST/{wildcards.project}.QF_{wildcards.quality_filter}.db_16S_ribosomal_RNA.BLAST.unfilter_top_hits.txt

		echo INITIATING BLASTn IN DATABASE: {params.db_GGC} ON FASTA FILE: {wildcards.project}.QF_{wildcards.quality_filter}.fasta
		{params.blastn} -db {params.db_GGC}.fasta -query OTU_{params.project}/BLAST/{wildcards.project}.QF_{wildcards.quality_filter}.fasta -num_threads {params.n_threads} -outfmt '6 qseqid evalue qcovhsp salltitles pident' -max_target_seqs {params.max_target_seqs} -evalue {params.evalue} > OTU_{params.project}/BLAST/{wildcards.project}.QF_{wildcards.quality_filter}.db_GGC.BLAST.unfilter_top_hits.txt

		"""

#------------------------------------------------------------------------------------------------------------------------
# Guppy_Aligner_NCBI: Compares the sequences against NCBI db_16S_ribosomal_RNA blast database (minimap2_based)
# sed -e 's/\s\+/_/g' 16S_ribosomal_RNA.fasta  to 16S_ribosomal_RNA.fasta.names_corrected
#------------------------------------------------------------------------------------------------------------------------
rule Guppy_Aligner:
	input:
		rules.ProwlerTrimmer.output,
	output:
		guppy_16S_ribosomal_RNA="OTU_{project}/GUPPY_ALIGNER_QF_{quality_filter}/{project}.ProwlerTrimmer.QF_{quality_filter}.bam",
		summary_Guppy_aligner="OTU_{project}/GUPPY_ALIGNER_QF_{quality_filter}/{project}.Guppy_aligner.QF_{quality_filter}.summary",
		summary_Guppy_aligner_GGC="OTU_{project}/GUPPY_ALIGNER_QF_{quality_filter}/{project}.Guppy_aligner.QF_{quality_filter}.db_GGC.summary",
	params:
		project=config["project"],
		n_threads=config["n_threads"],
		guppy_aligner=config["guppy_aligner"],	
		db_16S_ribosomal_RNA=config["db_16S_ribosomal_RNA"],
		db_GGC=config["db_ggc"],
	shell:
		"""
		echo STARTING ALIGNING OF FASTQ VIA GUPPY_BARCODER IN DATABASE: {params.db_16S_ribosomal_RNA} ON FASTA FILE: ProwlerTrimmer_QF_{wildcards.quality_filter}
		{params.guppy_aligner} -i OTU_{params.project}/ProwlerTrimmer_QF_{wildcards.quality_filter}/ -s OTU_{params.project}/GUPPY_ALIGNER_QF_{wildcards.quality_filter} --align_ref {params.db_16S_ribosomal_RNA}.fasta.names_corrected --align_type auto -t {params.n_threads} --bam_out
		cp OTU_{params.project}/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/alignment_summary.txt OTU_{params.project}/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/{wildcards.project}.Guppy_aligner.QF_{wildcards.quality_filter}.summary

		rm OTU_{params.project}/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/alignment_summary.txt
		
		echo STARTING ALIGNING OF FASTQ VIA GUPPY_BARCODER IN DATABASE: {params.db_GGC} ON FASTA FILE: ProwlerTrimmer_QF_{wildcards.quality_filter}
		{params.guppy_aligner} -i OTU_{params.project}/ProwlerTrimmer_QF_{wildcards.quality_filter}/ -s OTU_{params.project}/GUPPY_ALIGNER_QF_{wildcards.quality_filter} --align_ref {params.db_GGC}.fasta --align_type auto -t {params.n_threads}
		cp OTU_{params.project}/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/alignment_summary.txt OTU_{params.project}/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/{wildcards.project}.Guppy_aligner.QF_{wildcards.quality_filter}.db_GGC.summary
		"""

#---------------------------------------------------------------------------------------
# BLAST_16S_ppm: Transform blast and guppy_aligner results into a percentage content file
#----------------------------------------------------------------------------------------
rule BLAST_16S_ppm:
	input:
		BLAST_results_NCBI=rules.BLAST.output.db_16S_ribosomal_RNA,
		BLAST_results_GGC=rules.BLAST.output.db_GGC,
	output:
		NCBI_species="OTU_{project}/16S_ppm_results/Out16S_{project}.QF_{quality_filter}.db_16S_ribosomal_RNA.BLAST_species_counts.txt",
		GGC_species="OTU_{project}/16S_ppm_results/Out16S_{project}.QF_{quality_filter}.db_GGC.BLAST_species_counts.txt",
	params:
		project=config["project"],
		python_16S_ppm=config["16S_ppm"],
	shell:
		"""
		BASEDIR=$PWD
		echo STARTING 16S_ppm ON: {input.BLAST_results_NCBI}
		cp {input.BLAST_results_NCBI} {input.BLAST_results_NCBI}.bak
		python {params.python_16S_ppm} {input.BLAST_results_NCBI} OTU_{params.project}/BLAST/ OTU_{params.project}/16S_ppm_results
		mv {input.BLAST_results_NCBI}.bak {input.BLAST_results_NCBI}
		
		echo STARTING 16S_ppm ON: {input.BLAST_results_GGC}
		cp {input.BLAST_results_GGC} {input.BLAST_results_GGC}.bak
		Rscript --vanilla scripts/GGC_BLAST_name_correction.R $BASEDIR/{input.BLAST_results_GGC}
		python {params.python_16S_ppm} {input.BLAST_results_GGC} OTU_{params.project}/BLAST/ OTU_{params.project}/16S_ppm_results
		mv {input.BLAST_results_GGC}.bak {input.BLAST_results_GGC}
		"""
		
#-------------------------------------------------------------------------------------------------
# GUPPY_ALIGNER_16S_ppm: Transform blast and guppy_aligner results into a percentage content file
#-------------------------------------------------------------------------------------------------
rule GUPPY_ALIGNER_16S_ppm:
	input:
		GUPPY_ALIGNER_results=rules.Guppy_Aligner.output.summary_Guppy_aligner,
		GUPPY_ALIGNER_results_GGC=rules.Guppy_Aligner.output.summary_Guppy_aligner_GGC,
	output:
		NCBI_species="OTU_{project}/16S_ppm_results/Out16S_{project}.QF_{quality_filter}.db_16S_ribosomal_RNA.GUPPY_ALIGNER_species_counts.txt",
		GGC_species="OTU_{project}/16S_ppm_results/Out16S_{project}.QF_{quality_filter}.db_GGC.GUPPY_ALIGNER_species_counts.txt",

	params:
		project=config["project"],
		python_16S_ppm=config["16S_ppm"],
	shell:
		"""
		BASEDIR=$PWD
		echo STARTING 16S_ppm ON: {input.GUPPY_ALIGNER_results}
		Rscript --vanilla scripts/guppy_aligner_parser.R $BASEDIR/{input.GUPPY_ALIGNER_results} $BASEDIR/OTU_{params.project}/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/{wildcards.project}.QF_{wildcards.quality_filter}.db_16S_ribosomal_RNA.GUPPY_ALIGNER.unfilter_top_hits.txt
		python {params.python_16S_ppm} $BASEDIR/OTU_{params.project}/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/{wildcards.project}.QF_{wildcards.quality_filter}.db_16S_ribosomal_RNA.GUPPY_ALIGNER.unfilter_top_hits.txt OTU_{params.project}/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/ OTU_{params.project}/16S_ppm_results
		
		echo STARTING 16S_ppm ON: {input.GUPPY_ALIGNER_results_GGC}
		Rscript --vanilla scripts/guppy_aligner_parser_GGC.R $BASEDIR/{input.GUPPY_ALIGNER_results_GGC} $BASEDIR/OTU_{params.project}/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/{wildcards.project}.QF_{wildcards.quality_filter}.db_GGC.GUPPY_ALIGNER.unfilter_top_hits.txt
		Rscript --vanilla scripts/GGC_BLAST_name_correction.R $BASEDIR/OTU_{params.project}/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/{wildcards.project}.QF_{wildcards.quality_filter}.db_GGC.GUPPY_ALIGNER.unfilter_top_hits.txt
		python {params.python_16S_ppm} $BASEDIR/OTU_{params.project}/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/{wildcards.project}.QF_{wildcards.quality_filter}.db_GGC.GUPPY_ALIGNER.unfilter_top_hits.txt OTU_{params.project}/GUPPY_ALIGNER_QF_{wildcards.quality_filter}/ OTU_{params.project}/16S_ppm_results
		"""		
		
#-------------------------------------------------------------------------------------------------
# OTU_COMBINER: Combine all the Species_count results into one expanded file
#-------------------------------------------------------------------------------------------------
rule OTU_COMBINER_CLEANER:
	input:
		rules.GUPPY_ALIGNER_16S_ppm.output,
		rules.BLAST_16S_ppm.output,
	output:
		All_Combined_Species_count="OTU_{project}/OTU_Analysis/QF_{quality_filter}_All_Combined_species_counts.txt",
		Expanded_Species_count="OTU_{project}/OTU_Analysis/QF_{quality_filter}_All_Combined_species_counts.expanded",
	params:
		project=config["project"],
	shell:
		"""
		BASEDIR=$PWD
		cat $(find OTU_{params.project}/16S_ppm_results -type f -name "*QF_{wildcards.quality_filter}*species_counts.txt") > OTU_{params.project}/OTU_Analysis/QF_{wildcards.quality_filter}_All_Combined_species_counts.txt
		sed -i '/^#/d' OTU_{params.project}/OTU_Analysis/QF_{wildcards.quality_filter}_All_Combined_species_counts.txt
		Rscript --vanilla scripts/OTU_expander.R $BASEDIR/OTU_{params.project}/OTU_Analysis/QF_{wildcards.quality_filter}_All_Combined_species_counts.txt $BASEDIR/OTU_{params.project}/OTU_Analysis/QF_{wildcards.quality_filter}_All_Combined_species_counts.expanded {wildcards.quality_filter}
		"""				