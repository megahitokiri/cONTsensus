library (dplyr)

args = commandArgs(trailingOnly=TRUE)

input_file=args[1]
output_file=input_file

BLAST_results <- read.table(input_file, header = FALSE)

colnames(BLAST_results) <- c("qseqid", "evalue", "qcovhsp", "salltitles", "pident")

BLAST_results$salltitles <- gsub("_"," ",BLAST_results$salltitles)

write.table(BLAST_results,file=output_file, quote = FALSE, sep ="\t", row.names = FALSE, col.names = FALSE, na ="")
