library (dplyr)

args = commandArgs(trailingOnly=TRUE)

input_file=args[1]
output_file=args[2]

Guppy_summary_results <- read.table(input_file, header = TRUE)
Guppy_summary_results <- Guppy_summary_results[,c("read_id","alignment_score","alignment_coverage","alignment_genome","alignment_identity")]
colnames(Guppy_summary_results) <- c("qseqid", "evalue", "qcovhsp", "salltitles", "pident")

Guppy_summary_results <- filter(Guppy_summary_results,evalue >0)
Guppy_summary_results$evalue <- (1/(Guppy_summary_results$evalue^2))

Guppy_summary_results$salltitles <- gsub("NR_"," ",Guppy_summary_results$salltitles)
Guppy_summary_results$salltitles <- sub(".*?_"," ",Guppy_summary_results$salltitles)
Guppy_summary_results$salltitles <- gsub("_"," ",Guppy_summary_results$salltitles)

Guppy_summary_results$qcovhsp <- round((Guppy_summary_results$qcovhsp*100), digits = 0)
Guppy_summary_results$pident <- round((Guppy_summary_results$pident*100), digits = 2)

write.table(Guppy_summary_results,file=output_file, quote = FALSE, sep ="\t", row.names = FALSE, col.names = FALSE, na ="")
