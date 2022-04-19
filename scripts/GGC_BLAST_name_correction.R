library (dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

input_file=args[1]
output_file=input_file

BLAST_results <- read.table(input_file, header = FALSE)

colnames(BLAST_results) <- c("qseqid", "evalue", "qcovhsp", "salltitles", "pident")

BLAST_results$salltitles <- sub("_"," ",BLAST_results$salltitles)
BLAST_results$Word_count <- str_count(BLAST_results$salltitles,"\\w+")
BLAST_results$salltitles <-ifelse(BLAST_results$Word_count==1,
                                  paste0(BLAST_results$salltitles," ",BLAST_results$salltitles),
                                  BLAST_results$salltitles)

BLAST_results$word1 <- word(BLAST_results$salltitles,1)
BLAST_results$word2 <- word(BLAST_results$salltitles,2)

BLAST_results$salltitles <-ifelse(BLAST_results$word1=="Candidatus",
                                  gsub("_"," ",BLAST_results$salltitles),
                                  BLAST_results$salltitles)

BLAST_results$Word_count <- str_count(BLAST_results$salltitles,"\\w+")

BLAST_results$Candidatus_fix <-ifelse(BLAST_results$Word_count==2&BLAST_results$word1=="Candidatus",
                                  TRUE,
                                  FALSE)
BLAST_results$Candidatus_fix <-ifelse(BLAST_results$Candidatus_fix==TRUE,
                                      paste0(BLAST_results$salltitles," ",BLAST_results$word2),
                                      BLAST_results$salltitles)

BLAST_results$salltitles <- BLAST_results$Candidatus_fix
BLAST_results$salltitles <- gsub("_"," ",BLAST_results$salltitles)

BLAST_results$Word_count <-NULL
BLAST_results$word1 <- NULL
BLAST_results$word2 <- NULL
BLAST_results$Candidatus_fix <- NULL

write.table(BLAST_results,file=output_file, quote = FALSE, sep ="\t", row.names = FALSE, col.names = FALSE, na ="")
