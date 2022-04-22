library("ape")
library("seqinr")
library(Rsamtools)
library(bamsignals)
library(GenomicRanges)
library(Biostrings)
library("phyloseq")
library(dplyr)
library(tidyr)
library(stringr)
library("ggplot2")
library("ggtree") 
library("ggimage")
library("treeio")
library("tidytree")
library(ggtreeExtra)
library(ggstar)
library("stargazer")
#BiocManager::install("msa")
library(msa)

setwd("/DATA/home/jmlazaro/Projects/Microbiome/Ref")

Assembly_Microbiome <- readDNAStringSet("gg_13_5.fasta.gz")

setwd("/DATA/home/jmlazaro/Projects/cONTsensus/databases/Green_genes_consensus")

tax.clean <- readRDS(file = "tax.clean.rds")
tax.clean$gg_id <- rownames(tax.clean)

Genbank.info <- readRDS(file = "Microbiome_Assembly_Genbank.rds")
Genbank.info$gg_id <- NULL
Genbank.info$accession_type <- NULL
Genbank.info$accession <- NA
Genbank.info$accession <- sub(" .*","",Genbank.info$Description)

Genbank.info <- filter(Genbank.info,is.na(accession)==FALSE)

Microbiome_Assembly_Genbank <- read.table("gg_13_5_accessions.txt", header = TRUE)

Curated_Genbank_info <- merge(Genbank.info,Microbiome_Assembly_Genbank, all.x=TRUE, by = "accession" )
Curated_Genbank_info <- distinct(Curated_Genbank_info, .keep_all = FALSE)

Curated_Genbank_info <- group_by(Curated_Genbank_info,Species)
Curated_Genbank_info_Summary <- summarize(Curated_Genbank_info,Count=n())
Curated_Genbank_info_Summary <- filter(Curated_Genbank_info_Summary, Count < 1000)


gg_id_names_list <- c(as.vector(Curated_Genbank_info$gg_id))

#writeXStringSet(Filtered_Aseembly_Microbiome, "Microbiome04_curated.fas")

Green_Genes_Bacterial_list <- Curated_Genbank_info_Summary$Species

system("touch BACTERIAL_GROUPS/GreenGenesConsensus.fas")
###Group by Family to understand variation
counter=0;
for (Bacterial_group_lookup in Green_Genes_Bacterial_list)
{
counter = counter +1;
Consensus_list <- list()
Combined_Taxa_msa_filter <- filter(Curated_Genbank_info,grepl(Bacterial_group_lookup, Species, fixed = TRUE))

msa_analysis <- c(as.character(as.vector(Combined_Taxa_msa_filter$gg_id)))

Filtered_Assembly_msa <- Assembly_Microbiome[msa_analysis]
if (length(Filtered_Assembly_msa)>1)
{
Consensus_msa <- msaConsensusSequence(msa(Filtered_Assembly_msa))
Consensus_msa <- as.character(Consensus_msa)
}else {Consensus_msa <- as.character(Filtered_Assembly_msa) }


Bacterial_group_lookup <- gsub("[[:punct:]]", "_", Bacterial_group_lookup)

Consensus_msa <- gsub("[[:punct:]]", "N", Consensus_msa)

Consensus_list <- append(Consensus_list,c(paste0(">",Bacterial_group_lookup),Consensus_msa))

#writeXStringSet(Filtered_Assembly_msa,paste0("BACTERIAL_GROUPS/",Bacterial_group_lookup,".fas"))
sink("BACTERIAL_GROUPS/GreenGenesConsensus.fas", append = TRUE)
writeLines(unlist(lapply(Consensus_list, paste, collapse=" ")))
sink()

print(paste0(counter,":",Bacterial_group_lookup)) 

}

