library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

input_file=args[1]
output_file=args[2]

Species_Count <- read.csv(input_file, header = FALSE)
colnames(Species_Count) <- c("Genus","species","counts","relative_abundance")
Species_Count$relative_abundance <- round(Species_Count$relative_abundance, digits = 0) +1
Species_Count$Genus <- toupper(gsub("[[:punct:]]", "", Species_Count$Genus))
Species_Count$species <-toupper(gsub("[[:punct:]]", "", Species_Count$species))

Species_Count$Bacteria <- paste0(Species_Count$Genus," ",Species_Count$species)

Expanded_List <- list()
for (i in 1:nrow(Species_Count))
  {
  Bacteria_List<- list(replicate(Species_Count$relative_abundance[i], Species_Count$Bacteria[i]))
  Expanded_List <- append(Expanded_List,Bacteria_List)
}

sink(output_file, append = TRUE)
writeLines(unlist(Expanded_List))
sink()

