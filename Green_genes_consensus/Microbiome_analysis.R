#https://joey711.github.io/phyloseq/import-data.html
#https://www.yanh.org/2021/01/01/microbiome-r/
#zcat *.fastq.gz > Microbiome_barcode01.fastq
#BLAST comparison in Terminal v 2.11

#../ncbi-blast-2.11.0+/bin/blastn -query test1.test -subject PSC8_GFF_DNA.genes.fasta  -evalue 1e-50  -outfmt 6 > HAN412_vs_PSC8.blast.table


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

#BiocManager::install("phyloseq")
#BiocManager::install("bamsignals")

setwd("/DATA/home/jmlazaro/Projects/Microbiome/Ref")

Assembly_Microbiome <- readDNAStringSet("gg_13_5.fasta.gz")
Assembly_Microbiome_info <- as.data.frame(names(Assembly_Microbiome))
colnames(Assembly_Microbiome_info) <- "gg_id"
Assembly_Microbiome_info$width <- Assembly_Microbiome@ranges@width

Microbiome_Assembly_Genbank <- read.table("gg_13_5_accessions.txt", header = TRUE)
best_single_hits.blastn <- read.table("best_single_hits.blastn")
colnames(best_single_hits.blastn) <- c("qseqid","accession","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")


HMQCP_gg_13_merge <- merge(best_single_hits.blastn,Microbiome_Assembly_Genbank, all.x = TRUE, by = "accession")

Assembly_Microbiome_info_with_Genbank <-merge(Assembly_Microbiome_info,Microbiome_Assembly_Genbank, all.x = TRUE, by = "gg_id")

#taxonomy <- read.table(file = "gg_13_5_taxonomy.txt", sep = "\t", header = F)
#colnames(taxonomy) <- c("gg_id","Taxon")
#rownames(taxonomy) <- taxonomy$gg_id
#tax <- taxonomy %>%
#  select(Taxon) %>% 
#  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "; ")

# clean the taxonomy, Greengenes format
#tax.clean <- data.frame(row.names = row.names(tax),
#                        Kingdom = str_replace(tax[,1], "k__",""),
#                        Phylum = str_replace(tax[,2], "p__",""),
#                        Class = str_replace(tax[,3], "c__",""),
#                        Order = str_replace(tax[,4], "o__",""),
#                        Family = str_replace(tax[,5], "f__",""),
#                        Genus = str_replace(tax[,6], "g__",""),
#                        Species = str_replace(tax[,7], "s__",""),
#                        stringsAsFactors = FALSE)

#tax.clean[is.na(tax.clean)] <- ""
#tax.clean[tax.clean=="__"] <- ""

#saveRDS(tax.clean, file = "tax.clean.rds")

tax.clean <- readRDS(file = "tax.clean.rds")

tax.clean$gg_id <- rownames(tax.clean)


#Reading from Genbank
#seq_1_DNAbin <- read.GenBank("NR_025230.1")
#attr(seq_1_DNAbin, "species")
#attr(seq_1_DNAbin, "description")

Microbiome_Assembly_gr <- with(Assembly_Microbiome_info_with_Genbank, GRanges(gg_id, IRanges(1, width),"+", accession = accession))
Microbiome_Assembly_gr <- sortSeqlevels(Microbiome_Assembly_gr)
Microbiome_Assembly_gr <- sort(Microbiome_Assembly_gr, ignore.strand = TRUE)


setwd("/DATA/home/jmlazaro/Projects/Microbiome/Minimap_Aligned_gg_13_5")

Barcode_number="04"
#Reading Bam file
Bam_File_Name <- paste0("Microbiome_minimap2_barcode",Barcode_number,".sorted.bam")
#Bam_File <- Rsamtools::BamFile(Bam_File_Name)
#seqinfo(Bam_File)

#############################
#Filtering BAM
bam_coverage <- scanBam(Bam_File_Name)


#function for collapsing the list of lists into a single list
#as per the Rsamtools vignette
.unlist <- function (x){
  ## do.call(c, …) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#store names of BAM fields
bam_field <- names(bam_coverage[[1]])

#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam_coverage, "[[", y)))

#store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

dim(bam_df)

Microbiome_counts <- as.data.frame(bam_df[,c("rname","flag","strand","mapq","pos","qwidth","cigar","seq")])
Microbiome_counts$MAPQ_filter <- ifelse(Microbiome_counts$mapq>=20, "MAPQ20",
                                        ifelse(Microbiome_counts$mapq>=10,"MAPQ10","MAPQ_BELOW10"))


Microbiome_counts_filtered <-
  Microbiome_counts %>%
  group_by(rname) %>%
  summarise(
    counts = dplyr::n(),
    counts_mq10 = sum(MAPQ_filter == "MAPQ10"),
    counts_mq20 = sum(MAPQ_filter == "MAPQ20"),
    Average_mq =  round(mean(mapq), digits = 4)
  )

colnames(Microbiome_counts_filtered)[1] <- "gg_id"

Assembly_Microbiome_info_with_Genbank_counts <- merge(Microbiome_counts_filtered,Assembly_Microbiome_info_with_Genbank, all.x=TRUE, by = "gg_id")

#Assembly_Microbiome_info_with_Genbank_counts_filtered <- filter(Assembly_Microbiome_info_with_Genbank_counts,counts_mq10 >0)
#Assembly_Microbiome_info_with_Genbank_counts_filtered <- filter(Assembly_Microbiome_info_with_Genbank_counts_filtered,
#                                                                counts_mq10 >= 2 | counts_mq20 > 0)

Assembly_Microbiome_info_with_Genbank_counts_filtered <- filter(Assembly_Microbiome_info_with_Genbank_counts,counts >1)
#End Filtering BAM
##############################


#List_of_Accesions <- as.vector(Assembly_Microbiome_info_with_Genbank_counts_filtered$accession[1:100])
#Organisms_info <- read.GenBank(List_of_Accesions) 
Assembly_Microbiome_info_with_Genbank_counts_filtered <- Assembly_Microbiome_info_with_Genbank_counts_filtered[order(Assembly_Microbiome_info_with_Genbank_counts_filtered$counts_mq10, decreasing = TRUE), ]

Assembly_Microbiome_info_with_Genbank_counts_filtered$Species <- NA
Assembly_Microbiome_info_with_Genbank_counts_filtered$Description <- NA
#for (i in 1:nrow(Assembly_Microbiome_info_with_Genbank_counts_filtered))
#  {
#    seq_1_DNAbin <- read.GenBank(Assembly_Microbiome_info_with_Genbank_counts_filtered$accession[i])
#    Assembly_Microbiome_info_with_Genbank_counts_filtered$Species[i] <-  attr(seq_1_DNAbin, "species")
#    Assembly_Microbiome_info_with_Genbank_counts_filtered$Description[i] <- attr(seq_1_DNAbin, "description")
#    print(i)
#  }

No_accesion = nrow(Assembly_Microbiome_info_with_Genbank_counts_filtered)
count = 1
Multiplier = 100

while (No_accesion >= 0)
{
  start = (count *Multiplier)-Multiplier +1 
   end = start + (Multiplier-1)
 
   No_accesion = No_accesion - Multiplier
   end <-ifelse (No_accesion >= 0, end, (start + No_accesion + Multiplier - 1))
   
   print(paste0(start,":",end))
   print(paste0("Left Bacteria to get info: ",No_accesion))
   List_of_Accesions <- as.vector(Assembly_Microbiome_info_with_Genbank_counts_filtered$accession[start:end])
   Organisms_info <- read.GenBank(List_of_Accesions) 
   Assembly_Microbiome_info_with_Genbank_counts_filtered$Species[start:end] <-  attr(Organisms_info, "species")
   Assembly_Microbiome_info_with_Genbank_counts_filtered$Description[start:end] <- attr(Organisms_info, "description")
   
  count = count+1
}


Combined_Taxa_count <- merge(Assembly_Microbiome_info_with_Genbank_counts_filtered, tax.clean, all.x = TRUE, by = "gg_id")
Combined_Taxa_count$Species <- ifelse(Combined_Taxa_count$Species.x=="uncultured_bacterium",
                                      paste0("Unclassified ",Combined_Taxa_count$Genus," ",Combined_Taxa_count$Species.y),
                                      Combined_Taxa_count$Species.x)
Combined_Taxa_count$Species <- ifelse(is.na(Combined_Taxa_count$Species)==TRUE,
                                      paste0("Unclassified ",Combined_Taxa_count$Genus," ",Combined_Taxa_count$Species.y),
                                      Combined_Taxa_count$Species)

Combined_Taxa_count <- Combined_Taxa_count[order(Combined_Taxa_count$counts_mq10, decreasing = TRUE), ]

saveRDS(Combined_Taxa_count,"Microbiome04_Combined_Taxa_count.rds")

gg_id_names_list <- c(as.vector(Combined_Taxa_count$gg_id))

Filtered_Aseembly_Microbiome <- Assembly_Microbiome[gg_id_names_list,]
writeXStringSet(Filtered_Aseembly_Microbiome, "Microbiome04_curated.fas")

###Group by Family to understand variation
Bacterial_group_lookup="prausnitzii"
Combined_Taxa_msa_filter <- filter(Combined_Taxa_count,grepl(Bacterial_group_lookup, Species.y, fixed = TRUE))

msa_analysis <- c(as.vector(Combined_Taxa_msa_filter$gg_id))

Filtered_Aseembly_msa <- Filtered_Aseembly_Microbiome[msa_analysis,]
writeXStringSet(Filtered_Aseembly_msa,paste0(Bacterial_group_lookup,".fas"))




#Creating Otu matrix
##Plot relative abundace all read in bam
otumat <- as.matrix(Combined_Taxa_count[,c("counts","counts_mq10","counts_mq20")])
rownames(otumat) <- Combined_Taxa_count$accession
colnames(otumat) <- c("MAPQ_0","MAPQ10","MAPQ20")

taxmat = as.matrix(Combined_Taxa_count[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

physeq = phyloseq(OTU, TAX)

Plot_OTU_All <- plot_bar(physeq, fill = "Class")



###Plot relative abundace only above mapq10 and pass filter (at least 2 Mapq10 or 1 Mapq20)
otumat <- as.matrix(Combined_Taxa_count[,c("counts_mq10","counts_mq20")])
rownames(otumat) <- Combined_Taxa_count$accession
colnames(otumat) <- c("MAPQ10","MAPQ20")


taxmat = as.matrix(Combined_Taxa_count[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

physeq = phyloseq(OTU, TAX)

Plot_OTU_above_MapQ10 <- plot_bar(physeq, fill = "Class")

p1
p2


###Circular Layout tree with relative average quality
physeq_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
#Tables
tree_tibble <-as_tibble(physeq_tree)
AA_tibble <- tibble(label = Combined_Taxa_count$accession,
                         Bacterial_Class = Combined_Taxa_count$Class,
                          Average_MAPQ = Combined_Taxa_count$Average_mq)

tree_with_info <- as.treedata(full_join(tree_tibble, AA_tibble, by = 'label'))

tr <- ape::as.phylo(full_join(tree_tibble, AA_tibble, by = 'label'))


dt = data.frame(id=tr$tip.label, group=Combined_Taxa_count$Class)

picture_Bacterial_Class <- ggtree(tree_with_info, branch.length='none', layout='circular', ladderize = FALSE, size=2) +  aes(color=Average_MAPQ)+
  theme_tree() + geom_tiplab(offset=0.5, align=TRUE,size =2.5,aes(angle=angle, color=Average_MAPQ)) 


picture_Bacterial_Class_layer2 <- picture_Bacterial_Class +
  geom_fruit(
    data=dt,
    geom=geom_star,
    mapping=aes(y=id, fill=group),
    size=2.5,
    starstroke=0
  )


##Graph2
dd = data.frame(id=tr$tip.label, value=abs(Combined_Taxa_count$counts_mq10),Class =Combined_Taxa_count$Class )
dt = data.frame(id=tr$tip.label, group=Combined_Taxa_count$Species)
p <- ggtree(tr, branch.length='none', layout='circular', ladderize = FALSE, size=2.5) + geom_tiplab(offset=1, align=TRUE,size =2.5,aes(angle=angle))

p1 <- p + 
  geom_fruit(
    data=dt,
    geom=geom_star,
    mapping=aes(y=id, fill=group),
    size=2.5,
    starstroke=0
  )
Relative_Quantitites_layer2 <- p1 + 
  geom_fruit(
    data=dd, 
    geom=geom_bar, 
    mapping=aes(x=value, y=id,fill=Class),
    orientation="y",
    stat="identity",
    offset = 0.37
  )  

#Images to PDF
pdf(paste0("Functional_Microbiome_tree_barcode",Barcode_number,".pdf"),  width=25, height=15)
Plot_OTU_All
Plot_OTU_above_MapQ10
picture_Bacterial_Class_layer2
Relative_Quantitites_layer2
dev.off()

row.names(Combined_Taxa_count) <- 1:nrow(Combined_Taxa_count)

#Table to text file
stargazer(Combined_Taxa_count,                 # Export txt
          summary = FALSE,
          type = "text",
          out = paste0("Functional_Info_barcode",Barcode_number,".txt"))

#Table to text file
stargazer(Combined_Taxa_count,                 # Export txt
          summary = TRUE,
          type = "text",
          out = paste0("Sumary_Info_barcode",Barcode_number,".txt"))

#Table to text file
stargazer(Microbiome_counts,                 # Export txt
          summary = TRUE,
          type = "text",
          out = paste0("Sumary_Raw_bam_barcode",Barcode_number,".txt"))



###Extracting new FASTA

gg_new_fasta_id_list = as.character(Microbiome_counts_filtered$gg_id)

New_microbiome_match_list = as.data.frame(match(gg_new_fasta_id_list, Assembly_Microbiome@ranges@NAMES))
colnames(New_microbiome_match_list) <- "Position"
New_microbiome_match_list = filter(New_microbiome_match_list, is.na(New_microbiome_match_list$Position)==FALSE)
New_Microbiome_Fasta <- Assembly_Microbiome[New_microbiome_match_list$Position]
#New_Microbiome_Fasta[grepl("593129", New_Microbiome_Fasta@ranges@NAMES)]
writeXStringSet(New_Microbiome_Fasta, "New_Microbiome_Fasta_iteration1.fas")

