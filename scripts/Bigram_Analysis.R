library(dplyr)
library(stringr)
library(tidytext)
library(tidyr)
library(ggplot2)
library(stargazer)

#install.packages("tidytext")
args = commandArgs(trailingOnly=TRUE)

input_file=args[1]
Percentile_filter_arg=as.numeric(args[2])
output_folder=args[3]
project=args[4]

print(output_folder)

Species_Count <- read.csv(input_file, header = FALSE)
colnames(Species_Count) <- "Species"
Species_tibble <- tibble(txt = Species_Count$Species)

Species_bigrams <- Species_tibble %>%  unnest_tokens(bigram, txt, token = "ngrams", n = 2)

Summary_1 <- Species_bigrams %>%
  count(bigram, sort = TRUE)

bigrams_separated <- Species_bigrams %>%
  separate(bigram, c("word1", "word2"), sep = " ")

# new bigram counts:
bigram_counts <- bigrams_separated %>% 
  count(word1, word2, sort = TRUE)
Total_bigram_counts <- as.numeric(sum(bigram_counts$n))

bigram_counts$Relative_OTU <- round((bigram_counts$n/Total_bigram_counts)*100,digits=2)

bigram_counts_grouped <-bigram_counts %>%
                        group_by(word1) %>%
                        summarize(Total=sum(n)) %>%
                        filter(Total >= quantile(Total, probs =Percentile_filter_arg)) %>%
                        arrange(desc(Total))

Total_bigram_counts <- as.numeric(sum(bigram_counts_grouped$Total))
bigram_counts_grouped$relative_OTU_group <- round((bigram_counts_grouped$Total/Total_bigram_counts)*100,digits=2)

bigram_counts_grouped$word1 <- factor(bigram_counts_grouped$word1, levels = bigram_counts_grouped$word1[order(bigram_counts_grouped$relative_OTU_group)])


out_name=paste0(output_folder,"/OTU_groups_",project,".jpg")
jpeg(out_name, height = 1500, width = 1500 ,res=300)

ggplot(bigram_counts_grouped,aes(relative_OTU_group, factor(word1), fill = word1)) +
  geom_col(show.legend = FALSE)

dev.off()



###Percentile filter
Percentile_Filter <- quantile(bigram_counts$Relative_OTU, probs = Percentile_filter_arg)

bigram_counts_filtered <- filter(bigram_counts, Relative_OTU >= Percentile_Filter)
Total_bigram_counts <- as.numeric(sum(bigram_counts_filtered$n))
bigram_counts_filtered$Relative_OTU <- round((bigram_counts_filtered$n/Total_bigram_counts)*100,digits=2)

bigram_counts_grouped <-bigram_counts_filtered %>%
  group_by(word1) %>%
  summarize(Total=sum(n)) %>%
  filter(Total >= quantile(Total, probs =(Percentile_filter_arg-0.2))) %>%
  arrange(desc(Total))

Total_bigram_counts <- as.numeric(sum(bigram_counts_grouped$Total))
bigram_counts_grouped$relative_OTU_group <- round((bigram_counts_grouped$Total/Total_bigram_counts)*100,digits=2)

bigram_counts_grouped$word1 <- factor(bigram_counts_grouped$word1, levels = bigram_counts_grouped$word1[order(bigram_counts_grouped$relative_OTU_group)])

out_name=paste0(output_folder,"/OTU_groups_filtered_",project,"_percentile_",Percentile_filter_arg,".jpg")
jpeg(out_name, height = 1500, width = 1500 ,res=300)

ggplot(bigram_counts_grouped,aes(relative_OTU_group, factor(word1), fill = word1)) +
  geom_col(show.legend = FALSE)

dev.off()




###Graphical Network
library(ggraph)
library(igraph)
library(stats)
#devtools::install_github('thomasp85/ggraph')
#install.packages("igraph")

bigram_graph <- bigram_counts_filtered %>%
  graph_from_data_frame(directed = FALSE)

out_name=paste0(output_folder,"/Dendogram_",project,"_percentile_",Percentile_filter_arg,".jpg")
jpeg(out_name, width = 15000, height = 5000)


c1 <-  cluster_fast_greedy(bigram_graph)
hc <-  as.hclust(c1)
plot(hc, hang = -1, cex = 1.5)

dev.off()

####Write results
colnames(bigram_counts) <- c("Genus","Species","Relative_Counts","Relative_Abundance%")  
#Table to text file
stargazer(bigram_counts,                 # Export txt
          summary = FALSE,
          type = "text",
          out = paste0(output_folder,"/Bigram_Analysis_",project,".summary"))

colnames(bigram_counts_filtered) <- c("Genus","Species","Relative_Counts","Relative_Abundance%")  
#Table to text file
stargazer(bigram_counts_filtered,                 # Export txt
          summary = FALSE,
          type = "text",
          out = paste0(output_folder,"/Bigram_Analysis_",project,"_percentile_",Percentile_filter_arg,".summary"))
