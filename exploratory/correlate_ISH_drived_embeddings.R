#this code compares pairwise distances in the ProtT5 embedding to embeddings derived from gene expression images of the human brain (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0262717)
#uses the Mantel test - correlation of pairwise distances
#a follow-up test could test if removing the non-expressing genes help (as annotated in the brain data)
library(vegan)
library(magrittr)
library(here)
library(tidyverse)
triplet <- read_csv(here("data", "Abed et al brain image embeddings", "crtx_triplet_all_training_embeddings_gene_level_with_info.csv"))

triplet %<>% select(-entrez_id)

#this is split into two to handle github file size limits
protT5_1 <- read_csv(here("data", "processed", "gene_symbol_summarized_prottrans_t5_xl_u50.1.csv.zip"))
protT5_2 <- read_csv(here("data", "processed", "gene_symbol_summarized_prottrans_t5_xl_u50.2.csv.zip"))

protT5 <- bind_rows(protT5_1,protT5_2)

dim(protT5)

common_genes <- intersect(protT5 %>% pull(gene_symbol), triplet %>% pull(gene_symbol))
print(paste("Genes in both spaces:", length(common_genes)))

protT5 %<>% filter(gene_symbol %in% common_genes)
triplet %<>% filter(gene_symbol %in% common_genes)

protT5 %<>% arrange(gene_symbol)
triplet %<>% arrange(gene_symbol)

#permutations can be set to 1000 and it still won't get a better correlation
mantel(triplet %>% select(-gene_symbol) %>% dist(), protT5 %>% select(-gene_symbol) %>% dist(), permutations = 99)

mantel(triplet %>% select(-gene_symbol) %>% dist(), protT5 %>% select(-gene_symbol) %>% dist(),
       method = "spearman", permutations = 99)
