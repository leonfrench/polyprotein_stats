#this code compares pairwise distances in the ProtT5 embedding to embeddings derived from gene expression images of the human brain (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0262717)
#uses the Mantel test - correlation of pairwise distances
#this is just a rough exploratory notebook at this point
library(vegan)
library(magrittr)
library(here)
library(tidyverse)

#given two matrices with gene_symbol column, correlate gene-gene distances in the two spaces using the Mantel test
run_mantel_on_data_frames <- function(df_a, df_b, permutations) {
  common_genes <- intersect(df_a %>% pull(gene_symbol), df_b %>% pull(gene_symbol))
  print(paste("Genes in both spaces:", length(common_genes)))
  
  df_a %<>% filter(gene_symbol %in% common_genes)
  df_b %<>% filter(gene_symbol %in% common_genes)
  
  df_a %<>% arrange(gene_symbol)
  df_b %<>% arrange(gene_symbol)
  
  #permutations can be set to 1000 and it still won't get a better correlation
  result <- mantel(df_b %>% select(-gene_symbol) %>% dist(), df_a %>% select(-gene_symbol) %>% dist(), permutations = permutations)
  result
}

triplet <- read_csv(here("data", "Abed et al brain image embeddings", "crtx_triplet_all_training_embeddings_gene_level_with_info.csv"))

triplet %<>% select(-entrez_id)

#this is split into two to handle github file size limits
protT5_1 <- read_csv(here("data", "processed", "gene_symbol_summarized_prottrans_t5_xl_u50.1.csv.zip"))
protT5_2 <- read_csv(here("data", "processed", "gene_symbol_summarized_prottrans_t5_xl_u50.2.csv.zip"))

protT5 <- bind_rows(protT5_1,protT5_2)

dim(protT5)

CZI_image_spectra <- read_csv(here("data", "CZI_image_spectra","feature_spectra.csv"))
dim(CZI_image_spectra)
run_mantel_on_data_frames(protT5, CZI_image_spectra, 99)
run_mantel_on_data_frames(CZI_image_spectra, triplet, 999)
setdiff(CZI_image_spectra$gene_symbol, protT5$gene_symbol)

#main result - comparing triplet and protT5
run_mantel_on_data_frames(protT5, triplet, 99)

#how does this compare to using proportions data?
proportions <- read_csv(here("data", "processed", "gene_symbol_summarized_proportions.csv"))

run_mantel_on_data_frames(triplet, proportions, permutations = 99)

#load in the annotations from Zeng
zeng_annotations <- read_csv(here("data","Zeng et al.","Cleaned_Zeng_dataset.csv"))

zeng_annotations %>% group_by(Expression.level) %>% count()

subset <- triplet %>% filter(gene_symbol %in% (zeng_annotations %>% filter(Expression.level == "+++++") %>% pull(gene_symbol)))
run_mantel_on_data_frames(subset, protT5, permutations = 99)

subset <- triplet %>% filter(gene_symbol %in% (zeng_annotations %>% filter(Expression.level == "++++") %>% pull(gene_symbol)))
run_mantel_on_data_frames(subset, protT5, permutations = 99)

subset <- triplet %>% filter(gene_symbol %in% (zeng_annotations %>% filter(Expression.level == "+++") %>% pull(gene_symbol)))
run_mantel_on_data_frames(subset, protT5, permutations = 99)

subset <- triplet %>% filter(gene_symbol %in% (zeng_annotations %>% filter(Expression.level == "++") %>% pull(gene_symbol)))
run_mantel_on_data_frames(subset, protT5, permutations = 99)

subset <- triplet %>% filter(gene_symbol %in% (zeng_annotations %>% filter(Expression.level == "+") %>% pull(gene_symbol)))
run_mantel_on_data_frames(subset, protT5, permutations = 99)

subset <- triplet %>% filter(gene_symbol %in% (zeng_annotations %>% filter(Expression.level == "-") %>% pull(gene_symbol)))
run_mantel_on_data_frames(subset, protT5, permutations = 99)

#marker information
zeng_annotations %>% group_by(Cortical.marker..human.) %>% count() %>% arrange(n)
subset <- triplet %>% filter(gene_symbol %in% (zeng_annotations %>% filter(Cortical.marker..human. == "interneuron") %>% pull(gene_symbol)))
run_mantel_on_data_frames(subset, protT5, permutations = 99)

subset <- triplet %>% filter(gene_symbol %in% (zeng_annotations %>% filter(Cortical.marker..human. == "layer 6") %>% pull(gene_symbol)))
run_mantel_on_data_frames(subset, protT5, permutations = 999)

#not much difference between annotated and not
subset <- triplet %>% filter(gene_symbol %in% (zeng_annotations %>% filter(is.na(Cortical.marker..human.)) %>% pull(gene_symbol)))
run_mantel_on_data_frames(subset, protT5, permutations = 99)

subset <- triplet %>% filter(gene_symbol %in% (zeng_annotations %>% filter(!is.na(Cortical.marker..human.)) %>% pull(gene_symbol)))
run_mantel_on_data_frames(subset, protT5, permutations = 99)


#embedding matrix based on expression level as a single dimension
zeng_annotations 
zeng_annotations %<>% filter(!is.na(Expression.level))
zeng_annotations %>% group_by(Expression.level) %>% count()
zeng_annotations %<>% mutate(Expression_level_numeric=recode(Expression.level,
                         `-`=0,
                         `+`=1,
                         `++`=2,
                         `+++`=3,
                         `++++`=4,
                         `+++++`=5
                         ))
run_mantel_on_data_frames(zeng_annotations %>% select(gene_symbol, Expression_level_numeric) %>%group_by(gene_symbol) %>% summarize_all(mean), 
                          protT5, permutations = 999)

#compared to triplet - Zeng
run_mantel_on_data_frames(zeng_annotations %>% select(gene_symbol, Expression_level_numeric) %>%group_by(gene_symbol) %>% summarize_all(mean), 
                          triplet, permutations = 9)
run_mantel_on_data_frames(protT5, triplet, permutations = 9)

#test donor embedding matrix
donor_matrix <- read_csv(here("data","Zeng et al.","cortex_gene_donor_one_hot_coded.csv"))
donor_matrix %<>% select(-entrez_id, '...1')


run_mantel_on_data_frames(donor_matrix, triplet, permutations = 99)
run_mantel_on_data_frames(donor_matrix, protT5, permutations = 99)

#remove genes assayed in H08-0002 which had all the g-protein receptors
run_mantel_on_data_frames(triplet %>% filter(gene_symbol %in% (donor_matrix %>% filter(`H08-0002` == 0) %>% pull(gene_symbol))), protT5, permutations = 99)

