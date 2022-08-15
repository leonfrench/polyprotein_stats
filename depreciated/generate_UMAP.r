library(magrittr)
library(dplyr)
library(readr)
library(umap)
feats <- read_csv("/Users/leon/projects/polyprotein_stats/data/CZI_image_spectra/feature_spectra.csv")
z <- umap(feats %>% select(-gene_symbol))
z <- z$layout %>% as_tibble()
colnames(z) <- c("UMAP 0","UMAP 1")
z %<>% mutate(gene_symbol = feats$gene_symbol)
z %<>% select(gene_symbol, everything())

plot(z$`UMAP 0`, z$`UMAP 1`)

write_csv(z, "/Users/leon/projects/polyprotein_stats/data/CZI_image_spectra/gene_symbol_summarized_UMAP.csv")
