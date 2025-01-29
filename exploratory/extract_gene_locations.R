#TODO - add strand
library(tidyr)
library(readr)
library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(caret)

gene_locations <- read_csv("/Users/leon/projects/polyprotein_stats/data/processed/gene_symbol_summarized_proportions.csv") %>% select(gene_symbol)


# Step 1: Map gene symbols to Entrez IDs
gene_map <- AnnotationDbi::select(org.Hs.eg.db, 
                   keys = gene_locations$gene_symbol, 
                   keytype = "SYMBOL", 
                   columns = c("ENTREZID", "CHR", "CHRLOC"))

# Step 2: Process gene locations
gene_map <- gene_map %>%
  mutate(CHRLOC = abs(CHRLOC)) %>%  # Convert negative values to absolute positions
  distinct(ENTREZID, .keep_all = TRUE)  # Keep unique mappings

# Step 3: Merge with original tibble
gene_locations <- gene_locations %>%
  left_join(gene_map, by = c("gene_symbol" = "SYMBOL")) %>%
  select(gene_symbol, CHR, CHRLOC)  # Keep relevant columns

gene_locations %<>% filter(!is.na(CHR))
gene_locations %<>% filter(!is.na(CHRLOC))
gene_locations %>% group_by(gene_symbol) %>% count() %>% filter(n!=1)
# View output
print(gene_locations)

gene_locations <- gene_locations %>%
  group_by(CHR) %>%
  mutate(CHRLOC_normalized = (CHRLOC - min(CHRLOC)) / (max(CHRLOC) - min(CHRLOC))) %>%
  ungroup()
hist(gene_locations$CHRLOC_normalized)

gene_locations <- gene_locations %>%
  mutate(CHR = as.character(CHR)) %>% 
  pivot_wider(names_from = CHR, values_from = CHR, names_prefix = "CHR_") %>%
  mutate(across(starts_with("CHR_"), ~ ifelse(is.na(.), 0, 1)))  # Convert to binary

gene_locations %>% write_csv("/Users/leon/projects/polyprotein_stats/data/processed/gene_symbol_locations.csv")

