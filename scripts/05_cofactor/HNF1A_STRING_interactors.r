# examining HNF1A interacting partners from the STRING database (v12.0)
# produces the DNA-binding TF candidate list feeding the cofactor enrichment

# load libraries
library(readr)
library(dplyr)

# read the STRING interaction export
HNF1A_Interactions <- read_tsv("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/HNF1A_STRING_interactors.tsv", show_col_types = FALSE)

HNF1A_Partners <- HNF1A_Interactions %>%
  filter(`#node1` == "HNF1A" | node2 == "HNF1A") %>%
  mutate(
    partner = ifelse(`#node1` == "HNF1A", node2, `#node1`)) %>%
  arrange(desc(combined_score)) %>%
  select(partner, combined_score,
         coexpression,
         experimentally_determined_interaction,
         database_annotated,
         automated_textmining) %>%
  distinct()

# filter to DNA-binding TF families of interest
DNABINDING_TFs <- HNF1A_Partners %>%
  filter(grepl("HNF|FOXA|ONECUT|CEBP|PPAR|NR", partner))