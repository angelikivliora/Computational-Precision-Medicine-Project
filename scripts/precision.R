library(tidyverse)
library(readr)
library(dplyr)

count_matrix <- readr::read_delim("LIHC-counts.tsv", delim = '\t')
metadata <- readr::read_delim("TCGA-LIHC.clinical.tsv", delim = '\t')
gene_matrix <- readr::read_delim("gencode.v36.annotation.gtf.gene.probemap", delim = '\t')

gene_matrix <-  gene_matrix |>
  dplyr::select(c(id, gene)) |>
  rename(Ensembl_ID=id)

count_matrix <- count_matrix |>
  inner_join(gene_matrix, join_by(Ensembl_ID)) |>
  dplyr::select(-c(Ensembl_ID)) |>
  relocate(gene)
  
count_matrix <-  count_matrix |> 
  pivot_longer(cols = -1, names_to = "sample", values_to = "log2_value") |>
  mutate(count = 2^log2_value - 1) |>
  group_by(gene, sample) |>
  summarise(count = sum(count), .groups = "drop") |>
  pivot_wider(names_from = sample, values_from = count)

# Convert Ensembl IDs to gene symbols
rownames(count_matrix) <- count_matrix$gene

#write.table(count_matrix, file='count_TCGA.tsv', quote=FALSE, sep='\t')

# Step 1: Filter out normal samples
tumor_samples <- metadata |>
  subset(tissue_type.samples == "Tumor" 
         & primary_diagnosis.diagnoses == "Hepatocellular carcinoma, NOS" 
         & sample_type.samples == "Primary Tumor" 
         & name.tissue_source_site == "Mayo Clinic - Rochester"
         & race.demographic == "white") |>
  (\(df) df$sample)() |>            # Extract sample IDs
  intersect(colnames(count_matrix)) 

# Step 2: Subset Count Matrix to Retain Only Tumor Samples
filtered_count_matrix <- count_matrix |>
  subset(select = c("gene", tumor_samples))

# Step 3: Filter Metadata to Match Filtered Sample IDs
filtered_metadata <- metadata |>
  subset(sample %in% colnames(filtered_count_matrix)[-1])

