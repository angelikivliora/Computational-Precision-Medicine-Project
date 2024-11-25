# Load required libraries
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(gtable)
library(readr)


# ------------------------------------------
# Step 1: Load and Prepare the Data
# ------------------------------------------

# Read in the data. tsv file with raw counts. read.delim is used to read in tsv files
cancer <- read.delim("data/filtered_count_TCGA.tsv", row.names = 1) # cancer samples

# From no_cancer we can drop id column and Name column. We also rename column Description to gene and row names to gene names
no_cancer <- read.delim("data/healthy_controls.tsv", row.names = 1) # non-cancer samples
no_cancer <- no_cancer[, -c(1, 2)] # Drop id and Name columns
colnames(no_cancer)[1] <- "gene" # Rename Description to gene

# Check for duplicate gene names in cancer
sum(duplicated(rownames(cancer)))

# Check for duplicate gene names in no_cancer
sum(duplicated(no_cancer$gene))

# Inspect duplicated genes
duplicate_genes <- no_cancer[duplicated(no_cancer$gene) | duplicated(no_cancer$gene, fromLast = TRUE), ] # Get duplicate genes from no_cancer
head(duplicate_genes)

# Mean the duplicate genes
no_cancer <- aggregate(. ~ gene, data = no_cancer, mean) # Mean the duplicate genes. Unique gene names now
sum(duplicated(no_cancer$gene)) # Check for duplicates again. Should be 0
# Will introduce non-integers and we will take care of that later


# Count the amount of samples (amount of columns)
ncol(cancer)  # Should be 47
ncol(no_cancer) # Should be 47

# Merge by the 'gene' column
combined_data <- merge(cancer, no_cancer, by = "gene", all = TRUE)
# Set gene column as row names
rownames(combined_data) <- combined_data$gene
combined_data <- combined_data[, -1]  # Remove gene column
combined_data <- round(combined_data) # Round to integers

# Make sure all values are integers
all(sapply(combined_data, function(x) all.equal(x, as.integer(x))))


# Dynamically create metadata
metadata <- data.frame(
  Sample = colnames(combined_data),
  Group = rep(c("Cancer", "No_Cancer"), each = 47)
)

# Ensure Group is a factor with the correct levels
metadata$Group <- factor(metadata$Group, levels = c("No_Cancer", "Cancer"))
# Example: Adding batch information to metadata
metadata$Batch <- rep(c("Lab1", "Lab2"), each = 47) 
# Check for confounding between Batch and Group
table(metadata$Batch, metadata$Group)


# Verify alignment
nrow(metadata) == ncol(combined_data)  # Should return TRUE
# Drop rows with any NA values
combined_data <- na.omit(combined_data)


# ------------------------------------------
# Step 2: PCA + PCA Plot. Normalizing Raw Counts with DESeq2
# ------------------------------------------

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = combined_data, colData = metadata, design = ~ Group)
dds <- dds[rowSums(counts(dds)) > 10, ] # Remove noise. Extremely low counts

# Perform normalization
dds <- estimateSizeFactors(dds) # Estimate size factors
normalized_counts <- counts(dds, normalized = TRUE) # Get normalized counts

# Perform PCA on normalized counts
pca <- prcomp(t(log2(normalized_counts + 1)))  # log2 transform for PCA stability
pca_data <- as.data.frame(pca$x)
pca_data$Group <- metadata$Group

# Plot the PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  ggtitle("PCA Plot (Normalized Counts)") +
  xlab("PC1") +
  ylab("PC2") +
  theme_minimal()

# ------------------------------------------
# Step 3: Hierarchical Clustering + Heatmap
# ------------------------------------------

# Compute hierarchical clustering
dist_matrix <- dist(t(normalized_counts)) # Compute distance matrix
hc <- hclust(dist_matrix, method = "average")

# Generate heatmap
pheatmap(
  normalized_counts,
  scale = "row",  # Scale rows for better visualization
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  annotation_col = metadata
)

# ------------------------------------------
# Step 4: Differential Gene Expression Analysis
# ------------------------------------------

# Perform differential expression analysis
dds <- DESeq(dds) # Estimate dispersion and fit model
res <- results(dds) # Get results

# Order results by adjusted p-value
res <- res[order(res$padj), ] # Order by adjusted p-value

# Res into a data frame
res_df <- as.data.frame(res) # Convert to data frame

# Remove na values from res_df since they are no good
res_df <- na.omit(res_df)

# Find the smallest non-zero padj value
min_nonzero_padj <- min(res_df$padj[res_df$padj > 0], na.rm = TRUE)

# Replace padj values that are 0 with the smallest non-zero padj
res_df$padj[res_df$padj == 0] <- min_nonzero_padj

# Define thresholds for significance
logfc_threshold <- 3  # Fold-change threshold
padj_threshold <- 1e-6  # Very strict adjusted p-value

# Add a column for significant gene labels
res_df <- res_df %>%
  filter(!is.na(log2FoldChange), !is.na(padj)) %>%  # Remove invalid rows
  mutate(
    label_adj = case_when(
      padj < padj_threshold & log2FoldChange > logfc_threshold ~ "up",
      padj < padj_threshold & log2FoldChange < -logfc_threshold ~ "down",
      TRUE ~ "not_sig"
    )
  )

# Dynamically set y-axis limit
y_limit <- ceiling(-log10(min(res_df$padj, na.rm = TRUE)))

# Create the volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = label_adj), size = 1.2, alpha = 0.8) +
  scale_color_manual(
    values = c('up' = "#008c41", 'down' = "#EC2326", 'not_sig' = "#BFBFBF")
  ) +
  ylim(c(0, 310)) +  # Slightly above the max of -log10(padj)
  xlab(bquote(Log["2"]("Fold Change"))) +
  ylab(bquote(-Log["10"]("Adjusted P-value"))) +
  geom_hline(
    yintercept = -log10(padj_threshold), 
    linetype = "dashed", 
    color = "gray40", 
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = c(-logfc_threshold, logfc_threshold), 
    linetype = "dashed", 
    color = c("gray40", "gray40"), 
    linewidth = 1
  ) +
  theme_minimal() +
  ggtitle("Volcano Plot (Differential Gene Expression)")


# Display the plot
print(volcano_plot)

# Based on the results lets filter out the significant genes for fgsea analysis
ranked_genes <- res_df$log2FoldChange
names(ranked_genes) <- rownames(res_df)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Write the ranked genes to a file
write.table(ranked_genes, file = "data/ranked_genes.txt", sep = "\t", col.names = FALSE, quote = FALSE)
