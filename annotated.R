# ====Load Libraries====

library(tidyverse)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(readxl)
library(msigdbr)
library(fgsea)
library(pheatmap)
library(matrixStats)
library(tibble)
library(GEOquery)

# ==== Choosing data ====

### Choosing data 
gset <- getGEO("GSE281695", GSEMatrix = TRUE, getGPL = TRUE)
#view(gset)

pheno = gset[["GSE281695_series_matrix.txt.gz"]]@phenoData@data
#view(pheno)

pheno_sel = pheno |> 
  select("sample id:ch1","characteristics_ch1.2","characteristics_ch1.7")
view(pheno_sel) # sample identification by ID, ir, and time

# summary of chosen samples
sample_vec = c("NK6809", "NK5262", "NK4844", "NK6566", "NK5342", "NK2414",	"NK3736",	"NK5813")
our_samples = pheno_sel |>
  filter(`sample id:ch1` %in% sample_vec)
view(our_samples)



# ==== Data exploration ====

### Load data

raw_counts <- read_excel("GSE281695_raw_counts_all_samples.xlsx", skip = 1)

colnames(raw_counts)[1] <- "Gene"
raw_counts <- column_to_rownames(raw_counts, "Gene")
head(raw_counts)

dim(raw_counts)
colnames(raw_counts)
head(raw_counts)
str(raw_counts)
glimpse(raw_counts) 
summary(raw_counts) 
sapply(raw_counts, class)
anyNA(raw_counts) # TRUE if any NA
#is.na(raw_counts)
#colSums(is.na(raw_counts)) # count of NA per column

# counting positive expressions in each column
colSums(raw_counts == 0) 
colSums(raw_counts != 0) 


### Histogram

p_raw_counts <- raw_counts %>% 
  pivot_longer(cols = NK6809:NK5813, 
               names_to = "sample", 
               values_to = "expression") %>%
  filter(expression > 0 & expression < 100) %>%
  ggplot(aes(expression)) +
  geom_histogram(bins = 100) +
  labs(
    title = "histogram of raw expression data < 100",
    x = "expression",
    y = "count") +
  theme_minimal()

p_raw_counts


# ==== PCA principle component analysis ====

### Remove constant rows before PCA analysis - because it needs variance
constant_cols <- apply(raw_counts, 1, var) == 0 # finds all rows with no variance
raw_counts_pca <- raw_counts[!constant_cols, ] # deletes all the rows with no variance

### PCA plot

pca_count <- prcomp(
  t(raw_counts_pca), # t - rows to columns
  scale. = TRUE) # scaling results so the large values don't overcome


df <- data.frame(PC1=pca_count$x[,1], PC2=pca_count$x[,2]) 
df # PC1 and PC2 are the sample’s coordinates on the first two principal components.
df$Sample <- colnames(raw_counts) # Adds sample names as a column.
condition = c("IR 24h", "IR 24h", "Non IR 24h", "Non IR 24h", "IR 7d", "IR 7d", "Non IR 7d", "Non IR 7d")
df$Condition = condition

# ggplot
ggplot(df, aes(x=PC1, y=PC2, color=Sample, shape = Condition)) + 
  geom_point(size = 3) + 
  labs(x = "PC1", y = "PC2", title = "PCA Plot") +
  theme_minimal() 

# ==== DEA - requires raw counts ====

### Create metadata table for DESeq2 analysis
# define experimental groups in DEA
sample_metadata <- data.frame(
  sample_id = colnames(raw_counts),
  condition = c("IR 24h", "IR 24h", "Non IR 24h", "Non IR 24h", "IR 7d", "IR 7d", "Non IR 7d", "Non IR 7d")
)
rownames(sample_metadata) <- sample_metadata$sample_id # Sets sample names as rownames
view(sample_metadata) 

### DESeq2 dataset and analysis
dds <- DESeqDataSetFromMatrix( #  creates a DESeq2 dataset object
  countData = raw_counts, 
  colData = sample_metadata, 
  design = ~ condition) # defines the model — tests for differences between conditions
# ~ condition It's the simplest model for comparing groups
dds <- DESeq(dds) 
# Runs the full DESeq2 differential expression pipeline:
# - Normalization
# - Dispersion estimation
# - Model fitting
# - Statistical testing
dds
res <- results(dds) # extracts DEA results:
# - log2 fold changes - Measures how much a gene’s expression changes between two conditions
# -  p-values
# - adjusted p- values - padj significant for <0.05

### DESeq2 comparison function
# define the function
perform_DESeq2 <- function(dds, condition1, condition2) {
  res <- results(dds, contrast=c("condition", condition1, condition2))
  res_df <- as.data.frame(res) %>% 
    rownames_to_column("Gene") %>% 
    mutate(significant = padj < 0.05, # True if yes
           condition = ifelse(log2FoldChange > 0, # gene upregulated in condition 1
                              condition1, 
                              condition2)) 
  res_df <- na.omit(res_df) 
  
  # volcano plot
  ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant, shape=condition)) +
    geom_point(alpha=0.7, size=2) + # alpha - transparency, 0 fully transparent
    scale_color_manual(values = c("turquoise", "red")) +
    scale_shape_manual(values = c(16, 17)) +
    labs(title = paste("Volcano Plot:", condition1, "vs", condition2),
         x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
    theme_minimal() +
    theme(legend.position = "top")
}

# Plot volcano plots for each condition pair
plot1 <- perform_DESeq2(dds, "IR 24h", "Non IR 24h")
plot2 <- perform_DESeq2(dds, "IR 7d", "Non IR 7d")
plot3 <- perform_DESeq2(dds, "IR 24h", "IR 7d")
plot4 <- perform_DESeq2(dds, "Non IR 24h", "Non IR 7d")

library(gridExtra)

grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)



### data visualization - heatmap
vsd <- vst(dds) # Variance Stabilizing Transformation - transforms to log2-like scale
row_variances <- rowVars(assay(vsd)) # computes variance across the samples

### Generate heatmap of top 100 most variable genes
top_genes_indices <- order(row_variances, 
                           decreasing = TRUE)[1:10]
mat <- assay(vsd)[top_genes_indices, ] # extracts top 10
annotation_df <- sample_metadata[colnames(mat),
                                 , 
                                 drop = FALSE]  # results are a df

# changing the names of genes from Ensembl to names
top_genes <- rownames(vsd)[top_genes_indices] # gets top 10 genes
top_genes
gene_symbols <- c("Alb", "Myh4", "Ahsg", "Atp2a1", "Ckm", "Trim63", "Actn3", "Apoa2", "Serpina3k", "Mylk2")
rownames(mat) = gene_symbols


# actual heatmap
pheatmap(mat, 
         scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         annotation_col = annotation_df, 
         show_rownames = T, 
         show_colnames = F, 
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         main = "Top 10 Most Variable Genes")

# single gene comparison

my_gene_rows =as.vector(top_genes_indices)
raw_counts_matrix = as.matrix(raw_counts)
colnames(raw_counts_matrix) # lists sample names in matrix

condition_colors <- c("IR 24h" = "orange", 
                      "Non IR 24h" = "skyblue", 
                      "IR 7d" = "red", 
                      "Non IR 7d" = "green")

bar_colors <- condition_colors[condition]

for (g in my_gene_rows) {
  barplot(raw_counts_matrix[g,], 
          main = paste("Expression of", gene_symbols[which(my_gene_rows == g)]),
          col = bar_colors,
          las = 2) # rotates label
  legend("topright", legend = names(condition_colors), 
         fill = condition_colors, 
         cex = 0.8, # shrinks text
         bty = "n") # removes the legend boarder
}


# ==== GSEA - Gene Set Enrichment Analysis ====

### GSEA: IR 7d vs Non IR 7d

# Load MSigDB gene sets (GO terms)
mouse_gene_sets <- msigdbr(species = "Mus musculus", # loads GO gene sets
                           category = "C5") 
gene_sets_list <- split(mouse_gene_sets$ensembl_gene, # organizes gene set into lists
                        mouse_gene_sets$gs_name)

# DESeq2 contrast: IR 7d vs Non IR 7d
res <- results(dds, contrast = c("condition", "IR 7d", "Non IR 7d"))

# Rank genes by log2FC / p-value
genes_ranked <- res$log2FoldChange / res$pvalue # effect by confidence
names(genes_ranked) <- rownames(res) # assigns names to vector
genes_ranked <- genes_ranked[!is.na(genes_ranked)] # takes away NAs
genes_ranked
# Run GSEA
fgsea_results <- fgsea(pathways = gene_sets_list, 
                       stats = genes_ranked, 
                       minSize = 15, # to reduce noise
                       maxSize = 500)

# Select top 10 enriched pathways
topPathways <- fgsea_results %>%
  dplyr::arrange(padj) %>%
  dplyr::slice(1:10)

# Plot - (Normalized Enrichment Score)
set.seed(23)
ggplot(topPathways, aes(x = reorder(pathway, NES), y = NES, fill = NES)) +
  geom_col(show.legend = FALSE) +
  scale_fill_gradient(low = "blue", high = "lightblue") +
  coord_flip() +
  labs(x = "Pathway", 
       y = "Normalized Enrichment Score (NES)", 
       title = "Top 10 Enriched Pathways: IR 7d vs Non IR 7d") +
  theme_minimal()


