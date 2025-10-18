setwd("~/R/AI & Biotech/Bioinformatics")
getwd()

gc()  # Clear memory to free up space before analysis

#Installing and Loading Required Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Installing Bioconductor packages required for microarray analysis
BiocManager::install(c("limma", "AnnotationDbi", "hgu133plus2.db"))

#Installing CRAN packages for data manipulation and visualization
install.packages(c("dplyr", "tibble", "ggplot2", "pheatmap"))

# Load Bioconductor packages
library(AnnotationDbi)   # Handles annotation and probe–gene mapping
library(hgu133plus2.db)  # Annotation database for Affymetrix HG-U133 Plus 2.0 array
library(limma)           # Performs linear modeling and differential expression
library(dplyr)           # Simplifies data manipulation tasks
library(tibble)          
library(ggplot2)         # Used for plotting and visualization
library(pheatmap)        # Generates heatmaps for gene expression data

#Probe-to-Gene Mapping using AnnotationDbi

#Loading preprocessed expression and phenotype data
load("GSE16515.RData")

#checking annotation slot of the dataset
annotation(raw_data)

raw_data

#Displaying objects available in the annotation package
ls("package:hgu133plus2.db")

columns(hgu133plus2.db)
keytypes(hgu133plus2.db)

#Extracting probe IDs from processed microarray data
probe_ids <- rownames(processed_data)

#Mapping probe IDs to gene symbols using the platform annotation database
gene_symbols <- mapIds(
  hgu133plus2.db,          # Database used for mapping
  keys = probe_ids,        # Input probe IDs
  keytype = "PROBEID",     # Probe ID key type
  column = "SYMBOL",       # Desired annotation column (gene symbols)
  multiVals = "first"      # Return first match if multiple exist
)

#Converting mapping to a data frame and rename columns
gene_map_df <- gene_symbols %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::rename(SYMBOL = 2)

#Handling multiple probes mapping to a single gene

#Summarizing number of probes per gene symbol
duplicate_summary <- gene_map_df %>%
  group_by(SYMBOL) %>%
  summarise(probes_per_gene = n()) %>%
  arrange(desc(probes_per_gene))

#Identifying genes associated with multiple probes
duplicate_genes <- duplicate_summary %>%
  filter(probes_per_gene > 1)

sum(duplicate_genes$probes_per_gene)

#Merging annotation mapping with expression data

#Verifying if probe IDs in mapping correspond to expression data
all(gene_map_df$PROBEID == row.names(processed_data))

#Merging annotation (SYMBOL) with expression matrix
processed_data_df <- processed_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::mutate(SYMBOL = gene_symbols[PROBEID]) %>%
  dplyr::relocate(SYMBOL, .after = PROBEID)

#Removing probes without valid gene symbol annotation
processed_data_df <- processed_data_df %>%
  dplyr::filter(!is.na(SYMBOL))

#Selecting only numeric expression columns
expr_only <- processed_data_df %>%
  dplyr::select(-PROBEID, -SYMBOL)

#Collapsing multiple probes per gene using average expression
averaged_data <- limma::avereps(expr_only, ID = processed_data_df$SYMBOL)

dim(averaged_data)

#Converting averaged expression data to matrix format
data <- as.data.frame(averaged_data)
data <- data.matrix(data)
str(data)        # Structure check
is.numeric(data) # Confirm numeric matrix

#Differential Gene Expression Analysis

#Defining sample groups based on phenotype data
#Adjusting group labels according to dataset annotation

#Extracting condition info and create factor
groups <- ifelse(grepl("Tumor", phenotype_data$source_name_ch1, ignore.case = TRUE),
                 "cancer", "normal")

#Converting to factor with clear order
groups <- factor(groups, levels = c("normal", "cancer"))

class(groups)
levels(groups)

#Creating design matrix for linear modeling
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

#Fitting linear model to expression data
fit_1 <- lmFit(data, design)

#Defining contrast to compare cancer vs normal samples
contrast_matrix <- makeContrasts(cancer_vs_normal = cancer - normal,
                                 levels = design)

#Applying contrasts and computing moderated statistics
fit_contrast <- contrasts.fit(fit_1, contrast_matrix)

fit_2 <- eBayes(fit_contrast)

# Extracting list of differentially expressed genes (DEGs)
deg_results <- topTable(fit_2,
                        coef = "cancer_vs_normal",  # Specify contrast of interest
                        number = Inf,               # Return all genes
                        adjust.method = "BH")       # Benjamini-Hochberg correction

#Classifying DEGs into Upregulated, Downregulated, or Not Significant
deg_results$threshold <- as.factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated",
         "No")
))

#Subseting genes by regulation direction
upregulated <- subset(deg_results, threshold == "Upregulated")
downregulated <- subset(deg_results, threshold == "Downregulated")

#Combining both sets of DEGs
deg_updown <- rbind(upregulated, downregulated)

write.csv(deg_results, file = "Results/DEGs_Results.csv")
write.csv(upregulated, file = "Results/Upregulated_DEGs.csv")
write.csv(downregulated, file = "Results/Downregulated_DEGs.csv")
write.csv(deg_updown, file = "Results/Updown_DEGs.csv")

#Data Visualization

#Saving volcano plot as PNG
png("Result_Plots/volcano_plot.png", width = 2000, height = 1500, res = 300)

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "No" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10(P-value)",
       color = "Regulation")

dev.off()

#Heatmap of Top Differentially Expressed Genes

#Selecting top genes with smallest adjusted p-values
top_genes <- head(rownames(deg_updown[order(deg_updown$adj.P.Val), ]), 25)

#Subsetting averaged expression matrix for selected genes
heatmap_data <- data[top_genes, ]

#Generating unique column names per sample group for display
group_char <- as.character(groups)
heatmap_names <- ave(group_char, group_char, FUN = function(x) paste0(x, "_", seq_along(x)))

#Assigning formatted names to heatmap columns
colnames(heatmap_data) <- heatmap_names

#Saving heatmap as PNG
png("heatmap_top25_DEGs.png", width = 2000, height = 1500, res = 300)

#Generating heatmap without additional scaling
pheatmap(
  heatmap_data,
  scale = "none", # for already normalized data
  cluster_rows = FALSE,              # Disable row clustering
  cluster_cols = TRUE,               # Cluster samples
  show_rownames = TRUE,              # Display gene names
  show_colnames = TRUE,              # Display sample labels
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 6,
  fontsize_col = 8,
  main = "Top 25 Differentially Expressed Genes"
)

dev.off()
