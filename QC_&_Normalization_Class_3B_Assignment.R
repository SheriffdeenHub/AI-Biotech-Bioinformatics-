#                     Microarray Data Analysis
# =====================================================================

setwd("~/R/AI & Biotech/Bioinformatics")
getwd()
#Installing and Loading Required Packages

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

#Installing Bioconductor packages
BiocManager::install(c("GEOquery","affy","arrayQualityMetrics"))

#Installing CRAN packages for data manipulation
install.packages("dplyr")

#Loading Required Libraries
library(GEOquery)            
library(affy)                 
library(arrayQualityMetrics)  
library(dplyr)                

url("https://ftp.ncbi.nlm.nih.gov/")
#Downloading Series Matrix Files

gse_data <- getGEO("GSE16515", GSEMatrix = TRUE)

#Alternatively...
gse_data <- getGEO(filename = "GSE16515_series_matrix.txt.gz")
gse_data
#Extracting expression data matrix (genes/probes × samples)
#Rows corresponds to probes and columns corresponds to samples
expression_data <- exprs(gse_data$GSE16515_series_matrix.txt.gz)


#Extracting feature (probe annotation) data
#Rows corresponds to probes and columns corresponds to samples
feature_data <-  fData(gse_data$GSE16515_series_matrix.txt.gz)


#Extracting phenotype (sample metadata) data
#Rows corresponds to samples and columns corresponds to probes
phenotype_data <-  pData(gse_data$GSE16515_series_matrix.txt.gz)

#Checking missing values in sample annotation
sum(is.na(phenotype_data$source_name_ch1)) 


#To download Raw Data (CEL files) 

#Fetching GEO supplementry files
getGEOSuppFiles("GSE16515", baseDir = "Raw_Data", makeDirectory = TRUE)


#Untarring CEL files if compressed as .tar
untar("Data/GSE16515_RAW.tar", exdir = "Raw_Data/CEL_Files")


#Reading CEL files into R as an AffyBatch object
raw_data <- ReadAffy(celfile.path = "Raw_Data/CEL_Files")

raw_data 

#Quality Control (QC) Before Pre-processing 

arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)

#RMA (Robust Multi-array Average) Normalization 

normalized_data <- rma(raw_data)

#QC after data normalization 
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)

#Extractig normalized expression values into a data frame
processed_data <- as.data.frame(exprs(normalized_data))

dim(processed_data)   # Dimensions: number of probes × number of samples

#Filtering Low-Variance Transcripts (“soft” intensity based filtering)

#Calculating median intensity per probe across samples
row_median <- rowMedians(as.matrix(processed_data))

#Visualizing distribution of probe median intensities
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")

#Setting a threshold to remove low variance probes (dataset-specific, adjust accordingly)
threshold <- 3.5 
abline(v = threshold, col = "black", lwd = 2) 

#Selecting probes above threshold
indx <- row_median > threshold 
filtered_data <- processed_data[indx, ] 

#Renaming filtered expression data with sample metadata
colnames(filtered_data) <- rownames(phenotype_data)

#Overwriting processed data with filtered dataset
processed_data <- filtered_data 


#Phenotype Data Preparation

class(phenotype_data$source_name_ch1) 

#Defining experimental groups (normal vs cancer)
groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("pancreatic normal", "pancreatic tumor"),
                 label = c("normal", "cancer"))

class(groups)
levels(groups)
