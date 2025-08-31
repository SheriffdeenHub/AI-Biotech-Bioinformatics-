setwd("C:/Users/sheri/Documents/R/AI & Biotech/Bioinformatics")
getwd()

# -----------------------
# Typical loop workflow:
# -----------------------

# Define the input folder (where raw data files are stored) and the output folder (where results will be saved).
# Specify the list of files that need to be processed.
# Prepare an empty list in R to store results for later use.
# For each file in the list:
#          Import the data into R.
#          Check and handle missing values (NA).
#          Classify genes using the classify_gene function.
#          Save the processed results both as a CSV file and in the R results list.

# ----------------------------------------------------
# Gene Classification of two dataset within loop

# Defining input and output folders
input_dir <- "Data" 
output_dir <- "Results"

# Listing which files to process
files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv") 

# Preparing empty list to store results in R 
result_list <- list()

# Creating classify_gene() function
classify_gene <- function(logFC, padj) {
  ifelse(logFC > 1 & padj < 0.05, "Upregulated",
         ifelse(logFC < -1 & padj < 0.05, "Downregulated", "Not_Significant"))
}

# Creating the loop

for (file_names in files_to_process) {
  cat("\nProcessing:", file_names, "\n")
  
  input_file_path <- file.path(input_dir, file_names)
  
  # Import dataset
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported. Checking for missing values...\n")
  
  # handling missing values
  
  if("padj" %in% names(data)){
    missing_count <- sum(is.na(data$padj))
    
    cat("Missing values in 'padj':", missing_count, "\n")
    data$padj[is.na(data$padj)] <- 1
  }
  
  if("logFC" %in% names(data)){
    missing_count <- sum(is.na(data$logFC))
    
    cat("Missing values in 'logFC':", missing_count, "\n")
    data$logFC[is.na(data$logFC)] <- mean(data$logFC, na.rm = TRUE)
  }
  
  # Classifying genes
  data$status <- classify_gene(data$logFC, data$padj)
  cat("Gene has been classified successfully.\n")
  
  # Saving results in R
  result_list[[file_names]] <- data 
  
  # Saving results in Results folder
  output_file_path <- file.path(output_dir, paste0("Gene_classification_results", file_names))
  write.csv(data, output_file_path, row.names = FALSE)
  cat("Results saved to:", output_file_path, "\n")
}

# Repeating the loop to process all files.

results_1 <- result_list[[1]] 
results_2 <- result_list[[2]]

# Printing summary counts of significant, upregulated, and downregulated genes
cat("\nSummary for", file_names, ":\n")
print(table(data$status))






