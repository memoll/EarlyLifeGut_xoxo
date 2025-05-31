###############################################################
# Script to extract data from Excel sheets                    #
# Data: qPCR results - 16S                                    #
# Mona Parizadeh - May 2025                                   #
###############################################################

library(readxl)
library(dplyr)
library(stringr)

# Specify the directory where your Excel files are located
directory <- "~/Documents/xoxo_article/xoxo_qpcr_FINAL/16s"
# List all Excel files in the directory
excel_files <- list.files(directory, pattern = "\\.xlsx$|\\.xls$", full.names = TRUE); excel_files
# Create an empty list to store the data frames
result_sheets <- list()

# Loop through each Excel file
for (file in excel_files) {
  # Find the row index where "Well" is found in the first column of the "Results" sheet
  start_row <- which(readxl::read_excel(file, sheet = "Results")[, 1] == "Well")[1]
  
  # Read the "Results" sheet starting from the row with "Well" header
  result_sheet <- readxl::read_excel(file, sheet = "Results", skip = start_row - 1)
  
  # Keep only the "Sample Name" and "Quantity Mean" columns 
  result_sheet <- result_sheet[, c("Sample Name", "Quantity Mean")]
  
  # Add the result sheet to the list
  result_sheets[[file]] <- result_sheet
}
result_sheets[[1]]

# Combine all data frames into one
all_results <- bind_rows(result_sheets); dim(all_results)

# 
# Remove rows with "_2" or "_5" in the "Sample Name" column (concentration tests)
filtered_data <- all_results[!(grepl("_2|_5", all_results$`Sample Name`)), ]; dim(filtered_data)
# Remove "_1" from the "Sample Name" column
filtered_data$`Sample Name` <- gsub("_1", "", filtered_data$`Sample Name`); dim(filtered_data)

# Filter unique rows based on "Sample Name" (remove replicates since the mean quantities are the same for them)
unique_results <- filtered_data %>% 
  na.omit %>% #remove NAs
  filter(!duplicated(`Sample Name`)); dim(unique_results)

cleaned_results <- unique_results %>%
  mutate(
    `Sample Name` = stringr::str_replace(`Sample Name`, "-(\\d)$", "-0\\1"), # (\\d) matches a single digit (before or after the dash).
    `Sample Name` = paste0("XC", str_pad(`Sample Name`, width = 5, pad = "0", side = "left"))); dim(cleaned_results)

colnames(cleaned_results) <- c("SAMPLE_ID", "bac_QUANTITY_MEAN")
head(cleaned_results)
tail(cleaned_results)

# Sort the data frame based on the "New Sample Name" column
sorted_results <- cleaned_results %>% 
  arrange(`SAMPLE_ID`); dim(sorted_results)
head(sorted_results)
# Convert it to a dataframe
qpcr_data = as.data.frame(sorted_results); dim(qpcr_data)

# Now add the qPCR data to the metadata table:

# Import metadata
meta <- read.delim2("~/Documents/xoxo_article/files/metadata_xoxo_mi_18months_qpcr_euk.csv", header = T, sep = ",")
head(meta);dim(meta)

# First clean the metadata table
# remove mom samples
#meta_babies <- subset(meta, grepl("^XC", SAMPLE_ID)); dim(meta_babies)
meta_babies <- subset(meta, !grepl("^VM", SAMPLE_ID)); dim(meta_babies)

# Remove the second dash and everything after it
meta_babies$SAMPLE_ID <- gsub("^([^\\-]*\\-[^\\-]*)-.*$", "\\1", meta_babies$SAMPLE_ID)

# Add qpcr data to the meta dataframe
meta_qpcr <- merge(meta_babies, qpcr_data, by = "SAMPLE_ID", all.x = TRUE); dim(meta_qpcr)
meta_qpcr

sum(is.na(meta_qpcr$bac_QUANTITY_MEAN)) 

#check for the differences between sample IDs in both dataframes
sum(is.na(meta_qpcr$bac_QUANTITY_MEAN)) - (length(meta_babies$SAMPLE_ID) - length(qpcr_data$SAMPLE_ID))
setdiff(qpcr_data$SAMPLE_ID, meta_babies$SAMPLE_ID) #"XC04-19" 

# correction:

# Remove XC04-19 from qpcr_data$SAMPLE_ID (we don't have such a sample in the metadata)
qpcr_data_cor <- qpcr_data[qpcr_data$SAMPLE_ID != "XC04-19", ]; dim(qpcr_data_cor)

# merge again
meta_qpcr_cor <- merge(meta_babies, qpcr_data_cor, by = "SAMPLE_ID", all.x = TRUE); dim(meta_qpcr_cor)
head(meta_qpcr_cor)

sum(is.na(meta_qpcr_cor$bac_QUANTITY_MEAN)) 
sum(is.na(meta_qpcr_cor$bac_QUANTITY_MEAN)) - (length(meta_babies$SAMPLE_ID) - length(qpcr_data_cor$SAMPLE_ID))

# Write the data frame to a CSV file
write.csv(meta_qpcr_cor, file = "~/Documents/xoxo_article/files/metadata_xoxo_mi_18months_qpcr_euk_bac.csv", row.names = FALSE)

