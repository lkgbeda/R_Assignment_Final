# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)
library(purrr)  

# Load SNP genotype data and SNP position data from input files using read_tsv
fang_data <- read_tsv("~/ISU/SPRING 2025/BCB4560/UNIX_Assignment/fang_et_al_genotypes.txt")
snp_data <- read_tsv("~/ISU/SPRING 2025/BCB4560/UNIX_Assignment/snp_position.txt")

#Alternative means of loading the data
fang_data <- read_tsv("https://raw.githubusercontent.com/EEOB-BioData/BCB546_Spring2025/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt")
snp_data <- read_tsv("https://raw.githubusercontent.com/EEOB-BioData/BCB546_Spring2025/main/assignments/UNIX_Assignment/snp_position.txt")

# Data Inspection for Genotype data (fang_et_al_genotypes)

# Display the first few rows of the fang_data dataset to get an overview of the data.
head(fang_data)
# Display the last few rows of the fang_data dataset to inspect how the data ends.
tail(fang_data)
# Check the size of the dataset
cat("The fang_data object occupies approximately", 
    format(object.size(fang_data), units = "MB"), "of memory.\n")
# Calculate the number of rows in the fang_data dataset
num_lines_fang_data <- nrow(fang_data)
# Calculate the number of columns  in the fang_data dataset.
num_columns_fang_data <- ncol(fang_data)
cat("The fang_et_al_genotypes dataset contains", num_lines_fang_data, "rows (samples) and", 
    num_columns_fang_data, "columns .\n")

# Display the dimensions of the dataset
#dim(fang_data)

cat("The dimensions of the dataset are:", dim(fang_data)[1], "rows and", dim(fang_data)[2], "columns.\n")

num_missing_fang <- sum(fang_data == "?/?", na.rm = TRUE)
cat("The sum of the missing values is:", num_missing_fang)

# Count the number of missing values represented as "?/?" in the dataset
num_missing_fang <- sum(fang_data == "?/?", na.rm = TRUE)

# Print the total number of missing values
cat("The sum of the missing values in the fang_data:", num_missing_fang)


# Data Inspection for Genotype data (fang_et_al_genotypes)

# Display the structure of the dataset
cat("Here is the structure of the fang_data dataset:\n")
#str(fang_data)


# Data Inspection for SNP Position data (snp_position)


# Display the first few rows of the snp_data dataset to get an overview of the data.
head(snp_data)

# Display the last few rows of the snp_data dataset to inspect how the data ends.
tail(snp_data)

# Check the size of the dataset
cat("The snp_data object occupies approximately", 
    format(object.size(snp_data), units = "MB"), "of memory.\n")

# Calculate the number of rows in the snp_data dataset
num_lines_snp_data <- nrow(snp_data)

# Calculate the number of columns in the snp_data dataset.
num_columns_snp_data <- ncol(snp_data)

cat("The snp_position dataset contains", num_lines_snp_data, "rows (SNP positions) and", 
    num_columns_snp_data, "columns.\n")

# Display the dimensions of the dataset
cat("The dimensions of the dataset are:", dim(snp_data)[1], "rows and", dim(snp_data)[2], "columns.\n")

# Count the number of missing values represented as "?/?" in the dataset
num_missing_snp <- sum(snp_data == "?/?", na.rm = TRUE)

# Print the total number of missing values
cat("The sum of the missing values in snp_data is:", num_missing_snp)

# View the snp_position dataset (This will open it in RStudio's Data Viewer)
cat("Opening the snp_position in the Viewer...\n")
View(snp_data) 

# Display the structure of the dataset
cat("Here is the structure of the snp_position:\n")
str(snp_data)


# Data Processing


# Transpose the fang_data dataset to switch rows and columns
fang_transposed <- as.data.frame(t(fang_data), stringsAsFactors = FALSE)

# View the transposed data (optional, for checking)
view(fang_transposed)

# Convert the third row into column names (assuming row 3 contains SNP names)
colnames(fang_transposed) <- fang_transposed[3, ]  

# Add the original row names as a new first column to preserve them
fang_transposed <- cbind(Original_Colnames = rownames(fang_transposed), fang_transposed)

# Remove the first 3 rows as they are not needed after setting column names
fang_transposed <- fang_transposed[-c(1:3), ]

# Rename the first column to "SNP_ID" for clarity
colnames(fang_transposed)[1] <- "SNP_ID"

# View the cleaned transposed data (optional, for checking)
view(fang_transposed)

# Reset row names to ensure they are sequential and start from 1
rownames(fang_transposed) <- NULL

# Joining the transposed genotype dat and the SNP position data by SNP 

# Select relevant columns (SNP_ID, Chromosome, Position) from snp_data  
snp_info <- select(snp_data, SNP_ID, Chromosome, Position)  

# Ensure unique column names in fang_transposed to avoid duplication issues  
colnames(fang_transposed) <- make.unique(colnames(fang_transposed))  

# Perform a left join to merge snp_info with fang_transposed based on SNP_ID  
df_joined <- left_join(snp_info, fang_transposed, by = "SNP_ID")  

# View the resulting dataframe in RStudio  
view(df_joined)

# Save the joined dataset as a tab-separated text file  
write.table(df_joined, "joined_fang_snp_data.txt", sep = "\t", row.names = FALSE, quote = FALSE) 


# Extract relevant columns for maize, including SNP_ID, Chromosome, Position, and columns starting with ZMMIL, ZMMLR, or ZMMMR  
df_maize <- df_joined %>% select(SNP_ID, Chromosome = Chromosome, Position = Position,  
                                 starts_with("ZMMIL"), starts_with("ZMMLR"), starts_with("ZMMMR"))  

# View the maize-specific dataset in RStudio  
view(df_maize)    

write.table(df_joined, "Maize/maize_data.txt", sep = "\t", row.names = FALSE, quote = FALSE) 

# Extract rows with 'multiple' or 'unknown' in Chromosome column for Maize

# Extract rows with 'multiple' in the Chromosome column  
df_multiple <- df_maize %>% filter(Chromosome == "multiple")  

# Extract rows with 'unknown' in the Chromosome column  
df_unknown <- df_maize %>% filter(Chromosome == "unknown")  

# Save the 'multiple' chromosome data to a tab-separated text file  
write.table(df_multiple, "Maize/maize_chromo_multiple.txt",sep = "\t", row.names = FALSE, quote = FALSE)  

# Save the 'unknown' chromosome data to a tab-separated text file  
write.table (df_unknown, "Maize/maize_chromo_unknown.txt",sep = "\t", row.names = FALSE, quote = FALSE)  

# Sort by ascending order of SNP position

# Sort by SNP position in increasing order  
# Sort by SNP position               
df_maize_chr_filtered <- data.frame()  # Initialize an empty data frame to store sorted SNPs

for (chr_num in 1:10) {  # Loop through chromosome numbers 1 to 10
  df_chr <- df_maize %>%  # Filter the data for the current chromosome
    filter(Chromosome == chr_num)  # Ensure column name matches the dataset structure
  
  df_maize_chr_filtered <- bind_rows(df_maize_chr_filtered, df_chr)  # Append filtered data for each chromosome
}

# Convert Position column to numeric type and sort in increasing order
df_maize_sorted_incre <- df_maize_chr_filtered %>%
  mutate(Position = suppressWarnings(as.numeric(Position))) %>%  # Convert Position to numeric type for accurate sorting
  arrange(Position)  # Sort the data by Position in increasing order

# View the sorted dataset in RStudio  
view(df_maize_sorted_incre)  

# Save the sorted maize dataset as a tab-separated text file  
write.table(df_maize_sorted_incre, "Maize/maize_sorted_incre_data.txt", sep = "\t", row.names = FALSE, quote = FALSE)  


# Create separate chromosome files for Increasing order of SNP Position for maize 

#Increasing SNP Position (Maize) 
# Create separate chromosome files for chromosomes 1 to 10  
for (chr_num in 1:10) {  
  # Filter rows for each chromosome number  
  df_chr <- df_maize_sorted_incre %>% filter(Chromosome == chr_num)  
  
  # Save the filtered chromosome data to a tab-separated text file without row names or quotes  
  write.table(df_chr, paste0("Maize/Maize_Increasing_Chromo/Maize_incre_chromo", chr_num, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)  
}


# Create separate chromosome files for decreasing order of SNP Position for maize 

# Repalce "?/?" by "-/-" and Sort by decreasing order of SNP position
# maize decreasing SNP Positions 
# Create a copy of the maize dataset  
df_maize_chr_filtered_dash <- df_maize_chr_filtered  

# Replace all occurrences of "?/?" with "-/-" in the dataset  
df_maize_chr_filtered_dash[df_maize_chr_filtered_dash == "?/?"] <- "-/-"  

# View the modified dataset in RStudio  
view(df_maize_chr_filtered_dash)  

# Save the modified dataset to a tab-separated text file  
write.table(df_maize_chr_filtered_dash, "Maize/maize_chr_filtered_dash.txt", sep = "\t", row.names = FALSE, quote = FALSE)  

# Sort the modified dataset by Position in descending order  
#df_maize_sorted_decre_dash <- df_maize_dash %>% arrange(desc(Position))  
# Convert Position column to numeric type and sort in decreasing order
df_maize_sorted_decre_dash <- df_maize_chr_filtered_dash %>%
  mutate(Position = suppressWarnings(as.numeric(Position))) %>%  # Convert Position to numeric type for accurate sorting
  arrange(desc(Position))   # Sort the data by Position in decrasing order
# View the sorted dataset in RStudio  
view(df_maize_sorted_decre_dash)  

# Save the sorted dataset to a tab-separated text file  
write.table(df_maize_sorted_decre_dash, "Maize/maize_sorted_decre_dash.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Generate separate chromosome files for decreasing order of SNP Position for maize

# Create separate chromosome files for chromosomes 1 to 10  
for (chr_num in 1:10) {  
  # Filter rows for each chromosome number  
  df_chr <- df_maize_sorted_decre_dash %>% filter(Chromosome == chr_num)  
  
  # Save the filtered chromosome data to a tab-separated text file  
  write.table(df_chr, paste0("Maize/Maize_Decreasing_Chromo/Maize_decre_chromo", chr_num, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)  
}



# Teosinte
#Similar to the codes used in maize, the codes below do same in teosinte. The summary of what the codes is below with detailed explanation seen in maize 

# Extract relevant columns for Teosinte 

# Define the directory path for saving the output file
file_dir = "~/ISU/SPRING 2025/BCB4560/R_Assignment_Final/Teosinte/"

# Extract relevant columns for Teosinte from the joined dataframe
# Selecting SNP_ID, Chromosome, Position, and columns that start with "ZMPBA", "ZMPIL", or "ZMPJA"
df_teosinte <- df_joined %>% select(SNP_ID, Chromosome = Chromosome, Position = Position, 
                                    starts_with("ZMPBA"), starts_with("ZMPIL"), starts_with("ZMPJA"))

# View the extracted dataframe in RStudio
view(df_teosinte)

# Save the extracted Teosinte data to a text file with tab-separated values
write.table(df_teosinte, paste0(file_dir, "data_teosinte.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


# Extract rows with 'multiple' or 'unknown' in Chromosome column (Teosinte)


# Extract rows where the Chromosome column has the value "multiple"
df_multiple <- df_teosinte %>% filter(Chromosome == "multiple")

# Extract rows where the Chromosome column has the value "unknown"
df_unknown <- df_teosinte %>% filter(Chromosome == "unknown")

# Save the extracted data for "multiple" chromosomes to a text file with tab-separated values
write.table(df_multiple, paste0(file_dir, "Teosinte_chromo_multiple.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Save the extracted data for "unknown" chromosomes to a text file with tab-separated values
write.table(df_unknown, paste0(file_dir, "Teosinte_chromo_unknown.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Sort by ascending order of SNP position (Teosinte)

# Sort by SNP position               
df_teosinte_chr_filtered <- data.frame()  # Initialize an empty data frame to store sorted SNPs

for (chr_num in 1:10) {  # Loop through chromosome numbers 1 to 10
  df_chr <- df_teosinte %>%  # Filter the data for the current chromosome
    filter(Chromosome == chr_num)  # Ensure column name matches the dataset structure
  
  df_teosinte_chr_filtered <- bind_rows(df_teosinte_chr_filtered, df_chr)  # Append filtered data for each chromosome
}

# Convert Position column to numeric type and sort in increasing order
df_teosinte_sorted_incre <- df_teosinte_chr_filtered %>%
  mutate(Position = suppressWarnings(as.numeric(Position))) %>%  # Convert Position to numeric type for accurate sorting
  arrange(Position)  # Sort the data by Position in ascending order
# View the sorted dataframe in RStudio
view(df_teosinte_sorted_incre)

# Save the sorted Teosinte data to a text file with tab-separated values
write.table(df_teosinte_sorted_incre, paste0(file_dir, "data_teosinte_sorted_incre.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)


# Increasing SNP Position (Teosinte)
# Create separate chromosome files for Teosinte (Increasing SNP Positions)
# Create separate files for each chromosome (1 to 10) in Teosinte dataset, sorted by increasing SNP positions
for (chr_num in 1:10) {
  
  # Filter the dataframe to include only rows where Chromosome matches the current chromosome number
  df_chr <- df_teosinte_sorted_incre %>% filter(Chromosome == chr_num)
  
  # Save the filtered data to a text file, named according to the chromosome number
  write.table(df_chr, paste0(file_dir, "/Teosinte_Increasing_Chromo/Teosinte_incre_chromo", chr_num, ".txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Decreasing SNP Position (Teosinte)
# Repalce "?/?" by "-/-" and Sort by decreasing order of SNP position 

# Create a copy of the Teosinte dataframe to modify the genotype format
df_teosinte_chr_filtered_dash <- df_teosinte_chr_filtered
# Replace all occurrences of "?/?" with "-/-" in the dataframe
df_teosinte_chr_filtered_dash[df_teosinte_chr_filtered_dash == "?/?"] <- "-/-"

# View the modified dataframe in RStudio
view(df_teosinte_chr_filtered_dash)

# Save the modified Teosinte data to a text file with tab-separated values
write.table(df_teosinte_chr_filtered_dash, paste0(file_dir, "teosinte_chr_filtered_dash.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Sort the modified dataframe by SNP position in decreasing order
df_teosinte_sorted_decre_dash <- df_teosinte_chr_filtered_dash %>% arrange(desc(Position))

# Convert Position column to numeric type and sort in decreasing order
df_teosinte_sorted_decre_dash <- df_teosinte_chr_filtered_dash %>%
  mutate(Position = suppressWarnings(as.numeric(Position))) %>%  # Convert Position to numeric type for accurate sorting
  arrange(Position)  # Sort the data by Position in ascending order
# View the sorted dataframe in RStudio
view(df_teosinte_sorted_decre_dash)

# Save the sorted Teosinte data to a text file with tab-separated values
write.table(df_teosinte_sorted_decre_dash, paste0(file_dir, "data_teosinte_sorted_decre_dash.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Create Separate Files for Each Chromosome (Decreasing SNP Positions)
# Create separate chromosome files for Teosinte (Decreasing SNP Positions)
# Create separate files for each chromosome (1 to 10) in the Teosinte dataset, sorted by decreasing SNP positions
for (chr_num in 1:10) {
  
  # Filter the dataframe to include only rows where Chromosome matches the current chromosome number
  df_chr <- df_teosinte_sorted_decre_dash %>% filter(Chromosome == chr_num)
  
  # Save the filtered data to a text file, named according to the chromosome number, using tab-separated values
  write_tsv(df_chr, paste0(file_dir, "/Teosinte_Decreasing_Chromo/Teosinte_decre_chromo", chr_num, ".txt"))
}

#Part II Visualization 

# Bar plot of SNP distribution by Chromosome in Maize
ggplot(data = df_maize, aes(x = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown")), 
                            fill = factor(Chromosome))) + 
  geom_bar() +  
  ggtitle("Distribution of SNPs by Chromosome in Maize") +
  
  # Improve theme and text positioning
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  # Improve axis labels
  xlab("Chromosome") +
  ylab("SNP Count") +
  
  # Improve legend
  scale_fill_viridis_d(name = "Chromosome")  # Better color mapping

# Bar plot of SNP distribution by Chromosome in Teosinte
ggplot(data = df_teosinte, aes(x = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown")), 
                               fill = factor(Chromosome))) + 
  geom_bar() +  
  ggtitle("Distribution of SNPs by Chromosome in Teosinte") +
  
  # Improve theme and text positioning
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  # Improve axis labels
  xlab("Chromosome") +
  ylab("SNP Count") +
  
  # Improve legend
  scale_fill_viridis_d(name = "Chromosome")  # Better color mapping

# Count SNPs per chromosome for maize

# Count the number of SNPs per chromosome for the maize dataset
snp_count_maize <- df_maize %>%
  
  # Group the data by the Chromosome column
  group_by(Chromosome) %>%
  
  # Count the number of SNPs in each chromosome group
  summarise(SNP_Count = n()) %>%
  
  # Add a new column to indicate the dataset source as "Maize"
  mutate(Species = "Maize")

# Count SNPs per chromosome for teosinte

# Count the number of SNPs per chromosome for the Teosinte dataset
snp_count_teosinte <- df_teosinte %>%
  
  # Group the data by the Chromosome column
  group_by(Chromosome) %>%
  
  # Count the number of SNPs in each chromosome group
  summarise(SNP_Count = n()) %>%
  
  # Add a new column to indicate the dataset source as "Teosinte"
  mutate(Species = "Teosinte")


#Combining Data from Two Dataframes

# Combine SNP count data from both Maize and Teosinte datasets
snp_counts <- bind_rows(snp_count_maize, snp_count_teosinte)


#Creating the Plot

# Plot SNP count per chromosome for Maize and Teosinte
ggplot(snp_counts, aes(x = Chromosome, y = SNP_Count, fill = Species)) +
  
  # Create a bar plot with "Chromosome" on the x-axis and "SNP_Count" on the y-axis
  # The bars are grouped by "Group" (Maize vs. Teosinte)
  geom_bar(stat = "identity", position = "dodge") +
  
  # Add labels to the plot
  labs(title = "SNP Distribution Across Chromosomes",
       x = "Chromosome", y = "Number of SNPs") +
  
  # Set custom colors for the "Maize" and "Teosinte" groups
  scale_fill_manual(values = c("Maize" = "blue", "Teosinte" = "green")) +
  
  # Ensure chromosome labels appear in the correct order
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown"))+
  # Use a clean and simple theme for the plot
  theme_classic()



# Scatter plot of SNP count per chromosome for Maize and Teosinte
ggplot(snp_counts, aes(x = Chromosome, y = SNP_Count, color = Species)) +
  
  # Create a scatter plot with "Chromosome" on the x-axis and "SNP_Count" on the y-axis
  geom_point(size = 3, position = position_jitter(width = 0.2, height = 0)) +  # Add slight jitter for better visibility
  
  # Add labels to the plot
  labs(title = "SNP Distribution Across Chromosomes",
       x = "Chromosome", y = "Number of SNPs") +
  
  # Set custom colors for the "Maize" and "Teosinte" groups
  scale_color_manual(values = c("Maize" = "blue", "Teosinte" = "green")) +
  
  # Ensure chromosome labels appear in the correct order
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown")) +
  
  # Use a clean and simple theme for the plot
  theme_classic()


# Scatter plot of SNP distribution by Chromosome in Maize
ggplot(data = df_maize, aes(x = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown")), 
                            y = Position, 
                            color = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown")))) + 
  geom_point(alpha = 0.7, size = 1.5) +  # Adjust transparency and size
  
  ggtitle("Distribution of SNPs by Chromosome in Maize") +
  
  # Improve axis labels
  xlab("Chromosome") +
  ylab("Genomic Position") +
  
  # Improve theme and text positioning
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)) +
  
  # Explicitly define color mapping and ensure correct legend order
  scale_color_manual(name = "Chromosome",
                     values = c("1" = "#F8766D", "2" = "#B79F00", "3" = "#7CAE00", "4" = "#00BE67",
                                "5" = "#00BFC4", "6" = "#00A9FF", "7" = "#619CFF", "8" = "#C77CFF",
                                "9" = "#D89000", "10" = "#E76BF3", "multiple" = "#FF61CC", "unknown" = "#F564E3"),
                     breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown")) +
  
  guides(color = guide_legend(ncol = 1))  # Ensure legend is in one column


# Scatter plot of the distribution of SNPs by Chromosome in Teosinte
ggplot(data = df_teosinte, aes(x = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown")), 
                               y = Position, color = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown")))) + 
  geom_point(alpha = 0.7, size = 1.5) +  # Adjust point transparency and size for better visibility
  ggtitle("Distribution of SNPs by Chromosome in Teosinte") +
  
  # Improve axis labels
  xlab("Chromosome") +
  ylab("Genomic Position") +
  
  # Improve theme and text positioning
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)) +
  
  # Explicitly define color mapping and ensure correct legend order
  scale_color_manual(name = "Chromosome",
                     values = c("1" = "#F8766D", "2" = "#B79F00", "3" = "#7CAE00", "4" = "#00BE67",
                                "5" = "#00BFC4", "6" = "#00A9FF", "7" = "#619CFF", "8" = "#C77CFF",
                                "9" = "#D89000", "10" = "#E76BF3", "multiple" = "#FF61CC", "unknown" = "#F564E3"),
                     breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown")) +
  
  guides(color = guide_legend(ncol = 1))  # Ensure legend is in one column


# For optimized mapp

# Function to classify SNPs
classify_snp <- function(genotype) {
  switch(genotype,
         "?/?" = "Missing",
         ifelse(grepl("A/A|C/C|G/G|T/T", genotype), "Homozygous", "Heterozygous"))
}
df_maize_long <- df_maize  %>%
  pivot_longer(cols = -c(SNP_ID, Chromosome, Position),  # Exclude metadata columns
               names_to = "Group", values_to = "Genotype") %>%
  
  # Standardize Group names (remove .1, .2 suffixes)
  mutate(Group = str_remove(Group, "\\.\\d+$")) %>%
  
  filter(str_starts(Group, "ZMMIL") | str_starts(Group, "ZMMLR") | str_starts(Group, "ZMMMR")) %>%  
  mutate(Species = "Maize", SNP_Type = map_chr(Genotype, classify_snp))    # Assign species and classify SNPs
head(df_maize_long)

# Convert teosinte data to long format
df_teosinte_long <- df_teosinte %>%
  pivot_longer(cols = -c(SNP_ID, Chromosome, Position),  # Exclude metadata columns
               names_to = "Group", values_to = "Genotype") %>%
  
  # Standardize Group names (remove .1, .2 suffixes)
  mutate(Group = str_remove(Group, "\\.\\d+$")) %>%
  
  filter(str_starts(Group, "ZMPBA") | str_starts(Group, "ZMPIL") | str_starts(Group, "ZMPJA")) %>%
  mutate(Species = "Teosinte", SNP_Type = map_chr(Genotype, classify_snp))    # Assign species and classify SNPs
head(df_teosinte_long)

# Combine both datasets into one
df_combined <- bind_rows(df_maize_long, df_teosinte_long)
head(df_combined)
# Merge with filtered SNP metadata
#df_snp_filtered <- merge(df_combined, by = "SNP_ID")

# Summarize SNP types by species (Maize or Teosinte)
df_snp_summary <- df_combined %>%
  group_by(Species, SNP_Type, Group) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Proportion = Count / sum(Count))  # Calculate proportion of SNPs

# Print SNP summary table
print(df_snp_summary)



# Plot SNP distribution by chromosome for maize vs teosinte
by_species_plot <- ggplot(data = df_combined) +
  geom_bar(mapping = aes(x = Chromosome, fill = Species)) +
  xlab("Chromosome") + 
  ylab("SNPs") +
  ggtitle("Single Nucleotide Polymorphism per Sample") +
  theme_classic()

# Display the plot
print(by_species_plot)

### Distribution of Homozygous Genotypes by Group in Maize
# Select only homozygous genotypes
maize_homozygous_snp <- df_maize_long %>%
  filter(SNP_Type == "Homozygous")

# Convert "Group" into a factor with the correct order
maize_homozygous_snp$Group <- factor(
  maize_homozygous_snp$Group,
  levels = c("ZMMIL", "ZMMLR", "ZMMMR")
)

# Exclude NA groups from the plot
maize_homozygous_snp <- maize_homozygous_snp %>%
  filter(!is.na(Group))


# Plot: Homozygous Genotype Distribution by Group
maize_homozygous_plot <- ggplot(maize_homozygous_snp, aes(x = Group, fill = Genotype)) +
  geom_bar(position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = c("A/A" = "#1b9e77", "T/T" = "#d95f02", "C/C" = "#7570b3", "G/G" = "#e7298a"),
    name = "Genotype"
  ) +
  xlab("Group") +
  ylab("Homozygous SNP Count") +
  ggtitle("Distribution of Homozygous Genotypes by Group in Maize") +
  theme_classic()
maize_homozygous_plot


### Distribution of Homozygous Genotypes by Group in Teosinte
# Select only homozygous genotypes
teosinte_homozygous_snp <- df_teosinte_long  %>%
  filter(SNP_Type == "Homozygous")

# Convert "Group" into a factor with the correct order
teosinte_homozygous_snp$Group <- factor(
  teosinte_homozygous_snp$Group,
  levels = c("ZMPBA", "ZMPIL", "ZMPJA")
)

# Exclude NA groups from the plot
teosinte_homozygous_snp <- teosinte_homozygous_snp %>%
  filter(!is.na(Group))

# Plot: Homozygous Genotype Distribution by Group
teosinte_homozygous_plot <- ggplot(teosinte_homozygous_snp, aes(x = Group, fill = Genotype)) +
  geom_bar(position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = c("A/A" = "#1b9e77", "T/T" = "#d95f02", "C/C" = "#7570b3", "G/G" = "#e7298a"),
    name = "Genotype"
  ) +
  xlab("Group") +
  ylab("Homozygous SNP Count") +
  ggtitle("Distribution of Homozygous Genotypes by Group in Teosinte") +
  theme_classic()
teosinte_homozygous_plot

# Select only heterozygous genotypes in maize groups
maize_heterozygous_snp <- df_maize_long %>%
  filter(SNP_Type == "Heterozygous")

# Convert "Group" into a factor with the correct order
maize_heterozygous_snp$Group <- factor(
  maize_heterozygous_snp$Group,
  levels = c("ZMMIL", "ZMMLR", "ZMMMR")
)


# Exclude NA groups from the plot
maize_heterozygous_snp <- maize_heterozygous_snp %>%
  filter(!is.na(Group))

# Plot: Heterozygous Genotype Distribution by Group
maize_heterozygous_plot <- ggplot(maize_heterozygous_snp, aes(x = Group, fill = Genotype)) +
  geom_bar(position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = c("C/G" = "#1b9e77", "G/T" = "#d95f02", "A/G" = "#7570b3",
               "C/T" = "#e7298a", "A/T" = "#66a61e", "A/C" = "#e6ab02"),
    name = "Genotype"
  ) +
  xlab("Group") +
  ylab("Heterozygous SNP Count") +
  ggtitle("Distribution of Heterozygous Genotypes by Group in Maize") +
  theme_classic()

maize_heterozygous_plot

# Select only heterozygous genotypes
teosinte_heterozygous_snp <- df_combined %>%
  filter(SNP_Type == "Heterozygous")

# Convert "Group" into a factor with the correct order
teosinte_heterozygous_snp$Group <- factor(
  teosinte_heterozygous_snp$Group,
  levels = c("ZMPBA", "ZMPIL", "ZMPJA")
)

# Exclude NA groups from the plot
teosinte_heterozygous_snp <- teosinte_heterozygous_snp %>%
  filter(!is.na(Group))

# Plot: Heterozygous Genotype Distribution by Group
teosinte_heterozygous_plot <- ggplot(teosinte_heterozygous_snp, aes(x = Group, fill = Genotype)) +
  geom_bar(position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = c("C/G" = "#1b9e77", "G/T" = "#d95f02", "A/G" = "#7570b3",
               "C/T" = "#e7298a", "A/T" = "#66a61e", "A/C" = "#e6ab02"),
    name = "Genotype"
  ) +
  xlab("Group") +
  ylab("Heterozygous SNP Count") +
  ggtitle("Distribution of Heterozygous Genotypes by Group in Teosinte") +
  theme_classic()

teosinte_heterozygous_plot

# Count proportions of SNP types per group
# Summarize SNP types by group (Maize or Teosinte)
df_snp_summary <- df_combined %>%
  
  # Group the data by "Group" (Maize or Teosinte) and "SNP_Type" (Homozygous, Heterozygous, Missing)
  group_by(Species, SNP_Type) %>%
  
  # Count the number of occurrences for each SNP type in each group
  summarise(Count = n(), .groups = "drop") %>%
  
  # Calculate the proportion of each SNP type within each group
  mutate(Proportion = Count / sum(Count))

# View the summarized results (SNP counts and proportions)
print(df_snp_summary)


# Plot SNP Type distribution per Chromosome

ggplot(df_combined, aes(x = Chromosome, fill = SNP_Type)) +
  geom_bar(position = "fill") +  # "fill" makes proportions instead of counts
  facet_wrap(~ Species) +  # Separate plots for Maize & Teosinte
  labs(title = "Zygosity Distribution Across Chromosomes",
       x = "Chromosome", y = "Proportion",
       fill = "SNP Type") +
  theme_classic()


#Proportional Bar Plot (SNP Distribution by Group)

#Bar Plot – SNP Distribution by Group
ggplot(df_snp_summary, aes(x = Species, y = Proportion, fill = SNP_Type)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked proportional bar chart
  labs(title = "SNP Type Proportions in Maize vs Teosinte",
       x = "Group", y = "Proportion") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic()


# Stacked Bar Plot (Absolute SNP Counts):

#Stacked Bar Plot – Absolute SNP Counts
ggplot(df_snp_summary, aes(x = Species, y = Count, fill = SNP_Type)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar chart
  labs(title = "Total SNP Counts in Maize vs Teosinte",
       x = "Group", y = "SNP Count") +
  theme_classic()

#SNP Distribution Bar Plot per Sample:
# Bar Plot – SNP Distribution per Species
ggplot(df_combined, aes(x = SNP_Type, y = ..count.., fill = Species)) +
  geom_bar(position = "dodge") +  
  labs(title = "SNP Type Distribution in Maize vs Teosinte",
       x = "SNP Type", y = "Count") +
  theme_classic()

# Bar Plot – SNP Distribution per Group in Maize
ggplot(df_maize_long, aes(x = SNP_Type, y = ..count.., fill = Group)) +
  geom_bar(position = "dodge") +  
  labs(title = "SNP Type Distributionper Group in Maize",
       x = "SNP Type", y = "Count") +
  theme_classic()

# Bar Plot – SNP Distribution per Group in Teosinte
ggplot(df_teosinte_long, aes(x = SNP_Type, y = ..count.., fill = Group)) +
  geom_bar(position = "dodge") +  
  labs(title = "SNP Type Distribution per Group in Teosinte",
       x = "SNP Type", y = "Count") +
  theme_classic()


# Heatmap of SNP distribution across chromosomes for Maize and Teosinte
ggplot(snp_counts, aes(x = Chromosome, y = Species, fill = SNP_Count)) +
  
  # Create a heatmap using geom_tile()
  geom_tile(color = "white ") +
  
  # Add labels to the plot
  labs(title = "SNP Distribution Heatmap Across Chromosomes",
       x = "Chromosome", y = "Species", fill = "SNP Count") +
  
  # Use a gradient color scale to represent SNP counts
  scale_fill_gradient(low = "yellow", high = "blue") +
  
  # Ensure chromosome labels appear in the correct order
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown")) +
  
  # Use a clean and simple theme
  theme_classic()















