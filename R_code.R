library(ggplot2)
library(dplyr)



# Read the input files
fang_data <- read_tsv("~/ISU/SPRING 2025/BCB4560/UNIX_Assignment/fang_et_al_genotypes.txt")
snp_data <- read_tsv("~/ISU/SPRING 2025/BCB4560/UNIX_Assignment/snp_position.txt")
#fang_transposed <- read_tsv("~/ISU/SPRING 2025/BCB4560/UNIX_Assignment/transposed_genotypes.txt")

# Data Inspection
# Count lines, words, and characters in fang_et_al_genotypes.txt
num_lines <- nrow(fang_data)
num_columns <- ncol(fang_data)
size_bytes <- file.info("~/ISU/SPRING 2025/BCB4560/UNIX_Assignment/fang_et_al_genotypes.txt")$size

# Count occurrences of '?'
num_missing_fang <- sum(fang_data == "?", na.rm = TRUE)
num_missing_snp <- sum(snp_data == "?", na.rm = TRUE)

# Data Processing
# Remove first two lines from transposed data
view(fang_transposed)
# Transpose the data
fang_transposed <- as.data.frame(t(fang_data), stringsAsFactors = FALSE)
view(fang_transposed)
# Convert first row to column names
colnames(fang_transposed) <- fang_transposed[3, ]  

# Remove the first row as it's now the column names
#fang_transposed <- fang_transposed[-c(1,2) ]

# Add original column names as a new first column
fang_transposed <- cbind(Original_Colnames = rownames(fang_transposed), fang_transposed)
fang_transposed <- fang_transposed[-c(1:3), ]
colnames(fang_transposed)[1] <- "SNP_ID"
view(fang_transposed)
# Reset row names
rownames(fang_transposed) <- NULL 

snp_info <- select(snp_data, SNP_ID, Chromosome, Position)
colnames(fang_transposed) <- make.unique(colnames(fang_transposed))
df_joined <- left_join(snp_info, fang_transposed, by = "SNP_ID")
view(df_joined)

# Extract relevant columns for maize
df_maize <- df_joined %>% select(SNP_ID, Chromosome = Chromosome, Position = Position, 
                                 starts_with("ZMMIL"), starts_with("ZMMLR"), starts_with("ZMMMR"))
view(df_maize)

# Sort by SNP position
df_maize_sorted_incre <- df_maize %>% arrange(Position)

# Create separate chromosome files
for (chr_num in 1:10) {
  df_chr <- df_maize_sorted_incre %>% filter(chromosome == chr_num)
  write_tsv(df_chr, paste0("Maize_incre_chromo", chr_num, ".txt"))
}

# Extract rows with 'multiple' or 'unknown' in Chromosome column
df_multiple <- df_maize_sorted_incre %>% filter(Chromosome == "multiple")
df_unknown <- df_maize_sorted_incre %>% filter(Chromosome == "unknown")
write_tsv(df_multiple, "Maize_chromo_multiple.txt")
write_tsv(df_unknown, "Maize_chromo_unknown.txt")


# maize decreasing chromosomes
df_maize_dash <- df_maize
df_maize_dash[df_maize_dash == "?/?"] <- "-/-"
#write.csv(df_maize_dash,"data_maize_dash.tsv")
df_maize_sorted_decre_dash <- df_maize_dash %>% arrange(desc(Position))
view(df_maize_sorted_decre_dash)

# Create separate chromosome files
for (chr_num in 1:10) {
  df_chr <- df_maize_sorted_decre_dash %>% filter(chromosome == chr_num)
  write_tsv(df_chr, paste0("Maize_decre_chromo", chr_num, ".txt"))
}


# Extract relevant columns for Teosinte
file_dir = "~/ISU/SPRING 2025/BCB4560/R_Assignment_Final/Teosinte/"
df_teosinte <- df_joined %>% select(SNP_ID, Chromosome = Chromosome, Position = Position, 
                                    starts_with("ZMPBA"), starts_with("ZMPIL"), starts_with("ZMPJA"))
view(df_teosinte)
#write_tsv(df_teosinte,paste0(file_dir,"data_teosinte.txt"))
# Sort by SNP position
df_teosinte_sorted_incre <- df_teosinte %>% arrange(Position)
#write_tsv(df_teosinte_sorted_incre,paste0(file_dir,"data_teosinte_sorted_incre.txt"))
# Create separate chromosome files
for (chr_num in 1:10) {
  df_chr <- df_teosinte_sorted_incre %>% filter(Chromosome == chr_num)
  write_tsv(df_chr, paste0(file_dir,"/Teosinte_Increasing_Chromo/Teosinte_incre_chromo", chr_num, ".txt"))
}

# Extract rows with 'multiple' or 'unknown' in Chromosome column
df_multiple <- df_teosinte_sorted_incre %>% filter(Chromosome == "multiple")
df_unknown <- df_teosinte_sorted_incre %>% filter(Chromosome == "unknown")
write_tsv(df_multiple, paste0(file_dir,"Teosinte_chromo_multiple.txt"))
write_tsv(df_unknown, paste0(file_dir,"Teosinte_chromo_unknown.txt"))

# teosinte decreasing chromosomes
df_teosinte_dash <- df_teosinte
df_teosinte_dash[df_teosinte_dash == "?/?"] <- "-/-"
#write_tsv(df_teosinte_dash,paste0(file_dir,"data_teosinte_dash.txt"))
df_teosinte_sorted_decre_dash <- df_teosinte_dash %>% arrange(desc(Position))
view(df_teosinte_sorted_decre_dash)
#write_tsv(df_teosinte_sorted_decre_dash, paste0(file_dir,"data_teosinte_sorted_decre_dash.txt"))
# Create separate chromosome files
for (chr_num in 1:10) {
  df_chr <- df_teosinte_sorted_decre_dash %>% filter(Chromosome == chr_num)
  write_tsv(df_chr, paste0(file_dir,"/Teosinte_Decreasing_Chromo/Teosinte_decre_chromo", chr_num, ".txt"))
}

#Part II Visualization

# Count SNPs per chromosome for maize
snp_count_maize <- df_maize %>%
  group_by(Chromosome) %>%
  summarise(SNP_Count = n()) %>%
  mutate(Group = "Maize")

# Count SNPs per chromosome for teosinte
snp_count_teosinte <- df_teosinte %>%
  group_by(Chromosome) %>%
  summarise(SNP_Count = n()) %>%
  mutate(Group = "Teosinte")

# Combine data
snp_counts <- bind_rows(snp_count_maize, snp_count_teosinte)

# Plot SNP count per chromosome
ggplot(snp_counts, aes(x = Chromosome, y = SNP_Count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "SNP Distribution Across Chromosomes",
       x = "Chromosome", y = "Number of SNPs") +
  scale_fill_manual(values = c("Maize" = "blue", "Teosinte" = "green")) +
  theme_classic()


library(dplyr)
library(tidyr)
library(purrr)  

# For optimized mapping

# Function to classify SNPs
classify_snp <- function(genotype) {
  switch(genotype,
          "?/?" = "Missing",
         ifelse(grepl("A/A|C/C|G/G|T/T", genotype), "Homozygous", "Heterozygous"))
}

# Apply classification to maize data
df_maize_long <- df_maize %>%
  pivot_longer(cols = -c(SNP_ID, Chromosome, Position),  # Select all except metadata columns
               names_to = "Sample", values_to = "Genotype") %>%
  mutate(SNP_Type = map_chr(Genotype, classify_snp))
# Faster than sapply
head(df_maize_long)
# Apply classification to teosinte data
df_teosinte_long <- df_teosinte %>%
  pivot_longer(cols = -c(SNP_ID, Chromosome, Position),  # Select all except metadata columns
               names_to = "Sample", values_to = "Genotype") %>%
  mutate(SNP_Type = map_chr(Genotype, classify_snp))

head(df_teosinte_long)
# Combine both datasets
df_combined <- bind_rows(
  df_maize_long %>% mutate(Group = "Maize"),
  df_teosinte_long %>% mutate(Group = "Teosinte")
)

# Count proportions of SNP types per group
df_long <- bind_rows(maize_long, teosinte_long)
df_snp_summary <- df_combined %>%
  group_by(Group, SNP_Type) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Proportion = Count / sum(Count))
# View results
print(df_snp_summary)


# Plot SNP Type distribution per Chromosome
ggplot(df_combined, aes(x = Chromosome, fill = SNP_Type)) +
  geom_bar(position = "fill") +  # "fill" makes proportions instead of counts
  facet_wrap(~ Group) +  # Separate plots for Maize & Teosinte
  labs(title = "Zygosity Distribution Across Chromosomes",
       x = "Chromosome", y = "Proportion",
       fill = "SNP Type") +
  theme_classic()



#Bar Plot – SNP Distribution by Group
ggplot(df_snp_summary, aes(x = Group, y = Proportion, fill = SNP_Type)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked proportional bar chart
  labs(title = "SNP Type Proportions in Maize vs Teosinte",
       x = "Group", y = "Proportion") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic()

#Stacked Bar Plot – Absolute SNP Counts
ggplot(df_snp_summary, aes(x = Group, y = Count, fill = SNP_Type)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar chart
  labs(title = "Total SNP Counts in Maize vs Teosinte",
       x = "Group", y = "SNP Count") +
  theme_classic()


# Bar Plot – SNP Distribution per Sample
ggplot(df_combined, aes(x = SNP_Type, y = ..count.., fill = Group)) +
  geom_bar(position = "dodge") +  
  labs(title = "SNP Type Distribution in Maize vs Teosinte",
       x = "SNP Type", y = "Count") +
  theme_classic()

