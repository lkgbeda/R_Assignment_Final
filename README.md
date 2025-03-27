Data Inspection, Processing and Visualization of SNP Genotype Data in R

This analysis focuses on inspecting, processing, and visualizing SNP genotype data for maize, teosinte, and an outgroup species using R. Two datasets were utilized:

fang_et_al_genotypes.txt – Contains SNP genotype data, with SNPs as rows and individual samples as columns.

snp_position.txt – Provides chromosomal positions and nucleotide locations for each SNP.

1. Data Inspection
The datasets were loaded and examined for structure, dimensions, and memory usage.

Missing values and inconsistencies (e.g., unknown chromosome assignments) were identified and addressed.

2. Data Processing 
The genotype dataset was transposed to make SNPs columns for easier merging.

The datasets were merged using SNP_ID as the common key, ensuring each SNP had corresponding chromosomal data.

The merged dataset was structured with SNP_ID, Chromosome, and Position as the first three columns, followed by genotype data.

SNPs were sorted in both ascending and descending order by chromosome.

The dataset was split into chromosome-specific files (1–10) for both maize and teosinte.

Missing data was represented as "?" in ascending-order files and "-" in descending-order files.

A total of 44 output files (20 for maize, 20 for teosinte) were generated.
Also, two unsorted files one for unknown chromosomes and one for multiple chromosomes, were generated for both maize and teosinte.

3. Data Visualization
SNP distribution plots were created using ggplot2, illustrating SNP counts per chromosome.

Genotype composition analysis assessed proportions of homozygous, heterozygous, and missing data.

Conclusion
This analysis demonstrated efficient handling of SNP genotype data, from preprocessing and merging to structured file generation and visualization. 
The workflow enhances reproducibility and provides a framework for large-scale genomic data processing in R.


Summary

The data inspection phase of the analysis was essential for understanding the structure, content, and quality of the input datasets before proceeding with further processing.
The two primary datasets, fang_et_al_genotypes.txt and snp_position.txt, were loaded into R as data frames and examined using various functions to assess their dimensions, memory usage, and overall structure.
The genotype dataset (fang_et_al_genotypes.txt) contained SNP information for maize, teosinte, and an outgroup species, with SNPs as rows and individual samples as columns.
Initial inspection involved displaying the first and last few rows of the dataset to check its format and contents.
Additionally, the object size function was used to determine the dataset's memory usage, ensuring efficient data handling during processing.
The snp_position.txt dataset provided chromosomal positions and nucleotide locations for each SNP.
Inspection of this dataset involved checking the number of columns and rows, as well as verifying the presence of key variables like SNP_ID, Chromosome, and Position.
Since this file served as the reference for merging genotype data with chromosomal locations, it was crucial to confirm its integrity and ensure no missing or ambiguous values in key columns.
During this phase, missing values and inconsistencies were also examined. 
Ensuring that column names were consistent and that all SNPs had a corresponding position entry was a necessary step before merging the datasets.
This R assignment focused on replicating a previous UNIX-based analysis of SNP genotype data while incorporating additional processing and visualization techniques. 
The primary objective was to handle and manipulate two datasets: 
fang_et_al_genotypes.txt, a published SNP dataset containing genetic information for maize, teosinte, and an outgroup species, and snp_position.txt, which provided the chromosomal positions for each SNP.
The task required merging, restructuring, and formatting these datasets to generate structured output files for downstream analysis.
The analysis began with loading the genotype and SNP position data into R as data frames.
Since the genotype dataset contained SNPs as rows and individual samples as columns, a critical preprocessing step involved transposing the fang_et_al_genotypes.txt dataset.
This is to reformat it so that SNPs became columns, making it compatible for merging with the SNP position file. 
Ensuring that column names remained unique after transposition was necessary to prevent errors in subsequent processing.
Once transposed, the genotype data was joined with the SNP position dataset using SNP_ID as the common key.
This integration provided a unified dataset containing each SNP's chromosome, nucleotide position, and genotype data for all samples.
The merged dataset was then formatted so that the first three columns represented the SNP ID, chromosome, and position, while the remaining columns contained genotype data.
For maize and teosinte, SNPs were sorted in both increasing and decreasing order based on their chromosomal positions.
Each dataset was then split into separate files for chromosomes 1 to 10, ensuring that missing data was represented by "?" in the increasing-order files and "-" in the decreasing-order files. 
This resulted in 40 structured output files 20 for maize and 20 for teosinte, each containing SNPs sorted by position and formatted correctly for further analysis.
Also, two unsorted files one for unknown chromosomes and one for multiple chromosomes, were generated for both maize and teosinte.
The assignment also involved visualizing SNP distributions and genetic characteristics using ggplot2.
Bar plots were generated to illustrate SNP counts per chromosome for both maize and teosinte, revealing differences in genetic variation across species.
Additional analysis was conducted to examine the proportions of homozygous, heterozygous, and missing data within each sample, providing insights into genetic diversity and sequencing quality.
In conclusion, this project demonstrated the ability to transpose, merge, clean, and restructure SNP genotype data in R, efficiently preparing it for downstream genomic analyses. 
The structured workflow not only automated file generation but also enabled effective visualization of genetic variation, offering a robust framework for large-scale SNP data processing.

