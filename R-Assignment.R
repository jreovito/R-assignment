# R Assignment

# Data Inspection

setwd("/Users/jordyn/Desktop/BCB546_Spring2025/assignments/Data files")

library(dplyr)

# assigning txt files to fang or snp
fang <- fang_et_al_genotypes
snp <- snp_position

dim(fang) # dim 2783, 986. this lists colms and rows
dim(snp) # dim 983, 15. this lists colms and rows

ncol(fang) # number of columns is 986
ncol(snp) # number of columns is 15

nrow(fang) # 2783 rows
nrow(snp) # 983 rows

names(snp) # tells the names of column labels

head(fang, 10) # shows first 10 rows of data
head(snp, 10) # shows first 10 rows of data

tail(fang, 10) # shows last 10 rows of data
tail(snp, 10) # shows last 10 rows of data

# this forloop displays the contents of each variable in the snp file
for(i in 1:ncol(snp)) {
  print(names(snp)[i])
  print(table(snp[, i], useNA = "ifany"))
  cat("\n") # Adds a line before the next iteration
} 

# same thing but for fang file
for(i in 1:ncol(fang)) {
  print(names(fang)[i])
  print(table(fang[, i], useNA = "ifany"))
  cat("\n") # Adds a line before the next iteration
} 

summary(snp) # this shows the length, class, mode, min's/maxes, quantiles of the file

# Data Processing

# read the file 
genotype_data <- read.table(file = "fang_et_al_genotypes.txt", header = TRUE)

# For Maize: filter these specified groups, ZMMIL, ZMMLR, ZMMR and assign to genotypes_maize. Essentially, this is extracting the rows when the group is one of the three indicated groups
genotypes_maize <- filter(genotype_data, Group == "ZMMIL" | Group == "ZMMLR" | Group == "ZMMR")

print(genotypes_maize)

# Same thing but for teosinte data
genotypes_teosinte <- filter(genotype_data, Group == "ZMPBA" | Group == "ZMPIL" | Group == "ZMPJA")

print(genotypes_teosinte)

# Transpose the data for maize and teosinte
maize_transpose <- t(genotypes_maize)
teosinte_transpose <- t(genotypes_teosinte)

print(maize_transpose)
print(teosinte_transpose)

# remove the unwanted columns from the snp data as indicated by the snp[-c(2,5:15)]. This means to remove column 2 and 5 through 15.
SNP_clean <- snp[-c(2,5:15)]
SNP_clean # we can see we have the columns we want: SNP_ID, Chromosome, and Position

# label rownames as SNP_ID so we can then merge SNP_clean and the teosinte and maize files by row.names
rownames(SNP_clean) <- SNP_clean$SNP_ID
rownames(SNP_clean)

# merge maize and snp_clean
joined_maize <- merge(SNP_clean, maize_transpose)

# merge teosinte and snp_clean
joined_teosinte <- merge(SNP_clean, teosinte_transpose)

# create a for-loop that cycles through chromosomes 1-10 in the joined_maize data file in my global env and creates txt files by position in increasing order.
for (i in 1:10) {inc_maize <- subset(joined_maize, joined_maize$Chromosome == i)
inc_maize[is.na(inc_maize)] <- "?"
inc_maize <- inc_maize[order(inc_maize$Position), ]
assign(paste0("/Users/jordyn/Desktop/BCB546_Spring2025/assignments/R-Assignment/chr", i, "_increasing_maize"), inc_maize)
write.table(inc_maize, file = paste0("/Users/jordyn/Desktop/BCB546_Spring2025/assignments/R-Assignment/chr", i, "_increasing_maize.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

# same thing but for teosinte data
for (i in 1:10) {inc_teosinte <- subset(joined_teosinte, joined_teosinte$Chromosome == i)
inc_teosinte[is.na(inc_teosinte)] <- "?"
inc_teosinte <- inc_teosinte[order(inc_teosinte$Position), ]
assign(paste0("/Users/jordyn/Desktop/BCB546_Spring2025/assignments/R-Assignment/chr", i, "_increasing_teosinte"), inc_teosinte)
write.table(inc_teosinte, file = paste0("/Users/jordyn/Desktop/BCB546_Spring2025/assignments/R-Assignment/chr", i, "_increating_teosinte.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

# Maize in decreasing order replacing "?" with "-". Creating a single file for each chromosome and sorting by position in decreasing order.
for (i in 1:10) {
  dec_maize <- subset(joined_maize, joined_maize$Chromosome == i)
  dec_maize[dec_maize == "?/?"] <- "-/-"
  dec_maize <- dec_maize[order(dec_maize$Position, decreasing = TRUE), ]
  assign(paste0("chr", i, "_decreasing_maize"), dec_maize)
  write.table(dec_maize, file = paste0("/Users/jordyn/Desktop/BCB546_Spring2025/assignments/R-Assignment/chr", i, "_decreasing_maize.txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t") 
}

# same thing but for teosinte data
for (i in 1:10) {
  dec_teosinte <- subset(joined_teosinte, joined_teosinte$Chromosome == i)
  dec_teosinte[dec_teosinte == "?/?"] <- "-/-"
  dec_teosinte <- dec_teosinte[order(dec_teosinte$Position, decreasing = TRUE), ]
  assign(paste0("chr", i, "_decreasing_maize"), dec_maize)
  write.table(dec_teosinte, file = paste0("/Users/jordyn/Desktop/BCB546_Spring2025/assignments/R-Assignment/chr", i, "_decreasing_maize.txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

# unknown positions for maize, this is looking at rows where the third column has a "?"
write.table(subset(joined_maize, joined_maize$Position == "?"), 
            file = file.path("/Users/jordyn/Desktop/BCB546_Spring2025/assignments/R-Assignment", "unknown_positions_maize.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# unknown positions for teosinte, this is also looking at rows where the third column has a "?"
write.table(subset(joined_teosinte, joined_teosinte$Position == "?"), 
            file = file.path("/Users/jordyn/Desktop/BCB546_Spring2025/assignments/R-Assignment", "unknown_positions_teosinte.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Data Visualization: Part II

library(ggplot2)
library(dplyr)

# Combine maize and teosinte data with a new "Type" column which is either maize or teosinte
maize_type_data <- joined_maize %>% mutate(Type = "Maize")
teosinte_type_data <- joined_teosinte %>% mutate(Type = "Teosinte")

# Merge teosinte and maize datasets
combined_data <- bind_rows(maize_type_data, teosinte_type_data)

# This creates the plot that shows the distribution of SNPs across chromosomes in both species
ggplot(combined_data, aes(x = as.factor(combined_data$Chromosome), fill = Type)) +
  geom_bar(position = "dodge") +
  labs(title = "SNPs across Chromosomes",
       x = "Chromosome",
       y = "SNP #",
       fill = "Species") +
  theme_minimal() +
  scale_fill_manual(values = c("Maize" = "blue", "Teosinte" = "red"))
# From what I can tell based on the graph it seems as though teosinte and maize have the same number of SNPs in each chromosome

# Heterozygosity and Missing Data

library(dplyr)
library(ggplot2)
library(tidyr)

# I am not entirely sure if this is what the questions was asking but I tried my best to create the graphs that I thought answered it.

maize_type_data2 <- joined_maize %>% mutate(Type = "Maize")
teosinte_type_data2 <- joined_teosinte %>% mutate(Type = "Teosinte")

combined_data2 <- bind_rows(maize_type_data2, teosinte_type_data2)

# this is where i use the information in combined_data2 EXCEPT SNP_ID, Chromsome, Row.names, and Position
# so that I am only left with V1:1966 which are the alleles.
# i use the recommended pivot_longer command to keep the type column the same but re-format the other columns so that I have the samples and the genotypes associated
# next we assign alleles to homozygous, heterozygous and missing.
# then I group by the Type column which indicates maize or teosinte and then the genotypes
# then the proportion is calculated for each genotype and assigned to a new column which is count
snp_counts <- combined_data2 %>%
  select(-Row.names, -SNP_ID, -Chromosome, -Position) %>%
  pivot_longer(cols = -Type, names_to = "Sample", values_to = "Genotype") %>%
  mutate(Genotype = case_when(
    Genotype %in% c("A/A", "G/G", "C/C", "T/T") ~ "Homozygous",
    Genotype %in% c("A/T", "G/T", "C/T", "A/G") ~ "Heterozygous",
    Genotype == "?/?" ~ "Missing"
  )) %>%
  group_by(Type, Genotype) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count/sum(Count))

ggplot(snp_counts, aes(x = Type, y = Proportion, fill = Genotype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Genotypes in Maize and Teosinte",
       x = "Species",
       y = "Proportion",
       fill = "Genotype") +
  theme_minimal() +
  scale_fill_manual(values = c("Homozygous" = "blue", 
                               "Heterozygous" = "red", 
                               "Missing" = "gray"))
# we can see that maize has slightly more heterozygous alleles, much more homozygous alleles, and a few more missing alleles

# Part III: My own visualization

# Calculate the minimum SNP length for each species
min_snp_length <- combined_data %>%
  group_by(Type) %>%
  summarize(min_length = min(Position, na.rm = TRUE))

# Plot the minimum SNP lengths
ggplot(min_snp_length, aes(x = Type, y = min_length, fill = Type)) +
  geom_bar(stat = "identity") +
  labs(title = "Minimum SNP length",
       x = "Species",
       y = "Minimum SNP Length") +
  theme_minimal() +
  scale_fill_manual(values = c("Maize" = "blue", "Teosinte" = "red"))
# As we can see both maize and teosinte have the same minimum SNP length







