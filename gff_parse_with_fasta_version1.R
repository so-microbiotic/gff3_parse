library(tidyverse) ## used for data manipulation 
library(ape) ## used to import a gff file
library(foreach) ##used for parsing fasta file in the end

## read in gff3 file (directly downloaded from GenBank)
gff <- read.gff("sequence.gff3", na.strings = c(".", "?"), GFF3 = TRUE)

##select only necessary columns - source, start, and end
gff2 <- gff %>% select(source, start, end) %>% filter(source == "RefSeq")

## removed first row which was a summary of the entire sequence
gff2 = gff2[-1, ]

##editted start column because the first row isn't needed

start2 = c(gff2$start[2:nrow(gff2)], NA)

##created two new columns to obtain coordinates and intergenic region length

gff3 <- gff2 %>% mutate(intergenic_region = start2[1:nrow(gff2)] - gff2$end[1:nrow(gff2)], intergenic_coordinates = str_c(gff2$end[1:nrow(gff2)], start2[1:nrow(gff2)],sep= "-")) 


## output file with unfiltered information 
write.csv(gff3, "unfiltered_gff_w_intergenic_info.csv")

##filtered the previous gff dataframe with intergenic sequence lengths amenable to sequencing

gff4 <- gff3 %>% filter(gff3$intergenic_region > 250 & gff3$intergenic_region < 500)

##output file for filtered information

write.csv(gff4, "filtered_gff_w_intergenic_info.csv")

##read in fasta file as a string, directly downloaded from GenBank without the header

fasta <- read_file("sequence.fasta")

##trim and clean up fasta file by getting rid of white space and removing new lines

fasta_trim <- str_trim(fasta, side = c("both", "left", "right"))

fasta_clean <- gsub("\n", "", fasta)

#using coordinates from gff file to parse the fasta file

fasta_list <- list(foreach(i = gff4$start, j= gff4$end) %do% substr(fasta_clean, i, j))

##generate output for downstream analysis

dput(fasta_list, "fasta_list.txt")


