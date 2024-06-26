#analyze the WUSTL patients cohort
#'met' stands for metastasis samples
#'primary' stands for matching primary samples

#load required packages
library('MutationalPatterns') 
library('BSgenome')
library('dplyr')
library('stringr')
#load reference genome - GRCh38
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg38' 
library(ref_genome, character.only = TRUE)

#load vcf files
met_vcfs <- list.files(path = paste(input_dir, 'variants_calling/WUSTL', sep = '/'),
                       full.names = TRUE, pattern = '_M_')
primary_vcfs <- list.files(path = paste(input_dir, 'variants_calling/WUSTL', sep = '/'),
                           full.names = TRUE, pattern = '_P_')
all_vcfs <- c(primary_vcfs, met_vcfs)
#use the sample identifiers as sample names 
sample <- all_vcfs %>%
  str_extract('(?<=WUSTL/).+') %>%
  str_remove('\\.vcf.gz')

#load the VCF files into a GRangesList
grl <- read_vcfs_as_granges(all_vcfs, sample, ref_genome)

#number of different types of base substitution 
type_occurrences <- mut_type_occurrences(grl, ref_genome)
write.csv(type_occurrences, paste(output_dir, 'WUSTL_number_base_subs.csv', sep = '/'))
#trinucleotide context
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
write.csv(mut_mat, paste(output_dir, 'WUSTL_trinucleotide_matrix.csv', sep = '/'))

#fit to signature COSMIC v3.2
cosmic_3.2 <- get_known_signatures(muttype = 'snv', 'COSMIC_v3.2')
fit <- fit_to_signatures(mut_mat, cosmic_3.2)
write.csv(fit$contribution, paste(output_dir, 'WUSTL_sbs_contribution.csv', sep = '/'))