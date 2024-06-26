#analyze the SBS mutational signatures in OVCAR4

#load required packages
library('MutationalPatterns')
library('BSgenome')
library('dplyr')
library('stringr')
#load reference genome
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg38'
library(ref_genome, character.only = TRUE)
#load vcf files 
vcfs <- list.files(path = paste(input_dir, 'variants_calling/OVCAR4_WES', sep = '/'),
                   full.names = TRUE, pattern = '.vcf.gz')

sample <- vcfs %>%
  str_extract('(?<=OVCAR4_WES/).+') %>%
  str_remove('_pass.vcf.gz')
#load the VCF files into a GRangesList
grl <- read_vcfs_as_granges(vcfs, sample, ref_genome)

#number of different types of base substitution 
type_occurrences <- mut_type_occurrences(grl, ref_genome)
write.csv(type_occurrences, paste(output_dir, 'OVCAR4_number_base_subs.csv', sep = '/'))
#trinucleotide context
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
write.csv(mut_mat, paste(output_dir, 'OVCAR4_trinucleotide_matrix.csv', sep = '/'))

#fit to signature COSMIC v3.2
cosmic_3.2 <- get_known_signatures(muttype = 'snv', 'COSMIC_v3.2')
fit <- fit_to_signatures(mut_mat, cosmic_3.2)
write.csv(fit$contribution, paste(output_dir, 'OVCAR4_sbs_contribution.csv', sep = '/'))