
###################
## Run MR-PRESSO ## 
###################

## Preprocess Data ##

# Import functions from mr_presso_functions.R:
source('/scratch/aalab/jared/aging_suds_project/mr_presso/scripts/mr_presso_functions.r')

# Read in the LDSC results:
ldsc_results <- fread('/scratch/aalab/jared/aging_suds_project/ldsc/results/EUR/sud_aging_ldsc_results_betweencategory.txt', h=T)
colnames(ldsc_results)[1:2] <- c('sudPhen', 'agingPhen')

# Separate the FDR significant columns:
ldsc_results$FDR <- p.adjust(ldsc_results$p, method = "fdr")
ldsc_results_sig <- data.frame(ldsc_results[ldsc_results$FDR <= .05,])

# Read in 1KG EAF file (this will be used later): 
eaf_file <- fread("/scratch/aalab/jared/aging_suds_project/mr_presso/1KG_Reference_Files/g1000_eur_allele_freqs.afreq", h=T)
eaf_file <- eaf_file[,c('ID', 'REF', 'ALT', 'ALT_FREQS')]
colnames(eaf_file) <- c('SNP', 'A1', 'A2', 'MAF')

# Set main directories
exposure_input_directory <- "/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/raw_sumstats/EUR"
outcome_input_directory <- "/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/agingTraits/raw_sumstats/EUR"

# Set up output directories:
exposure_output_directory <- "/scratch/aalab/jared/aging_suds_project/mr_presso/sumstats/sudTraits/EUR/"
outcome_output_directory <- "/scratch/aalab/jared/aging_suds_project/mr_presso/sumstats/agingTraits/EUR/"

# Get a list of all .txt files in subdirectories (make sure these are all sumstats files):
exposure_file_paths <- list.files(path = exposure_input_directory, recursive = TRUE, full.names = TRUE)
outcome_file_paths <- list.files(path = outcome_input_directory, recursive = TRUE, full.names = TRUE)

# Make sure not to accidentally read in any README files:
exposure_file_paths <- exposure_file_paths[!grepl('README', exposure_file_paths)]
outcome_file_paths <- outcome_file_paths[!grepl('README', outcome_file_paths)]

# Use lapply to apply the function to each file in the list we created above:
exposure_vars_raw <- lapply(exposure_file_paths, read_and_assign)
outcome_vars_raw  <- lapply(outcome_file_paths, read_and_assign)

# Use mapply to iterate over both exposure_vars_raw and exposure_output_directory
modified_exposure_datasets <- mapply(process_variable, exposure_vars_raw, exposure_output_directory, SIMPLIFY = FALSE)
modified_outcome_datasets <- mapply(process_variable, outcome_vars_raw, outcome_output_directory, SIMPLIFY = FALSE)

## Run MR-PRESSO ##

# Import functions from mr_presso_functions.R:
source('/scratch/aalab/jared/aging_suds_project/mr_presso/scripts/mr_presso_functions.r')

# Set main directories
sud_directory <- "/scratch/aalab/jared/aging_suds_project/mr_presso/sumstats/sudTraits/EUR"
aging_directory <- "/scratch/aalab/jared/aging_suds_project/mr_presso/sumstats/agingTraits/EUR"

# Read in the LDSC results:
ldsc_results <- fread('/scratch/aalab/jared/aging_suds_project/ldsc/results/EUR/sud_aging_ldsc_results_betweencategory.txt', h=T)
colnames(ldsc_results)[1:2] <- c('sudPhen', 'agingPhen')

# Separate the FDR significant columns:
ldsc_results$FDR <- p.adjust(ldsc_results$p, method = "fdr")
ldsc_results_sig <- data.frame(ldsc_results[ldsc_results$FDR <= .05,])

# Get a list of all .txt files in subdirectories (make sure these are all sumstats files):
sud_exposure_paths <- list.files(path = sud_directory, pattern = "\\_mrpresso_exposure.txt$", recursive = TRUE, full.names = TRUE)
sud_outcome_paths <- list.files(path = sud_directory, pattern = "\\_mrpresso_outcome.txt$", recursive = TRUE, full.names = TRUE)
  
aging_exposure_paths <- list.files(path = aging_directory, pattern = "\\_mrpresso_exposure.txt$", recursive = TRUE, full.names = TRUE)
aging_outcome_paths <- list.files(path = aging_directory, pattern = "\\_mrpresso_outcome.txt$", recursive = TRUE, full.names = TRUE)

# If you aren't using *too* many files, this is preferable to reading them in each time you run MR:
sud_exposures <- lapply(sud_exposure_paths, read_and_assign)
aging_exposures <- lapply(aging_exposure_paths, read_and_assign)

sud_outcomes <- lapply(sud_outcome_paths, read_and_assign)
aging_outcomes <- lapply(aging_outcome_paths, read_and_assign)

results <- run_mr_presso_analyses(ldsc_results_sig, exposure_column = 'sudPhen', outcome_column = 'agingPhen')
#fwrite(results, '/scratch/aalab/jared/aging_suds_project/mr_presso/results/EUR/mrpresso_SudOnAging_results.txt', sep=',')

# Run it again with the directionality switched:
results2 <- run_mr_presso_analyses(ldsc_results_sig, exposure_column = 'agingPhen', outcome_column='sudPhen')
#fwrite(results2, '/scratch/aalab/jared/aging_suds_project/mr_presso/results/EUR/mrpresso_agingOnSud_results.txt', sep=',')

results_all <- rbind(results, results2)

# Apply the operations based on the calculated pthreshold
results_all$Raw_Bonferroni_Sig <- results_all$Raw_Pvalue <= (.05 / nrow(ldsc_results_sig))
results_all$Outlier_Bonferroni_Sig <- results_all$Outlier_Pvalue <= (.05 / nrow(ldsc_results_sig))
results_all$Global_Bonferroni_Sig <- results_all$GlobalTest_Pvalue <= (.05 / nrow(ldsc_results_sig))
results_all$Distortion_Bonferroni_Sig <- results_all$DistortionTest_Pvalue <= (.05 / nrow(ldsc_results_sig))

results_all$Raw_pvalue_FDR <- p.adjust(results_all$Raw_Pvalue, method = "fdr")
results_all$Outlier_pvalue_FDR <- p.adjust(results_all$Outlier_Pvalue, method = "fdr")
results_all$Global_pvalue_FDR <- p.adjust(results_all$GlobalTest_Pvalue, method = "fdr")
results_all$Distortion_pvalue_FDR <- p.adjust(results_all$DistortionTest_Pvalue, method = "fdr")

results_all$Raw_FDR_Sig <- results_all$Raw_pvalue_FDR <= .05
results_all$Outlier_FDR_Sig <- results_all$Outlier_pvalue_FDR <= .05
results_all$Global_FDR_Sig <- results_all$Global_pvalue_FDR <= .05
results_all$Distortion_FDR_Sig <- results_all$Distortion_pvalue_FDR <= .05

#fwrite(results_all, '/scratch/aalab/jared/aging_suds_project/mr_presso/results/EUR/mrpresso_sudAging_results_ALL.txt', sep=',')


sudPhen <- rep("Height", 9)
agingPhen <- c("Alzheimers", "Brain_Age_Gap", "Death_Any_Cause", "Grim_Age", "Shared_Aging_Factor", "Telomere_Length", "Parental_Lifespan", "Healthspan", "Frailty_Index")
ldsc_results_sig <- data.frame(sudPhen, agingPhen)
colnames(ldsc_results_sig) <- c("sudPhen", "agingPhen")
ldsc_results_sig

Height <- Height[,c("ref", "alt", "rsid", "AF", "n_complete_samples", "beta", "se", "pval CHRPOS")]
colnames(Height) <- c('A1', 'A2', 'SNP', 'EAF', 'N', 'BETA', 'SE', 'P CHRPOS')
Height <- separate(Height, 'P CHRPOS', into = c("P", "CHRPOS"), sep = " ")
Height$CHRPOS <- NULL

Height2 <- subset(Height, grepl("^rs", SNP))
fwrite(Height2, '/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/raw_sumstats/EUR/Height.txt', sep='\t')

df_reduced <- na.omit(df_reduced)
df_reduced <- df_reduced[!apply(df_reduced_complete == "", 1, all), ]
