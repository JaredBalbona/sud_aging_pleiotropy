
###################
## Run MR-PRESSO ## 
###################

# Import functions from mr_presso_functions.R:
source('/lts/aalab/arpana_data/jared/aging_suds_project/mr_presso/scripts/mr_presso_functions.r')

# Set main directories
sud_directory <- "/scratch/aalab/jared/aging_suds_project/reviewer_response/sumstats/sudTraits/EUR"
aging_directory <- "/scratch/aalab/jared/aging_suds_project/reviewer_response/sumstats/agingTraits/EUR"

# Read in the LDSC results:
ldsc_results <- fread('/lts/aalab/arpana_data/jared/aging_suds_project/ldsc/results/EUR/sud_aging_ldsc_results_betweencategory.txt', h=T)
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

## Run IVW Analyses ## 
    ivw_results  <- run_ivw_mr(sig_ldsc_results = ldsc_results_sig, exposure_column = 'sudPhen', outcome_column = 'agingPhen')
    ivw_results2 <- run_ivw_mr(sig_ldsc_results = ldsc_results_sig, exposure_column = 'agingPhen', outcome_column = 'sudPhen')

    ivw_results_all <- rbind(ivw_results, ivw_results2)

    # Apply the operations based on the calculated pthreshold
    ivw_results_all$P_FDR <- p.adjust(ivw_results_all$P, method = "fdr")

    ivw_results_all$Bonferroni_Sig <- ivw_results_all$P <= (.05 / nrow(ldsc_results_sig))
    ivw_results_all$FDR_Sig <- ivw_results_all$P_FDR <= .05

    # Save Results:
    fwrite(ivw_results_all, '/scratch/aalab/jared/aging_suds_project/reviewer_response/mr_ivw_results_ALL.txt', sep=',')

## Run MR-Egger Analyses ##
    egger_results  <- run_mr_egger(sig_ldsc_results = ldsc_results_sig, exposure_column = 'sudPhen', outcome_column = 'agingPhen')
    egger_results2 <- run_mr_egger(sig_ldsc_results = ldsc_results_sig, exposure_column = 'agingPhen', outcome_column = 'sudPhen')

    egger_results_all <- rbind(egger_results, egger_results2)

    # Apply the operations based on the calculated pthreshold
    egger_results_all$P_FDR <- p.adjust(egger_results_all$P, method = "fdr")

    egger_results_all$Bonferroni_Sig <- egger_results_all$P <= (.05 / nrow(ldsc_results_sig))
    egger_results_all$FDR_Sig <- egger_results_all$P_FDR <= .05

    # Save Results:
    fwrite(egger_results_all, '/scratch/aalab/jared/aging_suds_project/reviewer_response/mr_egger_results_ALL.txt', sep=',')

## Run Median-Weighted MR ##
    medweighted_results  <- run_medweighted_mr(sig_ldsc_results = ldsc_results_sig, exposure_column = 'sudPhen', outcome_column = 'agingPhen')
    medweighted_results2 <- run_medweighted_mr(sig_ldsc_results = ldsc_results_sig, exposure_column = 'agingPhen', outcome_column = 'sudPhen')

    medweight_results_all <- rbind(medweighted_results, medweighted_results2)

    # Apply the operations based on the calculated pthreshold
    medweight_results_all$P_FDR <- p.adjust(medweight_results_all$P, method = "fdr")

    medweight_results_all$Bonferroni_Sig <- medweight_results_all$P <= (.05 / nrow(ldsc_results_sig))
    medweight_results_all$FDR_Sig <- medweight_results_all$P_FDR <= .05

    # Save Results:
    fwrite(medweight_results_all, '/scratch/aalab/jared/aging_suds_project/reviewer_response/mr_medweighted_results_ALL.txt', sep=',')
