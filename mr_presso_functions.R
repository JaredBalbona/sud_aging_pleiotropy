library(TwoSampleMR)
library(MRPRESSO)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)
library(readr)
library(cause, lib.loc = "/ref/aalab/software/r-envs/aalab_R/4.1.3")
library(ieugwasr, lib.loc = "/ref/aalab/software/r-envs/aalab_R/4.1.3")


## Format each of these datasets ## 

# Define the original and new column lists
orig_column_list <- list(
  SNP_columns = unique(c("RSID", "MARKERNAME", "SNP", "SNP_ID", "SNPID", "VARIANT_ID")),
  effect_allele_columns = unique(c("ALLELE1", "EFFECT_ALLELE", "A1", "EA", "REF", "ALLELE 1")),
  other_allele_columns = unique(c("ALLELE2", "OTHER_ALLELE", "A2", "NEA", "ALT", "ALLELE 2")),
  beta_columns = unique(c("EFFECT", "BETA", "Beta")),
  pval_columns = unique(c("P", "P-VALUE", "P.VALUE", "PVALUE", "P_VALUE", "P VALUE")),
  eaf_columns = unique(c('EAF', 'FREQ1', 'EFFECT ALLELE FREQUENCY (EAF)')),
  se_columns = unique(c('STDERR', 'SE', 'STANDARD_ERROR')), 
  samplesize_columns = unique(c('N', 'SAMPLESIZE', 'EFFECTIVE_N', 'N_COMPLETE_SAMPLES', 'SAMPLE SIZE')),
  maf = 'MAF' # To prevent duplicate columns names when we merge
)

new_column_list <- list(
  SNP = "SNP",
  effect_allele = "effect_allele",
  other_allele = "other_allele",
  beta = "beta",
  pval = "pval",
  eaf = 'eaf',
  se = 'se', 
  samplesize = 'samplesize',
  maf = 'maf' # To prevent duplicate columns names when we merge
)

# Function to read a file and assign its content to a variable:
read_and_assign <- function(file_path) {
  # Extract the file name without extension to use as the variable name
  variable_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Provide update:
  cat("\nLoading", variable_name, "data:\n")

  # Read the file content into the variable
  assign(variable_name, fread(file_path, h=T), envir = .GlobalEnv)
  
  # Return the variable name for reference
  return(variable_name)
}

# Function to rename columns in a data frame based on the specified column lists
rename_columns <- function(df, orig_columns, new_columns) {
  for (i in seq_along(orig_columns)) {
    orig_col_names <- orig_columns[[i]]
    new_col_name <- new_columns[[i]]
    
    # Check if any of the specified columns exist in the data frame
    columns_to_rename_exist <- intersect(colnames(df), orig_col_names)
    
    # Rename the existing columns to the corresponding new name
    if (length(columns_to_rename_exist) > 0) {
      colnames(df)[colnames(df) %in% columns_to_rename_exist] <- new_col_name
    }
  }
  
  return(df)
}

# Function to process the data to be in a usable format for MR-PRESSO
process_variable <- function(variable_name, output_directory, threshold = 5e-8) {

    cat("\n\n\n-------------------------------------------------\n")
    cat('Processing', variable_name, 'data:\n')
    cat("-------------------------------------------------\n")

  # Get the data frame associated with the variable
  df <- get(variable_name)
  
  # Rename columns in the data frame
  df <- rename_columns(df, orig_column_list, new_column_list)

  # Calculate SE if necessary:
  if (!"se" %in% colnames(df)) {
    df$z <- qnorm(1 - df$pval / 2)
    df$se <- abs(df$beta / df$z)
  }

  # Capitalize the allele columns to be consistent with the eaf file:
  df <- mutate(df, effect_allele = toupper(effect_allele))
  df <- mutate(df, other_allele = toupper(other_allele))

  # Add eaf column if necessary:
  if (!"eaf" %in% colnames(df)) {
    cat('Merging with MAF data...\n')
    df <- df %>% 
      inner_join(eaf_file, by = "SNP") %>% 
      mutate(eaf = case_when( 
        A1 == effect_allele ~ MAF,
        A1 == other_allele ~ 1 - MAF,
        TRUE ~ NA_real_))  %>%  
      filter(!is.na(eaf))
  }    

  # Add sample size column if necessary:
  if (!"samplesize" %in% colnames(df)) {
    df$samplesize <- ""
  }

  # Return the modified data frame
  df$Phenotype = variable_name
  df <- df[, c('Phenotype', 'SNP', 'beta', 'se', 'eaf', 'effect_allele', 'other_allele', 'pval', 'samplesize')]

  # Set it up as the exposure:
  cat('Formatting for MR...\n')
  cat("Creating Exposure Files... \n")
  df_reduced <- suppressWarnings(format_data(df, type = "exposure"))
  df_reduced$exposure = variable_name

  df_reduced <- df_reduced[which(df_reduced$pval.exposure < threshold), ]
  df_reduced <- clump_data(df_reduced)

  # Set it up as the outcome:
  cat("Creating Outcome Files... \n")
  df <- suppressWarnings(format_data(df, type = "outcome"))
  df$outcome = variable_name

  # Save the modified data frame to a new file using fwrite 
  cat('Saving Results...\n')
  fwrite(df, file = paste0(output_directory, variable_name, "_mrpresso_outcome.txt"), row.names = FALSE)
  fwrite(df_reduced, file = paste0(output_directory, variable_name, "_mrpresso_exposure.txt"), row.names = FALSE)

  cat('Done.\n')
  cat('-------------------------------------\n\n')
  return(list(df = df))
}

run_mr_presso_analyses <- function(sig_ldsc_results, exposure_column, outcome_column) {
    
    results_colnames <- c(
            'Exposure_Name',
            'Outcome_Name',
            'Raw_CausalEstimate',
            'Raw_SD',
            'Raw_Tstat',
            'Raw_Pvalue',
            'Outliers_Detected',
            'Outlier_CausalEstimate',
            'Outlier_SD',
            'Outlier_Tstat',
            'Outlier_Pvalue',
            'GlobalTest_RSSobs',
            'GlobalTest_Pvalue',
            'Distortion_Coefficient',
            'DistortionTest_Pvalue',
            'Error')
    
    results <- data.frame(matrix(NA, nrow = nrow(sig_ldsc_results), ncol = length(results_colnames)))
    colnames(results) <- results_colnames

    # Run the Loop:
    for (i in 1:nrow(sig_ldsc_results)) {
        
        pair_temp <- sig_ldsc_results[i, c(exposure_column, outcome_column)]

        results[i,'Exposure_Name'] <- exposure_name <- as.character(pair_temp[1])
        results[i,'Outcome_Name']  <- outcome_name  <- as.character(pair_temp[2])

        # Print update:
        cat("\n\n--------------------------------------------------------------------")
        cat("\nConducting MR-PRESSO for:", exposure_name, "--➜", outcome_name, "\n")
        cat("(Analysis", i, "of", nrow(sig_ldsc_results),")\n")
        cat("--------------------------------------------------------------------\n")

        # Get data ready for analysis:
        exposure <- data.frame(get(paste0(exposure_name, '_mrpresso_exposure')))
        outcome  <- data.frame(get(paste0(outcome_name,  '_mrpresso_outcome')))

        outcome <- outcome[outcome$SNP %in% exposure$SNP,]

        if(nrow(outcome) > 0){
        harmonized_mr_data <- harmonise_data(
            exposure_dat = exposure,
            outcome_dat = outcome
        )
        }
        
        cat("\nRunning analysis...\n")

        tryCatch({
            # Run MR-PRESSO:
            raw_results <- run_mr_presso(harmonized_mr_data, NbDistribution = 1000, SignifThreshold = 0.05)

            cat("Collecting output...\n")

            # Check if there are outliers:
            Outliers_Detected <- !is.na(raw_results[[1]]$'Main MR results'$'Causal Estimate'[2])
            
            # Collect Relevant Results:
            results_temp <- data.frame(
                exposure_name,
                outcome_name,
                raw_results[[1]]$'Main MR results'[1,3:6],
                Outliers_Detected,
                raw_results[[1]]$'Main MR results'[2,3:6],
                raw_results[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs,
                raw_results[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
            )
            # If the distortion tests did not run, supply NAs for those columns:
            if(is.null(raw_results[[1]]$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`)){
                results_temp <- cbind(results_temp, NA, NA)
            }else{
                results_temp <- cbind(results_temp,
                    raw_results[[1]]$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`,
                    raw_results[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue)
            }

            results_temp$Error <- 'None'
            colnames(results_temp) <- results_colnames
            print(results_temp)
            results[i,] <- results_temp
            cat("...Done.")
            
            }, error = function(e) {
                # If there were not enough SNPs, write that into the results table:
                cat("Error occurred:", conditionMessage(e), "\n")
                cat("Skipping to the next iteration.\n\n")
        })
    }
    results[is.na(results$Error),]$Error <- 'Not enough instrumental variables'
    return(results)
}

run_cause_analyses  <- function(sig_ldsc_results, output_directory, r2_thresh = 0.01, pval_thresh = 1e-3) {
    results1 <- data.frame()
    results2 <- data.frame()

    # Run the Loop:
    for (i in 1:nrow(sig_ldsc_results)) {
        pair_temp <- sig_ldsc_results[i,1:2]

        Trait1_name <- as.character(pair_temp[1])
        Trait2_name <- as.character(pair_temp[2])
        
        cat("\n\n\n-------------------------------------------------\n")
        cat('Performing CAUSE for', Trait1_name, 'and', Trait2_name, "\n")
        cat("-------------------------------------------------\n")

        Trait1 <- data.frame(get(Trait1_name))
        Trait2 <- data.frame(get(Trait2_name))

        cat('\nReformatting columns...\n')
        Trait1 <- rename_columns(Trait1, orig_column_list, new_column_list)
        Trait2 <- rename_columns(Trait2, orig_column_list, new_column_list)
    
        # Calculate SE if necessary:
        if (!"se" %in% colnames(Trait1)) {
            Trait1$z <- qnorm(1 - Trait1$pval / 2)
            Trait1$se <- abs(Trait1$beta / Trait1$z)
        }

        if (!"se" %in% colnames(Trait2)) {
            Trait2$z <- qnorm(1 - Trait2$pval / 2)
            Trait2$se <- abs(Trait2$beta / Trait2$z)
        }
        
        cat('\n\nMerging Traits...\n')
        traits_merged <- gwas_merge(Trait1, Trait2,
                        snp_name_cols = c("SNP", "SNP")%>% rev(), 
                        beta_hat_cols = c("beta", "beta")%>% rev(),  
                        se_cols = c("se", "se")%>% rev(),  
                        A1_cols = c("effect_allele", "effect_allele")%>% rev(),  
                        A2_cols = c("other_allele", "other_allele")%>% rev())

        plink_bin = '/ref/aalab/software/plink1.9/plink'
        bfile = '/scratch/aalab/jared/aging_suds_project/cause/reference_files/g1000_eur'

        tryCatch({
            cat('\nClumping SNPs...\n')
            traits_clumped = traits_merged %>%
            rename(rsid = snp,
                    pval = p1) %>%
            ieugwasr::ld_clump(dat = .,
                                clump_r2 = r2_thresh,
                                clump_p = pval_thresh,
                                plink_bin = plink_bin,
                                bfile = bfile)

            top_vars = traits_clumped$rsid

            cat('\n\nRunning CAUSE (this may take a minute)...\n')
            
            varlist = with(traits_merged, sample(snp, size=1000000, replace=F))
            params = est_cause_params(traits_merged, varlist) # This is the step that takes the longest
            
            # cause_results = cause(traits_merged, variants = top_vars, param_ests = params)
            cause_results = cause(traits_merged, variants = top_vars, param_ests = params)
            
            # Gather the relevant Results:
                cat('\n\nGathering Results...\n')
                results1_temp <- cause_results$elpd
                results1_temp$p_value <- 2 * (1 - pnorm(abs(results1_temp$z)))
                results1_temp$Traits <- paste0(Trait1_name, '-', Trait2_name)

                results1 <- rbind(results1, results1_temp)

                results2_temp <- data.frame(summary(cause_results, ci_size=0.95)$tab)
                results2_temp$Traits <- paste0(Trait1_name, '-', Trait2_name)
                results2 <- rbind(results2, results2_temp)

                fwrite(results1_temp, paste0(output_directory, 'cause_ELPD_Results_',Trait1_name, '_', Trait2_name,'.txt'))
                fwrite(results2_temp, paste0(output_directory, 'cause_medianInterval_Results_',Trait1_name, '_', Trait2_name,'.txt'))

                cat('\n\nSaving Output..\n')
                save(cause_results, file = paste0(output_directory, 'rData/', Trait1_name, '_', Trait2_name, '.RData'))

                cat('\n...done.\n\n')
            
            }, error = function(e) {
                # If there were not enough SNPs, write that into the results table:
                cat("Error occurred:", conditionMessage(e), "\n")
                cat("Skipping to the next iteration.\n\n")
        })
    }
    fwrite(results1, paste0(output_directory, 'ELPD/EUR/cause_ELPD_Results_ALL.txt'))
    fwrite(results2, paste0(output_directory, 'medianInterval/EUR/cause_medianInterval_Results_ALL.txt'))
}
