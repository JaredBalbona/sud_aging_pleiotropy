#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G

## Before munging the sumstats ##
    # First, copy each SS file into my project directory so that I can use it:
    # Source directory
    source_dir="/lts/aalab/arpana_data/SummaryStats/"

    # Destination directories
    sud_destination_dir="/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/raw_sumstats/EUR/"
    aging_destination_dir="/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/agingTraits/raw_sumstats/EUR/"

    # List of source files:
    sud_files=(
        # TUD: 
        "Toikumo_2023_TUD/published_version/TUD_meta_eur_maf0001_withUKBB1_chr_pos.txt"
        # DPW:
        "GSCAN_2022_sumstats/GSCAN_DrnkWk_2022_GWAS_SUMMARY_STATS_EUR.txt.gz"
        # SI: 
        "GSCAN_2022_sumstats/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt.gz"
        # PAU: 
        "Zhou_2023_PAU/PAU_Zhou_NatMed/PAU_EUR_2023.txt"
        # CanEU: 
        "Pasman_sumstats/Cannabis_ICC_UKB.gz"
        # CUD: 
        "Johnsonetal_CUD_Sumstats/CUD_EUR_casecontrol_wFreq1.gz"
        # mvSUD:
        "Hatoum_2023_SUDs/Hatoum2023AddictionEuropean.txt.gz"
    )

    sud_names=(
        # New names corresponding to each file in sud_files (Make sure it's in the same order as the list above)
        "Tobacco_Use_Disorder.txt"
        "Drinks_Per_Week.txt.gz"
        "Smoking_Initiation.txt.gz"
        "Problematic_Alcohol_Use.txt"
        "Cannabis_Ever_Use.gz"
        "Cannabis_Use_Disorder.gz"
        "Shared_Substance_Use_Factor.txt.gz"
    )

    aging_files=(
        # Alzheimer's: 
        "Bellenguez_2022_AlZ/GCST90027158_buildGRCh38.tsv.gz"
        # Brain Age Gap: 
        "Leonardsen_2023_BrainAgeGap/meta_allfolds_BAG_noIndel.txt.gz"
        # Death Any Cause:
        "FinnGen_2022_Death/summary_stats_finngen_R10_DEATH.gz"
        # Telomere Length:
        "Codd_2021_TelomereLength/UKB_telomere_gwas_summarystats.tsv.gz"
        # Frailty:
        "Rosoff_2023_Aging/Atkins_2021_FrailtyIndex/Atkins_2021_FrailtyIndex_GWAS_EuropeansAged60to70years_GRCh37_MetaAnalysis.txt.gz"
        # Grim Age:
        "Rosoff_2023_Aging/McCartney_2020_GrimAge/GrimAge_EUR_summary_statistics.txt.gz"
        # Healthspan:
        "Rosoff_2023_Aging/Zenin_2019_Healthspan/healthspan_summary.csv.gz"
        # Longevity:
        "Rosoff_2023_Aging/Deelen_2019_Longevity/Results_90th_percentile.txt.gz"
        # Parental Lifespan:
        "Rosoff_2023_Aging/Timmers_2019_ParentalLifespan/lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz"
        # mvAge:
        "Rosoff_2023_Aging/Rosoff_2023_mvAge/mvAge.summary.EUR.txt.gz"
    )

    aging_names=(
        # New names corresponding to each file in aging_files (Make sure it's in the same order as the list above)
        "Alzheimers.tsv.gz"
        "Brain_Age_Gap.txt.gz"
        "Death_Any_Cause.gz"
        "Telomere_Length.tsv.gz"
        "Frailty_Index.txt.gz"
        "Grim_Age.txt.gz"
        "Healthspan.csv.gz"
        "Longevity.txt.gz"
        "Parental_Lifespan.tsv.gz"
        "Shared_Aging_Factor.txt.gz"
    )

    # Loop through each file in SUDs and copy it to the destination directory with the new name
    for ((i=0; i<${#sud_files[@]}; i++)); do
        scp "$source_dir${sud_files[$i]}" "$sud_destination_dir${sud_names[$i]}" &
    done

    # Loop through each file in Aging and copy it to the destination directory with the new name
    for ((i=0; i<${#aging_files[@]}; i++)); do
        scp "$source_dir${aging_files[$i]}" "$aging_destination_dir${aging_names[$i]}" &
    done

    # Have to do Kember OUD separately bc it's being borrowed from the Elliot's lab: 
    scp /lts/enlab/elliot_data_3/overdoseFilesForAlexanderEmma/Combined/OUDForMunge.txt.gz "$sud_destination_dir""Opioid_Use_Disorder.txt.gz"

    # Next, manually change column headers so that they don't cause errors when munging:
    # (I tried to automate this but there were too many idiosyncrasies between each file) 
    gunzip /scratch/aalab/jared/aging_suds_project/ldsc/sumstats/agingTraits/raw_sumstats/EUR/*.gz
    gunzip /scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/raw_sumstats/EUR/*.gz    

# PRE-MUNGE THE SUMSTATS BEFORE MUNGING # 

    # In R: 
        # Healthspan needs to be tab separated, not comma separated, and the p-val column needs to be converted from -log10:
            # library(data.table)
            # healthspan <- fread("/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/agingTraits/raw_sumstats/EUR/Healthspan.csv", h=T)
            # colnames(healthspan)[10] <- "P"
            # healthspan$P <- 10^(-healthspan$P)
            # fwrite(healthspan, "/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/agingTraits/raw_sumstats/EUR/Healthspan.txt", sep='\t')

        # Longevity needs rsids instead of CHR:BP
            # longevity <- fread("/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/agingTraits/raw_sumstats/EUR/Longevity.txt", h=T)

        # TUD needs betas instead of Z scores:
            # I validated this approach using the healthspan sumstats which have both Z and Beta, and the estimated Beta and actual beta were correlated at .99
            # TUD <- fread("/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/raw_sumstats/EUR/Tobacco_Use_Disorder.txt", h=T)
            
            # Sample_Size = 739895

            # TUD$Beta <- TUD$ZSCORE  / sqrt((2 * TUD$EAF) * (1 - TUD$EAF) * (Sample_Size + TUD$ZSCORE^2))
            # TUD$SE   <- 1           / sqrt((2 * TUD$EAF) * (1 - TUD$EAF) * (Sample_Size + TUD$ZSCORE^2))

            # zero_rows <-        which(sqrt((2 * TUD$EAF) * (1 - TUD$EAF) * (Sample_Size + TUD$ZSCORE^2)) == 0)
            # TUD$Beta[zero_rows] <- TUD$SE[zero_rows] <- 1e-100
            # colnames(TUD)[9] <- "ZSc"

            # fwrite(TUD, "/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/raw_sumstats/EUR/Tobacco_Use_Disorder.txt", sep='\t')


    # In Bash: 
        # Remove the .csv file we just changed:
            # rm /scratch/aalab/jared/aging_suds_project/ldsc/sumstats/agingTraits/raw_sumstats/EUR/Healthspan.csv

        # Tobacco_Use_Disorder has some SNPs in CHR:POS format, so those need to be removed:    
            # # Save the header row in a separate file
            # head -n 1 /scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/raw_sumstats/EUR/Tobacco_Use_Disorder.txt > /scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/raw_sumstats/EUR/temp.txt

            # # Filter rows starting with "rs" and append them to the header
            # grep "^rs" /scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/raw_sumstats/EUR/Tobacco_Use_Disorder.txt >> /scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/raw_sumstats/EUR/temp.txt

            # # Replace yourfile.txt with the filtered content
            # mv /scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/raw_sumstats/EUR/temp.txt /scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/raw_sumstats/EUR/Tobacco_Use_Disorder.txt


    # Fix the column headers so that they work with LDSC: 
    modify_sumstats() {
        local directory="$1"
        local old_snp_names=("VARIANT_ID" "RSIDS" "RSID" "SNPID" "SNP_ID", 'MARKERNAME')
        local old_a2_names=("REF" "ALLELE2" "RA" "NEA" "A0" "OTHER_ALLELE")
        local old_a1_names=("ALT" "ALLELE1" "EA" "EFFEC_ALLELE" "EFFECT_ALLELE")
        local old_n_names=("EFFECTIVE_N")
        local old_p_names=("P_VALUE" "P-VALUE" "PVALUE" "P.VALUE" "PVAL")
        local old_eaf_names=("FREQ1" "A1F" "A1_FREQUENCY" "AF_A1" "AF_1000G" "effect_allele_frequency")
        local old_beta_names=("BETA1" "EFFECT")

        for file in "$directory"/*; do
            if [ -f "$file" ]; then
                local first_line=$(head -n 1 "$file" | tr '[:lower:]' '[:upper:]')  # Capitalize the first line
                local modified_line="$first_line"

                for string in "${old_snp_names[@]}" "${old_a2_names[@]}" "${old_a1_names[@]}" "${old_n_names[@]}" "${old_p_names[@]}" "${old_eaf_names[@]}" "${old_beta_names[@]}"; do
                    if [[ $first_line == *"$string"* ]]; then
                        case "$string" in
                            VARIANT_ID|RSIDS|RSID|SNPID|SNP_ID) modified_line=${modified_line//$string/SNP} ;;
                            REF|ALLELE2|RA|NEA|A0|OTHER_ALLELE) modified_line=${modified_line//$string/A2} ;;
                            ALT|ALLELE1|EA|EFFEC_ALLELE|EFFECT_ALLELE) modified_line=${modified_line//$string/A1} ;;
                            EFFECTIVE_N) modified_line=${modified_line//$string/N} ;;
                            P_VALUE|P-VALUE|PVALUE|P.VALUE|PVAL) modified_line=${modified_line//$string/P} ;;
                            FREQ1|A1F|A1_FREQUENCY|AF_A1|AF_1000G|effect_allele_frequency) modified_line=${modified_line//$string/EAF} ;;
                            BETA1|EFFECT) modified_line=${modified_line//$string/Beta} ;;
                        esac
                    fi
                done

                local temp_file=$(mktemp)
                echo "$modified_line" > "$temp_file"
                tail -n +2 "$file" >> "$temp_file"
                mv "$temp_file" "$file"
            fi
        done
    }

    # Apply the above function: 
    modify_sumstats '/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/agingTraits/raw_sumstats/EUR/'
    modify_sumstats '/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/raw_sumstats/EUR/'

    # Manually addressing errors from the munge_sumstats output: 
        # Healthspan: Too many signed sumstat columns
            # Change Z to "ZSc" so it is ignored
            # Note, I could also use the "--ignore" flag, but eh it's just one column
        # Parental Lifespan: 2 columns named SNP
            # Change the CHR:POS SNP column to "CHR_POS"

    # gzip /scratch/aalab/jared/aging_suds_project/ldsc/sumstats/agingTraits/raw_sumstats/EUR/*
    # gzip /scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/raw_sumstats/EUR/*
