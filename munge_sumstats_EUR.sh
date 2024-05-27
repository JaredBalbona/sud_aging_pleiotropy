#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G


## Munge Sumtats in Parallel ##

eval $(spack load --sh python@3.9.12 )

run_munge_sumstats() {
    input_file="$1"
    output_prefix="$2"
    merge_alleles="/scratch/aalab/jared/aging_suds_project/ldsc/reference_panels/eur_w_ld_chr/w_hm3.snplist"
    chunksize="500000"
    N="${3:-}"  # Third argument is optional
    N_cas="${4:-}"  # Fourth argument is optional
    N_con="${5:-}"  # Fifth argument is optional

    /ref/aalab/software/ldsc/env/bin/munge_sumstats.py \
    --sumstats "$input_file" \
    --merge-alleles "$merge_alleles" \
    ${N:+--N "$N"} \
    ${N_cas:+--N-cas "$N_cas"} \
    ${N_con:+--N-con "$N_con"} \
    --chunksize "$chunksize" \
    --out "$output_prefix"
}

# Input directories
    sud_input_dir="/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/raw_sumstats/EUR/"
    aging_input_dir="/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/agingTraits/raw_sumstats/EUR/"
    
# Output directories
    sud_output_dir="/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/munged_sumstats/EUR/"
    aging_output_dir="/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/agingTraits/munged_sumstats/EUR/"


# Substance Use Traits:
    # TUD: 
        run_munge_sumstats $sud_input_dir"Tobacco_Use_Disorder.txt" $sud_output_dir"Tobacco_Use_Disorder_Munged" 739895 &
    # DPW:
        run_munge_sumstats $sud_input_dir"Drinks_Per_Week.txt"  $sud_output_dir"Drinks_Per_Week_Munged" &
    # SI: 
        run_munge_sumstats $sud_input_dir"Smoking_Initiation.txt" $sud_output_dir"Smoking_Initiation_Munged" &
    # PAU: 
        run_munge_sumstats $sud_input_dir"Problematic_Alcohol_Use.txt" $sud_output_dir"Problematic_Alcohol_Use_Munged" 903147 &
    # CanEU: 
         run_munge_sumstats $sud_input_dir"Cannabis_Ever_Use" $sud_output_dir"Cannabis_Ever_Use_Munged" &
    # CUD: 
        run_munge_sumstats $sud_input_dir"Cannabis_Use_Disorder" $sud_output_dir"Cannabis_Use_Disorder_Munged" &
    # OUD: 
       run_munge_sumstats $sud_input_dir"Opioid_Use_Disorder.txt" $sud_output_dir"Opioid_Use_Disorder_Munged" &
    # mvSUD:
        run_munge_sumstats $sud_input_dir"Shared_Substance_Use_Factor.txt" $sud_output_dir"Shared_Substance_Use_Factor_Munged" 1025550 &

# Aging Traits:
    # Alzheimer's Bellenguez: 
        run_munge_sumstats $aging_input_dir"Alzheimers.tsv" $aging_output_dir"Alzheimers_Munged" &
    # Brain Age Gap: 
        run_munge_sumstats $aging_input_dir"Brain_Age_Gap.txt" $aging_output_dir"Brain_Age_Gap_Munged" 28104 &
    # Death Any Cause: 
        run_munge_sumstats $aging_input_dir"Death_Any_Cause" $aging_output_dir"Death_Any_Cause_Munged" 412181 &
    # Telomere Length: 
        run_munge_sumstats $aging_input_dir"Telomere_Length.tsv" $aging_output_dir"Telomere_Length_Munged" 472174 &
    # Frailty INdex:  
        run_munge_sumstats $aging_input_dir"Frailty_Index.txt" $aging_output_dir"Frailty_Index_Munged" 175226 &
    # Grim Age: 
        run_munge_sumstats $aging_input_dir"Grim_Age.txt" $aging_output_dir"Grim_Age_Munged" &
    # Healthspan: 
        run_munge_sumstats $aging_input_dir"Healthspan.txt" $aging_output_dir"Healthspan_Munged" 300447 &
    # Longevity: 
#        run_munge_sumstats $aging_input_dir"Longevity.txt" $aging_output_dir"Longevity_Munged" &
    # Parental Lifespan: 
        run_munge_sumstats $aging_input_dir"Parental_Lifespan.tsv" $aging_output_dir"Parental_Lifespan_Munged" &
    # mvAge: 
        run_munge_sumstats $aging_input_dir"Shared_Aging_Factor.txt" $aging_output_dir"Shared_Aging_Factor_Munged" 1958774 &

# Wait for all background processes to finish
wait
