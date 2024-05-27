#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G

# Directories for our two categories of variables:
sud_files=$(find /scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/munged_sumstats/EUR/*.gz -maxdepth 1 -type f)
aging_files=$(find /scratch/aalab/jared/aging_suds_project/ldsc/sumstats/agingTraits/munged_sumstats/EUR/*.gz -maxdepth 1 -type f)

# Loop through sud and aging files to get all possible pairwise combinations: 
for sud_file in $sud_files; do
    for aging_file in $aging_files; do    
       sbatch /scratch/aalab/jared/aging_suds_project/ldsc/scripts/sud_aging_ldsc_EUR.sh "$sud_file" "$aging_file"
    done
done

# Process the output:
sleep 45 # To ensure that the next script isn't run until the previous one has finished
sbatch /scratch/aalab/jared/aging_suds_project/ldsc/scripts/collect_LDSC_Results_EUR.sh # Gather the results into a single file
mv sud_aging_ldsc_results.txt sud_aging_ldsc_results_betweencategory.txt


# Loop through sud files to get all possible pairwise combinations within sud category:
for sud_file1 in $sud_files; do
    for sud_file2 in $sud_files; do
            sbatch /scratch/aalab/jared/aging_suds_project/ldsc/scripts/sud_aging_ldsc_EUR.sh "$sud_file1" "$sud_file2"
    done
done

# Loop through aging files to get all possible pairwise combinations within aging category:
for aging_file1 in $aging_files; do
    for aging_file2 in $aging_files; do
            sbatch /scratch/aalab/jared/aging_suds_project/ldsc/scripts/sud_aging_ldsc_EUR.sh "$aging_file1" "$aging_file2"
    done
done

# Process the output:
sleep 45 # To ensure that the next script isn't run until the previous one has finished
sbatch /scratch/aalab/jared/aging_suds_project/ldsc/scripts/collect_LDSC_Results_EUR.sh # Gather the results into a single file
sleep 15 # To ensure that the next command isn't run until the previous script has finished
mv slurm-* /scratch/aalab/jared/aging_suds_project/ldsc/scripts/slurms/ # Move the 10,000 slurm output files into a new directory
