#!/bin/bash

# Header line
header="p1,p2,rg,se,z,p,h2_obs,h2_obs_se,h2_int,h2_int_se,gcov_int,gcov_int_se"

# Create a temporary file to store the content
temp_file=$(mktemp)

# Loop through each .log file
for logfile in /scratch/aalab/jared/aging_suds_project/ldsc/results/EUR/*/*.log; do
    # Extract content between the specified patterns excluding the header line and remove the first column
    content=$(awk '/p1\s+p2\s+rg\s+se\s+z\s+p\s+h2_obs\s+h2_obs_se\s+h2_int\s+h2_int_se\s+gcov_int\s+gcov_int_se/,/Analysis finished at/' "$logfile" | sed '1d;$d' | cut -d ',' -f 2-)
    # Reformat content to have comma as separator
    formatted_content=$(echo "$content" | awk '{$1=$1}1' OFS=",")
    # Print content to the temporary file
    echo "$formatted_content" >> "$temp_file"
done

# Prepend the header line to the temporary file
echo "$header" > /scratch/aalab/jared/aging_suds_project/ldsc/results/EUR/sud_aging_ldsc_results.txt
cat "$temp_file" >> /scratch/aalab/jared/aging_suds_project/ldsc/results/EUR/sud_aging_ldsc_results.txt

# Remove unwanted strings from the resultant file

# Define patterns to remove
patterns=(
    "/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/sudTraits/munged_sumstats/EUR/"
    "/scratch/aalab/jared/aging_suds_project/ldsc/sumstats/agingTraits/munged_sumstats/EUR/"
    "_Munged.sumstats.gz"
)

# Iterate over each pattern and remove it from the file using sed
for pattern in "${patterns[@]}"; do
    sed -i "s|$pattern||g" "/scratch/aalab/jared/aging_suds_project/ldsc/results/EUR/sud_aging_ldsc_results.txt"
done
# Remove the temporary file
rm "$temp_file"
