#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G

#eval $(spack load --sh python@2.7.18 )

trait1_filename=$(basename "$1" _Munged.sumstats.gz)
trait2_filename=$(basename "$2" _Munged.sumstats.gz)

mkdir "/scratch/aalab/jared/aging_suds_project/ldsc/results/EUR/${trait1_filename}"

# Run LDSC with each file pair
/ref/aalab/software/ldsc/env/bin/ldsc.py \
--rg "$1","$2" \
--ref-ld-chr "/scratch/aalab/jared/aging_suds_project/ldsc/reference_panels/eur_w_ld_chr/" \
--w-ld-chr "/scratch/aalab/jared/aging_suds_project/ldsc/reference_panels/eur_w_ld_chr/" \
--out "/scratch/aalab/jared/aging_suds_project/ldsc/results/EUR/${trait1_filename}/${trait1_filename}_${trait2_filename}_LDSC_results"
