#!/bin/bash
#SBATCH --job-name=shards    # Job name
#SBATCH --mem=5G      # Job memory request
#SBATCH --time=0-08:00:00
#SBATCH --output=shards_%j.log   # Standard output and error log
#SBATCH --partition=normal,long,bigmem

export DATA_DIR=/PHShome/jl2251/Ines/Merge_Data/merged_AMD_GSA_CommVar
#export REF_DIR=/data/mee-ogi/data/OGI_shared/data/Resources_WGS_WES/b38/
#export REF_DIR=/data/mee-ogi/data/OGI_shared/data/Resources_WGS_WES/1000Genomes/1000G_Phase3_GRCh37_ref/ALL.chr.MERGED.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
export REF_DIR=ALL.chr.phase3_shapeit2_mvncall_integrated_v5a.20220325.genotypes.vcf.gz
export CODE_DIR=shards.dev_version_012225.py
export plink=/apps/software/plink/plink
export plink2=/apps/lib/plink/2.00a5.11/plink2
export covar=batches.tsv

${plink2} --bfile "${DATA_DIR}" --double-id --export vcf id-paste=iid --out "merged_AMD_GSA_CommVar"

python ${CODE_DIR} --filename_study "merged_AMD_GSA_CommVar.vcf" --format_study 'vcf' --filename_reference ${REF_DIR} --format_reference 'vcf' --reference_populations "1000_genomes_ancestry_hg19.tsv" --genome_build "b37" --underscore_sample_id --verbose --plink1_path ${plink} --plink2_path ${plink2} --study_covariates ${covar}

