#!/bin/bash
#SBATCH --job-name=PCA    # Job name
#SBATCH --mem=20G      # Job memory request
#SBATCH --output=PCA_%j.log   # Standard output and error log
#SBATCH --partition=normal,mee-compute,long

#module load PLINK/1.9b-6.10
#module load plink3/alpha3

#OUT="/data/rossinlab/users/jlama/AMD/Post_imputation/Post-Impute/results"
#"${OUT}/chr${chr}.finalQC.pvar"

plink2=/data/rossinlab/users/jlama/AMD/Merge_Data/plink2/plink2
plink=/data/rossinlab/users/jlama/AMD/Merge_Data/plink1.9/plink
SEX=/data/rossinlab/users/jlama/AMD/pheno/AMD_selfReportedSex_updated.tsv
ances=/data/rossinlab/users/jlama/AMD/Shards/inferred_ancestry.tsv
export All=/data/rossinlab/users/jlama/AMD/Post_imputation/Post-Impute/Final_vcf/concat/allchrom.finalQC.vcf.gz
export OUT=/data/rossinlab/users/jlama/AMD/PCA/AMD_705_EUR_imputed_QC

${plink2} --vcf ${All} --update-sex ${SEX} --split-par hg38 --double-id --make-bed --out ${OUT}

python pca.py --study_in ${OUT} --format_in bfile --output_file study_pca --verbose --min_maf 0.01 --max_missingness 0.02 --max_r2 0.1 --window_size 200 --step_size 100 --plink1_path ${plink} --plink2_path ${plink2} --underscore_sample_id --file_inferred_ancestry ${ances} --flag_plot True --number_of_PCs 20

