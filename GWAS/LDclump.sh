#!/bin/bash
#SBATCH --job-name=LDclumping_071424    # Job name #name_bfile_date
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlama@meei.harvard.edu     # Where to send mail    
#SBATCH --mem=1G --nodes=1             # Job memory request
#SBATCH --output=sub.LDclumping_071424_%j.log   # Standard output and error log
#SBATCH --partition=short,medium,long


#you need the pgen files and output from gwas with only two columns (SNP and P)

module load plink3/alpha3
module load plink2/alpha5.11

export outFile=/data/Segre_Lab/users/jlama/GSA_new.All_040424/GWAS/FAME/Responder.UpdatedSteroidMedication/WithoutSham.PCupdated/tmp

for chr in {1..22} X
do

echo subseting pgen files for chr${chr} ...

plink2 --pfile ${outFile}/chr${chr}.finalQC.FAME.530 \
         --clump chr${chr}.ld \
         --clump-kb 500 \
         --clump-p1 1e-5 \
         --clump-r2 0.1 \
         --out chr${chr}.ld
done

