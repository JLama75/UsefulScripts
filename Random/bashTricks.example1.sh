#!/bin/bash
#SBATCH --job-name=concat.assoc_pgen   # Job name
#SBATCH --mem=5G            # Job memory request
#SBATCH --partition=short,medium,long,xnetwork

outFile=$1
trait=$2
shift 2             # Remove the first two arguments so "$@" contains only SNPs

# Store remaining arguments (SNPs) in an array
snpList=("$@")
echo "Trait: $trait"
echo "SNPs: ${snpList[@]}"

dir="/gpfs/fs1/data/Segre_Lab/users/jlama/GSA_new.All_040424/Regenie_GWAS/Raw/results/GWAS"
mkdir -p ${dir}

echo CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA > "${output}/${trait}.step2.Max_IOPRise.regenie"

for chr in {1..22} X 
do
awk 'NR > 1 {print $0}' "${outFile}/step2.${trait}.${analysis}.chr${chr}"_Max_IOPRise.regenie >> "${dir}/${trait}.step2.Max_IOPRise.regenie"
done

awk 'NR==1 {print $3, $5, $6, $10, $11 "P"; next} {P=10^(-$13); print $0, P}' "${dir}/${trait}.step2.Max_IOPRise.regenie" > "${dir}/${trait}.step2.Max_IOPRise.with_P.regenie" 


# Use grep with word boundaries to match SNPs exactly
echo ID ALLELE1 A1FREQ BETA SE P > "${dir}/${trait}.raw.sumStat.tsv"
grep -wF -f <(printf "%s\n" "${snpList[@]}") "${dir}/${trait}.step2.Max_IOPRise.with_P.regenie" >> "${dir}/${trait}.raw.sumStat.tsv"
#grep -F -f patterns.txt file where patterns.txt is temporary list with snpIds from <(printf "%s\n" "${snpList[@]}")
