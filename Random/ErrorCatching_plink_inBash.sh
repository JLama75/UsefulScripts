#!/bin/bash
#SBATCH --job-name=haplotype_1000G # Job name #name_bfile_date
#SBATCH --mem=70G 
#SBATCH --time=0-05:00:00
#SBATCH --output=haplotype_Rhodopsin_1000G_%j.log   # Standard output and error log
#SBATCH --partition=mee-compute,normal,long

module load BCFtools/1.21-GCC-13.3.0

cd /data/piercelab/users/jlama/Kinga/jl_reanalysis/imputation/Phasing_Imputation_Michigan_HRC
QC=/data/piercelab/users/jlama/Kinga/jl_reanalysis/imputation/Phasing_Imputation_Michigan_HRC/QC
wdir=/data/piercelab/users/jlama/Kinga/jl_reanalysis/imputation/Phasing_Imputation_Michigan_HRC/hg19_extract_RHO_SNPs/Manual_extracting_haplotypeFreq/AfterQC

plink2=/data/rossinlab/users/jlama/AMD/Merge_Data/plink2/plink2
plink=/data/rossinlab/users/jlama/AMD/Merge_Data/plink1.9/plink

SEX=/data/piercelab/users/jlama/Kinga/data/merged_130_samples_RHO_P23H.sex.tsv
data=/data/piercelab/users/jlama/Kinga/data/
file="chr3.dose.vcf.gz"
outfile="${wdir}/chr3.QCed_100kb_RHO.vcf.gz"

#set -e  

run_cmd() {
    "$@"
    local exit_code=$?
    if [ ${exit_code} -ne 0 ]; then
        echo "ERROR: Command failed with exit code ${exit_code}:" >&2
        echo "  $*" >&2
        exit ${exit_code}
    fi
}

run_cmd ${plink2} --vcf "${file}" dosage=HDS --double-id --exclude-if-info "R2<=0.4" --export vcf-4.2 vcf-dosage=HDS id-paste=iid --out "${QC}/chr3.dose.rsq0.4" --threads 12
echo -e "updating sex information"
#Assign predicted sex info into .fam file
run_cmd ${plink2} --vcf "${QC}/chr3.dose.rsq0.4" dosage=HDS --double-id --update-sex ${SEX} --make-pgen --out "${QC}/chr3.dose.rsq0.4.sexUpd" --threads 12
run_cmd ${plink2} --pfile  "${QC}/chr3.dose.rsq0.4.sexUpd" --set-hh-missing --make-pgen --out "${QC}/chr3.dose.rsq0.6.sexUpd.ind.hhmissing" --threads 12

run_cmd ${plink2} --pfile  "${QC}/chr3.dose.rsq0.6.sexUpd.ind.hhmissing" --make-pgen --geno 0.02 dosage --out "${QC}/chr3.dose.rsq0.6.sexUpd.ind.hhmissing.geno" --threads 12
run_cmd ${plink2} --pfile "${QC}/chr3.dose.rsq0.6.sexUpd.ind.hhmissing.geno" --maf 0.01 --hardy midp --keep "${data}/predicted_${ances}_samples.tsv" --out "${QC}/chr3_${ances}_6" --threads 12 

run_cmd awk '$10<0.0000001 {print $0}'  "${QC}/chr3}_${ances}_6.hardy" >>  "${INPUT}/chr3_step6_HWE_remove.txt"

echo -e "\n N hwe failed for autosomes chr3 \n" 
run_cmd awk '{print $2}' "${QC}/chr3_step6_HWE_remove.txt" | sort | uniq | wc -l > ${QC}/HWE.failed.log

run_cmd ${plink2} --pfile "${QC}/chr3.dose.rsq0.6.sexUpd.ind.hhmissing.geno" --exclude "${QC}/chr3_step6_HWE_remove.txt" --make-pgen --out "${QC}/chr3.finalQC" --threads 12
run_cmd ${plink2} --pfile "${QC}/chr3.dose.rsq0.6.sexUpd.ind.hhmissing.geno" --exclude "${QC}/chr3_step6_HWE_remove.txt" --export vcf-4.2 vcf-dosage=HDS id-paste=iid --out "${QC}/chr3.finalQC" --threads 12

run_cmd ${plink2} --pfile "${QC}/chr3.finalQC" --maf 0.01 --make-pgen --out "${QC}/chr3.finalQC.maf0.01" --threads 12

echo "Finished QC for chromosome 3"

##################################################################################################


#Check from 100kb from the Rodopsin's start and end site 
#Hg19
#chr3:129247577-129252561
#
CHR=3
START=$((129247577 - 100000))
END=$((129252561 + 100000))

echo $START
echo $END

# Check chromosome naming in your VCF
run_cmd bcftools view -h "${QC}/chr3.finalQC.vcf" | grep "^##contig" | head -5

run_cmd bcftools view \
        -r "${CHR}:${START}-${END}" \
        "${QC}/chr3.finalQC.vcf" \
        -o ${outfile} -Oz

bcftools index ${outfile}
