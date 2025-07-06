>>> python
>>> import math
>>> result = -math.log10(5e-8)
>>> print(result)
7.301029995663981
>>> quit()

#ssh jl2251@172.27.81.202  

##chmod u+x,g+x file.sh file.py #giving group and user file execution permission.
#source $HOME/miniconda3/bin/activate/hail
##source ~/.bash_profile 

#short-cut
alias WES="cd /data/Segre_Lab/users/jlama/WES_new.ALL_050824"
alias GSA="cd /data/Segre_Lab/users/jlama/GSA_new.All_040424"


alias GSA="cd /data/mee-ogi/data/Sobrin_Lab/users/jlama/GSA_new.All_040424"
alias WES="cd /data/mee-ogi/data/Sobrin_Lab/users/jlama/WES_new.ALL_050824"

/data/Segre_Lab/data/DRCR/analysis/colocalization/run_retina_013124
MAGMA: /data/Segre_Lab/data/DRCR/analysis/MAGMA

#Transfer folder/files from computer to cluster 
rsync -ravoP your_file user@ogi-mbc-login.meei.harvard.edu:/gpfs/fs1/home/jl2251/GWAS_tutorial
#rsync -ravoP PLINK_2019_math.etal jl2251@ogi-mbc-login.meei.harvard.edu:/gpfs/fs1/home/jl2251/GWAS_tutorial
rsync -ravoP Phenotype_data_steroid_IOP.xlsx jl2251@ogi-mbc-login.meei.harvard.edu:/data/Segre_Lab/users/jlama/WES/QC/sampleQcTool/Draft_data


#Copying from OGI cluster to computer:
rsync -ravoP user@ogi-mbc-login.meei.harvard.edu:/gpfs/fs1/data/Segre_Lab/data/DIR_OF_INTEREST/FILE_NAME.txt 
rsync -ravoP jl2251@ogi-mbc-login.meei.harvard.edu:/gpfs/fs1/data/Segre_Lab/data/DIR_OF_INTEREST/FILE_NAME.txt 
rsync -ravoP jl2251@eristwofs.partners.org:/data/mee-ogi/data/Sobrin_Lab/users/jlama/WES_new.ALL_050824/PCA/PCplot.png /Users/jl2251/Downloads 

#invoke sinfo to see the list of nodes and their states in SLURM hpcc
sinfo

#PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
#short        up   23:59:00      1   idle compute01
#medium       up 8-00:00:00      5   idle compute[02-03,05-07]
#medium       up 8-00:00:00      1   down compute04
#long*        up 107-00:00:      3   idle compute[09-10,12]
#xnetwork     up   infinite      1   idle compute11

#This command is showing detailed information about a compute node named compute02 in an HPCC system. 
scontrol show node compute02

#	NodeName=compute02 		Arch=x86_64 	CoresPerSocket=7 
#   CPUAlloc=0		 CPUTot=14		 CPULoad=0.01
#   OS=Linux 3.10.0-1160.25.1.el7.x86_64 #1 SMP Wed Apr 28 21:49:45 UTC 2021 
#   RealMemory=115000(115GB)	 AllocMem=0 	FreeMem=99340 (99GB) 	Sockets=2 	Boards=1
#   State=IDLE 	ThreadsPerCore=1 TmpDisk=0 Weight=1 Owner=N/A MCS_label=N/A
#   Partitions=medium 
#   BootTime=2023-03-02	T09:49:38 SlurmdStartTime=2023-06-21	T10:25:13
#   CfgTRES=cpu=14,mem=115000M,billing=14
#....

#HPCC might have 11 nodes with 14 cpus in one node.
#Nodes:
#Asking for 2 nodes with cpus-per-task=1 means that the job will run on two separate nodes, each with its own set of resources, including memory, disk space, and possibly other hardware components.
#Asking for 1 node with cpus-per-task=2 means that the job will run on a single node, utilizing the CPU cores available on that node.
#Resource Allocation:
#When requesting 2 nodes with cpus-per-task=1, the job's tasks will be spread across two different nodes, potentially allowing for better load balancing and avoiding resource contention if the nodes are heavily utilized.
#When requesting 1 node with cpus-per-task=2, all tasks of the job will run on the same node, which could potentially lead to resource contention if multiple tasks are competing for the same resources, such as memory bandwidth or disk I/O
#Communication Overhead:
#Running tasks on separate nodes may incur communication overhead if the tasks need to communicate frequently during execution, as data may need to be transferred between nodes over the network.
#Running tasks on the same node may reduce communication overhead since tasks can communicate more efficiently through shared memory.

#To extract the column number of the column with column name = pct_chimeras
awk 'NR==1 { for(i=1; i<=NF; i++) if($i=="pct_chimeras") { print i; exit } }' samples_steroid_IOP.tsv
awk 'BEGIN {FS="\t"; OFS="\t"} NR>1 {print $3, $1, $2}' $file > $snpLoc

#Iterate over the column names only
head -n 1 score_bam_QC.tsv | awk 'BEGIN {FS="\t"} {for (i=1; i<=NF; i++) print $i}' | more
#How many columns
head -n 1 plink_qt_bfile_022224.MaxIOPRise_rank.glm.linear.adjusted | awk -F'\t' '{print NF}'

#making sub vcf file
zcat input.vcf.gz | bcftools view -Oz -o sub.vcf.gz -s $(bcftools query -l input.vcf.gz | head -n 1000 | tr '\n' ',') 
zcat chr11.dose.vcf.gz | bcftools view -Oz -o sub.vcf.gz -s $(zcat chr11.dose.vcf.gz |bcftools query -l  | head -n 1000 | tr '\n' ',') --force-samples
bcftools view -i 'INFO/RRSQ>0.8' sub.vcf.gz -Oz -o filtered.vcf.gz

#ignoring header.
awk '!/^#/ { print $18 }' filename
grep -E "Highlighted.*" your_file.txt

#exact matching
grep '^rs7773324$' your_file.txt
grep -w 'rs7773324' your_file.txt
awk '!/^#/ { print $1 }'

gzip 1412_GSA_Callset_03292024.vcf
conda activate gsutil
gcloud auth login

gsutil -m cp -r 'gs://fc-2e94ced4-e58e-4804-8b4a-a407f668a7e1/Exome_VCF/Sobrin_MEEI_14_WES_VCF.vcf.bgz' Exome_VCF/Sobrin_MEEI_14_WES_VCF.vcf.bgz
vcftools --vcf $VCF --remove RemoveOldIDs.041124.txt --recode --recode-INFO-all --out 1412_GSA_Callset_03292024.Removed_oldID.vcf
/data/Segre_Lab/data/Sobrin_R01_GC_OHTN/genotype/WES/TM_SC_WES/Sobrin_MEEI_14_WES_VCF_04_22_24 


awk '{print $3}' 1412_GSA_Callset_042624_5.hh | sort | uniq -c | sort -k1,1nr > 1412_GSA_Callset_042624_5.uniqVariants.hh 

bcftools query -f'%CHROM\t%POS\t%ID\t%REF,%ALT,[ %GT]\n' 1412_GSA_Callset_03292024.vcf.gz > sample.vcf| 
grep "GSA-rs75891733" sample.vcf > sample.var


head -n 1 Sobrin_MEEI_14_WES_sample_qc.tsv | awk 'BEGIN {FS="\t"} {for (i=1; i<=NF; i++) print $i}' | less

module load vcftools
vcftools --vcf /data/Steroid_IOP_Sobrin_R01/GSA/04022024_MEE_GSA_merged_JL/1412_GSA_Callset_03292024.vcf.gz --remove RemoveOldIDs.041124.txt --recode --recode-INFO-all --out 1412_GSA_Callset_03292024.Removed_oldID.vcf #1377 samples
vcftools --gzvcf /data/Steroid_IOP_Sobrin_R01/GSA/04022024_MEE_GSA_merged_JL/1412_GSA_Callset_03292024.vcf.gz --snps /data/Segre_Lab/users/jlama/GSA_new.All_040424/qc2/1412_GSA_Callset_042624_11_Duplicates.list --recode --recode-INFO-all --out duplicateSubset

bcftools query -f'%CHROM\t%POS\t%ID\t%REF,%ALT,[ %GT]\n' duplicateSubset.recode.vcf  > duplicateSubset.recode.sample

vcftools --gzvcf /data/Steroid_IOP_Sobrin_R01/GSA/04022024_MEE_GSA_merged_JL/1412_GSA_Callset_03292024.vcf.gz --snp GSA-rs113445782 --snp rs113445782 --recode --recode-INFO-all --out duplicateSubset1
vcftools --gzvcf /data/Steroid_IOP_Sobrin_R01/GSA/04022024_MEE_GSA_merged_JL/1412_GSA_Callset_03292024.vcf.gz --snp GSA-rs36083022 --snp rs36083022 --recode --recode-INFO-all --out duplicateSubset2

bcftools query -f'%CHROM\t%POS\t%ID\t%REF,%ALT,[ %GT]\n' duplicateSubset1.recode.vcf  > duplicateSubset1.recode.sample
bcftools query -f'%CHROM\t%POS\t%ID\t%REF,%ALT,[ %GT]\n' duplicateSubset2.recode.vcf  > duplicateSubset2.recode.sample

bcftools query -f'%CHROM\t%POS\t%ID\t%REF,%ALT,[ %MAF;%R2]\n' duplicateSubset2.recode.vcf  > duplicateSubset2.recode.sample
bcftools query -f'%CHROM\t%POS\t%ID\t%REF,%ALT,[ %MAF]\n' duplicateSubset2.recode.vcf  > duplicateSubset2.recode.sample

OriginalVCR[c('chr', 'start', 'ID', 'ref')] <- str_split_fixed(OriginalVCR$V1, '\t', 4)
bcftools query -f'%CHROM\t%POS\t%ID\n'  
bcftools query -f'%CHROM\t%POS\t%ID\t%R2\n' fame.subset.530ID.fame.txt.recode.vcf
bcftools query -f'%ID\t%R2\n' fame.subset.530ID.fame.txt.recode.vcf > fame.subset.infoScore.txt

bcftools query -f'%CHROM\t%POS\t%ID\n' fame.subset.530ID.fame.txt.recode.vcf
awk {gsub(/^chr/, "", $1); print $0, $1 ":" $2}

bcftools view -R regions.bed.tabix vcf > subset.vcf
bcftools query -f'%ID\t%FILTER\n' SIOP.variant_qc.092624.annotated.vcf.bgz | less
 
bcftools view -f HWEEUR SIOP.variant_qc.093024.annotated.vcf.bgz  | bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/AF\n' > annotate_variants_forHWE.txt
bcftools view -f HWEEUR SIOP.variant_qc.093024.annotated.vcf.bgz | \
bcftools query -i 'INFO/AF < 0.01 || INFO/AF > 0.99' -f '%CHROM\t%POS\t%ID\t%INFO/AF\n' > annotate_variants_forHWE_MAF.01.threshold.txt


bcftools view -f HWEEUR SIOP.variant_qc.093024q.filtered.vcf.bgz | bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/AF\n' > filtered_variants_forHWE.txt
bcftools view -f HWEEUR SIOP.variant_qc.093024.filtered.vcf.bgz| \
bcftools query -i 'INFO/AF < 0.01 || INFO/AF > 0.99' -f '%CHROM\t%POS\t%ID\t%INFO/AF\n' > filtered_variants_forHWE_MAF.01.threshold.txt


bcftools query -i 'INFO/GAB || INFO/GGQ || INFO/GNP' -f '%CHROM\t%POS\t%ID\t%INFO/AF\t%INFO/AN\n' SIOP.variant_qc.093024.annotated.vcf.bgz  > missingness.others.txt

awk 'BEGIN{FS=OFS="\t"} {split($4, a, ","); $4=a[1]; $5=a[2]; print}' |  awk '$4<0.01 || $4>0.99' | wc -l

module load vcftools
vcftools -vcf gnomad.exomes.v4.1.sites.chr2.vcf --snps chr2.sig1e-5.527IDs.txt

sacct -j 585500 --format="JobID%20,JobName,User,Partition,NodeList,Elapsed,CPUTime,State,AllocTRES%32,MaxRSS,MaxVMSize"
sacct -j 912274 --format="JobID%20,JobName,User,Partition,NodeList,Elapsed,CPUTime,State,AllocTRES%32,MaxRSS,MaxVMSize"

total_size=$(find /data/Segre_Lab/users/jlama/GSA_new.All_040424/imputation -type f -name "chr*.dose.vcf.gz" -exec du -ch {} + | grep total$ | awk '{print $1}')
total_size=$(find /data/Segre_Lab/users/jlama/GSA_new.All_040424/imputation -type f -name "chr2*.dose.vcf.gz" -exec du -ch {} + | grep total$ | awk '{print $1}')

export LD_LIBRARY_PATH=$metal_env/lib:$LD_LIBRARY_PATH
ls -l $metal_env/lib/libstdc++.so.6

#######################################################################################################################################################################

Rscript my_r_script.R "$FILE_NAME"

#R script

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the filename argument is provided
if (length(args) == 0) {
  stop("No file name provided.")
}

# The first argument is the file name
file_name <- args[1]

# Print the file name (for testing purposes)
cat("File name provided:", file_name, "\n")

# Read the file (assuming it's a CSV file for this example)
data <- read_table(file_name, col_name=T)

# Print the first few lines of the file to verify it was read correctly
print(head(data))

# Further processing can be done here...

export out=/data/Segre_Lab/users/jlama/GSA_new.All_040424/GWAS/LD.clumping/FAME.shamRemoved.530IDs.UpdPCs/Genotype.population.scatter
bcftools view -i 'ID="rs13425173"' $OUT > ${out}/allchrom.finalQC.FAME.rs13425173.530ID.txt
bcftools query -f'%ID\t%[ %GT]\n' ${out}/allchrom.finalQC.FAME.rs13425173.530ID.txt > ${out}/allchrom.finalQC.FAME.rs13425173.530ID.final.txt
bcftools view --include ID==@snps.list 1000Genomes.Norm.bcf


(bcftools query -l allchrom.finalQC.FAME.rs13425173.530ID.txt | awk 'BEGIN { printf("ID"); } { printf("\t%s", $1); } END { printf("\n"); }' && bcftools query -f '%ID[\t%GT]\n' allchrom.finalQC.FAME.rs13425173.530ID.txt) > allchrom.finalExtract.txt
bcftools query -l allchrom.finalQC.FAME.rs13425173.530ID.txt | awk 'BEGIN{printf("ID")} {printf("\t%s_GT",$1} END{printf("\n")}'  &&  bcftools query -f '%FILTER[\t%GT]\n' allchrom.finalQC.FAME.rs13425173.530ID.txt | less


awk 'BEGIN {OFS="\t"} {split($1, MarkerName, ":"); print $0, MarkerName[1], MarkerName[2]}' metal.stderr.gc.FAMEv1.MEEv1.MaxIOPRise_rank.072924.add.1tbl > metal.stderr.gc.FAMEv1.MEEv1.MaxIOPRise_rank.072924.add_2.1tbl
awk 'BEGIN {OFS="\t"; print "MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P-value", "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal", "chromosome", "position"} NR>1 {split($1, snp, ":"); print $0, snp[1], snp[2]}' metal.stderr.gc.FAMEv1.MEEv1.MaxIOPRise_rank.072924.add.1tbl  >  metal.stderr.gc.FAMEv1.MEEv1.MaxIOPRise_rank.072924.add_2.1tbl
file1=metal.stderr.gc.FAMEv1.MEEv1.MaxIOPRise_rank.072924.add_2.1tbl
file2=qt.mee.pgen.577IDs.080524.MaxIOPriseRINT.glm.linear
awk 'BEGIN {OFS="\t"} NR==FNR {a[$3] = $7; next} FNR==1 {print $0, "A_freq"} FNR>1 {print $0, a[$1]}' qt.mee.pgen.577IDs.080524.MaxIOPriseRINT.glm.linear metal.stderr.gc.FAMEv1.MEEv1.MaxIOPRise_rank.072924.add_2.1tbl > metal.stderr.gc.FAMEv1.MEEv1.MaxIOPRise_rank.072924.add_4.1tbl
(head -n 1 metal.stderr.gc.FAMEv1.MEEv1.MaxIOPRise_rank.072924.add_4.1tbl && tail -n +2 metal.stderr.gc.FAMEv1.MEEv1.MaxIOPRise_rank.072924.add_4.1tbl | sort -k12,12 -k13,13) > metal.stderr.gc.FAMEv1.MEEv1.MaxIOPRise_rank.072924.add_5.sorted.1tbl
sacct -j 843221 --format=JobID,MaxRSS,ReqMem,State

while IFS= read -r input; do
# Extract the ID from the input
echo "$input" 

# Use bcftools to filter based on the ID and append to the output file
bcftools view -i "ID=\"$input\"" ${DIR}/tmp/allchrom.finalQC.FAME.530.vcf.gz >> ${out}/fame.subset.530ID.fame.txt
done < ID.list

metal.top <- metal.top[match(topGenes$MarkerName, metal.top$MarkerName), ]


bcftools query -l $vcf > samples.txt
shuf -n 50 samples.txt > selected_samples.txt
bcftools view -S selected_samples.txt -o subset.vcf.gz -Oz $vcf


df.columns=['SNP','Symbol','Annotation']
display(df)
df_combined = df.groupby('Symbol')['SNP'].agg(lambda x: ','.join(x)).reset_index()

 grep -v '^#' /data/Segre_Lab/data/SCORE/GSA/imputation/Post_imputation_QC_dosage/Merged_freeze_dosage.pvar | wc -l
 
 
 Replace space with tab
  sed -i 's/ \+/\t/g' npc2.muc22.variant.log
  
  #Check no. of columns in your file
  awk -F'\t' '{print NF; exit}' your_file.txt
  
  
#Gencode
zgrep "^chr14" gencode.v47.annotation.gtf.gz | grep -i "LTBP2" |awk '{if ($3~/exon/){print}}'| wc -l
119
-bash-4.2$ zgrep "^chr14" gencode.v47.annotation.gtf.gz | grep -i "LTBP2" |awk '{if ($3~/exon/){print}}'| less
zgrep "^chr14" gencode.v47.annotation.gtf.gz | grep -i "LTBP2" |awk '{if ($3~/gene|exon/){print}}'

##########

du -hs * | sort -hr
 df -h ./
Filesystem                                Size  Used Avail Use% Mounted on
panfs://10.129.86.180/hpc/groups/mee-ogi   91T   74T   18T  81% /data/mee-ogi

du -ch SCORE_removeSampleQC.* SCORE_removeSampleQC_newid_removeMT_withSex_indelremove_splitPAR_chrX_sethhmissing_mergeMaleFemale_mergePAR_mergeAuto_sort_rmdup_miss_mono_HWE_chr*
Find combined storage of all files in the folder except some selected
du -ch --max-depth=1 | grep -vE './(PCA|sumstat|Ancestry)$' | grep total
du -ch --max-depth=1 | grep -vE './(*sh| *py | *R | *ipynb | *tsv |*linear |*adjusted |*merged| *locuszoom|PCA|Ancestry)$' | grep total
or,
find . -maxdepth 1 -type f ! -name 'PCA' ! -name 'sumstat' ! -name 'Ancestry' -exec du -ch {} + | grep total

(echo -e "SNP\tGENE\tFeature\tFeature_type\tConsequences" && grep "1:17036:T:C" ./step1/SIOP_FAME_annotated.101524.vcf | awk 'BEGIN{FS=OFS="\t"} !/^#/ {print $1, $4, $5, $6, $7}' ) > 1_17036_T_C_transcriptID.txt
(echo -e "SNP\tGENE\tFeature\tFeature_type\tConsequences" && grep "1:17036:T:C" /gpfs/fs1/data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/FAME_updated/step1/SIOP_FAME_annotated.101524.vcf | awk 'BEGIN{FS=OFS="\t"} !/^#/ {print $1, $4, $5, $6, $7}' ) | less
df -h /data/sobrinlab/

########
list_snp=('chr5:82257902:G:A' 'chr10:58592770:A:G' 'chr7:110054496:T:C' 
          'chr16:63988826:A:G','chr14:29954101:A:G','chr16:63988826:A:G'
          'chrX:21320210:C:A')
# Loop through each SNP and append the sbatch line to the output file
for snp in "${list_snp[@]}"; do
echo "sbatch submit.violinPlot.sh ${snp} \"${trait}\" \"${pheno}\"" >> "$outfile"
done

#########
snp_list1=("chr12:30572747:C:A" "chr4:118413677:G:A"  "chr5:82257902:G:A" "chr10:58592770:A:G")
snp_list2=("chr1:96900773:G:A" "chr1:96863842:C:A" "chr1:96246678:T:C")
snp_list3=("chr14:96255017:G:T" "chr14:96265492:G:A" "chr14:96267038:AAGAG:A" "chr14:96267226:TG:T")
snp_list4=("chr9:134198147:A:G" "chr17:63045506:A:G")
snp_list5=("chr16:63988826:A:G")
#Round to three digits after decimal point
awk '{$3 = sprintf("%.3f", $3); print}' "$tmp_file" > ./MAF_cal/tmp_file_rounded.tsv

# Header
echo -n -e "predicted_ancestry" > "$output_file"
for snp in "${snp_list1[@]}" "${snp_list2[@]}" "${snp_list3[@]}" "${snp_list4[@]}" "${snp_list5[@]}"; do
    echo -n -e "\t$snp" >> "$output_file"
done
echo >> "$output_file"


topHits_fame <- topHits_fame %>% 
      separate(ID, into = c("chrom", "pos", "extra1", "extra2"), sep = ":", remove = FALSE) %>% 
      dplyr::select(-extra1, -extra2)
