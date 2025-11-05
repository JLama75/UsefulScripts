#You have 307,624,125 (~300 million) rows in info.tsv with 6 columns 
#and in gwas.linear you have 9 million rows and 16 columns. 
#You want to merge columns 1 and 7 named ["SNP ID", "R2"] from info.tsv with all columns of the gwas.linear, which has "SNP ID" column in column 3. 
#Essentially you want R2(info score) values for each SNPs in gwas.linear file. Write a efficient code to do this-

awk 'NR==FNR {r2[$1]=$7; next} 
     FNR==1 {print $0"\tR2"; next} 
     {print $0"\t" ( ($3 in r2) ? r2[$3] : "NA") }' info.tsv gwas.linear > gwas_with_r2.tsv
