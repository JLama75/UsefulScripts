#You have 307,624,125 (~300 million) rows in info.tsv with 6 columns 
#and in gwas.linear you have 9 million rows and 16 columns. 
#You want to merge columns 1 and 7 named ["SNP ID", "R2"] from info.tsv with all columns of the gwas.linear, which has "SNP ID" column in column 3. 
#Essentially you want R2(info score) values for each SNPs in gwas.linear file. Write a efficient code to do this-

awk 'NR==FNR {r2[$1]=$7; next} 
     FNR==1 {print $0"\tR2"; next} 
     {print $0"\t" ( ($3 in r2) ? r2[$3] : "NA") }' info.tsv gwas.linear > gwas_with_r2.tsv

#Building a in-memory lookup table r2 mapping the 'SNP ID' (column 1) of info.tsv --> R2 (column 7)
#Then read gwas.linear and for each line append the R2 value found for the SNP ID located in column 3; If no match, append NA
#Write the combined final output to gwas_with_r2.tsv

#awk processes input files in order of info.tsv first and then gwas.linear
#NR=total no. of input records in across all files
#FNR= no. of input recodrs seen so far in the current file
#NR==FNR is true only while awk is reading the first file (info.tsv). Once it moves to gwas.linera, NR continues to grow but FNR resets to 1 to the new file

#r2[$1] = $7 — create (or set) an element in the associative array r2 with key = field 1 of the current line ($1, assumed to be the SNP ID) and value = field 7 ($7, assumed to be R2).
#next — stop processing the current record and move to the next input line.

FNR==1 { print $0"\tinfo"; next }
#This block runs when starting the second file (gwas.linear) because FNR resets to 1 at the first line of each file, which is the header line
#next — skip remaining blocks for the header line.

{ print $0"\t" ( ($3 in r2) ? r2[$3] : "NA") }
# This block is for every non-header line
#($3 in r2) — checks whether the SNP ID found in column 3 of the current gwas.linear line is a key present in the r2 associative array built from info.tsv.
#The ? : is a ternary operator:
#If ($3 in r2) is true → return r2[$3] (the matching R2 value).
#If false → return "NA".
