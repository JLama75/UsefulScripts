#In bash, Merge two files based on values in column "ID" in file A and in master_File.tsv so that the column "rsID" of the master_File.tsv is added 
# and save as a new file with the name of the file A in a new directory /data/results/. 
# The fileA will be the 48 files in the directory /data 
#fileA: ID A1 A1.FREQ SE BETA P 
#master_File.tsv: ID rsID 
#new_File <- merge(fileA, master_File.tsv, by="ID") 
#new_File ID A1 A1.FREQ SE BETA P rsID

#!/bin/bash
#SBATCH --job-name=AddRsID   # Job name #name_bfile_date
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlama@meei.harvard.edu     # Where to send mail    
#SBATCH --mem=10G --nodes=1             # Job memory request
#SBATCH --output=AddRsID_%j.log   # Standard output and error log
#SBATCH --partition=short,medium,long

mkdir -p /data/results

sort -k1,1 master_File.tsv > master_File.sorted.tsv

for f in /data/*; do
  fname=$(basename "$f")

  {
    # header: original header + rsID
    echo -e "$(head -n 1 "$f")\trsID"

    # left join on ID
    join -t $'\t' \
      -a 1 \
      -e "NA" \
      -o auto \
      -1 1 -2 1 \
      <(tail -n +2 "$f" | sort -k1,1) \
      <(tail -n +2 master_File.sorted.tsv)
  } > /data/results/"$fname"

done

###############################################
############# Explanation of the above code....

#What happens
head -n 1 "$f"
→ grabs the header line from fileA
"$( ... )\trsID"
→ appends a new column name rsID
-e
→ allows \t to be interpreted as a tab

#Resulting header
ID   A1   A1.FREQ   SE   BETA   P   rsID

#Left join
join -t $'\t' \
Use tab as the field delimiter

# -a 1 = left join
Keep all rows from fileA, even if there is no match in master_File.tsv

# -e "NA"
When there is no matching ID, output NA
Prevents empty fields

# -o auto
Automatically prints:
all columns from fileA
plus the matched column(s) from master file (rsID)

#Join key specification
-1 1 -2 1
Join on column 1 in:
fileA (-1 1)
master file (-2 1)
That column is ID
#As both ID column is in column 1 for both the files.

<(tail -n +2 "$f" | sort -k1,1)
tail -n +2
→ skips the header (join can’t handle headers)
sort -k1,1
→ sorts by ID (required by join)
This is fileA data only.

<(tail -n +2 master_File.sorted.tsv)
Skips header of master file
Already sorted earlier
