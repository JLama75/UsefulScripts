#From external hard drive to cluster
#!/bin/bash

max_jobs=1       # number of parallel jobs
count=0          # counter for successful transfers

while IFS=$'\t' read -r LOCAL_DIR; do
    echo "Processing directory: $LOCAL_DIR"

    if [ -d "$LOCAL_DIR" ]; then
        job_count=0

        # Find only .crai and .md5 files in the folder
        #find "$LOCAL_DIR" -maxdepth 1 -type f \( -name "*.crai" -o -name "*.md5" \) | while read -r file; do
        find "$LOCAL_DIR" -maxdepth 1 -type f -name "*.cram" | while read -r file; do
            {
                rsync -ratlzv "$file" \
                  ogi-mbc-login.meei.harvard.edu:/gpfs/archive1/SIOP/CRAMS

                if [ $? -eq 0 ]; then
                    echo "Successfully transferred $file"
                    ((count++))
                else
                    echo "Error transferring $file" >&2
                fi
            } &  # background job

            ((job_count++))
            if (( job_count % max_jobs == 0 )); then
                wait  # wait for batch to finish
            fi
        done

        wait  # wait for any leftover jobs
    else
        echo "Directory does not exist: $LOCAL_DIR" >&2
    fi
done < SIOP_transfer_files.tsv

echo "All files processed. Total files successfully transferred: $count"
