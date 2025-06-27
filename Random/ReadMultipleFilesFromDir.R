# List all *.txt files in both directories
files_dir1 <- list.files(dir1, pattern = "\\.ld$", full.names = TRUE)
files_dir2 <- list.files(dir2, pattern = "\\.txt$", full.names = TRUE)

# Extract keys from dir1 filenames (chr:pos:ref:alt)
get_key <- function(file_path) {
  fname <- basename(file_path)
  str_remove(fname, "\\.ld\\.txt$|\\.txt$") #remove .ld.txt or .ld or .txt from filename
}

keys_dir1 <- sapply(files_dir1, get_key)
# [1] "chr1:123456:A:G" "chr2:234567:T:C"

# Loop through each file in dir2
for (file2 in files_dir2) {
  fname2 <- basename(file2)
  
  # Check if any key from dir1 is part of filename
  matched_key <- keys_dir1[str_detect(fname2, fixed(keys_dir1))]
  #The fixed() function in R, used with stringr functions like str_detect(), tells the function to treat the pattern as a literal string, not a regular expression (where ':' has different meaning)
  
  #If multiple key is match just choose one
  if (length(matched_key) > 0) {
    key <- matched_key[1]  # Assuming one match
    file1 <- files_dir1[keys_dir1 == key]
    
    # Read the two files
    df1 <- read.table(file1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    df2 <- read.table(file2, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Merge by common columns or add logic to merge as needed
    merge_df <- full_join(df1, df2, by = intersect(names(df1), names(df2)))
    
    # Write to output directory
    out_file <- file.path(out_dir, paste0("merge.", key, ".tsv"))
    write.table(merge_df, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    message("Merged and saved: ", out_file)
  }
}
