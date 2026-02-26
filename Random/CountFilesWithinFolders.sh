#!/bin/bash

# Set your target path
TARGET_DIR="/shared/ssd_14T/home/jyotilama/Choroidalizer/Scripts_from_Souvick/inputs"

# Loop through each directory inside TARGET_DIR
for dir in "$TARGET_DIR"/*/; do
    # Skip if no directories found
    [ -d "$dir" ] || continue

    # Count regular files inside the directory (non-recursive)
    file_count=$(find "$dir" -maxdepth 1 -type f | wc -l) #-maxdepth 1  --> only look at one level deep

    # Print folder name and count
    echo "$(basename "$dir") : $file_count files"
done
