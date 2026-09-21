#!/usr/bin/env bash


#Code that from the main directory finds any file that ends in *_blueHighlight.tsv within the subfolders and adds a column called Trait and adds in the filename except for the _blueHighlight.tsv. Then concatenate all the  *_blueHighlight.tsv into a file called subthresholdHits.tsv.
set -euo pipefail

# Run this from your main directory
MAIN_DIR="."
OUTPUT="subthresholdHits.tsv"

# Remove old output if it exists, so we don't accidentally append to stale data
rm -f "$OUTPUT"

first_file=true

# Find all *_blueHighlight.tsv files in subfolders (not the main dir itself)
while IFS= read -r -d '' file; do #"raw" mode: don't let backslashes be treated as escape characters. -d '' — sets the delimiter IFS (Internal Field Separator) normally trims leading/trailing whitespace when read splits input.
    fname=$(basename "$file")
    trait="${fname%_blueHighlight.tsv}"

    tmpfile=$(mktemp)

    # Add a "Trait" column: header gets "Trait", every data row gets the trait name
    awk -F'\t' -v OFS='\t' -v trait="$trait" '
        NR==1 { print $0, "Trait"; next }
        { print $0, trait }
    ' "$file" > "$tmpfile"

    # Overwrite the original file with the new column added
    mv "$tmpfile" "$file"

    # Append to the concatenated output — header only from the first file
    if $first_file; then
        cat "$file" >> "$OUTPUT"
        first_file=false
    else
        tail -n +2 "$file" >> "$OUTPUT"
    fi

    echo "Processed: $file (Trait=$trait)"

done < <(find "$MAIN_DIR" -mindepth 2 -type f -name "*_blueHighlight.tsv" -print0)
#-mindepth 2: skip anything directly inside $MAIN_DIR itself; only look at files at least 2 levels deep (i.e., inside subfolders, not the top level)
#<(...) — process substitution
echo "Done. Combined file written to: $OUTPUT"



