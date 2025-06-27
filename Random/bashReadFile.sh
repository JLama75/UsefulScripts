list_of_SNPs=()

while IFS= read -r snp; do
  list_of_SNPs+=("$snp")
done < text.txt

for snp in "${list_of_SNPs[@]}"; do
  chrom=$(echo "$snp" | cut -d':' -f1)
  echo "SNP: $snp → Chromosome: $chrom"
done
