#If you want to quickly overlap the positions in one dataframe (x) to based on intervals provided in another dataframe (y): use foverlaps() of R

#Plink bim file with columns - chr snp_id pos
#ChIP-seq bed intervals with column - chr start stop

ExtractSNPs <- function(bim, ChIP){
  # Convert to data.tables
  dt_bim <- as.data.table(bim)
  dt_co_ord <- as.data.table(ChIP)
  
  # Rename for clarity
  setnames(dt_bim, c("chr", "snp_id", "pos"))
  setnames(dt_co_ord, c("chr", "start", "stop"))
  
  # Create keys for fast overlap joins
  setkey(dt_bim, chr, pos, pos)
  setkey(dt_co_ord, chr, start, stop)
  
  # Perform overlap join
  result <- foverlaps(
    dt_bim[, .(chr, pos, snp_id, pos2 = pos)], 
    dt_co_ord, 
    by.x = c("chr", "pos", "pos2"),
    by.y = c("chr", "start", "stop"),
    nomatch = 0
  )
  return(result)
}


The equivalent for loop code (time consuming) will be:

ForLoop <- function(bim, ChIP){
  # Initialize an empty vector to store matched SNP IDs
  snps_in_regions <- c()

  for (i in 1:nrow(ChIP)) {
    snps <- bim$snp_id[
      bim$chr == ChIP$chr[i] &
      bim$pos >= ChIP$start[i] &
      bim$pos <= ChIP$stop[i]
    ]
    snps_in_regions <- c(snps_in_regions, snps)
    return(snps_in_regions)
}


  
}
