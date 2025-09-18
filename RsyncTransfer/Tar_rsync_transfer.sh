#!/bin/bash

# On your Mac, in the folder above Exome
cd "/Volumes/One Touch/cloud_backups/Sobrin_MEEI_ExomesGSA/07122023_MEEI_Exomes/RP-2799/"

# Create a compressed tarball (1 TB → will take time)
tar -cvf Exome.tar Exome/

# Transfer in one go
rsync -ratlzvP Exome.tar ogi-mbc-login.meei.harvard.edu:/gpfs/archive1/SIOP/CRAMS/

# On the cluster, untar
#ssh ogi-mbc-login.meei.harvard.edu
#cd /gpfs/archive1/SIOP/CRAMS
#tar -xvf Exome.tar
