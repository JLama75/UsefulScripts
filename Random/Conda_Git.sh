#Creating conda env
#conda remove --name ENV_NAME --all
#login terminal: ssh jl2251@172.27.81.202
Hi, Jyoti. Yes, I accepted a position at UMass Chan. My phone number is (347)363-9183. My personal email is yuyang.luo0411@gmail.com
source ~/.bash_profile

###ERIS2
bsub -Is -n 2 -t 4:0:0 -q interactive  /bin/bash
module load miniconda3/latest

#######
/data/Segre_Lab/data/DRCR/analysis/colocalization/run_retina_013124
MAGMA: /data/Segre_Lab/data/DRCR/analysis/MAGMA

#Install miniconda3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh

#Run miniconda3 within your home directory. Add this directory to your path
bash miniconda.sh -b -p $HOME/miniconda  
export PATH="$HOME/miniconda/bin:$PATH"
https://jupyterhub2.partners.org/
#create an environment.yml file with all dependencies

#example:
name: xmitgcm
dependencies:
  - numpy
  - scipy
  - xarray
  - netcdf4
  - dask
  - jupyter
  - matplotlib
  - pip:
    - pytest
    - xmitgcm

conda env create --file environment.yml
conda activate xmitgcm

###############################################

#create Renvironment.yml file

name: R
channels:
  - conda-forge
dependencies:
  - r-base  #4.3.2
  - r-devtools
  - r-optparse
  - r-dplyr
  - r-ggplot2
  - r-gridextra
  
conda env create -f environment.yml
conda activate R
R
install.packages("SamplyzeR")

#source $HOME/miniconda3/bin/activate/xmitgcm

##############################################################
conda create --name myPython python=3.9
conda activate myPython
conda install numpy pandas seaborn plotly matplotlib scikit-learn
conda install scipy statsmodels xlsxwriter
conda install pysnptools
pip install tqdm


###########################################################

conda create --name myR r-base=4.3.2
conda activate myR
conda install -c conda-forge r-devtools r-dplyr r-readr r-tidyr r-optparse -y


######
conda create --name myConda
conda activate myConda
conda install bioconda::p7zip
##################################################
conda create --name myPython2 python=3.9
conda activate myPython2
conda install numpy pandas seaborn plotly matplotlib scikit-learn scipy
conda install statsmodels xlsxwriter json
conda install bioconda::regenie

 
###########################################

conda create --name myPython3 python=3.8.13
conda activate myPython3
conda install numpy=1.24.4 pandas=2.0.3 seaborn=0.13.2 plotly matplotlib scikit-learn scipy


conda create --name myPython3.8 python=3.8 seaborn numpy pandas plotly matplotlib scikit-learn scipy statsmodels -c conda-forge -c bioconda
conda install xlsxwriter pysnptools

conda install seaborn=0.13.2
conda install statsmodels xlsxwriter json
conda install bioconda::regenie
conda env remove --name myPython3


conda create --name myPython3.7 python=3.7 seaborn=0.13.2 numpy=1.24.4 pandas=2.0.3  plotly matplotlib scikit-learn scipy -c conda-forge
################################################
#For Dask

mamba create -n pythonDask python=3.10 dask pandas=2.2
conda activate pythonDask
#Confirm it works
python -c "import dask.dataframe as dd; print(dd.__version__)" 
#You should no longer get the ImportError
##############################################
#Installing Hail on OGI cluster
#Installation Steps for HAIL:

conda create --name Hail python=3.8
conda activate Hail
module load GCC/5.4.0-2.26
conda install -c conda-forge openjdk
conda install -c anaconda lz4
conda install -c conda-forge blas
conda install -c conda-forge lapack
pip install pyspark==3.3.0
########
git clone https://github.com/hail-is/hail.git
cd hail/hail
make install-on-cluster HAIL_COMPILE_NATIVES=1 SCALA_VERSION=2.12.15 SPARK_VERSION=3.3.0 (make clean)
#conda install -c conda-forge gxx
#export PATH=$Hail/bin:$PATH
#which g++

#*Please use the corresponding SPARK and SCALA version

#For variant qc that requires pyvcf:
conda install -c bioconda pyvcf
###
done
#conda remove --Hail --all

conda install -c bioconda -c conda-forge bioperl
conda create -n vep_env -c conda-forge -c bioconda ensembl-vep=110.1
###########################################################################################
#To access Jupyter notebook on OGI cluster via SSH tunnel:
#1.    On your local computer, run ssh -L 2001:localhost:2001 username@ogi-mbc-login.meei.harvard.edu
#2.    Create a conda environment and install Jupyter notebook
#3.    Under the activated conda environment, run jupyter notebook --no-browser --port=2001
#4.    Copy the pop-up URL to the browser.

ssh -L 2001:localhost:2001 jl2251@ogi-mbc-login.meei.harvard.edu
ssh -L 2007:localhost:2007 jl2251@ogi-mbc-login.meei.harvard.edu
ssh -L 2036:localhost:2036 jl2251@ogi-mbc-login.meei.harvard.edu
ssh -L 2034:localhost:2034 jl2251@ogi-mbc-login.meei.harvard.edu

ssh -L 2037:localhost:2037 jl2251@ogi-mbc-login.meei.harvard.edu
#Installing jupyter notebook
conda create --name JupyterNotebook python=3.8
conda activate JupyterNotebook
conda install -c conda-forge notebook

jupyter notebook --no-browser --port=2037
jupyter notebook --no-browser --port=2034
module load jupyter_notebook/4.4.0

ssh -L 2019:localhost:2019 jl2251@ogi-mbc-login.meei.harvard.edu
jupyter notebook --no-browser --port=2019
#control-C: stop the server and shut down all kernels.
jupyter notebook --no-browser --port=2007
###########################################

find . -name Makefile
vim ~./.bash_profile
source $HOME/miniconda/bin/activate 

##########################################################
#Installation steps for Hail updated from Liyin and Sudeep 09/16/24
#ERIS
module load miniconda3/latest
conda create -n hailenv.092724 
conda activate hailenv.092724
conda install bioconda::hail=0.2.33-0
conda install bioconda::hail=0.2.61
#For variant qc that requires pyvcf:
conda install -c bioconda pyvcf

############################################################

#Installiation steps by meghana 02/05/25
#ERIS
conda create -n hail python=3.9 #HAIL VERSION==0.2.133
conda activate hail 
#java version = 11 
#conda install openjdk
python3.9 -m pip install hail

pip install numpy pandas seaborn plotly matplotlib tqdm scikit-learn


######
conda create --n My_python3.6 python=3.6 -y

##############################################
#Installation Steps for HAIL from Yan: sept 27 2024
conda create --name Hail python=3.9
conda activate Hail
#conda install anaconda::java-1.8.0-openjdk-devel-cos7-s390x
#conda install anaconda::openjdk #Should be Java 8 or 11. Java 8 is already in the cluster #Java -version
module load GCC/5.4.0-2.26 
module load LLVM/3.8.1-foss-2016b 
conda install -c conda-forge blas
conda install -c conda-forge lapack
python3.9 -m pip install hail
conda install -c conda-forge lz4
pip install pyspark==3.3.2 cond

Proof its loading 
python
Python 3.9.18 | packaged by conda-forge | (main, Dec 23 2023, 16:33:10) 
[GCC 12.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import hail as hl


#For variant qc that requires pyvcf:
conda install bioconda::pyvcf==0.6.8
#####################################################################
#From Sudeep sucessfully installed
mamba create -n Hail3 python=3.9
mamba activate Hail3
mamba install anaconda::java-1.8.0-openjdk-devel-cos7-s390x
mamba install anaconda::openjdk 
mamba install GCC
mamba install LLVM
mamba install -c conda-forge blas lapack lz4
pip install pyspark==3.3.2
pip install hail

Proof its loading 
python
Python 3.9.18 | packaged by conda-forge | (main, Dec 23 2023, 16:33:10) 
[GCC 12.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import hail as hl
>>> hl.init(default_reference='GRCh38', log='./hail.log')
mamba install bioconda::pyvcf==0.6.8
#####################################################################

#Installation Steps for HAIL:

conda create --name Hail3 python=3.9
conda activate Hail3
#conda install anaconda::openjdk=8 #Already in the cluster #Java -version
conda install anaconda::java-1.8.0-openjdk-devel-cos7-s390x
module load GCC/5.4.0-2.26 #module av GCC
module load LLVM/3.8.1-foss-2016b #module av LLVM
conda install -c conda-forge blas
conda install -c conda-forge lapack
python3.9 -m pip install hail
conda install -c conda-forge lz4
pip install pyspark==3.3.2 #3.5.0 #not supported
#pip uninstall pyspark
#pip cache purge
pyspark --version #3.5.0 
pyspark
quit()

#source ~/.bash_profile #with all the paths to conda environment and pyspark# NOT NEEDED!

#git clone https://github.com/hail-is/hail #NOT NEEDED! ALREADY DONE!
cd hail/hail
make install-on-cluster HAIL_COMPILE_NATIVES=1 SCALA_VERSION=2.12.15 SPARK_VERSION=3.3.2 #make clean
#2.12.18
#make install-on-cluster HAIL_COMPILE_NATIVES=1 SCALA_VERSION=2.12.18 SPARK_VERSION=3.5.0 #make clean
conda remove -n vep_env --all
#conda remove -n hail_env lz4
#conda remove -n hail_env openjdk
#conda remove -n hail_env python
#conda install python=3.9
#conda install -c conda-forge gxx
#export PATH=$Hail/bin:$PATH
#which g++
#modify bokeh == 3.2.2
#find . -name Makefile
conda remove -n vep_env --all

source $HOME/miniconda/bin/activate Hail
#*Please use the corresponding SPARK and SCALA version

#For variant qc that requires pyvcf:
conda install -c bioconda pyvcf
###
done


##############################################
#Installation Steps for HAIL:

conda create --name Hail3 python=3.9
conda activate Hail3
#conda install anaconda::java-1.8.0-openjdk-devel-cos7-s390x
#conda install anaconda::openjdk #Should be Java 8 or 11. Java 8 is already in the cluster #Java -version
module load GCC/5.4.0-2.26 
module load LLVM/3.8.1-foss-2016b 
conda install -c conda-forge blas
conda install -c conda-forge lapack
pip install pyspark==3.3.2 
pip install hail #version 0.2.132-678e1f52b999

#pip install lz4==4.0.0
#conda install -c conda-forge lz4
#conda install anaconda::lz4==4.3.2


pyspark --version #3.3.2
pyspark
quit()

cd hail/hail
make install-on-cluster HAIL_COMPILE_NATIVES=1 SCALA_VERSION=2.12.15 SPARK_VERSION=3.3.2 #make clean
#*Please use the corresponding SPARK and SCALA version

##############
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh

source ~/miniconda3/bin/activate
conda --version

##############
conda create -n aspera-env
conda activate aspera-env
#conda install rpetit3::aspera-connect
conda install hcc::aspera-cli

#!/bin/bash

export PATH="/Applications/IBM Aspera Connect.app/Contents/Resources:$PATH"
export ASPERA_SCP_FILEPASS=xHct6V6cvkJ0FCmw

KEYFILE="/Users/jl2251/Documents/transferToNewLaptop/OGI/Steroid_IOP/dbGAP/aspera.openssh "
UPLOAD_URL="subasp@upload.ncbi.nlm.nih.gov:uploads/upload_requests/fZdogSBo/"

echo keyfile path: ${KEYFILE}
echo upload url: ${UPLOAD_URL}

for file in /Users/jl2251/Documents/transferToNewLaptop/OGI/Steroid_IOP/dbGAP/GENO/* ; do
  echo uploading ${file}
  ascp --file-crypt=encrypt -i "${KEYFILE}" "${file}" "${UPLOAD_URL}"
done




keyfile path: /Users/aspera.openssh
upload url: subasp@upload.ncbi.nlm.nih.gov:uploads/upload_requests/jdfdo/
uploading /Users/GENO/FAME.532.vcf.gz
./uploadFiles.sh: line 15: ascp: command not found

echo keyfile path: ${KEYFILE}
echo upload url: ${UPLOAD_URL}

for file in /Users/jl2251/Documents/transferToNewLaptop/OGI/Steroid_IOP/dbGAP/GENO/* ; do
  echo uploading ${file}
  ascp --file-crypt=encrypt -i "${KEYFILE}" "${file}" "${UPLOAD_URL}"
done

##########################################
gcloud auth login

gsutil -m cp -r 'gs://fc-2e94ced4-e58e-4804-8b4a-a407f668a7e1/Exome_VCF/Sobrin_MEEI_14_WES_VCF.vcf.bgz' Exome_VCF/Sobrin_MEEI_14_WES_VCF.vcf.bgz
gsutil -m cp -r 'gs://fc-2e94ced4-e58e-4804-8b4a-a407f668a7e1/Exome_VCF/Sobrin_MEEI_14_WES_VCF.vcf.bgz.tbi' Exome_VCF/Sobrin_MEEI_14_WES_VCF.vcf.bgz.tbi
gsutil -m cp -r 'gs://fc-4362e0a2-530d-48ca-93f5-42602409b7c4/PO-62806_07122023_MEEI_Exomes_terra_manifest.tsv' PO-62806_07122023_MEEI_Exomes_terra_manifest.tsv




#pyspark --version #3.3.2
#pyspark
#quit()

#git clone https://github.com/hail-is/hail.git
#cd hail/hail
#make install-on-cluster HAIL_COMPILE_NATIVES=1 SCALA_VERSION=2.12.15 SPARK_VERSION=3.3.2 #make clean
#*Please use the corresponding SPARK and SCALA version
#############################################
#Sudeep's code
#Hail Version 0.2.130

mamba create -n hailenv python=3.9
mamba activate hailenv
mamba install anaconda::java-1.8.0-openjdk-devel-cos7-s390x
mamba install anaconda::openjdk 
mamba install GCC
mamba install LLVM
mamba install -c conda-forge blas lapack lz4
pip install pyspark==3.3.2
pip install hail

Proof its loading 
python
Python 3.9.18 | packaged by conda-forge | (main, Dec 23 2023, 16:33:10) 
[GCC 12.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import hail as hl


mamba create -n hailenv2 python=3.9
mamba activate hailenv2
mamba install anaconda::java-1.8.0-openjdk-devel-cos7-s390x
mamba install anaconda::openjdk 
mamba install GCC
mamba install LLVM
mamba install -c conda-forge blas lapack lz4
pip install pyspark==3.3.2
pip install hail

Proof its loading 
python
Python 3.9.18 | packaged by conda-forge | (main, Dec 23 2023, 16:33:10) 
[GCC 12.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import hail as hl

############################################

conda create -n gsutil 
conda activate gsutil
conda install conda-forge::gsutil


mamba create --name bcftools bcftools
mamba install -c bioconda bcftools #1.9-9


##############################
#Git repo

git clone https://pat@github.com/JLama75/RepositoryName.git

git clone https://pat@github.com/JLama75/RepositoryName.git
#Fine grained PAT with Git for authentication

git clone https://pat@github.com/JLama75/RepositoryName.git
Username: your-github-username
Password: your-fine-grained-token #github_pat_number
git push origin main
#####################################
To pull the directory from privat github
First get your Personal access tokens (classic) from Settings --> Developer Settings --> Personal access tokens

git clone https://pat@github.com/JLama75/RepositoryName.git

git remote set-url origin https://pat@github.com/JLama75/RepositoryName.git
#############################
To push the directory
git add .
git commit -m "Updated project files"
git push origin main

######################
#Push new directory to github

#1.If the directory is not yet a Git repository, initialize it:
git init
#2.If you haven't added a remote repository yet, link your local directory to your GitHub repo:
git remote add origin https://github.com/your-username/your-repository.git
#3.Add your files to Git: Stage all files in the directory:
git add .
git commit -m "Initial commit of scripts"
git branch -M main

#Get your personal access token
git remote set-url origin https://pat@github.com/JLama75/RepositoryName.git

git push -u origin main


