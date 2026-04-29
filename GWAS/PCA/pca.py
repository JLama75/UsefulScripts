# Author: Yuyang Luo, PhD; Ayellet Segre, PhD\
# Massacusetts Eye and Ear, Havard Medical School
# 2023

import re
import pdb
import time
import typing

import numpy as np
import pandas as pd
import seaborn as sns
import plotly.express as px
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from tqdm import tqdm
from typing import Tuple
from pathlib import Path
from random import random
from shutil import copyfile
from pandas import DataFrame, Series
from argparse import ArgumentParser, Namespace
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_score
from seaborn import heatmap, scatterplot, clustermap, lineplot
from subprocess import run, CompletedProcess, CalledProcessError

def get_args() -> Namespace:
    """
    Parse command line arguments using argparse.
    Returns the parsed argparse.Namespace object.
    """
    parser = ArgumentParser()
    parser.add_argument('--study_in', type=str, required=True,
                        help='Study subjects to for PCA calculation')
    parser.add_argument('--format_in', type=str, required=True,
                        help='Format of study input')
    parser.add_argument('--underscore_sample_id', action='store_true',
                        help='The delimiter of samples ID is underscore or not')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose operation')
    parser.add_argument('--min_maf', type=float, default=0.05,
                        help='Minimum MAF for defining common variants. '
                        'Default: 0.05')
    parser.add_argument('--max_missingness', type=float, default=0.05,
                        help='Maximum missingness of common variants. '
                        'Default: 0.05')
    parser.add_argument('--max_r2', type=float, default=0.1,
                        help='Maximum R^2 for LD pruning via Plink')
    parser.add_argument('--window_size', type=int, default=200,
                        help='Window size for LD pruning via Plink')
    parser.add_argument('--step_size', type=int, default=100,
                        help='Step size for LD pruning via Plink')
    parser.add_argument('--output_file', type=str, default='study_pca',
                        help='Output name for study')
    parser.add_argument('--plink1_path', type=str, default='plink',
                        help='Path to plink. Needed if plink is not '
                             'in $PATH')    
    parser.add_argument('--plink2_path', type=str, default='plink2',
                        help='Path to plink. Needed if plink is not in $PATH')
    parser.add_argument('--file_inferred_ancestry', type=str, default='',
                        help='The SHARDS output, including the column of "IID" and "predicted_ancestry" ')
    parser.add_argument('--flag_plot', type=bool, default=False,
                        help='PCA plots') 
    parser.add_argument('--number_of_PCs', type=int, default=20,
                        help='The number os PCs to save in the eigenvec file. All PCs will be calculated but this number determines how many PCs will be saved in the file of eigenvec.')

    args = parser.parse_args()
    return args

def validate_args(args: Namespace) -> bool:
    """
    Given the parsed arguments, validate them.
    """
    files_which_need_to_exist = [args.study_in, args.format_in]
    return all(
        [Path(f).is_file() for f in files_which_need_to_exist]
        )


def die(epitaph: str,
        exit_code: int = 1) -> None:
    """
    Print a message and exit.
    """
    print(epitaph)
    exit(exit_code)

def read_pcs_and_populations(pcs_file: str, inferred_ancestry_file: str) -> DataFrame:
    """
    Read in a plink.eigenvec file from PCA.
    Simultaneously, read in the populatins file from
    the reference panel.
    Return a DataFrame of PCs for the merged study and reference
    panel samples with population labels on the reference samples.
    """
    pcs_df = pd.read_csv(pcs_file, sep="\t").drop('#FID', axis=1)
    #ancestry_df = pd.read_csv(inferred_ancestry_file, sep=",")
    ancestry_df = pd.read_csv(inferred_ancestry_file, sep="\t")
    merged = pcs_df.merge(ancestry_df, on='IID', how='left')
    return merged

def set_var_ids(bim_file: str,
                chrom_pos_id: bool = False,
                bak_prefix: str = '.bak') -> None:
    """
    Reset the variant IDs in a bim file
    to chr:pos:ref:alt in place.
    If chrom_pos_id is set to True,
    set the IDs just to chr:pos in place.
    Keeps a backup of the original bim file.
    """
    ##Chromosome, Variant identifier, Position in morgans or centimorgans (safe to use dummy value of '0')
    ##Base-pair coordinate, Allele 1, Allele 2.
    backup_bim_file = bim_file + bak_prefix
    copyfile(bim_file, backup_bim_file)
    with open(backup_bim_file) as bim_in, open(bim_file, 'w') as bim_out:
        for line in bim_in:
            chrom, snp_id, cM, pos, ref, alt = line.rstrip().split('\t')
            new_id = ':'.join(
                ['chr' + chrom, pos] if chrom_pos_id
                else ['chr' + chrom, pos, ref, alt])  
            print('\t'.join([chrom, new_id, cM, pos, ref, alt]), file=bim_out)
    Path(backup_bim_file).unlink()


def bim_to_var_ids(bim_file: str,
                   var_ids_file: str) -> None:
    """
    Extract the variant IDs (second column)
    from a Plink BIM file and save it to a single column
    text file. Useful to generate a file that Plink can use
    to subset variants.
    """
    with open(bim_file) as inny, open(var_ids_file, 'w') as outy:
        for line in inny:
            fields = line.rstrip().split()
            print(f"{fields[1]}", file=outy)

def multi_allelic_ids_excluded(bim_file: str,
                      var_ids_file: str) -> None:
    """
    Extract the variant IDs (second column)
    from a Plink BIM file and save it to a single column
    text file. Useful to generate a file that Plink can use
    to subset variants.
    """
    with open(bim_file) as inny, open(var_ids_file, 'w') as outy:
        dicts = dict()
        for line in inny:
            fields = line.rstrip().split('\t')
            var_chr_pos = fields[1]
            if var_chr_pos not in dicts:
                dicts[var_chr_pos] = 1
            else:
                dicts[var_chr_pos] += 1
        #iterate dicts
        for key in dicts:
            if dicts[key] == 1:
                print(f"{key}", file=outy)

def ambiguous_site_ids_excluded(bim_file: str,
                                var_ids_file: str) -> None:
    """
    remove A/T or G/C; Extract the variant IDs (second column)
    from a Plink BIM file and save it to a single column
    text file. Useful to generate a file that Plink can use
    to subset variants.
    """
    with open(bim_file) as inny, open(var_ids_file, 'w') as outy:
        for line in inny:
            chrom, var_chr_pos, cM, pos, ref, alt = line.rstrip().split('\t')
            if ref == 'A' and alt == 'T' or ref == 'T' and alt == 'A' or ref == 'G' and alt == 'C' or ref == 'C' and alt == 'G':
                continue
            else:
                print(f"{var_chr_pos}", file=outy)

def run_plink(plink_command: str,
              verbose: bool = False,
              check: bool = True) -> CompletedProcess:
    """
    Given a string of a plink command,
    run it in a subprocess and return the CompletedProcess
    object.
    """
    if verbose:
        print(f'Calling plink:\n{plink_command}\n')
    try:
        completed_process = run(plink_command.split(),
                                capture_output=True,
                                check=check)
        return completed_process
    except CalledProcessError as called_process_error:
        die(f'Got an error from Plink:\n'
            f'{completed_process.stderr}',
            completed_process.returncode)

def order_chromosome(file_in: str,
		     format_in: str = 'bfile',
                     path_to_plink2: str = 'plink2',
		     verbose: bool = False,
		     ordered_file: str = 'ordered_study',
		     vcf_half_call: str = 'haploid',
		     underscore_sample_id: bool = False) -> CompletedProcess:
    """
    Given a VCF file, order the chromosome
    """
    if format_in == 'vcf':
       plink_command = (f'{path_to_plink2} --vcf {file_in} --vcf-half-call haploid ' + ("--double-id " if underscore_sample_id else "") +
                        f'--set-all-var-ids @_#_$r_$a --allow-extra-chr --new-id-max-allele-len 10000 --max-alleles 2 --make-pgen --sort-vars --out {ordered_file}')
       return run_plink(plink_command, verbose=verbose)
    elif format_in == 'bfile':
       plink_command = (f'{path_to_plink2} --bfile {file_in} ' + ("--double-id " if underscore_sample_id else "") +
                       f'--set-all-var-ids @_#_$r_$a --allow-extra-chr --new-id-max-allele-len 10000 --max-alleles 2 --make-pgen --sort-vars --out {ordered_file}')
       return run_plink(plink_command, verbose=verbose)


def filter_variants(plink_file: str,
                    trimmed_study_file: str,
                    min_maf: float = 0.05,
                    max_missingness: float = 0.01,
                    path_to_plink2: str = 'plink2',
                    verbose: bool = False,
                    underscore_sample_id: bool = False) -> CompletedProcess:
    """
    Given a VCF file, subset common, biallelic SNVs with
    minimal missingness to a new Plink fileset.
    """
    plink_command = (f'{path_to_plink2} --bfile {plink_file} --maf {min_maf} ' + ("--double-id " if underscore_sample_id else "") +
                     f'--geno {max_missingness} --snps-only --max-alleles 2 --allow-extra-chr '
                     f'--make-bed --out {trimmed_study_file}')
    return run_plink(plink_command, verbose=verbose)

def pfile_to_bfile(pfile_file: str,
                     path_to_plink2: str = 'plink2',
                     verbose: bool = False,
                     ordered_file: str = 'ordered_study',
                     underscore_sample_id: bool = False) -> CompletedProcess:
    """
    Given a VCF file, order the chromosome
    """
    plink_command = (f'{path_to_plink2} --pfile {pfile_file} ' + ("--double-id " if underscore_sample_id else "") +
                     f'--make-bed --out {ordered_file}')
    return run_plink(plink_command, verbose=verbose)


def ld_prune_variants(plink_file: str,
                      pruned_file: str,
                      r2: float = 0.2,
                      window_size: int = 200,
                      step_size: int = 100,
                      path_to_plink2: str = 'plink2',
                      verbose: bool = False,
                      underscore_sample_id: bool = False) -> CompletedProcess:
    """
    Call Plink to LD-prune variants in a study.
    Doesn't actually subset variants, just generates
    a prune.in file for later use.
    """
    plink_call = (f'{path_to_plink2} --bfile {plink_file} ' + ("--double-id " if underscore_sample_id else "") +
                  f'--indep-pairwise {window_size} {step_size} {r2} --allow-extra-chr --rm-dup '
                  f'--out {pruned_file}')
    return run_plink(plink_call, verbose=verbose)

def vcf_to_plink(vcf_file: str,
                 plink_file: str,
                 path_to_plink2: str = 'plink',
                 verbose: bool = False,
                 underscore_sample_id: bool = False) -> CompletedProcess:
    """
    Call Plink to convert a VCF file to a Plink file.
    """
    plink_call = (f'{path_to_plink2} --vcf {vcf_file} ' + ("--double-id " if underscore_sample_id else "") +
                  f'--max-alleles 2 --make-bed --out {plink_file}')
    return run_plink(plink_call, verbose=verbose)

def subset_variants(file_to_subset: str,
                    subsetted_file: str,
                    variants_to_keep: str,
                    path_to_plink2: str = 'plink2',
                    verbose: bool = False,
                    underscore_sample_id: bool = False) -> CompletedProcess:
    """
    Call Plink to subset variants in a Plink file.
    """
    plink_call = (f'{path_to_plink2} --bfile {file_to_subset} '
                  f'--extract {variants_to_keep} --snps-only ' + ("--double-id " if underscore_sample_id else "") +
                  f'--max-alleles 2 --allow-extra-chr --make-bed --out {subsetted_file}')
    return run_plink(plink_call, verbose=verbose)

def flip_variants(plink_file: str,
                  to_flip: str,
                  path_to_plink1: str = 'plink',
                  verbose: bool = False,
                  underscore_sample_id: bool = False) -> CompletedProcess:
    """
    Call Plink to flip variants specified in a text file.
    Variants are flipped in place (mutates the input file).
    """
    plink_call = (f'{path_to_plink1} --bfile {plink_file} ' + ("--double-id " if underscore_sample_id else "") + 
                  f'--flip {to_flip} --make-bed --out '
                  f'{plink_file}')
    return run_plink(plink_call, verbose=verbose)

def exclude_variants(plink_file: str,
                     to_exclude: str,
                     path_to_plink2: str = 'plink2',
                     verbose: bool = False,
                     underscore_sample_id: bool = False) -> CompletedProcess:
    """
    Call plink to exclude variants specified in a text file.
    Variants are excluded in place (mutates the input file).
    """
    plink_call = (f'{path_to_plink2} --bfile {plink_file} ' + ("--double-id " if underscore_sample_id else "") +
                  f'--exclude {to_exclude} --make-bed --out '
                  f'{plink_file}')
    return run_plink(plink_call, verbose=verbose)

def run_pca(plink_file: str,
            num_of_samples: int,
            verbose: bool = False,
            file_out: str = 'study_pca',
            path_to_plink2: str = 'plink2') -> CompletedProcess:
    """
    Call Plink to compute PCA.
    """

    """
    if num_of_samples >= 8000:
       num_of_samples = 8000
    if num_of_samples > 5000:
       plink_call = (f'{path_to_plink2} --bfile {plink_file} --pca 100 approx --out {file_out}')
       return run_plink(plink_call, verbose=verbose)
    else:
       plink_call = (f'{path_to_plink2} --bfile {plink_file} --pca {num_of_samples} --out {file_out}')
       return run_plink(plink_call, verbose=verbose)
    """
    if num_of_samples >= 8000:
       num_of_samples = 8000
    plink_call = (f'{path_to_plink2} --bfile {plink_file} --pca {num_of_samples} --out {file_out}')
    return run_plink(plink_call, verbose=verbose)

def main() -> None:
    args = get_args()
    print(args.flag_plot)
   # Validate arguments.
   # if not validate_args(args):
   #     die('Invalid arguments: one or more required files '
   #         'do not exist.', 255)

    #order chromosome
    ordered_file = 'ordered_study'
    order_chromosome_completed_process = order_chromosome(file_in=args.study_in,
	                                                  format_in=args.format_in,
                                                          path_to_plink2=args.plink2_path,
                                                          verbose=args.verbose,
	                                                  ordered_file=ordered_file,
	                                                  underscore_sample_id=args.underscore_sample_id)

    pfile_to_bfile_completed_process = pfile_to_bfile(ordered_file,
                                                      path_to_plink2=args.plink2_path,
                                                      verbose=args.verbose,
                                                      ordered_file=ordered_file,
                                                      underscore_sample_id=args.underscore_sample_id)

    # Subset common, unlinked variants from study:
    ##This step will also conovert vcf format to bim format
    trimmed_study_file = 'trimmed_study'
    filter_variants_completed_process = filter_variants(ordered_file,
                                                        trimmed_study_file,
                                                        min_maf=args.min_maf,
                                                        max_missingness=args.max_missingness,
                                                        path_to_plink2=args.plink2_path,
                                                        verbose=args.verbose,
                                                        underscore_sample_id=args.underscore_sample_id)

    #remove ambiguous
    study_remove_ambiguous_file = 'study_remove_ambiguous'
    ambiguous_site_ids_excluded(bim_file='trimmed_study.bim',
                                var_ids_file='ambiguous.study.in')

    remove_ambiguous_study_completed_process = subset_variants(file_to_subset=trimmed_study_file,
                                                               subsetted_file=study_remove_ambiguous_file,
                                                               variants_to_keep='ambiguous.study.in',
                                                               path_to_plink2=args.plink2_path,
                                                               verbose=args.verbose)

    # Recode the variant IDs to "{chrom}:{pos}" for merging.
    set_var_ids('study_remove_ambiguous.bim', chrom_pos_id=True)
    
    #add one new added function for subsetted_study
    #remove multiallelic
    study_remove_multiallelic_file = 'study_remove_multiallelic'
    multi_allelic_ids_excluded(bim_file='study_remove_ambiguous.bim',
                               var_ids_file='multiallelic.study.in')

    remove_multiallelic_study_completed_process = subset_variants(file_to_subset=study_remove_ambiguous_file,
                                                                  subsetted_file=study_remove_multiallelic_file,
                                                                  variants_to_keep='multiallelic.study.in',
                                                                  path_to_plink2=args.plink2_path,
                                                                  verbose=args.verbose)

    # Check the number of sample of study_remove_multiallelic_file.fam
    with open(r"study_remove_multiallelic.fam","r") as fp:
        num_sample_study=len(fp.readlines())
        
    if num_sample_study<50:
        print('sample size too small to do ld pruning!')
    
    # LD prune the study's variants:
    ld_pruned_study_file = 'ld_pruned_study'
    ld_prune_variants_completed_process = ld_prune_variants(plink_file=study_remove_multiallelic_file,
                                                            pruned_file=ld_pruned_study_file,
                                                            r2=args.max_r2,
                                                            window_size=args.window_size,
                                                            path_to_plink2=args.plink2_path,
                                                            step_size=args.step_size,
                                                            verbose=args.verbose,
                                                            underscore_sample_id=args.underscore_sample_id)

    # Not all of the variants in the subsetted reference file
    # will also exist in the LD-pruned study, so now
    # reciprocally subset the LD-pruned study:
    subsetted_study_file = 'subsetted_study'
   
    subset_study_completed_process = subset_variants(file_to_subset=study_remove_ambiguous_file,
                                                     subsetted_file=subsetted_study_file,
                                                     variants_to_keep='ld_pruned_study.prune.in',
                                                     path_to_plink2=args.plink2_path,
                                                     verbose=args.verbose)
    
    #Determine the number os samples
    pca_completed_process = run_pca(plink_file=subsetted_study_file,
                                    num_of_samples=num_sample_study,
                                    verbose=args.verbose,
                                    file_out=args.output_file,
                                    path_to_plink2=args.plink2_path)

    #The previous steps will save all PCs into the file of eigenvec. But users can specify how many PCs they want.
    df_eigenvec=pd.read_csv(args.output_file+'.eigenvec', sep='\t')
    df_eigenvec_subset=df_eigenvec.iloc[:, :args.number_of_PCs+2]
    df_eigenvec_subset.to_csv(args.output_file+'.eigenvec', sep='\t', index=None)

    # % variance explained by PCA components
    eigenval=pd.read_csv(args.output_file+'.eigenval', sep='\s+', header=None)
    eigenval=eigenval.iloc[:,0].to_list()
    eigenval_fractional=eigenval/np.sum(eigenval)

    plt.figure(figsize=(6, 6))
    plt.plot((eigenval_fractional)[:20] * 100, '.-')
    plt.xticks(range(0, 20), list(range(1, 21)), fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('% variance explained', fontsize=16)
    plt.xlabel('number of PC', fontsize=16)
    plt.title('Percent variance explained', fontsize=18)
    plt.tick_params(axis='both', labelsize=12)
    plt.tight_layout()
    plt.savefig('Percent_variance_explained.png', dpi=300)

    if args.flag_plot == True:
        print('Plotting!')
        #Read the file from SHARDS if available
        file_inferred_ancestry = args.file_inferred_ancestry
        pcs_df = read_pcs_and_populations(pcs_file=args.output_file+'.eigenvec',
                                          inferred_ancestry_file=file_inferred_ancestry)
    
#        # % variance explained by PCA components
#        eigenval=pd.read_csv(args.output_file+'.eigenval', sep='\s+', header=None)
#        eigenval=eigenval.iloc[:,0].to_list()
#        eigenval_fractional=eigenval/np.sum(eigenval)

#        plt.figure(figsize=(6,6))
#        plt.plot((eigenval_fractional)[:20]*100, '.-')
#        plt.xticks(range(0,20), list(range(1,21)))
#        plt.ylabel('% variance explained')
#        plt.xlabel('number of PC')
#        plt.savefig('Percent_variance_explained.png')

        #PCA plots
        sns.set_context('poster', font_scale=2.0)
        fig, axes = plt.subplots(ncols=1, nrows=3, figsize=[14,36])

        scatterplot(x=pcs_df['PC1'],
                    y=pcs_df['PC2'],
                    hue=pcs_df['predicted_ancestry'],
                    linewidth=0.1,
                    ax=axes[0]
                    )
        axes[0].set_xlabel('PC1', fontsize=28)
        axes[0].set_ylabel('PC2', fontsize=28)
        axes[0].tick_params(axis='both', labelsize=24)
        axes[0].legend(fontsize=20, title_fontsize=22)
        scatterplot(x=pcs_df['PC1'],
                    y=pcs_df['PC3'],
                    hue=pcs_df['predicted_ancestry'],
                    linewidth=0.1,
                    ax=axes[1]
                   )
        axes[1].set_xlabel('PC1', fontsize=28)
        axes[1].set_ylabel('PC3', fontsize=28)
        axes[1].tick_params(axis='both', labelsize=24)
        axes[1].legend(fontsize=20, title_fontsize=22)
        scatterplot(x=pcs_df['PC2'],
                    y=pcs_df['PC3'],
                    hue=pcs_df['predicted_ancestry'],
                    linewidth=0.1,
                    ax=axes[2]
                   )
        axes[2].set_xlabel('PC2', fontsize=28)
        axes[2].set_ylabel('PC3', fontsize=28)
        axes[2].tick_params(axis='both', labelsize=24)
        axes[2].legend(fontsize=20, title_fontsize=22)
        #plt.legend(loc='best', bbox_to_anchor=(1,-0.05))
        plt.savefig("PCplot.pdf", format="pdf", bbox_inches="tight")
        plt.savefig('PCplot.png')

if __name__ == "__main__":
    st=time.time()
    main()
    et=time.time()
    elapsed_time=et-st
    print('Execution time:', elapsed_time, 'seconds')
