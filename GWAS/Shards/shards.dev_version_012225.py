# Author: Brian S. Cole, PhD; Yuyang Luo, PhD; Ayellet Segre, PhD\
# Massacusetts Eye and Ear, Havard Medical School
# 2020

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
    parser.add_argument('--filename_study', type=str, required=True,
                        help='Data from study subjects to infer ancestry, can be pfile, bfile, bgen or vcf, should include the extension of your file if vcf, or just filename if pfile, bfile or bgen. For vcf, the filename should be xxx.vcf; for pfile, bfile or bgen, would be just filename without suffix')
    parser.add_argument('--format_study', type=str, required=True,
                        help=' The format of data from study subjects, can be either pfile, bfile, bgen or vcf. ')
    parser.add_argument('--filename_reference', type=str, required=True,
                        help='filename from reference panel')
    parser.add_argument('--format_reference', type=str, required=True,
                        help='The format of data from reference panel, can be either pfile, bfile, bgen or vcf.')
    parser.add_argument('--reference_populations', type=str, required=True,
                        help='TSV of populations from reference panel. '
                             'Needs first column as sample ID, '
                             'second column as population')
    parser.add_argument('--output_file', type=str,
                        default='inferred_ancestry.tsv',
                        help='Output TSV of inferred ancestry for study')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose operation')
    parser.add_argument('--min_maf', type=float, default=0.05,
                        help='Minimum MAF for defining common variants. '
                        'Default: 0.05')
    parser.add_argument('--num_cross_validation', type=int, default=5,
                        help='The number of cross validation. '
                        'Default: 5')
    parser.add_argument('--max_missingness', type=float, default=0.01,
                        help='Maximum missingness of common variants. '
                        'Default: 0.01')
    parser.add_argument('--max_r2', type=float, default=0.2,
                        help='Maximum R^2 for LD pruning via Plink')
    parser.add_argument('--window_size', type=int, default=200,
                        help='Window size for LD pruning via Plink')
    parser.add_argument('--step_size', type=int, default=100,
                        help='Step size for LD pruning via Plink')
    parser.add_argument('--plink1_path', type=str, default='plink',
                        help='Path to plink. Needed if plink is not '
                             'in $PATH')
    parser.add_argument('--plink2_path', type=str, default='plink2',
                        help='Path to plink. Needed if plink is not '
                             'in $PATH')
    parser.add_argument('--genome_build', type=str, default='b37',
                        help='Genome build version')
    parser.add_argument('--delete_intermediate_files', action='store_true',
                        help='Delete intermediate files at end')
    parser.add_argument('--underscore_sample_id', action='store_true',
                        help='The delimiter of samples ID is underscore or not')
    parser.add_argument('--number_of_principal_components', type=int, default=20,
                        help='The maximum number of principal components for cross validation')
    parser.add_argument('--number_of_neighbors', type=int, default=50,
                        help='The maximum number of neighbors for cross validation')
    parser.add_argument('--vcf_half_call', type=str, default='haploid',
                        help='The current VCF standard does not specify how "0/." or "./1" '
                        'and similar GT values should be interpreted. Using "haploid" to treat '
                        'half-calls as haploid/homozygous; Using "missing" to treat as missing '
                        'Using "reference" to treat as reference. ' 
                        'Default: haploid')
    parser.add_argument('--study_covariates', type=str, help='TSV of study samples for cases/controls. '
                        'For the PCAplot purpose only, if there is a file indicating the cases and controls in the study, '
                        'the PCAplot for the study data will distinguish them. ')
    args = parser.parse_args()
    return args


def die(epitaph: str,
        exit_code: int = 1) -> None:
    """
    Print a message and exit.
    """
    print(epitaph)
    exit(exit_code)


def validate_args(args: Namespace) -> bool:
    """
    Given the parsed arguments, validate them.
    """
    files_which_need_to_exist = check_reference_data(args)
    if args.format_study == 'vcf':
        file_study_in=args.filename_study
        files_which_need_to_exist.append(file_study_in)

    elif args.format_study == 'bfile':
        file_study_in1=args.filename_study+'.bim'
        file_study_in2=args.filename_study+'.bed'
        file_study_in3=args.filename_study+'.fam'
        files_which_need_to_exist.append(file_study_in1)
        files_which_need_to_exist.append(file_study_in2)
        files_which_need_to_exist.append(file_study_in3)

    elif args.format_study == 'pfile':
        file_study_in1=args.filename_study+'.pvar'
        file_study_in2=args.filename_study+'.psam'
        file_study_in3=args.filename_study+'.pgen'
        files_which_need_to_exist.append(file_study_in1)
        files_which_need_to_exist.append(file_study_in2)
        files_which_need_to_exist.append(file_study_in3)

    elif args.format_study == 'bgen':#emphasize this in README file
        file_study_in1=args.filename_study+'.bgen'
        file_study_in2=args.filename_study+'.sample'
        files_which_need_to_exist.append(file_study_in1)
        files_which_need_to_exist.append(file_study_in2)
        
    # Check each file and print which ones are missing
    all_valid = True
    for f in files_which_need_to_exist:
        if not Path(f).is_file():
            print(f"  Missing file: '{f}'")
            all_valid = False
 
    #return all([Path(f).is_file() for f in files_which_need_to_exist])
    return all_valid

def check_reference_data(args: Namespace) -> list:
    reference_need_to_exist = []
    if args.format_reference == 'vcf':
        reference_need_to_exist.append(args.filename_reference)
    elif args.format_reference == 'bfile':
        reference_need_to_exist.append(args.filename_reference+'.bim')
        reference_need_to_exist.append(args.filename_reference+'.bed')
        reference_need_to_exist.append(args.filename_reference+'.fam')
    elif args.format_reference == 'pfile':
        reference_need_to_exist.append(args.filename_reference+'.pvar')
        reference_need_to_exist.append(args.filename_reference+'.pgen')
        reference_need_to_exist.append(args.filename_reference+'.psam')
    elif args.format_reference == 'bgen':
        reference_need_to_exist.append(args.filename_reference+'.bgen')
        reference_need_to_exist.append(args.filename_reference+'.sample')
    reference_need_to_exist.append(args.reference_populations)
    return reference_need_to_exist

def check_bgen_argument(**kwargs):
    if 'bgen_file' not in kwargs or 'sample_file' not in kwargs:
        raise ValueError("In 'bgen' mode, 'bgen_file' and 'sample_file' are required.")
    print(f"bgen mode with .bgen and .sample files")

def generate_sex_file_vcf_study(vcf_file: str,
                                path_to_plink1: str = 'plink',
                                underscore_sample_id: bool = False,
                                verbose: bool = False):
    plink_command = (f'{path_to_plink1} --vcf {vcf_file} ' + ("--double-id " if underscore_sample_id else "") +
                     f'--make-just-fam --out study_sex')
    return run_plink(plink_command, verbose=verbose)

def generate_sex_file_vcf_reference(vcf_file: str,
                                path_to_plink1: str = 'plink',
                                underscore_sample_id: bool = False,
                                verbose: bool = False):
    plink_command = (f'{path_to_plink1} --vcf {vcf_file} ' + ("--double-id " if underscore_sample_id else "") +
                     f'--make-just-fam --out reference_sex')
    return run_plink(plink_command, verbose=verbose)

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

def order_chromosome(study_in: str,
                     format_in: str,
                     path_to_plink2: str = 'plink2',
                     path_to_plink1: str = 'plink',
                     verbose: bool = False,
                     ordered_file: str = 'ordered_study',
                     vcf_half_call: str = 'haploid',
                     genome_build: str = 'b37',
                     underscore_sample_id: bool = False) -> CompletedProcess:
    """
    Given a VCF file, order the chromosome
    """
    if format_in == 'bfile':
        plink_command = (f'{path_to_plink2} --bfile {study_in} ' + ("--double-id " if underscore_sample_id else "") +
                         f'--set-all-var-ids @_#_$r_$a --allow-extra-chr --new-id-max-allele-len 10000 --max-alleles 2 --make-pgen --sort-vars --out {ordered_file}')
    elif format_in == 'vcf':
        plink_command = (f'{path_to_plink1} --vcf {study_in} --make-bed --out temp_file_bfile '+ ("--double-id " if underscore_sample_id else ""))
        run_plink(plink_command, verbose=verbose)
        plink_command = (f'{path_to_plink2} --bfile temp_file_bfile --make-pgen --out temp_file '+ ("--double-id " if underscore_sample_id else ""))
        run_plink(plink_command, verbose=verbose)
        plink_command = (f'{path_to_plink2} --pfile temp_file ' + ("--double-id " if underscore_sample_id else "") +
                     f'--set-all-var-ids @_#_$r_$a --allow-extra-chr --new-id-max-allele-len 10000 --max-alleles 2 --split-par {genome_build} --make-pgen --sort-vars --out {ordered_file}')
    elif format_in == 'pfile':
        plink_command = (f'{path_to_plink2} --pfile {study_in} ' + ("--double-id " if underscore_sample_id else "") +
                         f'--set-all-var-ids @_#_$r_$a --allow-extra-chr --new-id-max-allele-len 10000 --max-alleles 2 --make-pgen --sort-vars --out {ordered_file}')
    return run_plink(plink_command, verbose=verbose)

def order_chromosome_bgen(bgen_file: str,
                          sample_file: str,
                          path_to_plink2: str = 'plink2',
                          verbose: bool = False,
                          ordered_file: str = 'ordered_study',
                          underscore_sample_id: bool = False) -> CompletedProcess:
    """
    Given a VCF file, order the chromosome
    """
    plink_command = (f'{path_to_plink2} --bgen {bgen_file} --sample {sample_file}' + ("--double-id " if underscore_sample_id else "") +
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
                   updated_file: str,
                   path_to_plink2: str = 'plink2',
                   verbose: bool = False) -> CompletedProcess:
    """
    Given a VCF file, order the chromosome
    """
    plink_command = (f'{path_to_plink2} --pfile {pfile_file} --make-bed --out {updated_file}')
    return run_plink(plink_command, verbose=verbose)

def bfile_to_bfile(bfile_file: str,
                   updated_file: str,
                   path_to_plink2: str = 'plink2',
                   verbose: bool = False,
                   underscore_sample_id: bool = False) -> CompletedProcess:
    """
    """
    plink_command = (f'{path_to_plink2} --bfile {bfile_file} ' + ("--double-id " if underscore_sample_id else "") +
                     f'--make-bed --out {updated_file}')
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

def vcf_to_bfile(vcf_file: str,
                 plink_file: str,
                 path_to_plink1: str = 'plink',
                 path_to_plink2: str = 'plink2',
                 verbose: bool = False,
                 genome_build: str = 'b38',      
                 underscore_sample_id: bool = False) -> CompletedProcess:
    """
    Call Plink to convert a VCF file to a Plink file.
    """
    generate_sex_file_vcf_reference(vcf_file, path_to_plink1, underscore_sample_id, verbose)
    plink_call = (f'{path_to_plink2} --vcf {vcf_file} --fam reference_sex.fam ' + ("--double-id " if underscore_sample_id else "") +
                  f'--max-alleles 2 --split-par {genome_build} --make-bed --out {plink_file}')
    return run_plink(plink_call, verbose=verbose)

def bgen_to_bfile(bgen_file: str,
                  sample_file: str, 
                  plink_file: str,
                  path_to_plink2: str = 'plink',
                  verbose: bool = False,
                  genome_build: str = 'b38',
                  underscore_sample_id: bool = False) -> CompletedProcess:
    """
    Call Plink to convert a bgen file to a Plink file.
    """
    plink_call = (f'{path_to_plink2} --bgen {bgen_file} --sample {sample_file} ' + ("--double-id " if underscore_sample_id else "") +
                  f'--max-alleles 2 --split-par {genome_build} --make-bed --out {plink_file}')
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

def merge_plink(first_file: str,
                second_file: str,
                merged_file: str,
                path_to_plink1: str = 'plink',
                path_to_plink2: str = 'plink2',
                verbose: bool = False,
                underscore_sample_id: bool = False) -> CompletedProcess:
    """
    Call Plink to merge two Plink files.
    Don't use the 'check' option because this will fail
    with variants that need flipping or exclusion.
    """
    plink_call = (f'{path_to_plink1} --bfile {first_file} ' + ("--double-id " if underscore_sample_id else "") +
                  f'--bmerge {second_file} --allow-extra-chr --make-bed --out '
                  f'{merged_file}')
    return run_plink(plink_call, verbose=verbose,
                     check=False)

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
                  f'--flip {to_flip} --allow-extra-chr --make-bed --out '
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

def five_pass_merge(first_file: str,
                    second_file: str,
                    merged_file: str,
                    path_to_plink2: str = 'plink2',
                    path_to_plink1: str = 'plink',
                    verbose: bool = False) -> CompletedProcess:
    """
    Merge two Plink files in (up to) five passes.
    1) Try a merge. If it worked, return the merged file.
    2) Flip alleles for variants in the second_file that can't be merged.
    3) Repeat step 1.
    4) Exclude variants in the second_file that can't be merged.
    5) Repeat step 1.
    This step currently runs only with plink 1.9, which merges based on minor/major alleles.
    plink2.0 merges based on REF/ALT in pfile, but currently is under development, and is not supported for bfiles.
    """
    merge_completed_process = merge_plink(first_file=first_file,
                                          second_file=second_file,
                                          merged_file=merged_file,
                                          path_to_plink1=path_to_plink1,
                                          path_to_plink2=path_to_plink2,
                                          verbose=verbose)
    if merge_completed_process.returncode == 0:
        return merge_completed_process
    else:
        if verbose:
            print('Flipping variants in merged-merge.missnp')
        flip_completed_process = flip_variants(plink_file=second_file,
                                               to_flip='merged-merge.missnp',
                                               path_to_plink1=path_to_plink1,
                                               verbose=verbose)
        merge_completed_process = merge_plink(first_file=first_file,
                                              second_file=second_file,
                                              merged_file=merged_file,
                                              path_to_plink1=path_to_plink1,
                                              path_to_plink2=path_to_plink2,
                                              verbose=verbose)
        if merge_completed_process.returncode == 0:
            return merge_completed_process
        else:
            print('Excluding variants in merged-merge.missnp')
            exclude_proc = exclude_variants(plink_file=second_file,
                                            to_exclude='merged-merge.missnp',
                                            path_to_plink2=path_to_plink2,
                                            verbose=verbose)
            exclude2_proc = exclude_variants(plink_file=first_file,
                                             to_exclude='merged-merge.missnp',
                                             path_to_plink2=path_to_plink2,
                                             verbose=verbose)
            return merge_plink(first_file=first_file,
                               second_file=second_file,
                               merged_file=merged_file,
                               path_to_plink1=path_to_plink1,
                               verbose=verbose)

def run_pca(plink_file: str,
            num_of_samples: int,
            path_to_plink2: str = 'plink2',
            verbose: bool = False) -> CompletedProcess:
    """
    Call Plink to compute PCA.
    """
    if num_of_samples > 5000:
       plink_call = (f'{path_to_plink2} --bfile {plink_file} --pca approx')
#                  f'--pca approx header tabs var-wts')
       return run_plink(plink_call, verbose=verbose)
    else:
       plink_call = (f'{path_to_plink2} --bfile {plink_file} --pca')
       return run_plink(plink_call, verbose=verbose)

def read_pcs_and_populations(pcs_file: str,
                             populations_file: str) -> DataFrame:
    """
    Read in a plink.eigenvec file from PCA.
    Simultaneously, read in the populatins file from
    the reference panel.
    Return a DataFrame of PCs for the merged study and reference
    panel samples with population labels on the reference samples.
    """
    pcs_df = pd.read_csv(pcs_file, sep="\t").drop('#FID', axis=1)
    populations_df = pd.read_csv(populations_file, sep="\t")
    merged = pcs_df.merge(populations_df, on='IID', how='outer')
    return merged

def split_pcs(merged_df: DataFrame) -> Tuple[DataFrame, DataFrame]:
    """
    Given a DataFrame containing both reference and study samples,
    return two DataFrames: one of the study samples and one of the
    reference samples.
    """
    reference_df = merged_df[merged_df['Ancestry'].notna()]
    samples_df = merged_df[merged_df['Ancestry'].isna()]
    return samples_df, reference_df

def update_max_category(data: DataFrame) -> DataFrame:
    # Get columns with probabilities
    prob_columns = data.columns[2:]
    list_inferred_ancestry=[]
    for index, row in data.iterrows():
        # Find the categories with the max probability
        max_categories = [cat for cat in prob_columns if row[cat] == row['MAX']]
    
        # If multiple categories have the same max probability, concatenate them with a delimiter
        if len(max_categories) > 1:
            list_inferred_ancestry.append(','.join(max_categories))
        else:
            list_inferred_ancestry.append(max_categories)
    data['predicted_ancestry']=list_inferred_ancestry
    #clean the dataframe
    data['predicted_ancestry']=data['predicted_ancestry'].astype(str)
    data['predicted_ancestry'] = data['predicted_ancestry'].str.replace(r'[\[\]\'"]', '', regex=True)
    return data

def optimize_classifier(df: DataFrame,
                        max_n_pcs: int = 20,
                        max_n_neighbors: int = 50,
                        num_cross_validation: int = 5,
                        verbose: bool = False) -> Tuple[int, int, list, list, list]:
    """
    Sweep over the number of PCs and the number of neighbors
    and optimize a KNeighborsClassifier.
    Return the optimal number of PCs and the optimal number of neighbors.
    """
    trimmed_df = df.drop(['IID', 'Population'], axis=1)
    optimal_n_pcs, optimal_n_neighbors = 0, 0
    optimal_accuracy = 0
    y = trimmed_df['Ancestry']
    list_for_pcs=[]
    list_for_neighbors=[]
    list_for_accuracy=[]
    for n_pcs in range(2, max_n_pcs+1):
        X = trimmed_df.drop('Ancestry', axis=1).iloc[:, 0:n_pcs]
        for n_neighbors in range(5, max_n_neighbors+1):
            knc = KNeighborsClassifier(n_neighbors=n_neighbors)

	    #debugging 
            if X.isna().any().any():
                print(f"NaNs detected for n_pcs={n_pcs}")
                print(X[X.isna().any(axis=1)].head())
	    #debugging ends here 
	
            knc_scores = cross_val_score(knc, X, y, cv=num_cross_validation,
                                         scoring='balanced_accuracy')
            this_mean_accuracy = knc_scores.mean()
            list_for_pcs.append(n_pcs)
            list_for_neighbors.append(n_neighbors)
            list_for_accuracy.append(this_mean_accuracy)
            if this_mean_accuracy > optimal_accuracy:
                optimal_n_pcs = n_pcs
                optimal_n_neighbors = n_neighbors
                optimal_accuracy = this_mean_accuracy
        if verbose:
            print(f"Done with {n_pcs} principal components.")
    print(f"Optimal number of PCs: {optimal_n_pcs}")
    print(f"Optimal number of neighbors: {optimal_n_neighbors}")
    print(f"Optimal balanced accuracy: {optimal_accuracy:.3f}")
    return optimal_n_pcs, optimal_n_neighbors, list_for_pcs, list_for_neighbors, list_for_accuracy


def main() -> None:
    args = get_args()

    #Validate arguments.
    if not validate_args(args):
        die('Invalid arguments: one or more required files '
            'do not exist.', 255)

    files_to_delete = []  # Append files to clean up at end.

    if args.format_study == 'vcf':
        generate_sex_completed_process = generate_sex_file_vcf_study(args.filename_study, args.plink1_path, args.underscore_sample_id)
    #order chromosome
    ordered_file = 'ordered_study'
    if args.format_study == 'vcf' or args.format_study == 'bfile' or args.format_study == 'pfile':
        order_chromosome_completed_process = order_chromosome(
            study_in=args.filename_study,
            format_in=args.format_study,
            path_to_plink2=args.plink2_path,
            path_to_plink1=args.plink1_path,
            verbose=args.verbose,
	    ordered_file=ordered_file,
            vcf_half_call=args.vcf_half_call,
            genome_build=args.genome_build,
            underscore_sample_id=args.underscore_sample_id)
        pfile_to_bfile_completed_process = pfile_to_bfile(
            ordered_file,
            ordered_file,
            path_to_plink2=args.plink2_path,
            verbose=args.verbose)
        if args.delete_intermediate_files:
            # Delete the filtered variants Plink fileset later.
           files_to_delete.extend(
               ['ordered_study.bim', 'ordered_study.bed',
                'ordered_study.fam', 'ordered_study.log',
                'ordered_study.nosex'])

    #order chromosome for .bgen/.sample files
    if args.format_study == 'bgen':
        check_bgen_argument(args.format_study)
        order_chromosome_completed_process = order_chromosome_bgen(
            bgen_file=args.filename_study+'.bgen',
            sample_file=args.filename_study+'.sample',
            path_to_plink2=args.plink2_path,
            verbose=args.verbose,
            ordered_file=ordered_file,
            underscore_sample_id=args.underscore_sample_id)
        pfile_to_bfile_completed_process = pfile_to_bfile(
            ordered_file,
            ordered_file=ordered_file,
            path_to_plink2=args.plink2_path,
            verbose=args.verbose)
        if args.delete_intermediate_files:
            # Delete the filtered variants Plink fileset later.
           files_to_delete.extend(
               ['ordered_study.bim', 'ordered_study.bed',
                'ordered_study.fam', 'ordered_study.log',
                'ordered_study.nosex'])

    # Subset common, unlinked variants from study:
    ##This step will also conovert vcf format to bim format
    trimmed_study_file = 'trimmed_study'
    filter_variants_completed_process = filter_variants(
        ordered_file,
        trimmed_study_file,
        min_maf=args.min_maf,
        max_missingness=args.max_missingness,
        path_to_plink2=args.plink2_path,
        verbose=args.verbose,
        underscore_sample_id=args.underscore_sample_id)
    if args.delete_intermediate_files:
        # Delete the filtered variants Plink fileset later.
        files_to_delete.extend(
            ['trimmed_study.bim', 'trimmed_study.bed',
             'trimmed_study.fam', 'trimmed_study.log',
             'trimmed_study.nosex'])

    #remove ambiguous
    study_remove_ambiguous_file = 'study_remove_ambiguous'
    ambiguous_site_ids_excluded(bim_file='trimmed_study.bim',
                                var_ids_file='ambiguous.study.in')
    if args.delete_intermediate_files:
        files_to_delete.append('ambiguous.study.in')
    remove_ambiguous_study_completed_process = subset_variants(
        file_to_subset=trimmed_study_file,
        subsetted_file=study_remove_ambiguous_file,
        variants_to_keep='ambiguous.study.in',
        path_to_plink2=args.plink2_path,
        verbose=args.verbose)
    if args.delete_intermediate_files:
        files_to_delete.extend(
            ['study_remove_ambiguous.bed',
             'study_remove_ambiguous.fam',
             'study_remove_ambiguous.bim',
             'study_remove_ambiguous.nosex',
             'study_remove_ambiguous.log'])

    # Recode the variant IDs to "{chrom}:{pos}" for merging.
    set_var_ids('study_remove_ambiguous.bim', chrom_pos_id=True)
    
    #add one new added function for subsetted_study
    #remove multiallelic
    study_remove_multiallelic_file = 'study_remove_multiallelic'
    multi_allelic_ids_excluded(bim_file='study_remove_ambiguous.bim',
                               var_ids_file='multiallelic.study.in')
    if args.delete_intermediate_files:
        files_to_delete.append('multiallelic.study.in')
    remove_multiallelic_study_completed_process = subset_variants(
        file_to_subset=study_remove_ambiguous_file,
        subsetted_file=study_remove_multiallelic_file,
        variants_to_keep='multiallelic.study.in',
        path_to_plink2=args.plink2_path,
        verbose=args.verbose)
    if args.delete_intermediate_files:
        files_to_delete.extend(
            ['study_remove_multiallelic.bed',
             'study_remove_multiallelic.fam',
             'study_remove_multiallelic.bim',
             'study_remove_multiallelic.nosex',
             'study_remove_multiallelic.log'])

    # Check the number of sample of study_remove_multiallelic_file.fam
    with open(r"study_remove_multiallelic.fam","r") as fp:
        num_sample_study=len(fp.readlines())
    
    #if sample size is smaller than 50, gentoype data merged with reference panel for the LD pruning step    
    if num_sample_study<50:
        print('sample size too small to do ld pruning!')
        print('Trying merge with reference data')
         
        # Subset the LD-pruned variants from the reference VCF.
        # This requires the reference VCF to have similarly formatted IDs
        # and be in the same genome build.
        if args.format_reference == 'vcf':
            reference_plink_file = 'reference'
            vcf_to_bfile_completed_process = vcf_to_bfile(
                vcf_file=args.filename_reference,
                plink_file=reference_plink_file,
                path_to_plink1=args.plink1_path,
                path_to_plink2=args.plink2_path,
                verbose=args.verbose,
                underscore_sample_id=args.underscore_sample_id)
            set_var_ids('reference.bim', chrom_pos_id=True)
        elif args.format_reference == 'bgen':
            reference_plink_file = 'reference'
            bgen_to_bfile_completed_process = bgen_to_bfile(
                bgen_file=args.filename_reference+'.bgen',
                sample_file=args.filename_reference+'.sample',
                plink_file=reference_plink_file,
                path_to_plink2=args.plink2_path,
                verbose=args.verbose,
                underscore_sample_id=args.underscore_sample_id)
            set_var_ids('reference.bim', chrom_pos_id=True)
        elif args.format_reference == 'pfile':
            reference_plink_file = 'reference'
            pfile_to_bfile_completed_process = pfile_to_bfile(
                pfile_file=args.filename_reference,
                updated_file=reference_plink_file,
                path_to_plink2=args.plink2_path,
                verbose=args.verbose)
            set_var_ids('reference.bim', chrom_pos_id=True)
        elif args.format_reference == 'bfile':
            reference_plink_file = 'reference'
            bfile_to_bfile_completed_process = bfile_to_bfile(
                bfile_file=args.filename_reference,
                updated_file=reference_plink_file,
                path_to_plink2=args.plink2_path,
                verbose=args.verbose,
                underscore_sample_id=args.underscore_sample_id)
            set_var_ids('reference.bim', chrom_pos_id=True)  

            if args.delete_intermediate_files:
                files_to_delete.extend(
                    ['reference.bed', 'reference.fam',
                     'reference.bim', 'reference.log',
                     'reference.nosex'])

        ##add two new added functions for subsetted_reference
        #remove multiallelic variants
        reference_remove_multiallelic_file = 'reference_remove_multiallelic'
        multi_allelic_ids_excluded(bim_file='reference.bim',
                                   var_ids_file='multiallelic.reference.in')
        if args.delete_intermediate_files:
            files_to_delete.append('multiallelic.reference.in')
        remove_multiallelic_refernce_completed_process = subset_variants(
            file_to_subset=reference_plink_file,
            subsetted_file=reference_remove_multiallelic_file,
            variants_to_keep='multiallelic.reference.in',
            path_to_plink2=args.plink2_path,
            verbose=args.verbose)
        if args.delete_intermediate_files:
            files_to_delete.extend(
                ['reference_remove_multiallelic.bed',
                 'reference_remove_multiallelic.fam',
                 'reference_remove_multiallelic.bim',
                 'reference_remove_multiallelic.nosex',
                 'reference_remove_multiallelic.log'])

        #remove strand ambiguous variants
        reference_remove_ambiguous_file = 'reference_remove_ambiguous'
        ambiguous_site_ids_excluded(bim_file='reference_remove_multiallelic.bim',
                                    var_ids_file='ambiguous.reference.in')
        if args.delete_intermediate_files:
            files_to_delete.append('ambiguous.reference.in')
        remove_ambiguous_refernce_completed_process = subset_variants(
            file_to_subset=reference_remove_multiallelic_file,
            subsetted_file=reference_remove_ambiguous_file,
            variants_to_keep='ambiguous.reference.in',
            path_to_plink2=args.plink2_path,
            verbose=args.verbose)
        if args.delete_intermediate_files:
            files_to_delete.extend(
                ['reference_remove_ambiguous.bed',
                 'reference_remove_ambiguous.fam',
                 'reference_remove_ambiguous.bim',
                 'reference_remove_ambiguous.nosex',
                 'reference_remove_ambiguous.log'])

        # Not all of the variants in the subsetted reference file
        # will also exist in the LD-pruned study, so now
        # reciprocally subset the LD-pruned study:
        subsetted_reference_file = 'subsetted_reference'
        bim_to_var_ids(bim_file='study_remove_multiallelic.bim',
                       var_ids_file='intersection.prune.in')
        if args.delete_intermediate_files:
            files_to_delete.append('intersection.prune.in')
        subset_study_completed_process = subset_variants(
            file_to_subset=reference_remove_ambiguous_file,
            subsetted_file=subsetted_reference_file,
            variants_to_keep='intersection.prune.in',
            path_to_plink2=args.plink2_path,
            verbose=args.verbose)
        if args.delete_intermediate_files:
            files_to_delete.extend(
                ['subsetted_reference.bed',
                 'subsetted_reference.fam',
                 'subsetted_reference.bim',
                 'subsetted_reference.nosex',
                 'subsetted_reference.log'])

        merged_file = 'merged'
        merge_completed_process = five_pass_merge(
        first_file=study_remove_multiallelic_file,
        second_file=subsetted_reference_file,
        merged_file=merged_file,
        path_to_plink2=args.plink2_path,
        path_to_plink1=args.plink1_path,
        verbose=args.verbose)

        # LD prune the study's variants:
        ld_pruned_study_file = 'ld_pruned_merged'
        ld_prune_variants_completed_process = ld_prune_variants(
            plink_file=merged_file,
            pruned_file=ld_pruned_study_file,
            r2=args.max_r2,
            window_size=args.window_size,
            path_to_plink2=args.plink2_path,
            step_size=args.step_size,
            verbose=args.verbose,
            underscore_sample_id=args.underscore_sample_id)
        if args.delete_intermediate_files:
            files_to_delete.extend(
                ['ld_pruned_merged.prune.in',
                 'ld_pruned_merged.prune.out',
                 'ld_pruned_merged.log',
                 'ld_pruned_merged.nosex'])
        
        subsetted_merged_file = 'subsetted_merged'
        subset_reference_completed_process = subset_variants(
            file_to_subset=merged_file,
            subsetted_file=subsetted_merged_file,
            variants_to_keep='ld_pruned_merged.prune.in',
            path_to_plink2=args.plink2_path,
            verbose=args.verbose,
            underscore_sample_id=args.underscore_sample_id)
        if args.delete_intermediate_files:
            files_to_delete.extend(
                ['subsetted_merged.bed',
                 'subsetted_merged.fam',
                 'subsetted_merged.bim',
                 'subsetted_merged.nosex',
                 'subsetted_merged.log'])

   
        # Now run PCA on the merged file:
        pca_completed_process = run_pca(plink_file=subsetted_merged_file,
                                        num_of_samples=num_sample_study,
                                        path_to_plink2=args.plink2_path,
                                        verbose=args.verbose)
        
        if args.delete_intermediate_files:
            files_to_delete.extend(
                ['plink2.log', 'plink2.nosex',
                 'plink2.eigenvec', 'plink2.eigenval',
                 'plink2.eigenvec.var'])

    elif num_sample_study>=50:
        # Not all of the variants in the subsetted reference file
        # will also exist in the LD-pruned study, so now
        # reciprocally subset the LD-pruned study:
        # LD prune the study's variants:
        ld_pruned_study_file = 'ld_pruned_study'
        ld_prune_variants_completed_process = ld_prune_variants(
            plink_file=study_remove_multiallelic_file,
            pruned_file=ld_pruned_study_file,
            r2=args.max_r2,
            window_size=args.window_size,
            path_to_plink2=args.plink2_path,
            step_size=args.step_size,
            verbose=args.verbose,
            underscore_sample_id=args.underscore_sample_id)
        if args.delete_intermediate_files:
            files_to_delete.extend(
                ['ld_pruned_study.prune.in',
                 'ld_pruned_study.prune.out',
                 'ld_pruned_study.log',
                 'ld_pruned_study.nosex'])

        # Subset the LD-pruned varaints from the reference VCF.
        # This requires the reference VCF to have similarly formatted IDs
        # and be in the same genome build.
        if args.format_reference == 'vcf':
            reference_plink_file = 'reference'
            vcf_to_bfile_completed_process = vcf_to_bfile(
                vcf_file=args.filename_reference,
                plink_file=reference_plink_file,
                path_to_plink1=args.plink1_path,
                path_to_plink2=args.plink2_path,
                verbose=args.verbose,
                underscore_sample_id=args.underscore_sample_id)
            set_var_ids('reference.bim', chrom_pos_id=True)
        elif args.format_reference == 'bgen':
            reference_plink_file = 'reference'
            bgen_to_bfile_completed_process = bgen_to_bfile(
                bgen_file=args.filename_reference+'.bgen',
                sample_file=args.filename_reference+'.sample',
                plink_file=reference_plink_file,
                path_to_plink2=args.plink2_path,
                verbose=args.verbose,
                underscore_sample_id=args.underscore_sample_id)
            set_var_ids('reference.bim', chrom_pos_id=True)
        elif args.format_reference == 'pfile':
            reference_plink_file = 'reference'
            pfile_to_bfile_completed_process = pfile_to_bfile(
                pfile_file=args.filename_reference,
                updated_file=reference_plink_file,
                path_to_plink2=args.plink2_path,
                verbose=args.verbose)
            set_var_ids('reference.bim', chrom_pos_id=True)
        elif args.format_reference == 'bfile':
            reference_plink_file = 'reference'
            bfile_to_bfile_completed_process = bfile_to_bfile(
                bfile_file=args.filename_reference,
                updated_file=reference_plink_file,
                path_to_plink2=args.plink2_path,
                verbose=args.verbose,
                underscore_sample_id=args.underscore_sample_id)
            set_var_ids('reference.bim', chrom_pos_id=True)        

            if args.delete_intermediate_files:
                files_to_delete.extend(
                    ['reference.bed', 'reference.fam',
                     'reference.bim', 'reference.log',
                     'reference.nosex'])

        subsetted_reference_file = 'subsetted_reference'
        subset_reference_completed_process = subset_variants(
            file_to_subset=reference_plink_file,
            subsetted_file=subsetted_reference_file,
            variants_to_keep='ld_pruned_study.prune.in',
            path_to_plink2=args.plink2_path,
            verbose=args.verbose,
            underscore_sample_id=args.underscore_sample_id)
        if args.delete_intermediate_files:
            files_to_delete.extend(
                ['subsetted_reference.bed',
                 'subsetted_reference.fam',
                 'subsetted_reference.bim',
                 'subsetted_reference.nosex',
                 'subsetted_reference.log'])

        ##add two new added functions for subsetted_reference
        #remove multiallelic variants
        subsetted_reference_remove_multiallelic_file = 'subsetted_reference_remove_multiallelic'
        multi_allelic_ids_excluded(bim_file='subsetted_reference.bim',
                                   var_ids_file='multiallelic.reference.in')
        if args.delete_intermediate_files:
            files_to_delete.append('multiallelic.reference.in')
        remove_multiallelic_refernce_completed_process = subset_variants(
            file_to_subset=subsetted_reference_file,
            subsetted_file=subsetted_reference_remove_multiallelic_file,
            variants_to_keep='multiallelic.reference.in',
            path_to_plink2=args.plink2_path,
            verbose=args.verbose)
        if args.delete_intermediate_files:
            files_to_delete.extend(
                ['subsetted_reference_remove_multiallelic.bed',
                 'subsetted_reference_remove_multiallelic.fam',
                 'subsetted_reference_remove_multiallelic.bim',
                 'subsetted_reference_remove_multiallelic.nosex',
                 'subsetted_reference_remove_multiallelic.log'])    
    
        #remove strand ambiguous variants
        subsetted_reference_remove_ambiguous_file = 'subsetted_reference_remove_ambiguous'
        ambiguous_site_ids_excluded(bim_file='subsetted_reference_remove_multiallelic.bim',
                               var_ids_file='ambiguous.reference.in')
        if args.delete_intermediate_files:
            files_to_delete.append('ambiguous.reference.in')
        remove_ambiguous_refernce_completed_process = subset_variants(
            file_to_subset=subsetted_reference_remove_multiallelic_file,
            subsetted_file=subsetted_reference_remove_ambiguous_file,
            variants_to_keep='ambiguous.reference.in',
            path_to_plink2=args.plink2_path,
            verbose=args.verbose)
        if args.delete_intermediate_files:
            files_to_delete.extend(
                ['subsetted_reference_remove_ambiguous.bed',
                 'subsetted_reference_remove_ambiguous.fam',
                 'subsetted_reference_remove_ambiguous.bim',
                 'subsetted_reference_remove_ambiguous.nosex',
                 'subsetted_reference_remove_ambiguous.log'])

        # Not all of the variants in the subsetted reference file
        # will also exist in the LD-pruned study, so now
        # reciprocally subset the LD-pruned study:
        subsetted_study_file = 'subsetted_study'
        bim_to_var_ids(bim_file='subsetted_reference_remove_ambiguous.bim',
                       var_ids_file='intersection.prune.in')
        if args.delete_intermediate_files:
            files_to_delete.append('intersection.prune.in')
        subset_study_completed_process = subset_variants(
            file_to_subset=study_remove_ambiguous_file,
            subsetted_file=subsetted_study_file,
            variants_to_keep='intersection.prune.in',
            path_to_plink2=args.plink2_path,
            verbose=args.verbose)
        if args.delete_intermediate_files:
            files_to_delete.extend(
                ['subsetted_study.bed',
                 'subsetted_study.fam',
                 'subsetted_study.bim',
                 'subsetted_study.nosex',
                 'subsetted_study.log'])
    
        # Now we have subsetted reference and study files
        # that contain the same variants. Merge them.
        merged_file = 'merged'
        merge_completed_process = five_pass_merge(
            first_file=subsetted_study_file,
            second_file=subsetted_reference_remove_ambiguous_file,
            merged_file=merged_file,
            path_to_plink2=args.plink2_path,
            path_to_plink1=args.plink1_path,
            verbose=args.verbose)
        if args.delete_intermediate_files:
            files_to_delete.extend(
                ['merged.bed', 'merged.bim',
                 'merged.fam', 'merged.log',
                 'merged.nosex',
                 'merged-merge.missnp',
                 'subsetted_reference_remove_ambiguous.bed~',
                 'subsetted_reference_remove_ambiguous.bim~',
                 'subsetted_reference_remove_ambiguous.fam~',
                 'subsetted_study.bed~',
                 'subsetted_study.bim~',
                 'subsetted_study.fam~'])
    
        # Now run PCA on the merged file:
        pca_completed_process = run_pca(plink_file=merged_file,
                                        num_of_samples=num_sample_study,
                                        path_to_plink2=args.plink2_path,
                                        verbose=args.verbose)

        if args.delete_intermediate_files:
            files_to_delete.extend(
                ['plink2.log', 'plink2.nosex',
                 'plink2.eigenvec', 'plink2.eigenval',
                 'plink2.eigenvec.var'])

    populations_file = args.reference_populations
    pcs_df = read_pcs_and_populations(pcs_file='plink2.eigenvec',
                                      populations_file=populations_file)
    study_df, reference_df = split_pcs(pcs_df)
    # Plot PC1 vs PC2 from reference.
    # You should see spatial separation of the different ancestries.
    
    # % variance explained by PCA components
    print("Plot % variance explained by PCA components ")
    eigenval=pd.read_csv('plink2.eigenval', sep='\s+', header=None)
    eigenval=eigenval.iloc[:,0].to_list()
    eigenval_fractional=eigenval/np.sum(eigenval)

    plt.figure(figsize=(6,6))
    plt.plot((eigenval_fractional)[:20]*100, '.-')
    plt.xticks(range(0,20), list(range(1,21)))
    plt.ylabel('% variance explained')
    plt.xlabel('number of PC')
    plt.savefig('Percent_variance_explained.png')


    sns.set_context('poster')

    fig, axes = plt.subplots(ncols=1, nrows=3, figsize=[14,36])

    scatterplot(x=reference_df['PC1'],
                y=reference_df['PC2'],
                hue=reference_df['Ancestry'],
                linewidth=0.1,
                ax=axes[0],
                alpha=0.5
                )
    scatterplot(x=reference_df['PC1'],
                y=reference_df['PC3'],
                hue=reference_df['Ancestry'],
                linewidth=0.1,
                ax=axes[1],
                alpha=0.5
                )
    scatterplot(x=reference_df['PC2'],
                y=reference_df['PC3'],
                hue=reference_df['Ancestry'],
                linewidth=0.1,
                ax=axes[2],
                alpha=0.5
                )
    #plt.legend(loc='best',bbox_to_anchor=(1, -0.05))
    plt.savefig('PCplot_reference.png')
    plt.savefig("PCplot_reference.pdf", format="pdf", bbox_inches="tight")

    # Optimize a classifier on the reference data using CV.
    verbose = args.verbose
    optimal_n_pcs, optimal_n_neighbors, list_for_pcs, list_for_neighbors, list_for_accuracy = optimize_classifier(df=reference_df,
                                                                                                                  max_n_pcs=args.number_of_principal_components,
                                                                                                                  max_n_neighbors=args.number_of_neighbors,
                                                                                                                  num_cross_validation = args.num_cross_validation,
                                                                                                                  verbose=verbose)
    print(list_for_accuracy)

    # Plot validation accuracy figure as function of number of PC and N
    data_validation_accuracy = {'PC': list_for_pcs,
                                'K': list_for_neighbors,
                                'Accuracy': list_for_accuracy}
    df_validation_accuracy = pd.DataFrame(data_validation_accuracy)
    plt.figure(figsize=(10, 8))
    plt.scatter(df_validation_accuracy['PC'], df_validation_accuracy['K'], 
                c=df_validation_accuracy['Accuracy'], cmap='viridis', s=10)
    plt.colorbar(label='Validation Accuracy')
    plt.xlabel('Number of principal components')
    plt.ylabel('Number of neighbors')
    plt.show()
    plt.savefig('Accuracy_K_PC_training_validation.png')
    plt.savefig("Accuracy_K_PC_training_validation.pdf", format="pdf", bbox_inches="tight")
    # Fit the optimized classifier and use it generate predictions.
    ref_X = reference_df.drop(['IID', 'Population', 'Ancestry'], axis=1).iloc[:, 0:optimal_n_pcs]
    ref_y = reference_df['Ancestry']
    knc = KNeighborsClassifier(n_neighbors=optimal_n_neighbors)
    knc.fit(ref_X, ref_y)

    study_X = study_df.drop(['IID', 'Population', 'Ancestry'], axis=1).iloc[:, 0:optimal_n_pcs]
    study_y_pred_series = Series(knc.predict(study_X), index=study_df['IID'])

    study_y_pred_PCAplots = knc.predict(study_X)
    study_y_pred_proba = knc.predict_proba(study_X)
    study_y_pred_proba_df = pd.DataFrame(study_y_pred_proba,
                                         columns=knc.classes_,
                                         index=study_df['IID'])

    study_y_pred_proba_max = np.max(study_y_pred_proba, axis=1)
    # Make output table.
    output_df = pd.DataFrame(study_y_pred_proba,
                             columns=knc.classes_,
                             index=study_df['IID'])
    output_df.insert(0, 'predicted_ancestry', study_y_pred_series)
    output_df.insert(1, 'MAX', study_y_pred_proba_max)
    #add steps to identify the equal probability and update the predicted_ancestry column
    #by identifying two or more categories have the same prob and = MAX
    output_df_update = update_max_category(output_df)
    output_df_update.to_csv(args.output_file, sep='\t')
    print(f"Finished generating {args.output_file}.")

    # How many predictions of each Superpopulation are in our study?
    print(pd.Series(study_y_pred_series).value_counts())

    df_study_y_pred_iid_pc=study_df.loc[:,['IID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']]
    df_study_y_pred_ancestry=pd.DataFrame(study_y_pred_PCAplots)
    df_study_y_pred_iid_pc.reset_index(drop=True, inplace=True)
    df_study_y_pred_ancestry.reset_index(drop=True, inplace=True)
    df_study_y_pred=pd.concat([df_study_y_pred_iid_pc,df_study_y_pred_ancestry],ignore_index=True,axis="columns")
    df_study_PCAplot=df_study_y_pred
    df_study_PCAplot['Type']='study'
    df_study_PCAplot.columns=['IID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','Ancestry','Type']    
    
    # Visualize PCA for the study:
    if not args.study_covariates:
        sns.set_context('poster')
        fig, axes = plt.subplots(ncols=1, nrows=3, figsize=[14,36])
        scatterplot(x=df_study_PCAplot['PC1'],
                    y=df_study_PCAplot['PC2'],
                    hue=df_study_PCAplot['Ancestry'],
                    linewidth=0.1,
                    ax=axes[0],
                    alpha=0.5
                    )
        scatterplot(x=df_study_PCAplot['PC1'],
                    y=df_study_PCAplot['PC3'],
                    hue=df_study_PCAplot['Ancestry'],
                    linewidth=0.1,
                    ax=axes[1],
                    alpha=0.5
                    )
        scatterplot(x=df_study_PCAplot['PC2'],
                    y=df_study_PCAplot['PC3'],
                    hue=df_study_PCAplot['Ancestry'],
                    linewidth=0.1,
                    ax=axes[2],
                    alpha=0.5
                    )
        #plt.legend(loc='best', bbox_to_anchor=(1,-0.05))
        plt.savefig('PCplot_study.png')
        plt.savefig("PCplot_study.pdf", format="pdf", bbox_inches="tight")
    elif args.study_covariates:
        df_study_covariates = pd.read_csv(args.study_covariates, sep=None, engine="python")
        #df_study_covariates = pd.read_csv(args.study_covariates, sep="\t")
        df_study_covariates.columns=['IID','Covariate']
        df_study_PCAplot_covariates = pd.merge(df_study_PCAplot, df_study_covariates, on=['IID'], how='inner')
        print(df_study_PCAplot_covariates)
        
        sns.set_context('poster', font_scale=2.0)
        fig, axes = plt.subplots(ncols=1, nrows=3, figsize=[14,36])
        scatterplot(x=df_study_PCAplot_covariates['PC1'],
                    y=df_study_PCAplot_covariates['PC2'],
                    hue=df_study_PCAplot_covariates['Ancestry'],
                    style=df_study_PCAplot_covariates['Covariate'],
                    linewidth=0.1,
                    ax=axes[0],
                    alpha=0.5
                    )
                    
        axes[0].set_xlabel('PC1', fontsize=28)
        axes[0].set_ylabel('PC2', fontsize=28)
        axes[0].tick_params(axis='both', labelsize=24)
        axes[0].legend(fontsize=20, title_fontsize=22)
        
        scatterplot(x=df_study_PCAplot_covariates['PC1'],
                    y=df_study_PCAplot_covariates['PC3'],
                    hue=df_study_PCAplot_covariates['Ancestry'],
                    style=df_study_PCAplot_covariates['Covariate'],
                    linewidth=0.1,
                    ax=axes[1],
                    alpha=0.5
                    )
                    
        axes[1].set_xlabel('PC1', fontsize=28)
        axes[1].set_ylabel('PC3', fontsize=28)
        axes[1].tick_params(axis='both', labelsize=24)
        axes[1].legend(fontsize=20, title_fontsize=22)

        scatterplot(x=df_study_PCAplot_covariates['PC2'],
                    y=df_study_PCAplot_covariates['PC3'],
                    hue=df_study_PCAplot_covariates['Ancestry'],
                    style=df_study_PCAplot_covariates['Covariate'],
                    linewidth=0.1,
                    ax=axes[2],
                    alpha=0.5
                    )
        axes[2].set_xlabel('PC2', fontsize=28)
        axes[2].set_ylabel('PC3', fontsize=28)
        axes[2].tick_params(axis='both', labelsize=24)
        axes[2].legend(fontsize=20, title_fontsize=22)

        #plt.legend(loc='best', bbox_to_anchor=(1,-0.05))
        plt.savefig('PCplot_study.png')
        plt.savefig("PCplot_study.pdf", format="pdf", bbox_inches="tight")

    # Overlay PCA for the reference and study in one plot:
    df_reference_PCAplot=reference_df.loc[:,['IID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','Ancestry']]
    df_reference_PCAplot['Type']='reference'
    
    df_study_PCAplot_overlay=df_study_PCAplot
    df_study_PCAplot_overlay.columns=['IID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','Ancestry','Type']
    df_ref_study_PCAplot=pd.concat([df_reference_PCAplot,df_study_PCAplot_overlay],ignore_index=True,axis="rows")
    print(df_ref_study_PCAplot)
    sns.set_context('poster', font_scale=2.0)
    fig, axes = plt.subplots(ncols=1, nrows=3, figsize=[14,36])
    scatterplot(x=df_ref_study_PCAplot['PC1'],
                y=df_ref_study_PCAplot['PC2'],
                hue=df_ref_study_PCAplot['Ancestry'],
                style=df_ref_study_PCAplot['Type'],
                linewidth=0.1,
                ax=axes[0],
                alpha=0.5
                )
    axes[0].set_xlabel('PC1', fontsize=28)
    axes[0].set_ylabel('PC2', fontsize=28)
    axes[0].tick_params(axis='both', labelsize=24)
    axes[0].legend(fontsize=20, title_fontsize=22)
    scatterplot(x=df_ref_study_PCAplot['PC1'],
                y=df_ref_study_PCAplot['PC3'],
                hue=df_ref_study_PCAplot['Ancestry'],
                style=df_ref_study_PCAplot['Type'],
                linewidth=0.1,
                ax=axes[1],
                alpha=0.5
                )
    axes[1].set_xlabel('PC1', fontsize=28)
    axes[1].set_ylabel('PC3', fontsize=28)
    axes[1].tick_params(axis='both', labelsize=24)
    axes[1].legend(fontsize=20, title_fontsize=22)
    scatterplot(x=df_ref_study_PCAplot['PC2'],
                y=df_ref_study_PCAplot['PC3'],
                hue=df_ref_study_PCAplot['Ancestry'],
                style=df_ref_study_PCAplot['Type'],
                linewidth=0.1,
                ax=axes[2],
                alpha=0.5
                )
                
    axes[2].set_xlabel('PC2', fontsize=28)
    axes[2].set_ylabel('PC3', fontsize=28)
    axes[2].tick_params(axis='both', labelsize=24)
    axes[2].legend(fontsize=20, title_fontsize=22)
    #plt.legend(loc='best', bbox_to_anchor=(1,-0.05))
    plt.savefig('PCplot_overlay.png')
    plt.savefig("PCplot_overlay.pdf", format="pdf", bbox_inches="tight")

    if args.delete_intermediate_files:
        paths_to_delete = [Path(f) for f in files_to_delete]
        [f.unlink() for f in paths_to_delete if f.is_file()]

if __name__ == "__main__":
    st=time.time()
    main()
    et=time.time()
    elapsed_time=et-st
    print('Execution time:', elapsed_time, 'seconds')
