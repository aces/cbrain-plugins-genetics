#!/usr/bin/python3

# This script is a wrapper for the scRNAbox command line interface.
# It is used to run the scRNAbox pipeline from the command line.

# Based on a command line that contain all the parameters available in the scRNAbox pipeline,
# this script will generate the step1_*.txt file and the step2_*.txt file that will be used to run the pipeline.

# The script will also generate the command line that will be used to run the pipeline.

import argparse
import subprocess
import sys
import os
import glob
import time

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Utility methods                  #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

def print_log(message):
    subprocess.run(f"echo $(date +\"%Y-%m-%d %H:%M:%S\"): {message}", shell=True)

def is_dir(path):
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} is not a valid path")

def is_file(path):
    if os.path.isfile(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} is not a valid file")



def main():
    version = "1.0.0"
    print_log(f"Starting scRNAbox wrapper script")
    print_log(f"By Natacha Beck \(nbeck@mcin.ca\)")
    print_log(f"Version: {version}")

    parser = argparse.ArgumentParser(description='scRNAbox command line wrapper')

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
    #                                Pipeline parameters                             #
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

    ############################
    # Configuration parameters #
    ############################

    configuration_options = parser.add_argument_group('Configuration parameters')

    # -R  (--R_path)  =  Path to R executable
    configuration_options.add_argument('--R_path', type=is_dir, help='Path to R executable')

    # job_mode
    configuration_options.add_argument('--job_mode', type=str, help='Job mode (local, slurm)')

    # output_dir
    configuration_options.add_argument('--output_dir', type=str, help='Output directory')

    # Reference data directory
    configuration_options.add_argument('--ref_data_dir', type=is_dir, help='Reference data directory')

    ######################
    # Parameters #
    ######################

    general_options = parser.add_argument_group('General parameters')

    # id: General
    ################

    # -d  (--dir)  = Working directory (where all the outputs will be printed) (give full path)
    general_options.add_argument('-d', '--dir', type=str, help='Working directory (where all the outputs will be printed)')

    # --method  = Select your preferred method: HTO and SCRNA for hashtag, and Standard scRNA, respectively.
    general_options.add_argument('--method', type=str, help='Select the appropriate analysis track. For the Standard Track, select (SCRNA). For the HTO Track, select (HTO). The Standard Track is designed for experiments where each sample is captured and sequenced separately, while the HTO Track is designed for multiplexed experiments where samples are tagged with sample-specific oligonucleotide tagged Hashtag antibodies (HTO), pooled, and sequenced together.')

    # --steps  =  Specify what steps, e.g., 2 to run step 2. 2-6, run steps 2 through 6
    general_options.add_argument('--steps', type=str, help='Specify which steps to execute in this run. E.g., 2 will execute Step 2; 2-6 will execute Steps 2 through 6.')

    # id: Step 1
    ###############

    group1_options = parser.add_argument_group('Step 1')

    # id: Step 1 SCRNA Track: FASTQ to expression matrix
    #######################################################

    ## Do you want to perform automated library prep?
    #par_automated_library_prep= "no"
    group1_options.add_argument('--par_automated_library_prep_scrna', action='store_true', help='Whether or not to perform automated library prep. Alternatively, you may set this parameter to "no" and manually prepare the libraries. For most use cases, it is recommended to perform automated library prep.')

    ## Path to the directory containing the FASTQ files for the experiment. This folder should only contain the FASTQ files for the experiment.
    #par_fastq_directory= "/path/to/fastqs/directory"
    group1_options.add_argument('--par_fastq_directory_scrna', type=is_dir, help='Path to directory containing the FASTQ files. This directory should only contain FASTQ files for the experiment.')

    ## List the sample names used in the FASTQ nomenclature
    #par_sample_names= c("Sample1", "Sample2", "Sample3")
    group1_options.add_argument('--par_sample_names_step1', type=str, help='The sample names used to name the FASTQ files according to CellRanger nomenclature. E.g., if the FASTQ file is named ‘PD1_S1_L001_R1_001.fastq’, the sample name would be "PD1".')

    ## If you want to rename the samples, set par_rename_samples to yes.
    #par_rename_samples= "yes"
    group1_options.add_argument('--par_rename_samples_scrna', action='store_true', help='Whether or not you want to rename your samples. These names will be used to identify cells in the Seurat objects.')

    ## If you want to renames the samples (i.e. par_rename_samples= "yes"), list the new sample names in the same order as the old labels are listed in  par_sample_names.
    #par_new_sample_names= c("NewSample1", "NewSample2", "NewSample3")
    group1_options.add_argument('--par_new_sample_names', type=str, help='New sample names. Make sure they are defined in the same order as \'par_sample_names\'. E.g., you may wish to rename sample "PD1" to "Parkinson1".')

    ## If your sequencing is paired-end, set the following to TRUE. Otherwise set it as FALSE.
    #par_paired_end_seq=TRUE
    group1_options.add_argument('--par_paired_end_seq_scrna', action='store_true', help='Whether or not paired-end sequencing was performed.')

    # CellRanger counts pipeline parameters.
    ## Path to reference genome
    #par_ref_dir_grch='/path/to/CellRanger/reference/genome'
    group1_options.add_argument('--par_ref_dir_grch_scrna', type=str, help='Path to reference genome for FASTQ alignment.')

    ## Minimum number of bases to retain for R1 sequence of gene expression assay. If you want to use this parameter uncomment the line below and define your par_r1_length.
    #par_r1_length=20
    group1_options.add_argument('--par_r1_length_scrna', type=int, help='Minimum number of bases to retain for R1 sequence of gene expression. For most use-cases, this can be kept at 0 (default).')

    ## Minimum number of bases to retain for R2 sequence of gene expression assay. If you want to use this parameter uncomment the line below and define your par_r2_length.
    #par_r2_length=20
    group1_options.add_argument('--par_r2_length_scrna', type=int, help='Minimum number of bases to retain for R2 sequence of gene expression. For most use-cases, this can be kept at 0 (default).')

    ## If you want CellRnager to include introns when producing the gene expression matrices set the following parameter to "yes", otherwise keep the default as "no".
    #par_include_introns="no"
    group1_options.add_argument('--par_include_introns_scrna', action='store_true', help='Whether or not to include intronic reads in the gene expression matrix.')

    ## If you want to turn off CellRanger's target UMI filtering subpipeline uncomment the parameter below.
    #par_no_target_umi_filter="no"
    group1_options.add_argument('--par_no_target_umi_filter_scrna', action='store_true', help="Whether or not to turn off CellRanger's target UMI filtering sub-pipeline. For most use-cases, this can remain disabled (default).")

    ## If you want to specify the number of expected cells, uncomment the parameter below and enter the value. By default, CellRanger's auto-estimate algorithm will be used.
    #par_expect_cells=6000
    group1_options.add_argument('--par_expect_cells_scrna', type=int, help="Expected number of cells. By default, CellRanger's auto-estimate algorithm will be used. For most use-cases, this can be kept at 0 (default).")

    ## If you want to force the CellRanger count pipeline to use a certain number of cells, uncomment the parameter below and enter the number of cells
    #par_force_cells=6000
    group1_options.add_argument('--par_force_cells_scrna', type=int, help='Force the CellRanger count pipeline to use N cells. For most use-cases, this can be kept at 0 (default).')

    ## If you want to skip the bam file generation, uncomment the parameter below.
    #par_no_bam="no"
    group1_options.add_argument('--par_no_bam_scrna', action='store_true', help='Whether or not to skip the bam file generation in the CellRanger pipeline. For most use-cases, this can remain disabled (default).')

    # id: Step 1 HTO Track: FASTQ to expression matrix
    #####################################################

    ## Do you want to perform automated library prep?
    #par_automated_library_prep= "no"
    group1_options.add_argument('--par_automated_library_prep_hto', action='store_true', help='Whether or not to perform automated library prep. Alternatively, you may set this parameter to "no" and manually prepare the libraries. For most use cases, it is recommended to perform automated library prep.')

    ## Path to the directory containing the FASTQ files for the experiment. This folder should only contain the FASTQ files for the experiment.
    #par_fastq_directory= "/path/to/fastqs/directory"
    group1_options.add_argument('--par_fastq_directory_hto', type=is_dir, help='Path to directory containing the FASTQ files. This directory should only contain FASTQ files for the experiment.')

    # par_RNA_run_names	NULL	The names of the sequencing runs for the RNA assay
    group1_options.add_argument('--par_RNA_run_names', type=str, help='The names of the sequencing runs for the RNA assay.')

    # par_HTO_run_names	NULL The names of the sequencing runs for the HTO assay
    group1_options.add_argument('--par_HTO_run_names', type=str, help='The names of the sequencing runs for the HTO assay.')

    # par_seq_run_names	NULL	The user-selected name for the sequencing run. These names will be used to identify cells in the Seurat objects
    group1_options.add_argument('--par_seq_run_names', type=str, help='The user-defined name for the sequencing run. These names will be used to identify cells in the Seurat objects.')

    # id	NULL	Barcode ID which will be used to track the feature counts
    group1_options.add_argument('--id', type=str, help='Barcode ID that will be used to track the feature counts.')

    # name	NULL	The user-selected name for the barcode identifier
    group1_options.add_argument('--name', type=str, help='The user-defined name for the barcode identifier.')

    # read	R2	Which RNA sequencing read contains the barcode sequence. This value Will be either R1 or R2.
    group1_options.add_argument('--read', type=str, help='Which RNA sequencing read contains the barcode sequence. This value will be either R1 or R2.')

    # pattern	NULL	The pattern of the barcode identifiers
    group1_options.add_argument('--pattern', type=str, help='The pattern of the barcode identifiers.')

    # sequence	NULL	The nucleotide sequence associated with the barcode identifier
    group1_options.add_argument('--sequence', type=str, help='The nucleotide sequence for the barcode identifiers.')

    # CellRanger counts pipeline parameters.
    ## Path to reference genome
    #par_ref_dir_grch='/path/to/CellRanger/reference/genome'
    group1_options.add_argument('--par_ref_dir_grch_hto', type=str, help='Path to reference genome for FASTQ alignment.')

    ## If your sequencing is paired-end, set the following to TRUE. Otherwise set it as FALSE.
    #par_paired_end_seq=TRUE
    group1_options.add_argument('--par_paired_end_seq_hto', action='store_true', help='Whether or not paired-end sequencing was performed.')

    ## Minimum number of bases to retain for R1 sequence of gene expression assay. If you want to use this parameter uncomment the line below and define your par_r1_length.
    #par_r1_length=20
    group1_options.add_argument('--par_r1_length_hto', type=int, help='Minimum number of bases to retain for R1 sequence of gene expression. For most use-cases, this can be kept at 0 (default).')

    ## Minimum number of bases to retain for R2 sequence of gene expression assay. If you want to use this parameter uncomment the line below and define your par_r2_length.
    #par_r2_length=20
    group1_options.add_argument('--par_r2_length_hto', type=int, help='Minimum number of bases to retain for R2 sequence of gene expression. For most use-cases, this can be kept at 0 (default).')

    ## If you want CellRnager to include introns when producing the gene expression matrices set the following parameter to "yes", otherwise keep the default as "no".
    #par_include_introns="no"
    group1_options.add_argument('--par_include_introns_hto', action='store_true', help='Whether or not to include intronic reads in the gene expression matrix. ')

    ## If you want to turn off CellRanger's target UMI filtering subpipeline uncomment the parameter below.
    #par_no_target_umi_filter="no"
    group1_options.add_argument('--par_no_target_umi_filter_hto', action='store_true', help="Whether or not to turn off CellRanger's target UMI filtering sub-pipeline. For most use-cases, this can remain disabled (default).")

    ## If you want to specify the number of expected cells, uncomment the parameter below and enter the value. By default, CellRanger's auto-estimate algorithm will be used.
    #par_expect_cells=6000
    group1_options.add_argument('--par_expect_cells_hto', type=int, help="Expected number of cells. By default, CellRanger's auto-estimate algorithm will be used. For most use-cases, this can be kept at 0 (default)..")

    ## If you want to force the CellRanger count pipeline to use a certain number of cells, uncomment the parameter below and enter the number of cells
    #par_force_cells=6000
    group1_options.add_argument('--par_force_cells_hto', type=int, help='Force the CellRanger count pipeline to use N cells. For most use-cases, this can be kept at 0 (default).')

    ## If you want to skip the bam file generation, uncomment the parameter below.
    #par_no_bam="no"
    group1_options.add_argument('--par_no_bam_hto', action='store_true', help='Whether or not to skip the bam file generation in the CellRanger pipeline. For most use-cases, this can remain disabled (default).')

    # id: Step 2
    # Step 2: Create Seurat object and remove ambient RNA
    ########################################################

    group2_options = parser.add_argument_group('Step 2')

    # If you want to save an RNA expression matrix and metadata dataframe set the following to "yes"
    ###################################################################################################
    #step2_par_save_RNA= "yes"
    group2_options.add_argument('--par_save_RNA_step2', action='store_true', help='Whether or not to export an RNA expression matrix.')

    # If you want to save a metadata dataframe set the following to "yes"
    #step2_par_save_metadata= "yes"
    group2_options.add_argument('--par_save_metadata_step2', action='store_true', help='Whether or not to export a metadata dataframe.')

    ## If you want to remove ambient RNA from the expression matrix, keep the default as “yes”. If you do not want values changed to remove ambient RNA, change to “no”.
    # par_ambient_RNA= "yes"
    group2_options.add_argument('--par_ambient_RNA', action='store_true', help='Whether or not to correct the feature-barcode expression matrices for ambient RNA contamination.')

    ## Only retain genes that are present in at least a specified number of cells.
    # par_min.cells_L= 3
    group2_options.add_argument('--par_min_cells_L', type=int, help='Only retain genes expressed in a minimum number of cells.')

    ## Normalization method
    # par_normalization.method= "LogNormalize"
    group2_options.add_argument('--par_normalization_method_step2', type=str, help='Method to use for normalization (e.g., LogNormalize).')

    ## Scale factor
    # par_scale.factor= 10000
    group2_options.add_argument('--par_scale_factor_step2', type=int, help='Scale factor for scaling the data (e.g., 10000).')

    ## Method for choosing the top variable features. vst, mean.var.plot (mvp), dispersion (disp).
    # par_selection.method= "vst"
    group2_options.add_argument('--par_selection_method_step2', type=str, help='Method for selecting highly variable genes (e.g., vst).')

    ## Number of features to select as top variable features
    # par_nfeatures= 2500
    group2_options.add_argument('--par_nfeatures_step2', type=int, help='Number of highly variable genes to select (e.g., 2500).')

    # id: Step 3
    # Step 3: Quality control and filtering
    ##########################################

    group3_options = parser.add_argument_group('Step 3')

    # If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
    ############################################################################
    # par_save_RNA= "yes"
    group3_options.add_argument('--par_save_RNA_step3', action='store_true', help='Whether or not to export an RNA expression matrix.')
    # par_save_metadata= "yes"
    group3_options.add_argument('--par_save_metadata_step3', action='store_true', help='Whether or not to export a metadata dataframe.')

    # If you already have a processed Seurat RDS object, and did not perform Step 2 of scRNAbox,
    # use parameter this to add the path to the directory containing you Seurat object(s).
    # Make sure that no other files/objects are present in the directory besides Seurat RDS objects.
    # Uncomment the line to activate the parameter
    ############################################################################
    #par_seurat_object= "/path/to/directory/containing/seurat/object"
    group3_options.add_argument('--par_seurat_object_step3', type=is_dir, help='If you performed Step 2 of the pipeline, keep this parameter empty. Otherwise, if users have an existing Seurat object(s), they may provide the path to a directory that contains the Seurat object(s) to initiate the pipeline at Step 3.')

    # Quality control parameters
    # Uncomment the line to activate the parameter and add the desired value. Cells will be filtered out accordingly.
    # L = lower bound threshold
    # R = upper bound threshold
    ############################################################################
    ## Minimum number of unique RNA transcripts
    # par_nFeature_RNA_L= 300
    group3_options.add_argument('--par_nFeature_RNA_L', type=int, help='Only retain cells expressing a minimum number of unique RNA transcripts.')

    ## Maximum number of unique RNA transcripts
    # par_nFeature_RNA_U= 10000
    group3_options.add_argument('--par_nFeature_RNA_U', type=int, help='Only retain cells expressing a maximum number of unique RNA transcripts.')

    ## Minimum number of total RNA transcripts
    # par_nCount_RNA_L= 300
    group3_options.add_argument('--par_nCount_RNA_L', type=int, help='Only retain cells expressing a minimum number of total RNA transcripts.')

    ## Maximum number of total RNA transcripts
    # par_nCount_RNA_U= 20000
    group3_options.add_argument('--par_nCount_RNA_U', type=int, help='Only retain cells expressing a maximum number of total RNA transcripts.')

    ## Minimum mitochondrial RNA percentage
    # par_mitochondria_percent_L= 0
    group3_options.add_argument('--par_mitochondria_percent_L', type=int, help='Only retain cells with a minimum percentage of mitochondrial-encoded genes.')

    ## Maximum mitochondrial RNA percentage
    # par_mitochondria_percent_U= 20
    group3_options.add_argument('--par_mitochondria_percent_U', type=int, help='Only retain cells with a maximum percentage of mitochondrial-encoded genes.')

    ## Minimum ribosomal RNA percentage
    # par_ribosomal_percent_L= 0
    group3_options.add_argument('--par_ribosomal_percent_L', type=int, help='Only retain cells with a minimum percentage of ribosome genes.')

    ## Maximum ribosomal RNA percentage
    # par_ribosomal_percent_U= 100
    group3_options.add_argument('--par_ribosomal_percent_U', type=int, help='Only retain cells with a maximum percentage of ribosome genes.')

    # Parameters to filter out genes
    ############################################################################
    ## If you want to filter out mitochondrial and ribosomal genes set the following parameters to "yes". If you do not want to remove them keep the default as "no".
    # par_remove_mitochondrial_genes= "no"
    group3_options.add_argument('--par_remove_mitochondrial_genes', action='store_true', help='Whether or not to remove mitochondrial genes.')

    # par_remove_ribosomal_genes= "no"
    group3_options.add_argument('--par_remove_ribosomal_genes', action='store_true', help='Whether or not to remove ribosomal genes.')

    ## If you have specific genes that you want to remove, enter a vector of the genes. Uncomment the line to activate the parameter.
    #par_remove_genes= c("gene1", "gene2")
    group3_options.add_argument('--par_remove_genes', type=str, help='If users want to remove specific genes from their data, they may define a list of gene identifiers.')

    # Regress genes
    ############################################################################
    ## If you want to regress cell cycle genes, set the following parameters to "yes". If you do not want to regress them, keep the default as "no". Note: if you are using your own Seurat object (i.e. not from Step 2), you can only regress cell cycle genes if your Seurat object has the cell cycle score computed.
    # par_regress_cell_cycle_genes= "no"
    group3_options.add_argument('--par_regress_cell_cycle_genes', action='store_true', help='Whether or not to regress cell cycle genes.')

    ## If you want to regress a custom list of genes, set the following parameters to "yes". If you do not want to regress a custom list, keep the default as "no".
    # par_regress_custom_genes= "no"
    group3_options.add_argument('--par_regress_custom_genes', action='store_true', help='Whether or not to regress a custom list of genes.')

    ## Enter the genes that you want to regress in the list below.
    # par_regress_genes= c("gene1", "gene2")
    group3_options.add_argument('--par_regress_genes', type=str, help='List of custom genes to regress.')

    ## Normalization method
    # par_normalization.method= "LogNormalize"
    group3_options.add_argument('--par_normalization_method_step3', type=str, help='Method to use for normalization (e.g., LogNormalize).')

    ## Scale factor
    # par_scale.factor= 10000
    group3_options.add_argument('--par_scale_factor_step3', type=int, help='Scale factor for scaling the data (e.g., 10000).')

    ## Method for choosing the top variable features. vst, mean.var.plot (mvp), dispersion (disp).
    # par_selection.method= "vst"
    group3_options.add_argument('--par_selection_method_step3', type=str, help='Method for selecting highly variable genes (e.g., vst).')

    ## Number of features to select as top variable features
    # par_nfeatures= 2500
    group3_options.add_argument('--par_nfeatures_step3', type=int, help='Number of highly variable genes to select (e.g., 2500).')

    # Parameters for normalization and scaling after quality control
    ############################################################################

    ## Number of most variable features to be reported in csv file
    # par_top= 10
    group3_options.add_argument('--par_top', type=int, help='Number of highly variable genes to be reported in the exported csv file (e.g., 100).')

    ## Total Number of PCs to compute and store for RunPCA
    # par_npcs_pca= 30
    group3_options.add_argument('--par_npcs_pca_step3', type=int, help='Total number of principal components to compute and store for principal component analysis (e.g., 50).')

    # id: Step 4
    ###############

    # Step 4 SCRNA Track: Doublet removal
    ########################################

    group4_options = parser.add_argument_group('Step 4')

    # If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
    ############################################################################
    # par_save_RNA= "yes"
    group4_options.add_argument('--par_save_RNA_step4_scrna', action='store_true', help='Whether or not to export an RNA expression matrix.')
    # par_save_metadata= "yes"
    group4_options.add_argument('--par_save_metadata_step4_scrna', action='store_true', help='Whether or not to save a metadata dataframe.')

    # If you already have a processed Seurat RDS object, and did not perform Step 3 of scRNAbox use this to add the path to the directory containing you Seurat object(s).
    # Make sure that no other files/objects are present in the directory besides Seurat RDS objects.
    # Uncomment the line to activate the parameter
    ############################################################################
    #par_seurat_object= "/path/to/directory/containing/seurat/object"
    group4_options.add_argument('--par_seurat_object_step4_scrna', type=is_dir, help='If you performed Step 3 of the pipeline, keep this parameter empty. Otherwise, if users have an existing Seurat object(s), they may provide the path to a directory that contains the Seurat object(s) to initiate the pipeline at Step 4.')

    # Parameters for UMAP dimensional reduction
    ############################################################################
    ## Number of dimensions to use as input into UMAP
    # par_RunUMAP_dims= 25
    group4_options.add_argument('--par_RunUMAP_dims_step4_scrna', type=int, help='Number of dimensions to use as input features for uniform manifold approximation and projection (e.g., 25).')

    ## Number of neighbouring points to use in local approximation of manifold structure
    # par_RunUMAP_n.neighbors= 45
    group4_options.add_argument('--par_RunUMAP_n_neighbors_step4', type=int, help='Number of neighboring points used in local approximations of manifold structure (e.g., 45).')

    # Parameters for doublet detection and removal (optional)
    ############################################################################
    ## If you want to remove predicted doublets from downstream analyses set the following to "yes"
    ## If you want to keep predicted doublets for further analysis set the following to "no"
    # par_dropDN= "yes"
    group4_options.add_argument('--par_dropDN_scrna', action='store_true', help='Whether or not to remove predicted doublets from downstream analyses.')

    ## Number of principal components to use as input doublet analysis.
    ## This can be determined by the bend in the by elbow plot produced in Step 3
    # par_PCs= 25
    group4_options.add_argument('--par_PCs', type=int, help='The number of principal components to use for doublet detection by DoubletFinder (e.g., 25). Can be informed by the elbow plot produced in Step 3.')

    ## The number of artificial doublets to generate. DoubletFinder is largely invariant to this parameter. We suggest keeping 0.25
    # par_pN= 0.25
    group4_options.add_argument('--par_pN', type=float, help='The number of artificial doublets to generate for doublet detection by DoubletFinder (e.g., 0.25). DoubletFinder is largely invariant to this parameter. We suggest keeping 0.25.')

    ## Logical representing whether SCTransform was used during original Seurat object pre-processing
    # par_sct= FALSE
    group4_options.add_argument('--par_sct', action='store_true', help='Whether SCTransform was used while processing the Seurat object. For most use-cases, this parameter can remain deactivated. ')

    ##rate_nExp: the doublet rate according to the number of cells
    # par_rate_nExp=0.076
    group4_options.add_argument('--par_rate_nExp', type=float, help="Global expected doublet rate (e.g., for a 5% expected doublet rate, write 0.05). Expected doublet rates can be informed by the number of recovered cells. For more information, see the 10X documentation. ")

    ## Expected doublet rate for each sample. First list sample IDs, then list the corresponding expected doublet rate for each sample depending on the number of recovered or loaded cells. Sample names should be the same ones used in the library.csv file used for Step 1.
    # par_sample_names= c("Control1","Parkinson1")
    group4_options.add_argument('--par_sample_names_step4', type=str, help='A list of sample names for each sample in the experiment. Sample names should be the same as those defined in --par_rename_samples in Step 1 and must correspond to the sample-specific expected doublet rates defined in the parameter below. ')
    # par_expected_doublet_rate= c(0.05,0.05)
    group4_options.add_argument('--par_expected_doublet_rate', type=str, help="A vector of expected doublet rates for each sample (e.g., for a 5% sample-specific expected doublet rate, write 0.05). The expected doublet rates for each sample should be listed in the same order as the sample names in the above parameter. Make sure to have as many expected doublet rates listed as you have samples.")


    # Step 4 HTO Track: Demultiplexing
    #####################################

    # --msd  = You can get the hashtag labels by running the following code (HTO Step 4).
    group4_options.add_argument('--msd', action='store_true', help="Before demultiplexing, activate this parameter and leave all other parameters empty to retrieve the barcode labels used in the analysis. If the barcode labels are already known, leave this parameter deactivated and proceed with the remaining parameters.")

    # par_save_RNA= "yes"
    group4_options.add_argument('--par_save_RNA_step4_hto', action='store_true', help='Whether or not to export an RNA expression matrix.')
    # par_save_metadata= "yes"
    group4_options.add_argument('--par_save_metadata_step4_hto', action='store_true', help='Whether or not to save a metadata dataframe.')

    # If you already have a processed Seurat RDS object, and did not perform Step 3 of scRNAbox use this to add the path to the directory containing you Seurat object(s).
    # Make sure that no other files/objects are present in the directory besides Seurat RDS objects.
    # Uncomment the line to activate the parameter
    ############################################################################
    #par_seurat_object= "/path/to/directory/containing/seurat/object"
    group4_options.add_argument('--par_seurat_object_step4_hto', type=is_dir, help='If you performed Step 3 of the pipeline, keep this parameter empty. Otherwise, if users have an existing Seurat object(s), they may provide the path to a directory that contains the Seurat object(s) to initiate the pipeline at Step 4.')

    ## Normalization method
    # par_normalization.method= "LogNormalize"
    group4_options.add_argument('--par_normalization_method_step4', type=str, help="Method for normalizing the HTO assay (e.g., CLR).")

    ## Scale factor
    # par_scale.factor= 10000
    group4_options.add_argument('--par_scale_factor_step4', type=int, help='Scale factor for scaling the HTO assay (e.g., 10000).')


    # Parameters for integration
    ############################################################################
    ## Method for detecting top variable features. vst, mean.var.plot (mvp), dispersion (disp)
    # par_selection.method= "vst"
    group4_options.add_argument('--par_selection_method_step4', type=str, help='Method for selecting highly variable features in the HTO assay (e.g., vst).')

    ## Number of features to select as top variable features for integration
    # par_nfeatures= 2500
    group4_options.add_argument('--par_nfeatures_step4', type=int, help='Number of features to select as top variable features for the HTO assay (e.g., 5). This value is dependent on the number of sample specific barcodes used in the experiment.')

    group4_options.add_argument('--par_dims_umap', type=int, help='Number of dimensions to use as input features for uniform manifold approximation and projection of HTO assay (e.g., 5).')

    group4_options.add_argument('--par_n_neighbor', type=int, help='Number of neighboring points used in local approximations of manifold structure (e.g., 65).')

    group4_options.add_argument('--par_dimensionality_reduction', action='store_true', help='Whether or not to perform linear dimensionality reduction on the HTO assay.')

    ## Total Number of PCs to compute and store for RunPCA
    # par_npcs_pca= 30
    group4_options.add_argument('--par_npcs_pca_step4', type=int, help='The number of principal components to compute and store for principal component analysis of HTO assay (e.g., 30).')

    # Parameters for doublet detection and removal (optional)
    ############################################################################
    ## If you want to remove predicted doublets from downstream analyses set the following to "yes"
    ## If you want to keep predicted doublets for further analysis set the following to "no"
    # par_dropDN= "yes"
    group4_options.add_argument('--par_dropDN_hto', action='store_true', help='Whether or not to remove predicted doublets and negatives from downstream analyses.')

    group4_options.add_argument('--par_label_dropDN', type=str, help='Labels used to identify doublet and negative droplets (e.g., Doublet, Negative).')

    group4_options.add_argument('--par_quantile', type=float, help='The quantile to use for droplet classification using MULTIseqDemux (e.g., 0.9).')

    group4_options.add_argument('--par_autoThresh', action='store_true', help='Whether or not to perform automated threshold finding to define the best quantile for droplet classification using MULTIseqDemux. For most use cases, this parameter should be activated.')

    group4_options.add_argument('--par_maxiter', type=int, help='Maximum number of iterations to use if --par_autoThresh is activated (e.g., 5).')

    group4_options.add_argument('--par_RidgePlot_ncol', type=int, help='Number of columns used to display RidgePlots, which visualizes the enrichment of barcode labels across samples (e.g., 3).')

    group4_options.add_argument('--par_old_antibody_label', type=str, help='If you wish to rename the barcode labels, list the new labels corresponding to the old labels listed in the parameter above.')

    # id: Step 5
    # Step 5: Sample integration
    ###############################

    group5_options = parser.add_argument_group('Step 5')

    # If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
    ############################################################################
    # par_save_RNA= "yes"
    group5_options.add_argument('--par_save_RNA_step5', action='store_true', help='Whether or not to export an RNA expression matrix.')
    # par_save_metadata= "yes"
    group5_options.add_argument('--par_save_metadata_step5', action='store_true', help='Whether or not to save a metadata dataframe.')

    # If you already have a processed Seurat RDS object(s), and did not perform Step 4 of scRNAbox use this to add the path to the directory containing you Seurat object(s).
    # Make sure that no other files/objects are present in the directory besides Seurat RDS objects.
    # Uncomment the line to activate the parameter
    ############################################################################
    #par_seurat_object= "/path/to/directory/containing/seurat/objects"
    group5_options.add_argument('--par_seurat_object_step5', type=is_dir, help='If you performed Step 4 of the pipeline, keep this parameter empty. Otherwise, if users have an existing Seurat object(s), they may provide the path to a directory that contains the Seurat object(s) to initiate the pipeline at Step 5.')

    # If you only have one Seurat object and want to skip integration set the following to "yes"
    ############################################################################
    # par_one_seurat= "no"
    group5_options.add_argument('--par_one_seurat', action='store_true', help='Whether or not the experiment has only one sample. If this parameter is activated (only one sample), par_integrate_seurat and par_merge_seurat should be deactivated.')

    # If you have multiple Seurat objects, choose whether you want to integrate or merge the objects
    ############################################################################
    ## Integrate Seurat objects
    # par_integrate_seurat= "yes"
    group5_options.add_argument('--par_integrate_seurat', action='store_true', help='Whether or not to integrate the samples. If activated, par_merge_seurat must be deactivated.')

    ## Merge Seurat objects
    # par_merge_seurat= "no"
    group5_options.add_argument('--par_merge_seurat', action='store_true', help='Whether or not to merge the samples. If activated, par_integrate_seurat must be deactivated.')

    # Parameters for normalization and scaling
    # Even if you opt to skip integration, adjust the following parameters
    ############################################################################
    ## Assay to perform normalization and scaling (prior to integration). For most use cases this will be RNA
    # par_DefaultAssay= "RNA"
    group5_options.add_argument('--par_DefaultAssay', type=str, help='The assay to perform normalization, scaling, and linear dimensional reduction on (e.g., RNA)')

    ## Normalization method
    # par_normalization.method= "LogNormalize"
    group5_options.add_argument('--par_normalization_method_step5', type=str, help='Method to use for normalization (e.g., LogNormalize).')

    ## Scale factor
    # par_scale.factor= 10000
    group5_options.add_argument('--par_scale_factor_step5', type=int, help='Scale factor for scaling the data (e.g., 10000).')

    # Parameters for integration
    ############################################################################
    ## Method for detecting top variable features. vst, mean.var.plot (mvp), dispersion (disp)
    # par_selection.method= "vst"
    group5_options.add_argument('--par_selection_method_step5', type=str, help='Method for selecting highly variable genes (e.g., vst).')

    ## Number of features to select as top variable features for integration
    # par_nfeatures= 2500
    group5_options.add_argument('--par_nfeatures_step5', type=int, help='Number of highly variable genes to select (e.g., 2500).')

    ## Which dimensions to use from the CCA to specify the neighbour search space
    # par_FindIntegrationAnchors_dim= 25
    group5_options.add_argument('--par_FindIntegrationAnchors_dim', type=int, help='Which dimensions to use from the canonical correlation analysis (CCA) to specify the neighbor search space (e.g., 25).')

    # Parameters for linear dimensional reduction
    # even if you opt to skip integration, adjust the following parameters
    ############################################################################
    ## Total Number of PCs to compute and store for RunPCA
    # par_RunPCA_npcs= 30
    group5_options.add_argument('--par_RunPCA_npcs', type=int, help='Total number of principal components to compute and store for principal component analysis (e.g., 30).')

    ## Which dimensions to use as input features for RunUMAP
    # par_RunUMAP_dims= 25
    group5_options.add_argument('--par_RunUMAP_dims_step5', type=int, help='Number of dimensions to use as input features for uniform manifold approximation and projection (e.g., 25).')

    ## The number of neighbouring points used in local approximations of manifold structure.
    # par_RunUMAP_n.neighbors= 45
    group5_options.add_argument('--par_RunUMAP_n_neighbors_step5', type=int, help='Number of neighboring points used in local approximations of manifold structure (e.g., 45).')

    ## Whether or not to perform JackStraw computation. This computation takes a long time.
    # par_compute_jackstraw= "no"
    group5_options.add_argument('--par_compute_jackstraw', action='store_true', help='Whether or not to perform JackStraw computation to identify significant principal components.')

    # id: Step 6
    # Step 6: Clustering
    #######################

    group6_options = parser.add_argument_group('Step 6')

    # If you want to save an RNA expression matrix and metadata dataframe set the following to "yes"
    ############################################################################
    # par_save_RNA= "yes"
    group6_options.add_argument('--par_save_RNA_step6', action='store_true', help='Whether or not to export an RNA expression matrix.')
    # par_save_metadata= "yes"
    group6_options.add_argument('--par_save_metadata_step6', action='store_true', help='Whether or not to export a metadata dataframe.')

    # If you already have a processed Seurat RDS object, and did not perform Step 5 of scRNAbox
    # use this parameter to add the path to the directory containing your Seurat object.
    # Uncomment the line to activate the parameter. Note you can only have one Seurat object at this point.
    ############################################################################
    #par_seurat_object= "/path/to/seurat.rds"
    group6_options.add_argument('--par_seurat_object_step6', type=is_dir, help='If you performed Step 5 of the pipeline, keep this parameter empty. Otherwise, if users have an existing Seurat object(s), they may provide the path to a directory that contains the Seurat object(s) to initiate the pipeline at Step 6.')

    # If you skipped integration in step 5, set the following to "yes".
    # If you performed integration, keep the default as no.
    # This will keep the default assay as "RNA"
    # If you wish to cluster on with "RNA" and not the integrated data you can do so by setting this to "yes".
    ############################################################################
    # par_skip_integration= "no"
    group6_options.add_argument('--par_skip_integration', action='store_true', help='Whether or not the user skipped integration in Step 5.')

    # Clustering parameters
    ############################################################################
    ## Number of PC to use as input to find neighbours. Can be informed by the elbow and jackstraw plots produced in Step 5.
    # par_FindNeighbors_dims= 25
    group6_options.add_argument('--par_FindNeighbors_dims', type=int, help='Number of principal components used as input to identify neighbours (e.g., 25). Can be informed by the elbow and Jackstraw plots produced in Step 5.')

    ## Number of PCs to use as input for UMAP. Can be informed by the elbow and jackstraw plots produced in Step 5.
    # par_RunUMAP_dims= 25
    group6_options.add_argument('--par_RunUMAP_dims_step6', type=int, help='Number of dimensions to use as input features for uniform manifold approximation and projection (e.g., 25).')

    ## Defines k for the k-nearest neighbour algorithm. The number of neighbours to include when constructing the SNN graph.
    # par_FindNeighbors_k.param= 45
    group6_options.add_argument('--par_FindNeighbors_k_param', type=int, help='Defines k for the k-nearest neighbor algorithm (e.g., 45).')

    ## Sets the cutoff for acceptable Jaccard index when computing the neighbourhood overlap for the SNN construction
    # par_FindNeighbors_prune.SNN= 1/15
    group6_options.add_argument('--par_FindNeighbors_prune_SNN', type=float, help='Sets the cut off for acceptable Jaccard index when computing the neighborhood overlap for the shared nearest-neighbour construction (e.g., 0.067).')

    ## Value of the clustering resolution parameter. You may provide multiple resolution values. Use a value above 1.0 if you want to obtain a larger number of smaller communities.
    # par_FindClusters_resolution= c(0, 0.05, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
    group6_options.add_argument('--par_FindClusters_resolution', type=str, help='Value of the clustering resolution parameter. It is recommended to cluster at various resolutions. (e.g., 0, 0.05, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0).')

    # Adjusted Rand Index for clustering resolutions
    ############################################################################
    ## If you want to compute ARI, set the following to "yes". This takes a long time.
    # par_compute_ARI= "yes"
    group6_options.add_argument('--par_compute_ARI', action='store_true', help='Whether or not you want to compute the Adjusted Rand Index (ARI) to evaluate the stability of clusters across user-defined clustering resolutions.')

    ## Number of repetitions of generated clusters to calculate ARI.
    # par_RI_reps= 25
    group6_options.add_argument('--par_RI_reps', type=int, help='Number of iterations for clustering the data at a given resolution to compute the Adjusted Rand Index (e.g., 25).')

    # id: Step 7
    # Step 7: Cell type annotation
    #################################

    group7_options = parser.add_argument_group('Step 7')

    # --markergsea  = Identify marker genes for each cluster and run marker gene set enrichment analysis (GSEA) using EnrichR libraries (Step 7).
    group7_options.add_argument('--markergsea', action='store_true', help='Identify marker genes for each cluster and GSEA using EnrichR libraries.')

    # --knownmarkers  = Profile the individual or aggregated expression of known marker genes.
    group7_options.add_argument('--knownmarkers', action='store_true', help='Profile the individual or aggregated expression of known marker genes.')

    # --referenceannotation  = Generate annotation predictions based on the annotations of a reference Seurat object (Step 7).
    group7_options.add_argument('--referenceannotation', action='store_true', help='Generate annotation predictions based on the annotations of a reference Seurat object.')

    # --annotate  = Add clustering annotations to Seurat object metadata (Step 7).
    group7_options.add_argument('--annotate', action='store_true', help='Add clustering annotations to Seurat object metadata.')


    # If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
    ############################################################################
    # par_save_RNA= "yes"
    group7_options.add_argument('--par_save_RNA_step7', action='store_true', help='Whether or not to export an RNA expression matrix.')
    # par_save_metadata= "yes"
    group7_options.add_argument('--par_save_metadata_step7', action='store_true', help='Whether or not to export a metadata dataframe.')

    # If you already have a processed Seurat RDS object and did not perform Step 6 of scRNAbox use this parameter to add the path to the directory containing your Seurat object.
    # Uncomment the line to activate the parameter
    # Your Seurat object must already have clusters
    ############################################################################
    #par_seurat_object= "/path/to/seurat.rds"
    group7_options.add_argument('--par_seurat_object_step7', type=is_dir, help='If you performed Step 6 of the pipeline, keep this parameter empty. Otherwise, if users have an existing Seurat object(s), they may provide the path to a directory that contains the Seurat object(s) to initiate the pipeline at Step 7.')

    # General parameters for cluster annotation
    ############################################################################
    ## The cluster resolution that you want to use. If you skipped integration, use par_level_cluster="RNA_snn_res.0.75", for example, if you want to proceed with a clustering resolution of 0.75
    ## This parameter can also be set to your annotated cell type names if you want to check expression levels or find markers within the annotated groups.
    # par_level_cluster= "integrated_snn_res.0.75"
    group7_options.add_argument('--par_level_cluster', type=str, help='The cluster resolution that you want to annotate (e.g., integrated_snn_res.0.75).')

    # Tool 1: Marker GSEA
    ############################################################################
    ## Identify cluster specific markers
    # par_run_find_marker= "yes"
    group7_options.add_argument('--par_run_find_marker', action='store_true', help='Whether or not to identify marker genes for each cluster.')

    ## Number of top markers based on avg_log2FC
    ## This is the number of markers to include on a heatmap for visualization.
    # par_top_sel= 5
    group7_options.add_argument('--par_top_sel', type=int, help='Number of top markers to display in heatmap based on avg_log2FC (e.g., 5).')

    ## Run EnrichR GSEA on cluster-specific markers. This step should follow the identification of cluster-specific markers.  Additionally, this step can only be run if your HPC allows internet access.
    ## Note that code is provided to run this locally if your HPC cannot access the internet.
    # par_run_enrichR= "no"
    group7_options.add_argument('--par_run_enrichR', action='store_true', help='Whether or not to run gene set enrichment analysis (GSEA) on the marker genes for each cluster using the EnrichR tools. Note that the HPC must have access to the internet to run GSEA.')

    ## Character vector of EnrichR databases to search for enrichment
    ## This is only needed if you have internet connection and are running GSEA
    # par_db= c("Descartes_Cell_Types_and_Tissue_2021","CellMarker_Augmented_2021","Azimuth_Cell_Types_2021")
    group7_options.add_argument('--par_db', type=str, help='Character vector of EnrichR databases that define cell types. The top marker genes for each cluster will be tested for enrichment across these databases (e.g., Descartes_Cell_Types_and_Tissue_2021, CellMarker_Augmented_2021, Azimuth_Cell_Types_2021).')

    # Tool 2: Visualize the expression of known marker genes
    ######################################################################
    ## Perform module score computation
    # par_run_module_score= "yes"
    group7_options.add_argument('--par_run_module_score', action='store_true', help='Whether or not to compute module score to assess the aggregated expression of known cell type marker gene sets. ')

    ## Define the path to a csv file containing the genes sets for module score
    # par_module_score= "/path/to/gene_sets.csv"
    group7_options.add_argument('--par_module_score', type=is_file, help='Path to the csv file containing the gene sets for the module score (see ScRNAbox documentation for an example .csv file).')

    ## Visualize markers (dot plot, violin plot, feature plot)
    # par_run_visualize_markers= "yes"
    group7_options.add_argument('--par_run_visualize_markers', action='store_true', help='IWhether or not to visualize the individual expression of known marker genes.')

    ## List of markers that you want to visualize (dot plot, violin plot, feature plot)
    ## Be sure to use the official gene names
    # par_select_features_list= c("gene1", "gene2", "gene3")
    group7_options.add_argument('--par_select_features_list', type=str, help='List of marker genes whose expression will be visualized individually.')

    ## If you want to define multiple lists of markers to visualize, you can do so with a csv file. The header should contain the list names and all features belonging to the same list should be in the same column. Uncomment the below parameter and enter the location of the csv file. This can be the same csv file used for module score.
    #par_select_features_csv= "/path/to/visualize_features.csv"
    group7_options.add_argument('--par_select_features_csv', type=is_file, help='To define multiple feature lists for individual visualization, provide the path to a CSV file where each column represents a list. The column names should describe the list names, and all features belonging to the same list should be placed in the corresponding column.')

    # Tool 3: Reference-based annotation parameters
    ######################################################################
    ## Seurat RDS object to use as the reference
    # par_reference= "/path/to/reference_seurat_object.rds"
    group7_options.add_argument('--par_reference', type=is_file, help='Path to the reference Seurat object.')

    ## Define an arbitrary name for the reference object. This will be used to name the metadata slot.
    # par_reference_name= "reference"
    group7_options.add_argument('--par_reference_name', type=str, help='An arbitrary name for the reference object (e.g., Reference1). This will be used to name the metadata slot.')

    ## Name of a metadata column in the reference Seurat object that contains cell type annotations
    # par_level_celltype= "Cell_Type"
    group7_options.add_argument('--par_level_celltype', type=str, help='The name of the metadata column in the reference Seurat object that defines cell type labels to be transferred to the query object (e.g., Cell_Type).')

    ## How many dimensions to use to find transfer anchors between query and reference dataset
    # par_FindTransferAnchors_dim= 50
    group7_options.add_argument('--par_FindTransferAnchors_dim', type=int, help='Number of dimensions used to find transfer anchors between the reference and query Seurat objects (e.g., 50).')

    ## This will increase your RAM usage so set this number mindfully
    # par_futureglobalsmaxSize= 60000 * 1024^2
    group7_options.add_argument('--par_futureglobalsmaxSize', type=int, help='The maximum allowed size of global objects that can be exported to parallel workers (e.g., 62914560000). This will increase your RAM usage so set this number mindfully.')

    # Annotate parameters
    # Annotations from each iteration will be added to the Step 7 Seurat object
    ######################################################################
    ## The clustering resolution to annotate
    # par_annotate_resolution= "integrated_snn_res.0.75"
    group7_options.add_argument('--par_annotate_resolution', type=str, help='Which clustering resolution to annotate (e.g., integrated_snn_res.0.75).')

    ## the name of the metadata slot under which the cluster labels will be stored.
    # par_name_metadata= "Celltypes1"
    group7_options.add_argument('--par_name_metadata', type=str, help='The name of the metadata slot that will contain the annotations (e.g., Cell_Types1).')

    ## A list of cluster labels. Make sure you have as the same number of labels as clusters at the defined clustering resolution. Please do not use "_" when naming cell types.
    # par_annotate_labels= c("Annot1", "Annot2", "Annot3")
    group7_options.add_argument('--par_annotate_labels', type=str, help='A list of cluster labels, with one label per cluster at the specified clustering resolution. Avoid using underscores ("_") in the labels.')

    # id: Step 8
    # Step 8: Differential gene expression analysis
    ##################################################

    group8_options = parser.add_argument_group('Step 8')

    # --addmeta  = Add metadata columns to the Seurat object (Step 8).
    group8_options.add_argument('--addmeta', action='store_true', help='Add metadata information to the Seurat object to facilitate DGE contrasts.')

    # --rundge  = Perform differential gene expression contrasts (Step 8).
    group8_options.add_argument('--rundge', action='store_true', help='Perform DGE contrasts.')

    # If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
    ############################################################################
    # par_save_RNA= "yes"
    group8_options.add_argument('--par_save_RNA_step8', action='store_true', help='Whether or not to export an RNA expression matrix.')
    # par_save_metadata= "yes"
    group8_options.add_argument('--par_save_metadata_step8', action='store_true', help='Whether or not to save a metadata dataframe.')

    # If you already have a processed Seurat RDS object, and did not perform Step 7 of scRNAbox use this parameter to add the path to the Seurat object.
    # Uncomment the line to activate the parameter
    ############################################################################
    #par_seurat_object= "/path/to/seurat.rds"
    group8_options.add_argument('--par_seurat_object_step8', type=is_dir, help='If you performed Step 7 of the pipeline, keep this parameter empty. Otherwise, if users have existing Seurat object(s), they may provide the path to a directory that contains the Seurat object(s) to initiate the pipeline at Step 8.')

    # Add metadata parameters
    ############################################################################
    ## Enter the path to a csv file describing new metadata that should be added to the Seurat object to facilitate DEG analysis.
    ## The rows should contain the data to add in the order of the levels of "Sample_ID" or the metadata slot you will use to define your samples. The column names should be the desired name of the metadata slot to add.
    # par_metadata= "path/to/metadata.csv"
    group8_options.add_argument('--par_metadata', type=is_file, help='Path to the csv file containing the new metadata to be added to the Seurat object.')

    ## Define the column from the Seurat object that you want to use to add the new metadata
    # par_merge_meta= "Sample_ID"
    group8_options.add_argument('--par_merge_meta', type=str, help='The column from the Seurat metadata that will be used to merge the new metadata (e.g., orig.ident). This column must also exist in the csv file containing new metadata defined in the above parameter.')

    # Run DGE parameters
    # Choose which differential gene expression (DGE) methods you want to use for this submission
    # Be sure adjust the appropriate txt design files.
    ############################################################################
    ## Perform cell-base DGE with all cells
    # par_run_cell_based_all_cells= "yes"
    group8_options.add_argument('--par_run_cell_based_all_cells', action='store_true', help='Whether or not to compute cell-based DGE across all cell types together.')

    group8_options.add_argument('--contrast_cell_based_all_cells', type=is_file, help='Path to the contrast matrix for cell-based DGE with all cells. For an example contrast matrix see the scRNAbox documentation (https://neurobioinfo.github.io/scrnabox/site/Step8/).')

    ## Perform cell-based on each cell type group
    # par_run_cell_based_celltype_groups= "yes"
    group8_options.add_argument('--par_run_cell_based_celltype_groups', action='store_true', help='Whether or not to compute cell-based DGE within each cell type group independently. ')

    group8_options.add_argument('--contrast_cell_based_celltype_groups', type=is_file, help='Path to the contrast matrix for cell-based DGE with cell type groups. For an example contrast matrix see the scRNAbox documentation (https://neurobioinfo.github.io/scrnabox/site/Step8/).')

    # Cell-replicate DGE parameters
    ############################################################################
    ## Which statistical method to use when computing DGE using individual cells as replicates
    # par_statistical_method= "MAST"
    group8_options.add_argument('--par_statistical_method', type=str, help='Which statistical framework to use for computing cell-based DGE (e.g., MAST).')

    ## Perform Sample-based DGE with all cells (pseudobulk)
    # par_run_sample_based_all_cells= "yes"
    group8_options.add_argument('--par_run_sample_based_all_cells', action='store_true', help='Whether or not to compute cell-based DGE across all cell types together.')

    # contrast_sample_based_all_cells
    group8_options.add_argument('--contrast_sample_based_all_cells', type=is_file, help='Path to the contrast matrix for sample-based DGE with all cells. For an example contrast matrix see the scRNAbox documentation (https://neurobioinfo.github.io/scrnabox/site/Step8/).')

    ## Perform Sample-based DGE on each cell type group (pseudobulk)
    # par_run_sample_based_celltype_groups= "yes"
    group8_options.add_argument('--par_run_sample_based_celltype_groups', action='store_true', help='Whether or not to compute sample-based DGE within each cell type group independently.')

    # contrast_sample_based_celltype_groups
    group8_options.add_argument('--contrast_sample_based_celltype_groups', type=is_file, help='Path to the contrast matrix for sample-based DGE with cell type groups. For an example contrast matrix see the scRNAbox documentation (https://neurobioinfo.github.io/scrnabox/site/Step8/).')

    # General options
    ####################

    # --seulist  = You can directly call the list of Seurat objects to the pipeline.
    general_options.add_argument('--seulist', action='store_true', help='You can directly call the list of Seurat objects to the pipeline')

    # --rcheck  = You can identify which libraries are not installed.
    general_options.add_argument('--rcheck', action='store_true', help='You can identify which libraries are not installed')


    #################################
    # Validate the input parameters #
    #################################

    args = parser.parse_args()

    # check if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    if args.dir == None or args.steps == None or args.output_dir == None:
        print_log("Please provide the work directory, the steps to execute, and the output directory.")
        sys.exit(2)

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
    # Generate the step*.txt files based on the user input and the default parameters #
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

    # Create working directory
    # Check if the directory already exists
    current_dir = os.getcwd()
    dir_fullpath = current_dir + "/" + args.dir
    if os.path.exists(dir_fullpath):
        print_log("Directory already exists, use this directory to continue the analysis.")
    else:
        os.makedirs(dir_fullpath)
        # Create job_info directory
        os.makedirs(dir_fullpath + "/job_info")
        # Create configs directory
        os.makedirs(dir_fullpath + "/job_info/configs")
        # Create parameters directory
        os.makedirs(dir_fullpath + "/job_info/parameters")
        # Create .tmp directory
        os.makedirs(dir_fullpath + "/job_info/.tmp")
        # Create logs directory
        os.makedirs(dir_fullpath + "/job_info/logs")
        # touch summary_report.txt
        open(dir_fullpath + "/job_info/summary_report.txt", "w").close()

    #################################
    # Generate scrnabox_config.ini  #
    #################################

    config_file_path = dir_fullpath + "/job_info/configs/scrnabox_config.ini"
    print_log("Generating scrnabox_config.ini file: " + config_file_path)
    if os.path.exists(config_file_path):
        # Copy file with timestamp with standar mv command
        subprocess.run(["mv", config_file_path, config_file_path + "." + str(int(time.time()))])

    config_file = open(config_file_path, "w")

    # if args.R_path:
    #     config_file.write("R_LIB_PATH=" + args.R_path + "\n")
    # else:
    #     # Check for R_LIB_PATH in the environment variables
    #     if "R_LIB_PATH" in os.environ:
    #         config_file.write("R_LIB_PATH=" + os.environ["R_LIB_PATH"] + "\n")
    #     else:
    #         print_log("Please provide the path to the R library directory using the --R_path argument or set the R_LIB_PATH environment variable.")
    #         sys.exit(2)

    if args.method:
        config_file.write("SCRNA_METHOD=" + args.method + "\n")

    if args.job_mode == "slurm" or args.job_mode == "pbs":
        config_file.write("JOB_MODE=" + args.job_mode + "\n")
    else:
        config_file.write("JOB_MODE=local\n")

    # Add THREADS_ARRAY info for each step
    # THREADS_ARRAY["step_2"]=10
    # THREADS_ARRAY["step_3"]=10
    # THREADS_ARRAY["step_4"]=10
    # THREADS_ARRAY["step_5"]=10
    # THREADS_ARRAY["step_6"]=10
    # THREADS_ARRAY["step_7_markergsea"]=4
    # THREADS_ARRAY["step_7_knownmarkers"]=4
    # THREADS_ARRAY["step_7_referenceannotation"]=10
    # THREADS_ARRAY["step_7_annotate"]=4
    # THREADS_ARRAY["step_8_addmeta"]=4
    # THREADS_ARRAY["step_8_rundge"]=4
    # THREADS_ARRAY["integrate"]=10
    threads_array_steps = {
        "2": 10,
        "3": 10,
        "4": 10,
        "5": 10,
        "6": 10,
        "7_markergsea": 4,
        "7_knownmarkers": 4,
        "7_referenceannotation": 10,
        "7_annotate": 4,
        "8_addmeta": 4,
        "8_rundge": 4,
        "integrate": 10
    }


    for step in args.steps.split(','):

        step = step.strip()
        print_log("Configuring THREADS_ARRAY for step " + step)
        if step == '5' and args.par_integrate_seurat:
            config_file.write("THREADS_ARRAY[\"integrate\"]=" + str(threads_array_steps["integrate"]) + "\n")
        elif step == '7':
            if args.markergsea:
                config_file.write("THREADS_ARRAY[\"step_7_markergsea\"]="          + str(threads_array_steps["7_markergsea"])          + "\n")
            if args.knownmarkers:
                config_file.write("THREADS_ARRAY[\"step_7_knownmarkers\"]="        + str(threads_array_steps["7_knownmarkers"])        + "\n")
            if args.referenceannotation:
                config_file.write("THREADS_ARRAY[\"step_7_referenceannotation\"]=" + str(threads_array_steps["7_referenceannotation"]) + "\n")
            if args.annotate:
                config_file.write("THREADS_ARRAY[\"step_7_annotate\"]="            + str(threads_array_steps["7_annotate"])            + "\n")
        elif step == '8':
            if args.addmeta:
                config_file.write("THREADS_ARRAY[\"step_8_addmeta\"]="             + str(threads_array_steps["8_addmeta"])             + "\n")
            if args.rundge:
                config_file.write("THREADS_ARRAY[\"step_8_rundge\"]="              + str(threads_array_steps["8_rundge"])              + "\n")
        else:
            if step in threads_array_steps:
                config_file.write("THREADS_ARRAY[\"step_" + step + "\"]="          + str(threads_array_steps[step])                         + "\n")
            else:
                print_log("Warning: Step " + step + " is not recognized. Skipping THREADS_ARRAY entry for this step.")

    # close the file
    config_file.close()

    # Check if step 7 is selected but none of the required options are set
    if '7' in args.steps.split(','):
        if not any([args.markergsea, args.knownmarkers, args.referenceannotation, args.annotate]):
            print_log("Error: When selecting step 7, at least one of these options must be set: --markergsea, --knownmarkers, --referenceannotation, --annotate")
            sys.exit(2)

    #################################
    # Generate the step1_*.txt file #
    #################################

    step1_file_path = dir_fullpath + "/job_info/parameters/step1_par.txt"
    step1_file = open(step1_file_path, "w")

    if args.par_automated_library_prep_scrna or args.par_automated_library_prep_hto:
        step1_file.write("par_automated_library_prep=\"yes\"\n")
    else:
        step1_file.write("par_automated_library_prep=\"no\"\n")

    if args.par_fastq_directory_scrna or args.par_fastq_directory_hto:
        fastq_dir = args.par_fastq_directory_scrna if args.par_fastq_directory_scrna else args.par_fastq_directory_hto
        full_path_fastq = current_dir + "/" + fastq_dir
        step1_file.write("par_fastq_directory=\"" + full_path_fastq + "\"\n")

    if args.par_sample_names_step1:
        # split the list of sample names ',' and add double quotes join the list with
        split_par_sample_names_step1 = args.par_sample_names_step1.split(',')
        par_sample_names_step1 = ",".join([f"'{item}'" for item in split_par_sample_names_step1])
        step1_file.write("par_sample_names= c(" + par_sample_names_step1 + ")\n")

    if args.par_rename_samples_scrna:
        step1_file.write("par_rename_samples=\"yes\"\n")
    else:
        step1_file.write("par_rename_samples=\"no\"\n")

    if args.par_new_sample_names:
        split_par_new_sample_names = args.par_new_sample_names.split(',')
        par_new_sample_names = ",".join([f"'{item}'" for item in split_par_new_sample_names])
        step1_file.write("par_new_sample_names= c(" + par_new_sample_names + ")\n")

    if args.par_paired_end_seq_scrna or args.par_paired_end_seq_hto:
        step1_file.write("par_paired_end_seq=TRUE\n")
    else:
        step1_file.write("par_paired_end_seq=FALSE\n")

    ref_dir = args.par_ref_dir_grch_scrna if args.par_ref_dir_grch_scrna else args.par_ref_dir_grch_hto
    par_ref_dir_grch = args.par_ref_dir_grch_scrna if args.par_ref_dir_grch_scrna else args.par_ref_dir_grch_hto
    if args.ref_data_dir and ref_dir:
        step1_file.write("par_ref_dir_grch=\"" + args.ref_data_dir + '/' + ref_dir + "\"\n")
    elif par_ref_dir_grch:
        step1_file.write("par_ref_dir_grch=\"" + ref_dir + "\"\n")

    if args.par_r1_length_scrna or args.par_r1_length_hto:
        r1_length = args.par_r1_length_scrna if args.par_r1_length_scrna else args.par_r1_length_hto
        step1_file.write("par_r1_length=\"" + str(r1_length) + "\"\n")

    if args.par_r2_length_scrna or args.par_r2_length_hto:
        r2_length = args.par_r2_length_scrna if args.par_r2_length_scrna else args.par_r2_length_hto
        step1_file.write("par_r2_length=\"" + str(r2_length) + "\"\n")

    step1_file.write("par_mempercode=30\n")

    if args.par_include_introns_scrna or args.par_include_introns_hto:
        step1_file.write("par_include_introns=TRUE\n")
    else:
        step1_file.write("par_include_introns=FALSE\n")

    if args.par_no_target_umi_filter_scrna or args.par_no_target_umi_filter_hto:
        step1_file.write("par_no_target_umi_filter=TRUE\n")
    else:
        step1_file.write("par_no_target_umi_filter=FALSE\n")

    if args.par_expect_cells_scrna or args.par_expect_cells_hto:
        expect_cells = args.par_expect_cells_scrna if args.par_expect_cells_scrna else args.par_expect_cells_hto
        step1_file.write("par_expect_cells=\"" + str(expect_cells) + "\"\n")

    if args.par_force_cells_scrna or args.par_force_cells_hto:
        force_cells = args.par_force_cells_scrna if args.par_force_cells_scrna else args.par_force_cells_hto
        step1_file.write("par_force_cells=\"" + str(force_cells) + "\"\n")

    if args.par_no_bam_scrna or args.par_no_bam_hto:
        step1_file.write("par_no_bam=TRUE\n")
    else:
        step1_file.write("par_no_bam=FALSE\n")

    if args.par_RNA_run_names:
        step1_file.write("par_RNA_run_names=\"" + args.par_RNA_run_names + "\"\n")

    if args.par_HTO_run_names:
        step1_file.write("par_HTO_run_names=\"" + args.par_HTO_run_names + "\"\n")

    if args.par_seq_run_names:
        step1_file.write("par_seq_run_names=\"" + args.par_seq_run_names + "\"\n")

    if args.id:
        step1_file.write("id=\"" + args.id + "\"\n")

    if args.name:
        step1_file.write("name=\"" + args.name + "\"\n")

    if args.read:
        step1_file.write("read=\"" + args.read + "\"\n")

    if args.pattern:
        step1_file.write("pattern=\"" + args.pattern + "\"\n")

    if args.sequence:
        step1_file.write("sequence=\"" + args.sequence + "\"\n")

    step1_file.close()

    #################################
    # Generate the step2_*.txt file #
    #################################

    step2_file_path = dir_fullpath + "/job_info/parameters/step2_par.txt"
    step2_file = open(step2_file_path, "w")

    if args.par_save_RNA_step2:
        step2_file.write("par_save_RNA=\"yes\"\n")
    else:
        step2_file.write("par_save_RNA=\"no\"\n")

    if args.par_save_metadata_step2:
        step2_file.write("par_save_metadata=\"yes\"\n")
    else:
        step2_file.write("par_save_metadata=\"no\"\n")

    if args.par_ambient_RNA:
        step2_file.write("par_ambient_RNA=\"yes\"\n")
    else:
        step2_file.write("par_ambient_RNA=\"no\"\n")

    if args.par_min_cells_L:
        step2_file.write("par_min.cells_L=" + str(args.par_min_cells_L) + "\n")

    if args.par_normalization_method_step2:
        step2_file.write("par_normalization.method=\"" + args.par_normalization_method_step2 + "\"\n")

    if args.par_scale_factor_step2:
        step2_file.write("par_scale.factor=" + str(args.par_scale_factor_step2) + "\n")

    if args.par_selection_method_step2:
        step2_file.write("par_selection.method=\"" + args.par_selection_method_step2 + "\"\n")

    if args.par_nfeatures_step2:
        step2_file.write("par_nfeatures=" + str(args.par_nfeatures_step2) + "\n")

    step2_file.close()

    #################################
    # Generate the step_3.txt file #
    #################################

    step3_file_path = dir_fullpath + "/job_info/parameters/step3_par.txt"
    step3_file = open(step3_file_path, "w")

    if args.par_save_RNA_step3:
        step3_file.write("par_save_RNA=\"yes\"\n")
    else:
        step3_file.write("par_save_RNA=\"no\"\n")

    if args.par_save_metadata_step3:
        step3_file.write("par_save_metadata=\"yes\"\n")
    else:
        step3_file.write("par_save_metadata=\"no\"\n")

    if args.par_seurat_object_step3:
        full_path_seurat_object = current_dir + "/" + args.par_seurat_object_step3
        step3_file.write("par_seurat_object=\"" + full_path_seurat_object + "\"\n")

    if args.par_nFeature_RNA_L:
        step3_file.write("par_nFeature_RNA_L=" + str(args.par_nFeature_RNA_L) + "\n")

    if args.par_nFeature_RNA_U:
        step3_file.write("par_nFeature_RNA_U=" + str(args.par_nFeature_RNA_U) + "\n")

    if args.par_nCount_RNA_L:
        step3_file.write("par_nCount_RNA_L=" + str(args.par_nCount_RNA_L) + "\n")

    if args.par_nCount_RNA_U:
        step3_file.write("par_nCount_RNA_U=" + str(args.par_nCount_RNA_U) + "\n")

    if args.par_mitochondria_percent_L:
        step3_file.write("par_mitochondria_percent_L=" + str(args.par_mitochondria_percent_L) + "\n")

    if args.par_mitochondria_percent_U:
        step3_file.write("par_mitochondria_percent_U=" + str(args.par_mitochondria_percent_U) + "\n")

    if args.par_ribosomal_percent_L:
        step3_file.write("par_ribosomal_percent_L=" + str(args.par_ribosomal_percent_L) + "\n")

    if args.par_ribosomal_percent_U:
        step3_file.write("par_ribosomal_percent_U=" + str(args.par_ribosomal_percent_U) + "\n")

    if args.par_remove_mitochondrial_genes:
        step3_file.write("par_remove_mitochondrial_genes=\"yes\"\n")
    else:
        step3_file.write("par_remove_mitochondrial_genes=\"no\"\n")

    if args.par_remove_ribosomal_genes:
        step3_file.write("par_remove_ribosomal_genes=\"yes\"\n")
    else:
        step3_file.write("par_remove_ribosomal_genes=\"no\"\n")

    if args.par_remove_genes:
        step3_file.write("par_remove_genes= c(" + args.par_remove_genes + ")\n")

    if args.par_regress_cell_cycle_genes:
        step3_file.write("par_regress_cell_cycle_genes=\"yes\"\n")
    else:
        step3_file.write("par_regress_cell_cycle_genes=\"no\"\n")

    if args.par_regress_custom_genes:
        step3_file.write("par_regress_custom_genes=\"yes\"\n")
    else:
        step3_file.write("par_regress_custom_genes=\"no\"\n")

    if args.par_regress_genes:
        step3_file.write("par_regress_genes= c(" + args.par_regress_genes + ")\n")

    if args.par_normalization_method_step3:
        step3_file.write("par_normalization.method=\"" + args.par_normalization_method_step3 + "\"\n")

    if args.par_scale_factor_step3:
        step3_file.write("par_scale.factor=" + str(args.par_scale_factor_step3) + "\n")

    if args.par_selection_method_step3:
        step3_file.write("par_selection.method=\"" + args.par_selection_method_step3 + "\"\n")

    if args.par_nfeatures_step3:
        step3_file.write("par_nfeatures=" + str(args.par_nfeatures_step3) + "\n")

    if args.par_top:
        step3_file.write("par_top=" + str(args.par_top) + "\n")

    if args.par_npcs_pca_step3:
        step3_file.write("par_npcs_pca=" + str(args.par_npcs_pca_step3) + "\n")

    step3_file.close()

    #################################
    # Generate the step_4.txt file #
    #################################

    step4_file_path = dir_fullpath + "/job_info/parameters/step4_par.txt"
    step4_file = open(step4_file_path, "w")

    if args.par_save_RNA_step4_scrna or args.par_save_RNA_step4_hto:
        step4_file.write("par_save_RNA=\"yes\"\n")
    else:
        step4_file.write("par_save_RNA=\"no\"\n")

    if args.par_save_metadata_step4_scrna or args.par_save_metadata_step4_hto:
        step4_file.write("par_save_metadata=\"yes\"\n")
    else:
        step4_file.write("par_save_metadata=\"no\"\n")

    if args.par_seurat_object_step4_scrna or args.par_seurat_object_step4_hto:
        seurat_object = args.par_seurat_object_step4_scrna if args.par_seurat_object_step4_scrna else args.par_seurat_object_step4_hto
        full_path_seurat_object = current_dir + "/" + seurat_object
        step4_file.write("par_seurat_object=\"" + full_path_seurat_object + "\"\n")

    if args.par_RunUMAP_dims_step4_scrna:
        run_umap_dims = args.par_RunUMAP_dims_step4_scrna
        step4_file.write("par_RunUMAP_dims=" + str(run_umap_dims) + "\n")

    if args.par_RunUMAP_n_neighbors_step4:
        step4_file.write("par_RunUMAP_n.neighbors=" + str(args.par_RunUMAP_n_neighbors_step4) + "\n")

    if args.par_dropDN_scrna or args.par_dropDN_hto:
        step4_file.write("par_dropDN=\"yes\"\n")
    else:
        step4_file.write("par_dropDN=\"no\"\n")

    if args.par_PCs:
        step4_file.write("par_PCs=" + str(args.par_PCs) + "\n")

    if args.par_pN:
        step4_file.write("par_pN=" + str(args.par_pN) + "\n")

    if args.par_sct:
        step4_file.write("par_sct=TRUE\n")
    else:
        step4_file.write("par_sct=FALSE\n")

    if args.par_rate_nExp:
        step4_file.write("par_rate_nExp=" + str(args.par_rate_nExp) + "\n")

    if args.par_sample_names_step4:
        # split the list of sample names ',' and add double quotes join the list with
        split_par_sample_names_step4 = args.par_sample_names_step4.split(',')
        par_sample_names_step4 = ",".join([f"'{item}'" for item in split_par_sample_names_step4])
        step4_file.write("par_sample_names= c(" + par_sample_names_step4 + ")\n")

    if args.par_expected_doublet_rate:
        step4_file.write("par_expected_doublet_rate= c(" + args.par_expected_doublet_rate + ")\n")

    if args.par_normalization_method_step4:
        step4_file.write("par_normalization.method=\"" + args.par_normalization_method_step4 + "\"\n")

    if args.par_scale_factor_step4:
        step4_file.write("par_scale.factor=" + str(args.par_scale_factor_step4) + "\n")

    if args.par_selection_method_step4:
        step4_file.write("par_selection.method=\"" + args.par_selection_method_step4 + "\"\n")

    if args.par_nfeatures_step4:
        step4_file.write("par_nfeatures=" + str(args.par_nfeatures_step4) + "\n")

    if args.par_dims_umap:
        step4_file.write("par_dims_umap=" + str(args.par_dims_umap) + "\n")

    if args.par_dimensionality_reduction:
        step4_file.write("par_dimensionality_reduction=\"" + args.par_dimensionality_reduction + "\"\n")

    if args.par_npcs_pca_step4:
        step4_file.write("par_npcs_pca=" + str(args.par_npcs_pca_step4) + "\n")

    if args.par_label_dropDN:
        step4_file.write("par_label_dropDN=" + str(args.par_label_dropDN) + "\n")

    if args.par_quantile:
        step4_file.write("par_quantile=" + str(args.par_quantile) + "\n")

    if args.par_autoThresh:
        step4_file.write("par_autoThresh=\"yes\"\n")
    else:
        step4_file.write("par_autoThresh=\"no\"\n")

    if args.par_maxiter:
        step4_file.write("par_maxiter=" + str(args.par_maxiter) + "\n")

    if args.par_RidgePlot_ncol:
        step4_file.write("par_RidgePlot_ncol=" + str(args.par_RidgePlot_ncol) + "\n")

    if args.par_old_antibody_label:
        # split the list of sample names ',' and add double quotes join the list with
        split_par_old_antibody_label = args.par_old_antibody_label.split(',')
        par_old_antibody_label = ",".join([f"'{item}'" for item in split_par_old_antibody_label])
        step4_file.write("par_old_antibody_label= c(" + par_old_antibody_label + ")\n")

    step4_file.close()

    #################################
    # Generate the step_5.txt file #
    #################################

    step5_file_path = dir_fullpath + "/job_info/parameters/step5_par.txt"
    step5_file = open(step5_file_path, "w")

    if args.par_save_RNA_step5:
        step5_file.write("par_save_RNA=\"yes\"\n")
    else:
        step5_file.write("par_save_RNA=\"no\"\n")

    if args.par_save_metadata_step5:
        step5_file.write("par_save_metadata=\"yes\"\n")
    else:
        step5_file.write("par_save_metadata=\"no\"\n")

    if args.par_seurat_object_step5:
        full_path_seurat_object = current_dir + "/" + args.par_seurat_object_step5
        step5_file.write("par_seurat_object=\"" + full_path_seurat_object + "\"\n")

    if args.par_one_seurat:
        step5_file.write("par_one_seurat=\"yes\"\n")
    else:
        step5_file.write("par_one_seurat=\"no\"\n")

    if args.par_integrate_seurat:
        step5_file.write("par_integrate_seurat=\"yes\"\n")
    else:
        step5_file.write("par_integrate_seurat=\"no\"\n")


    if args.par_merge_seurat:
        step5_file.write("par_merge_seurat=\"yes\"\n")
    else:
        step5_file.write("par_merge_seurat=\"no\"\n")

    if args.par_DefaultAssay:
        step5_file.write("par_DefaultAssay=\'" + args.par_DefaultAssay + "\'\n")

    if args.par_normalization_method_step5:
        step5_file.write("par_normalization.method=\"" + args.par_normalization_method_step5 + "\"\n")

    if args.par_scale_factor_step5:
        step5_file.write("par_scale.factor=" + str(args.par_scale_factor_step5) + "\n")

    if args.par_selection_method_step5:
        step5_file.write("par_selection.method=\"" + args.par_selection_method_step5 + "\"\n")

    if args.par_nfeatures_step5:
        step5_file.write("par_nfeatures=" + str(args.par_nfeatures_step5) + "\n")

    if args.par_FindIntegrationAnchors_dim:
        step5_file.write("par_FindIntegrationAnchors_dim=" + str(args.par_FindIntegrationAnchors_dim) + "\n")

    if args.par_RunPCA_npcs:
        step5_file.write("par_RunPCA_npcs=" + str(args.par_RunPCA_npcs) + "\n")

    if args.par_RunUMAP_dims_step5:
        step5_file.write("par_RunUMAP_dims=" + str(args.par_RunUMAP_dims_step5) + "\n")

    if args.par_RunUMAP_n_neighbors_step5:
        step5_file.write("par_RunUMAP_n.neighbors=" + str(args.par_RunUMAP_n_neighbors_step5) + "\n")

    if args.par_compute_jackstraw:
        step5_file.write("par_compute_jackstraw=\"yes\"\n")
    else:
        step5_file.write("par_compute_jackstraw=\"no\"\n")

    step5_file.close()

    #################################
    # Generate the step_6.txt file #
    #################################

    step6_file_path = dir_fullpath + "/job_info/parameters/step6_par.txt"
    step6_file = open(step6_file_path, "w")

    if args.par_save_RNA_step6:
        step6_file.write("par_save_RNA=\"yes\"\n")
    else:
        step6_file.write("par_save_RNA=\"no\"\n")

    if args.par_save_metadata_step6:
        step6_file.write("par_save_metadata=\"yes\"\n")
    else:
        step6_file.write("par_save_metadata=\"no\"\n")

    if args.par_seurat_object_step6:
        full_path_seurat_object = current_dir + "/" + args.par_seurat_object_step6
        step6_file.write("par_seurat_object=\"" + args.par_seurat_object_step6 + "\"\n")

    if args.par_skip_integration:
        step6_file.write("par_skip_integration=\"yes\"\n")
    else:
        step6_file.write("par_skip_integration=\"no\"\n")

    if args.par_FindNeighbors_dims:
        step6_file.write("par_FindNeighbors_dims=" + str(args.par_FindNeighbors_dims) + "\n")

    if args.par_RunUMAP_dims_step6:
        step6_file.write("par_RunUMAP_dims=" + str(args.par_RunUMAP_dims_step6) + "\n")

    if args.par_FindNeighbors_k_param:
        step6_file.write("par_FindNeighbors_k.param=" + str(args.par_FindNeighbors_k_param) + "\n")

    if args.par_FindNeighbors_prune_SNN:
        step6_file.write("par_FindNeighbors_prune.SNN=" + str(args.par_FindNeighbors_prune_SNN) + "\n")

    if args.par_FindClusters_resolution:
        step6_file.write("par_FindClusters_resolution= c(" + args.par_FindClusters_resolution + ")\n")

    if args.par_compute_ARI:
        step6_file.write("par_compute_ARI=\"yes\"\n")
    else:
        step6_file.write("par_compute_ARI=\"no\"\n")

    if args.par_RI_reps:
        step6_file.write("par_RI_reps=" + str(args.par_RI_reps) + "\n")

    step6_file.close()

    #################################
    # Generate the step_7.txt file #
    #################################

    step7_file_path = dir_fullpath + "/job_info/parameters/step7_par.txt"
    step7_file = open(step7_file_path, "w")

    if args.par_save_RNA_step7:
        step7_file.write("par_save_RNA=\"yes\"\n")
    else:
        step7_file.write("par_save_RNA=\"no\"\n")

    if args.par_save_metadata_step7:
        step7_file.write("par_save_metadata=\"yes\"\n")
    else:
        step7_file.write("par_save_metadata=\"no\"\n")

    if args.par_seurat_object_step7:
        full_path_seurat_object = current_dir + "/" + args.par_seurat_object_step7
        step7_file.write("par_seurat_object=\"" + args.par_seurat_object_step7 + "\"\n")

    if args.par_level_cluster:
        step7_file.write("par_level_cluster=\"" + str(args.par_level_cluster) + "\"\n")

    if args.par_run_find_marker:
        step7_file.write("par_run_find_marker=\"yes\"\n")
    else:
        step7_file.write("par_run_find_marker=\"no\"\n")


    if args.par_top_sel:
        step7_file.write("par_top_sel=" + str(args.par_top_sel) + "\n")

    if args.par_run_enrichR:
        step7_file.write("par_run_enrichR=\"yes\"\n")
    else:
        step7_file.write("par_run_enrichR=\"no\"\n")

    if args.par_db:
        # split the list of databases ',' and add double quotes join the list with
        split_par_db = args.par_db.split(',')
        par_db = ",".join([f"'{item}'" for item in split_par_db])
        step7_file.write("par_db= c(" + par_db + ")\n")

    if args.par_run_module_score:
        step7_file.write("par_run_module_score=\"yes\"\n")
    else:
        step7_file.write("par_run_module_score=\"no\"\n")

    if args.par_module_score:
        full_path_module_score = current_dir + "/" + args.par_module_score
        step7_file.write("par_module_score=\"" + full_path_module_score + "\"\n")

    if args.par_run_visualize_markers:
        step7_file.write("par_run_visualize_markers=\"yes\"\n")
    else:
        step7_file.write("par_run_visualize_markers=\"no\"\n")

    if args.par_select_features_list:
        # split the list of sample names ',' and add double quotes join the list with
        split_par_select_features_list = args.par_select_features_list.split(',')
        par_select_features_list = ",".join([f"'{item}'" for item in split_par_select_features_list])
        step7_file.write("par_select_features_list= c(" + par_select_features_list + ")\n")

    if args.par_select_features_csv:
        full_path_select_features_csv = current_dir + "/" + args.par_select_features_csv
        step7_file.write("par_select_features_csv=\"" + full_path_select_features_csv + "\"\n")

    if args.par_reference:
        full_path_reference = current_dir + "/" + args.par_reference
        step7_file.write("par_reference=\"" + full_path_reference + "\"\n")

    if args.par_reference_name:
        step7_file.write("par_reference_name=\"" + args.par_reference_name + "\"\n")

    if args.par_level_celltype:
        step7_file.write("par_level_celltype=\"" + args.par_level_celltype + "\"\n")

    if args.par_FindTransferAnchors_dim:
        step7_file.write("par_FindTransferAnchors_dim=" + str(args.par_FindTransferAnchors_dim) + "\n")

    if args.par_futureglobalsmaxSize:
        step7_file.write("par_futureglobalsmaxSize=" + str(args.par_futureglobalsmaxSize) + "\n")

    if args.par_annotate_resolution:
        step7_file.write("par_annotate_resolution=\"" + args.par_annotate_resolution + "\"\n")

    if args.par_name_metadata:
        step7_file.write("par_name_metadata=\"" + args.par_name_metadata + "\"\n")

    if args.par_annotate_labels:
        # split the list of sample names ',' and add double quotes join the list with
        split_par_annotate_labels = args.par_annotate_labels.split(',')
        par_annotate_labels = ",".join([f"'{item}'" for item in split_par_annotate_labels])
        step7_file.write("par_annotate_labels= c(" + par_annotate_labels + ")\n")

    step7_file.close()

    #################################
    # Generate the step_8.txt file #
    #################################

    step8_file_path = dir_fullpath + "/job_info/parameters/step8_par.txt"
    step8_file = open(step8_file_path, "w")

    if args.par_save_RNA_step8:
        step8_file.write("par_save_RNA=\"yes\"\n")
    else:
        step8_file.write("par_save_RNA=\"no\"\n")

    if args.par_save_metadata_step8:
        step8_file.write("par_save_metadata=\"yes\"\n")
    else:
        step8_file.write("par_save_metadata=\"no\"\n")

    if args.par_seurat_object_step8:
        full_path_seurat_object = current_dir + "/" + args.par_seurat_object_step8
        step8_file.write("par_seurat_object=\"" + args.par_seurat_object_step8 + "\"\n")

    if args.par_metadata:
        full_path_metadata = current_dir + "/" + args.par_metadata
        step8_file.write("par_metadata=\"" + full_path_metadata + "\"\n")

    if args.par_merge_meta:
        step8_file.write("par_merge_meta=\"" + args.par_merge_meta + "\"\n")

    if args.par_run_cell_based_all_cells:
        step8_file.write("par_run_cell_based_all_cells=\"yes\"\n")
    else:
        step8_file.write("par_run_cell_based_all_cells=\"no\"\n")

    if args.par_run_cell_based_celltype_groups:
        step8_file.write("par_run_cell_based_celltype_groups=\"yes\"\n")
    else:
        step8_file.write("par_run_cell_based_celltype_groups=\"no\"\n")

    if args.par_run_sample_based_all_cells:
        step8_file.write("par_run_sample_based_all_cells=\"yes\"\n")
    else:
        step8_file.write("par_run_sample_based_all_cells=\"no\"\n")

    if args.par_run_sample_based_celltype_groups:
        step8_file.write("par_run_sample_based_celltype_groups=\"yes\"\n")
    else:
        step8_file.write("par_run_sample_based_celltype_groups=\"no\"\n")

    if args.par_statistical_method:
        step8_file.write("par_statistical_method=\"" + args.par_statistical_method + "\"\n")

    step8_file.close()

    if args.contrast_cell_based_all_cells:
        contrast_cell_based_all_cells_path = dir_fullpath + "/job_info/parameters/step8_contrast_cell_based_all_cells.txt"
        with open(contrast_cell_based_all_cells_path, "w") as contrast_file:
            contrast_file.write(args.contrast_cell_based_all_cells)

    if args.contrast_cell_based_celltype_groups:
        contrast_cell_based_celltype_groups_path = dir_fullpath + "/job_info/parameters/step8_contrast_cell_based_celltype_groups.txt"
        with open(contrast_cell_based_celltype_groups_path, "w") as contrast_file:
            contrast_file.write(args.contrast_cell_based_celltype_groups)

    if args.contrast_sample_based_all_cells:
        contrast_sample_based_all_cells_path = dir_fullpath + "/job_info/parameters/step8_contrast_sample_based_all_cells.txt"
        with open(contrast_sample_based_all_cells_path, "w") as contrast_file:
            contrast_file.write(args.contrast_sample_based_all_cells)

    if args.contrast_sample_based_celltype_groups:
        contrast_sample_based_celltype_groups_path = dir_fullpath + "/job_info/parameters/step8_contrast_sample_based_celltype_groups.txt"
        with open(contrast_sample_based_celltype_groups_path, "w") as contrast_file:
            contrast_file.write(args.contrast_sample_based_celltype_groups)

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
    # Generate the command line         #
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

    # Execute launch_scrnabox.sh script
    command = [ "launch_scrnabox.sh", "-d", dir_fullpath, "--steps", args.steps ]

    # Add msd if provided
    if args.msd:
        command.append("--msd")

    # Add markergsea if provided
    if args.markergsea:
        command.append("--markergsea T")

    # Add knownmarkers if provided
    if args.knownmarkers:
        command.append("--knownmarkers T")

    # Add referenceannotation if provided
    if args.referenceannotation:
        command.append("--referenceannotation T")

    # Add annotate if provided
    if args.annotate:
        command.append("--annotate T")

    # Add addmeta if provided
    if args.addmeta:
        command.append("--addmeta T")

    # Add rundge if provided
    if args.rundge:
        command.append("--rundge T")

    # Add seulist if provided
    if args.seulist:
        command.append("--seulist T")

    # Add rcheck if provided
    if args.rcheck:
        command.append("--rcheck T")

    print_log("Executing command:\n" + " ".join(command) + "\n")
    cmd    = " ".join(command)
    status = os.system(cmd)

    if status != 0:
        print_log("Error while executing launch_scrnabox.sh script")
        sys.exit(status)

    # Verify if each step* file in dir_fullpath exists
    execution_failed = False
    for step in args.steps.split(','):
        step_dir_path = dir_fullpath + f"/step{step}"
        # Check if the step directory exists
        if not os.path.exists(step_dir_path):
            print_log(f"Error: {step_dir_path} does not exist. Please check the steps you provided.")
            execution_failed = True
        else:
        # Verify if at least one file exists in the step director
            file_count = sum(len(files) for _, _, files in os.walk(step_dir_path))
            print_log(f"Number of files in {step_dir_path}: {file_count}")
            if file_count == 0:
                print_log(f"Error: {step_dir_path} is empty. Please check the steps you provided.")
                execution_failed = True

        log_files_pattern = dir_fullpath + f"/job_info/logs/step_{step}*"
        # list all log files that start with step_{step}
        last_log_files = sorted(glob.glob(log_files_pattern), key=os.path.getmtime, reverse=True)
        if len(last_log_files)  > 0:
            last_log_file = last_log_files[0]
            tail_cmd      = f"tail -n 50 {last_log_file}"
            tail_content  = subprocess.run(tail_cmd, shell=True, capture_output=True, text=True)
            print(f"Last 50 lines of log file: {last_log_file}")
            if tail_content.returncode == 0:
                print(tail_content.stdout)
            else:
                print_log(f"Error while reading log file {last_log_file}: {tail_content.stderr}")
        else:
            print_log("No log files found.")
            sys.exit(2)

        if execution_failed:
            print_log("Execution failed. Please check the steps you provided and the log files.")
            sys.exit(2)

        print_log("Execution completed successfully.")

        # if step1 was runned
        if '1' in args.steps.split(','):
            # find all bam files in `./step1/*/*/SC_RNA_COUNTER_CS/*"`
            bam_files = glob.glob(dir_fullpath + "/step1/*/*/SC_RNA_COUNTER_CS/*/*.bam")
            if bam_files:
                # remove the bam files
                print_log("Removing bam files from step1.")
                for bam_file in bam_files:
                    try:
                        os.remove(bam_file)
                        print_log(f"Removed {bam_file}")
                    except Exception as e:
                        print_log(f"Error removing {bam_file}: {e}")

        # copy the working directory to the output directory
        print_log("Copying the working directory to the output directory.")
        subprocess.run(["cp", "-r", dir_fullpath, args.output_dir])

    sys.exit(0)

if __name__ == "__main__":
    main()
