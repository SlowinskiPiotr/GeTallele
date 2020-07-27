# GeTallele

Repository with the toolbox (and code used for analysis and to produce figures) described in: "GeTallele: a method for analysis of DNA and RNA allele frequencies distributions". Preprint availible at https://www.biorxiv.org/content/10.1101/491209v4

## GeTallele descritpion
### Functionality
The toolbox that provides a suite of functions for analysis, statistical assessment and visualization of **Ge**nome and **T**ranscriptome **allele** frequencies distributions. The main functionality of GeTallele is estimation of **variant probability (vpr)** in chromosomal segments (continuous multi-SNV genomic regions). Variant probability allows for high-level description of variant allele frequencies (VAF) distributions in chromosomal segments. 

### Requirments
GeTallele requires basic knowledge of Matlab software environment. The Matlab skills requird to use GeTallele include: knowledge of basic Matlab syntax (e.g. *for* loop), importing data, working with Matlab scripts, working with sections in Matlab scripts, using Matlab help. 

GeTallel was created and tested with Matlab 2020a on macOS 10.14.6. 

GeTallele requires: 'Signal Processing Toolbox' and 'Statistics and Machine Learning Toolbox'.

### Input
The input to the GeTallele are variant read counts, and reference read counts. The variant and reference read counts can be generated by any available tool. As an input GeTallele requires a matrix that contains four columns that contain: 
1. chromosome number, 
1. coordinates in base pairs on the chromosome, 
1. variant read counts, 
1. reference read counts. 

Two first columns are required for interpretability and plotting. Vpr estimation requires only variant and reference read counts that can be obtaine by any availible method. Matrix with data for analysis can be imported into Matlab software environment from a wide range of text and binary file formats using buitl-in Matlab functionality.

Since VAF estimations can be affected by allele mapping bias (https://pubmed.ncbi.nlm.nih.gov/19808877) which can lead to overestimation of the reference allele count (https://pubmed.ncbi.nlm.nih.gov/25787242), we suggest that GetAllele input is based on SNV-aware alignments (e.g. using https://pubmed.ncbi.nlm.nih.gov/26366987).

To compare vpr values in different signals GeTallele requires matched sequencing sets. Namely, different signals should have reads at the same base-pair postions. Such matched signals can be for example obtained using RNA2DNAlign package (https://github.com/HorvathLab/RNA2DNAlign, https://pubmed.ncbi.nlm.nih.gov/27576531/)

### Example
Folder \example contains a script 'simple_example.m' that presents basic steps of using GeTallele. The script is the best starting point to get started with GeTallele. The script should be run from within the \example folder. The script shows how to:
* add the GeTallele toolbox folder to Matlab search path [*lines 1-2*]
* import data from an example file (which contains synthetic data) [*lines 4-10*]
* generate histogram bins based on a Farey sequence [*lines 12-14*]
* generate a set of synthetic VAF samples (and their experimental cumulative distribution functions CDFs) that are later used to estimate vpr values in the data [*lines 16-26*]
* convert read counts to VAF values [*line 29*] 
* segment and estimate vpr values in all the chromosomes in the example dataset [*lines 30-44*]
* plot results of segmentation for each chromosome [*lines 46-55*]
* generate a circos plot with results of segmentation of all chromosomes [*lines 57-94*]. Additionally, this part of the example provides a basic guidance on installing and using cirocs from whitin Matlab on macOS. The instructions provided in the guidance worked in (July 2020).

## Directory structure of the repository:
* toolbox_v1 - __the main GeTallele folder__, contains matlab functions that provide core and extended functionalites of the toolbox
  * external - folder with matlab tools by other developers (downloaded from https://mathworks.com/matlabcentral/fileexchange/)  
  * circos - folder with circos .conf (\circos_conf_files subfolder) and functions that convert GeTallele results into circos input text files
* example - folder with matlab script files containg walkthrough examples and files with example synthetic data
* paper_code - folder with the matlab code that has been used to produce results and figures presented in the paper. In comparison with the latest version of the toolbox the code is commented very sparingly. This part of the code is no longer being developed and is include for archival and reproducibilty purposes. If you have questions concering any part of the code please contact me at p.m.slowinski@exeter.ac.uk.
  * toolbox_v0 - the initial version of the toolbox that has been build around and optimised for analysis of datasets containing 4 matched signals (normal exome, normal transcriptome, tumor exome, tumor transcriptome)
  * figures - functions and scripts that have beene used to generate all the figures included in the paper 
  * data - some auxiliary datasets that have been used for the analysis presented in the paper, open access to the full dataset cannot be provided due to GDPR and U. S. Privacy Act regulations.
  
## List of functions of the GeTallele toolboox:
* __Core functions__
** farey_bins.m - this function generates bin edges and bin centers based on a Farey sequennce.
** convert_reads_to_VAF.m - this function converts read counts into VAF values.
** generate_VAF_cdfs.m - this function generates synthetic VAF samples (and their experimental cumulative distribution functions CDFs) with vpr from 0.5 to 1 with 0.01 step.
** find_segments_and_fit_vpr.m - this function segments data into chromosomal segments and estimates vpr in each segment.
** plot_VAF_and_vpr.m - this functions plots output of the function find_segments_and_fit_vpr.m or fit_vpr_in_segments.m.  
* __Functions used by the core functions__
** farey_sequence.m - used by the function farey_bins.m to generate farey_sequence.
** generate_VAF_sample.m - used by the function generate_VAF_cdfs.m to generate synthethic VAF samples.
* __Additional functions__
** __functions that show how variant probability eestimates can be used to estimate purity__ example of usage of this function can be found in script example_matched_signals_and_purity.m that can be found in the folder example).
*** fit_vpr_in_segments.m - this function uses segments generated by the find_segments_and_fit_vpr.m to estimate vpr in the same segments in a different dataset 
*** gen_all_possible_vprs.m - fuction that generates a full set of proportions of all the population in the mixture with step (increment of) 0.01 and compute all the possible vpr values that each of the mixtures could produce. Function used to produce results in section 4.4 of the paper.
*** VBP_estimator.m - function that estimates variant based purity using vpr values estimated in the data and output of the function gen_all_possible_vprs.m
*** plot_adm_mix.m - function that generates ternary plot showing admissible mixtures (for interpretation see the paper)
** generate_var_ref_sample.m - function that can be used to generate synthethic samples of variant and reference read counts (based on generate_VAF_sample.m). Was used to generate the synthethic datasets.

## To do:
* function to statistically (bootstrap based methodolgy) compare vpr values in two segments (that don't have to be matched)
* example of how to generate a synthetic dataset (method used to genereate files: example_data_Nex.tsv, example_data_Tex.tsv, example_data_Ntr.tsv, example_data_Ttr.tsv)
* R version of the toolbox

## Contact:
Any questions/ comments/ bugs please get in touch at p.m.slowinski@exeter.ac.uk

## Acknowledgments
Development of the GeTallel was generously supported by:
* McCormick Genomic and Proteomic Center (MGPC), The George Washington University; [MGPC_PG2018]. 
* Wellcome Trust Institutional Strategic Support Award [204909/Z/16/Z]. 
* EPSRC [EP/N014391/1].
