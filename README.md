# GeTallele

Repository with the toolbox (and code used for analysis and to produce figures) described in: "GeTallele: a method for analysis of DNA and RNA allele frequencies distributions". Preprint availible at https://www.biorxiv.org/content/10.1101/491209v4

## GeTallele descritpion
### Functionality
The toolbox that provides a suite of functions for analysis, statistical assessment and visualization of **Ge**nome and **T**ranscriptome **allele** frequencies distributions. The main functionality of GeTallele is estimation of **variant probability (v_{PR})** in chromosomal segments (continuous multi-SNV genomic regions). Variant p

### Input

### Example

## Directory structure of the repository:
* toolbox_v1 - folder with matlab functions that provide core and extended functionalites of the toolbox
  * external - folder with matlab tools by other developers (download from https://mathworks.com/matlabcentral/fileexchange/)  
  * circos - folder with circos .conf and functions that convert GeTallele results into circos input files
* example - folder with matlab script files containg walkthrough examples and files with example synthetic data
* paper_code - folder with the matlab code that has been used to produce results and figures presented in the paper. In comparison with the latest version of the toolbox the code is commented very sparingly. This part of the code is no longer being developed and is include for archival and reproducibilty purposes. If you have questions concering any part of the code please contact me at p.m.slowinski@exeter.ac.uk.
  * toolbox_v0 - the initial version of the toolbox that has been build around and optimised for analysis of datasets containing 4 matched signals (normal exome, normal transcriptome, tumor exome, tumor transcriptome)
  * figures - functions and scripts that have beene used to generate all the figures included in the paper 
  * data - some auxiliary datasets that have been used for the analysis presented in the paper, open access to the full dataset cannot be provided due to GDPR and U. S. Privacy Act regulations.
  
## List of functions of the GeTallele toolboox:
### Core functions
*
*
### Additional functions
*
*

## To do:
* function to statistically (bootstrap based methodolgy) compare vpr values in two segments
* example of how to generate a synthetic dataset (method used to genereate files: example_data_Nex.tsv, example_data_Tex.tsv, example_data_Ntr.tsv, example_data_Ttr.tsv)
* R version of the toolbox

## Contact:
Any questions/ comments/ bugs please get in touch at p.m.slowinski@exeter.ac.uk

## Acknowledgments


