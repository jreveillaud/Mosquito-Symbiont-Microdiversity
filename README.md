This project describes the bioinformatic reproductible workflow associated to the paper "The mosquito microbiome includes habitat-specific but rare symbionts". 

## Downloading the fastq files
Before everything else, you need to download the data and place it in a Data folder in the main directory. 
The raw illumina paired-end reads of 136 samples are stored in ENA and are accesible with the accession number PRJEB43079.

## Complete reproductible workflow

In the project main directory, you will find a R project file to open in Rstudio : metabarcoding.Rproj. 

You will have access to multiple sub-folders numbered from 0 to 2 that need to be run one after the other. Each numbered folder corresponds to a step of the global analysis : 
- Folder 0_Merge corresponding to the Merge step
- Folder 1_MED corresponding to the MED analysis
- Folder 2_Oligotyping corresponding to the Oligotyping analyses

And the folders 1_MED and 2_Oligotyping include multiple additional sub-folders starting with a letter ranging from A to D and corresponding to a sub-step of analysis. 
For example, folder 1_MED contains : 
- A_Run_MED corresponding to the MED run step
- B_Taxonomy corresponding to the taxonomy step
- C_Phyloseq corresponding to the creation of a phyloseq object
- D_Decontam corresponding to the detection and removing of contaminants
- E_Analysis corresponding to the analysis done on the dataset without contaminants. This folder contains a sub-folder for each analysis (distribution and rarefaction, alpha diversity, beta diversity, etc)

So, to do the complete workflow, you need to run the scripts by number and letter order (0-A, 0-B,....1-A,1-B, etc).

The folder "source" includes all functions and the main R packages used in the workflow. 

## Metadata and additionnal files

You will also find a metadata folder that contains the Supplementary Table 1 avaible in the MS and a additional_files folder that contains files used for the Beta diversity analysis. 


## Needed ressources
The workflow uses mainly two softwares, oligotyping version 2.1 and anvi'o version 7.1, using a docker container and a conda environment, respectively. 

To install oligotyping-2.1 in a docker container, you can follow this great tutorial : https://merenlab.org/2014/09/02/virtualbox/

To install anvi'o-7.1, you can follow this one : https://merenlab.org/2016/06/26/installation-v2/


You also need to download the silva database version 138 and put files into a Bank folder in the main directory :
- silva_nr99_v138_train_set.fa.gz
- silva_species_assignment_v138.fa files 

Available just here: https://zenodo.org/record/3986799#.Ya3eHy3pN71. 




