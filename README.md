This project describes the bioinformatic reproductible workflow associated to the paper "Mosquito symbiont microdiversity as revealed by oligotyping". 

## Downloading the fastq files
Before everything else, you need to download the data. 
The raw illumina paired-end reads of 136 samples are stored in ENA and accesible with the accession number PRJEB43079. 

## Complete reproductible workflow

In the project main directory, you will find a R project file to open in Rstudio : metabarcoding.Rproj. 

You will have access to multiple sub-folders numbered from 0 to 2 (each corresponding to a step of the global analysis) that need to be run one after the other. Each numbered folder corresponds to a step of the global analysis : 
- Folder 0_Merge corresponding to the Merge step
- Folder 1_MED corresponding to the MED analysis
- Folder 2_Oligotyping corresponding to the Oligotyping analysis

And the folders 1_MED and 2_Oligotyping include multiple additional sub-folders starting with a letter from A to D and corresponding to a sub-step of analysis. 
For example, the folder 1_MED contains : 
- A_Run_MED that corresponding to the MED run step
- B_Taxonomy that corresponding to the taxonomy step
- C_Phyloseq that corresponding to the creation of a phyloseq object
- D_Decontam that corresponding to the detection and removing of contaminants
- E_Analysis that corresponding to the analysis done with the dataset without contaminants. This folder contains a sub-folder for each analysis (Distribution and rarefaction, alpha diversity, beta diversity, etc)

So, to do the complete workflow, you thus need to run the scripts by number and letter order (0-A, 0-B,....1-A,1-B, etc).

By the way, the folder "source" including all functions and the main R packages used in the workflow. 

## Needed ressources
The workflow uses mainly two softwares, oligotyping version 2.1 and anvi'o version 7.1, using a docker container and a conda environment respectively. 
To install oligotyping-2.1 in a docker container, you can follow this great tutorial : https://merenlab.org/2014/09/02/virtualbox/
To install anvi'o-7.1, you can follow this one : https://merenlab.org/2016/06/26/installation-v2/




