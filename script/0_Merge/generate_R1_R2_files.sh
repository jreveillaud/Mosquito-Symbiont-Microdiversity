#!/bin/bash

PATH_MAIN="/Volumes/MY_PASSPORT/draft_final"
PATH_DATA=$PATH_MAIN/data
PATH_OUTPUT=$PATH_MAIN/output/0_useful_files

# export list of R1 and R2 sequences
cd $PATH_DATA
ls | grep "R1" > $PATH_OUTPUT/0_R1.txt
ls | grep "R2" > $PATH_OUTPUT/0_R2.txt