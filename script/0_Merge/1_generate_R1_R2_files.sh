#!/bin/bash

# load paths
source /Volumes/MY_PASSPORT/draft_final/script/source/paths_for_bash.sh

# export list of R1 and R2 sequences
cd $PATH_DATA
ls | grep "R1" > $PATH_R1_R2_FILES/0_R1.txt
ls | grep "R2" > $PATH_R1_R2_FILES/0_R2.txt
