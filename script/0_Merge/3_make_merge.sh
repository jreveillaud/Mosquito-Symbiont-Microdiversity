#!/bin/bash

# load paths
source /Volumes/MY_PASSPORT/draft_final/script/source/paths_for_bash.sh

# run config
cd $PATH_R1_R2_FILES
mkdir ini
source /Users/hschrieke/miniconda2/bin/activate illumina-utils
iu-gen-configs $PATH_R1_R2_FILES/0_sample_file.txt -o ini

# merge
cd $PATH_R1_R2_FILESini

for ini in *.ini
do
iu-merge-pairs --min-overlap-size 30 --max-num-mismatches 0 --enforce-Q30-check $ini
echo $ini "done" > 0_merge_progession.txt
done

# cp merged sequence in sequence_merged folder
mv $(ls | grep "MERGED" ) $PATH_R1_R2_FILES
mv $(ls | grep "MERGED" ) $PATH_MERGE_SCRIPT
