#!/bin/bash

PATH_MAIN="/Volumes/MY_PASSPORT/draft_final"
PATH_SCRIPT=$PATH_MAIN/script/0_Merge
PATH_OUTPUT=$PATH_MAIN/output/0_useful_files

# run config
cd $PATH_OUTPUT
mkdir ini
source /Users/hschrieke/miniconda2/bin/activate illumina-utils
iu-gen-configs $PATH_OUTPUT/0_sample_file.txt -o ini

# merge
cd $PATH_OUTPUT/ini

for ini in *.ini
do
iu-merge-pairs --min-overlap-size 30 --max-num-mismatches 0 --enforce-Q30-check $ini
echo $ini "done" > 0_merge_progession.txt
done

# cp merged sequence in sequence_merged folder
mv $(ls | grep "MERGED" ) $PATH_OUTPUT

mv $(ls | grep "MERGED" ) /Volumes/MY_PASSPORT/draft_final/output/0_merged