# CREATE SAMPLE FILE TO SETUP THE MERGE STEP

## library and path
require(tidyverse);packageVersion("tidyverse")
path_main <- "/Volumes/MY_PASSPORT/draft_final"
path <- paste0(path_main, "/output/0_useful_files")
path_data <- paste0(path_main, "/data/")


## import R1 and R2 files
setwd(path)
R1 <- read.table("0_R1.txt")
R2 <- read.table("0_R2.txt")
# R1 <- paste0(path_data, df)
# R2 <- paste0(path_data, R2)

## create file
df <- data.frame(R1, R2)
colnames(df) <- c("r1", "r2")
df$sample <- sub("\\_.*", "", df$r1)
df <- df %>% select(c(sample, r1, r2))
df$r1 <- paste0(path_data, df$r1)
df$r2 <- paste0(path_data, df$r2)

## print df
head(df)

## save file
write.table(df,'0_sample_file.txt', sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


