# CREATE SAMPLE FILE TO SETUP THE MERGE STEP

## Load main paths and package
source("../source/paths.R")
require(tidyverse);packageVersion("tidyverse")

## import R1 and R2 files
setwd(path_output)
if (!dir.exists("0_useful_files")) {dir.create("0_useful_files")}

setwd(paste0(path_output, "/0_useful_files"))
R1 <- read.table("0_R1.txt")
R2 <- read.table("0_R2.txt")

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


