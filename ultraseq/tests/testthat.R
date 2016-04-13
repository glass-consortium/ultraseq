library(testthat)
library(ultraseq)

test_check("ultraseq")


# test preprocess

# test mutect

# test bam_mutect

# test with dict file and WO dict file


# test flowr run, 
# - with flow in same folder, 
# - path to a file

# this document provides and overview of the tests included

# # testing
# Rscript -e 'library(flowr);run(x="bam_mutect", tumor_bam="tumor.bam", normal_bam="normal.bam", paired_name="tumor_normal", platform="local", execute=FALSE)'