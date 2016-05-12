#library(testthat)


context("Test fq sheet creation")

test_that("Creating a df with, from fq names works well", {
  
  # CASAVA 1.8, no index
  fq18 = c("ABCD_00647_NoIndex_L001_R1_001.fastq.gz", 
           "ABCD_00835_TGACCA_L001_R1_004.fastq.gz",
           "ABCD-FC112-MS11-Cap854-3-ID09_GATCAG_L008_R1_001.fastq.gz")

  fmt18 = detect_fq_format2(fq18)
  df = split_names_fastq2(fq18, fmt18)
  expect_that(class(df), is_equivalent_to("data.frame"))
  expect_that(ncol(df), is_equivalent_to(6))
  
  # CASAVA 2.0
  fq20 = c("ABCD_00914_S19_L008_R1_001.fastq.gz", 
          "AB-8-20m-DOX-serbp-1_S43_L008_R1_001.fastq.gz",
          "ABCD-ABC_S1_L001_R1_001.fastq.gz")

  fmt20 = detect_fq_format2(fq20)
  df = split_names_fastq2(fq20, fmt20)
  expect_that(class(df), is_equivalent_to("data.frame"))
  expect_that(ncol(df), is_equivalent_to(5))
  
  # ------ testing both together ----
  # expect warning
  fq1820 = c(fq18, fq20)
  fmt1820 = detect_fq_format2(fq1820)
  expect_warning(split_names_fastq2(fq1820, fmt1820), "there was a issue parsing this filename")
  
  # expect error:
  expect_warning(expect_error(split_names_fastq2(fq1820, fmt1820, strict_format_checking = TRUE), "Some file names do not have the correct format"))
  
  
  # ----- add a wrong format ----
  fq20.e = c("ABCD_00914_AB19_L008_R1_001.fastq.gz")
  expect_error(detect_fq_format2(fq20.e), "Looks like we could not understand pattern in names ")
  expect_warning(split_names_fastq2(fq20.e,  "{{samplename}}_{{index}}_L00{{lane}}_R{{read}}_{{num}}.fastq.gz"))


})