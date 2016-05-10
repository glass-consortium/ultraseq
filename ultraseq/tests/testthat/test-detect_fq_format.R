context("Detect FQ naming format")

test_that("naming format for CASAVA", {
  
  
  # CASAVA 1.8, no index
  fq = "ABCD_00647_NoIndex_L001_R1_001.fastq.gz"
  expect_message(detect_fq_format2(fq), "Using CASAVA 1.8 naming format")
  
  fq = "ABCD_00835_TGACCA_L001_R1_004.fastq.gz"
  expect_message(detect_fq_format2(fq), "Using CASAVA 1.8 naming format")

  fq = "ABCD-FC112-MS11-Cap854-3-ID09_GATCAG_L008_R1_001.fastq.gz"
  expect_message(detect_fq_format2(fq), "Using CASAVA 1.8 naming format")

  # CASAVA 2.0
  fq = "CCCT_00914_S19_L008_R1_001.fastq.gz"
  expect_message(detect_fq_format2(fq), "Using MiSeq/bcl2fastq 2.0 naming format")

  fq = "AB-8-20m-DOX-serbp-1_S43_L008_R1_001.fastq.gz"
  expect_message(detect_fq_format2(fq), "Using MiSeq/bcl2fastq 2.0 naming format")
  
  # miseq:
  fq = "ABCD-ABC_S1_L001_R1_001.fastq.gz"
  expect_message(detect_fq_format2(fq), "Using MiSeq/bcl2fastq 2.0 naming format")
  

})