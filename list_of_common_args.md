## A list of common names

**arguments**

- `tumor_bam` tumor bam file
- `normal_bam` normal bam file

- `tumor_sample`  tumor samplename
- `normal_sample`  normal reference

- `out_prefix`      a combination of tumor and normal name
                 used as out-prefix in several commands

- `fqs_tum_1`
- `fqs_tum_2`
- `fqs_ref_1`
- `fqs_ref_2`

**function names**

- `preprocess`
- `mutect`
- `bwa`
- `samtools_sort`
- `samtools_merge`
- `samtools_index`


**pipeline names**

- `bam_preprocess`
- `bam_mutect`: take two bams, preprocess and perform mutect
