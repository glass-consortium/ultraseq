
This repo is under active development, please raise a issue if you have trouble
using/installing it.


We have various pipelines available as part of ultraseq.

```
├── bin
│   ├── ultraseq
│   └── ultraseq.conf
└── pipes
    ├── fastq_bam_bwa.R
    ├── fastq_bam_bwa.def
    ├── fastq_mutect.R
    ├── fastq_mutect.def
    ├── fastq_mutect.md
    ├── fastq_mutect.pdf
    ├── fastq_mutect_mem.R
    ├── fastq_mutect_mem.def
    ├── fastq_mutect_mem.pdf
    ├── fastq_mutect_old.R
    ├── fastq_mutect_singlenode.def
    ├── fastq_preprocess.R
    ├── fastq_preprocess.def
    └── fastq_preprocess.pdf
```

## installation

```
# install dependencies
Rscript -e 'install.packages("flowr", repos = "http://cran.rstudio.com")'
Rscript -e 'install.packages("ngsflows", , repos = c(CRAN="http://cran.rstudio.com", DRAT="http://sahilseth.github.io/drat"))'

# get somatic variant calling workflows
git clone https://github.com/sahilseth/ultraseq.git

cd ultraseq
chmod u+x ultraseq
bin/ultraseq
```
