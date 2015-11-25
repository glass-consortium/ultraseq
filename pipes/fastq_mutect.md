## get data from chub:

Tumor Sample: 
Sequence alignment/map of sample CCLE-HCC1143-DNA-08 with 22 read groups providing alignments to HG19_Broad_variant (Human reference GRCh37) 
89324e86-5b7a-4f69-92c1-3b67293f8748

Normal Sample (non-tumor): 
Sequence alignment/map of sample CCLE-HCC1143 BL-DNA-08 with 16 read groups providing alignments to HG19_Broad_variant (Human reference GRCh37) 
ebdb53ae-6386-4bc4-90b1-4f249ff9fcdf

## Convert bam files to fastq

Issue with converting sam to fastq

http://sourceforge.net/p/samtools/mailman/message/29626514/

```
> SAM	154,555,012
> R1 	81,045,486
> R2 	73,509,526
> R1+R2	154,555,012
>
> When I count the  first and second reads using samtools flags
> 0x0040	1	the read is the first read in a pair - 73509526
> 0x0080	2	the read is the second read in a pair - 73509526
>
> I see that the number match, and also match to R2
>
> What is the reason of Picard extracting fastq differently and where
> extra reads in R1 come from.
>
```

Have included appropriate parameters as part of this function:

```
bam=402285b1-01fb-4ae6-8cdc-aab1d479f31b/TCGA-A6-6141-01A-11D-1771-10_wgs_Illumina.bam
flowr ngsflows::picard_bam_fastq x=$bam \
 samplename=$bam execute=TRUE
java -jar /risapps/noarch/picard/1.138/picard.jar SamToFastq \
 INPUT=C835.HCC1143.2.bam \ 
 FASTQ=C835.HCC1143.2_1.fastq \
 SECOND_END_FASTQ=C835.HCC1143.2_2.fastq
```

## read were split using:
```
script=~/Dropbox/public/github_ngsflows/inst/scripts/split_fq.sh
$script -n 4000000 -f C835.HCC1143_BL.4_1.fastq &
$script -n 4000000 -f C835.HCC1143_BL.4_2.fastq
```

```
cd /rsrch2/iacs/iacs_dep/sseth/flowr/wex/ebdb53ae-6386-4bc4-90b1-4f249ff9fcdf


```

Details on the tumor sample:

```
samtools flagstat 89324e86-5b7a-4f69-92c1-3b67293f8748/C835.HCC1143.2.bam
87572043 + 9184798 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
16337477 + 2161424 duplicates
86787319 + 7493523 mapped (99.10%:81.59%)
78931502 + 8337050 paired in sequencing
39465751 + 4168525 read1
39465751 + 4168525 read2
77280418 + 6544800 properly paired (97.91%:78.50%)
77593988 + 6587948 with itself and mate mapped
596328 + 256806 singletons (0.76%:3.08%)
273026 + 39042 with mate mapped to a different chr
234383 + 32595 with mate mapped to a different chr (mapQ>=5)
```

```
samtools flagstat ebdb53ae-6386-4bc4-90b1-4f249ff9fcdf/C835.HCC1143_BL.4.bam
68629600 + 6562518 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
10054468 + 1517557 duplicates
67842739 + 5593779 mapped (98.85%:85.24%)
68629600 + 6562518 paired in sequencing
34314800 + 3281259 read1
34314800 + 3281259 read2
66854380 + 5314442 properly paired (97.41%:80.98%)
67196156 + 5353240 with itself and mate mapped
646583 + 240539 singletons (0.94%:3.67%)
301196 + 35316 with mate mapped to a different chr
260127 + 29485 with mate mapped to a different chr (mapQ>=5)
```




## using WGS files:
```
script=/rsrch2/iacs/iacs_dep/sseth/Dropbox/public/github_ngsflows/inst/scripts/split_fq.sh
python=/risapps/rhel6/python/2.7.6/anaconda/bin/python2.7
bamtofastq=/risapps/src6/speedseq/bin/bamtofastq.py
cd=/rsrch2/iacs/iacs_dep/sseth/flowr/wgs/402285b1-01fb-4ae6-8cdc-aab1d479f31b
bam=TCGA-A6-6141-01A-11D-1771-10_wgs_Illumina.bam
$python $bamtofastq -r 0,0.1,0.2  -i $bam | $script -f /dev/stdin -n 4000000 -o TCGA-A6-6141-01A-11D-1771-10_wgs_Illumina_ -z
```


## split up data for evaluation:
```
sftp://hms00/rsrch2/iacs/iacs_dep/sseth/flowr/wgs/402285b1-01fb-4ae6-8cdc-aab1d479f31b/split
```