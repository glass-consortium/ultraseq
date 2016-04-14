# Contributing

Each file in the `ultraseq/R` folder represents wrapper to a tool OR module. The file has a function by the same name, and returns a set of commands to run the function. For example, `mutect.R`, has a function `mutect()`.

Multiple modules, can be stitched together to form pipelines - a few examples are in [pipelines](https://github.com/flow-r/ultraseq/tree/master/pipelines). And a few experimental pipelines in [pipelines/extra](https://github.com/flow-r/ultraseq/tree/master/pipelines/extra).

We use the following overarching ideas for modules - enabling consistent return expections, for each of them.

## A few suggested specifications, for modules:


1. **inputs**:
A module function should accept minimum of **two inputs**, 
    - samplename (is used to append a column to the flowmat)
    - first argument, should preferably be an input file. For example fastq, bam etc. This really depends on the module. 
    - One choose argument names from [list of common arguments](https://github.com/flow-r/ultraseq/blob/master/list_of_common_args.md) (under active development - the list will change!)
  ```
  samplename = opts_flow$get("samplename")
  
  # a few example inputs
  bam
  fqs
  
  tumor_bam
  normal_bam
  ```

2. **return**
A module function should always return a list arguments:
    - **flowmat** (required)   : contains all the commands to run
    - **outfiles** (recommended): could be used as an input to other tools

  ```
  return(list(outfiles = bam, flowmat = flowmat))
  ```

3. **conf file**:
Preferably, default values for all other arguments should be read from the `conf` file.
For example a pipeline called `bam_mutect`, will have a conf file `bam_mutect.conf`, in the same folder. All params from that file would be read, and can be fetched using `opts_flow$get()`. For example `opts_flow$get('bwa_exe')`

 ```
 ## Example conf file:
 bwa_exe	/apps/bwa/bin/bwa
 ```

4. **checks**

To make sure all arguments have values, and none of them are NULL, a module should use `check_args()`.

For example we have a function:

```
fastq_bwa <- function(..., bwa_exe = opts_flow$get("bwa_exe"){

	# asserts, that none of the arguments are NULL
	check_args()
	# ignore a few:
	check_args(ignore = "bwa_exe")
.
.
.

}

```

## Example

Let, use an example function `picard_merge` which merged bam files (`bams`), and creates a single `mergedbam`.

```{r picard_merge, echo=TRUE, comment=""}
picard_merge <- function(bams,
        samplename = opts_flow$get("samplename"),
        mergedbam,
        java_exe = opts_flow$get("java_exe"),
        java_mem = opts_flow$get("java_mem"),
        java_tmp = opts_flow$get("java_tmp"),
        picard_jar = opts_flow$get("picard_jar")){
	
  ## Make sure all args have a value (not null)
  ## If a variable was not defined in a conf. file opts_flow$get, will return NULL
  check_args()  
  
  # create a vector of bams
  bam_list = paste("INPUT=", x, sep = "", collapse = " ")
  
  # create a named list of commands
  cmds = list(merge = sprintf("%s %s -Djava.io.tmpdir=%s -jar %s MergeSamFiles %s OUTPUT=%s ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true USE_THREADING=true",
  java_exe, java_mem, java_tmp, picard_jar, bam_list, mergedbam))
  
  # Create a flowmat
  flowmat = to_flowmat(cmds, samplename)
  
  # return a list, flowmat AND outfiles
  return(list(outfiles = mergedbam, flowmat = flowmat))
}
```


<!-- Here are a few things to note regarding naming a function and what it should do:

- it is *recommended* to have a lower case function name, seperated by `_`. And in general we try to follow Advanced R's [style guide](http://adv-r.had.co.nz/Style.html).
- Each function has:
  - a few input files,
  - a few paths (to files and tools)
  - and default parameters
- Each function return a list, with elements:
  - outfiles: a list/character vector of output file names
  - flowmat: a data.frame, with a few extra attributes -->
  
