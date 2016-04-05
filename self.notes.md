


Adding this to DRAT repo

```{r}

# insert package
pkg = devtools::build(".")
out = drat::insertPackage(pkg, "~/Dropbox/public/github_drat", commit = TRUE)

# push repo
system("cd ~/Dropbox/public/github_drat;git push")


```


## renaming ngsflows to ultraseq

flowr run




flowr run
errors saying a patterns is required


funr, expectations.

funr: should fail if
1. no function by that name
2. function does not run properly


# error regarding no valid argument. x=ultraseq
flowr run ultraseq

# error regarding no valid argument
# no bams!
flowr run x=ultraseq
