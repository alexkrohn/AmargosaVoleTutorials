Bootstrapping 100 Fst Calculations From ANGSD
================

Fst Bootstrapping
-----------------

One way to determine whether K=1 is significantly different from K=2 is to compare the Fst values of the proposed K=2 scenario to a null expectation of Fst values. To generate this null distribution of Fst values you can randomly combine your individuals into two populations in the same proportions as the K=2 scenario and calculate Fst. Repeat this 100 times (with replacement), and you will be using parametric bootstrapping to generate a null distribution of Fst values without any assumptions about the distribution of Fst values in your population.

In other words, in the example from the paper the proposed K=2 scenario from ngsAdmix was made up of one population with 32 individuals and another with 20 individuals. Using the code below, I will generate Fst values for populations of 32 and 20 randomly selected individuals, and then compare those Fst values to the actual K=2 value. If the K=2 value is greater than 95% of those values, I could say that at a alpha = 0.05 level that the K=2 value is significant.

Generating Bam Lists
--------------------

First, generate 100 bam lists of 32 and 20 individuals. Start in a directory that contains the bam files that you'd like to include. In this example, the directory contains 52 bam files. bamlistA1 through bamlist100 will contain 32 individuals, bamlistB1 through bamlistB100 will contain 20 individuals.

``` r
bams <- list.files(getwd())

for (i in 1:100) {
    bamlistA <- sample(bams, 32, replace = F)
    bamlistB <- setdiff(bams, bamlistA)
    
    write.table(bamlistA, file = paste("bamlistfolder/bamlistA", i, sep = ""), quote = FALSE, 
        row.names = FALSE, col.names = FALSE)
    write.table(bamlistB, file = paste("bamlistfolder/bamlistB", i, sep = ""), quote = FALSE, 
        row.names = FALSE, col.names = FALSE)
    
}
```

Calculate Fst
-------------

The rest of these steps are basically from the [ANGSD webiste](http://popgen.dk/angsd/index.php/Fst). You may not be able to run all of the commands at once. Start within the folder containing the bamlists (bamlistfolder/ above).

``` r
bamlist <- list.files(getwd())

for (i in 1:length(bamlist)) {
    system(paste("angsd -anc kmer60-min500-scaffolds.fa -sites 90percent_keepfile.keep -rf 90percentrf.rf -only_proper_pairs 1 -dosaf 1 -gl 1 -bam ", 
        bamlist[i], " -out ", bamlist[i], "_sfs", sep = ""))
}

files <- list.files(getwd())
idxlist <- files[grep("*.saf.idx", files)]

for (i in 1:100) {
    idx <- idxlist[grep(paste("[A-B]", i, "_", sep = ""), idxlist)]
    system(paste("realSFS ", idx[1], " ", idx[2], " > ", idx[1], "_", idx[2], ".ml", 
        sep = ""))
}

for (i in 1:100) {
    idx <- idxlist[grep(paste("[A-B]", i, "_sfs.saf.idx$", sep = ""), idxlist)]
    system(paste("realSFS fst index ", idx[1], " ", idx[2], " -sfs ", idx[1], "_", 
        idx[2], ".ml -fstout ", idx[1], "_", idx[2], sep = ""))
}
```

Then, in bash, create a list of the 100 Fst values.

``` bash
(for pop in `ls *.fst.idx`
do
realSFS fst stats $pop
done) > 100random_fst_values
```

Find the Null Distribution of Fst Values
----------------------------------------

Load in the 100 values and run a t-test to find the 95% confidence interval (the range of Fst values that contains 95% of the values). If your K=2 value falls outside of this range, it's significantly different from random (i.e. K=1).

``` r
fstvalues <- read.table("100random_fst_values")

t.test(fstvalues$V2)
```
