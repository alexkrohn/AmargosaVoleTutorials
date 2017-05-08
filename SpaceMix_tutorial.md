SpaceMix Analyses from ANGSD Outputs
================

This tutorial will take you through the steps necessary to run a SpaceMix analysis from data run through ANGSD. Although ANGSD is not often used to call SNPs, it can do so quite easily with a variety of filters. To begin this tutorial, you will need to call SNPs using the ANGSD flag **-doGeno 2**. There is also an excellent [SpaceMix vignette](https://github.com/gbradburd/SpaceMix/blob/master/vignettes/spacemix_vignette.Rmd) from the author available online to guide you through using SpaceMix

Importing and Formatting ANGSD Data for SpaceMix
------------------------------------------------

Begin by importing your unzipped -doGeno 2 output. This is a very large matrix with each row being a locus, and each column being the genotype of one individual.

``` r
gen.matrix <- read.table("~/Documents/Berkeley/VoleProject/bestrad2016/Analysis/aligned_radtools_noclones/k90_outputs/spacemix/all_wild_inds_dogeno2_minInd.geno")
```

SpaceMix will not accept loci labels or position numbers (the first two columns of the -doGeno 2 output). Let's remove those. SpaceMix also calls for each row to be an individual and each column to be a locus, so let's transpose the matrix. I have 54 individuals in my analysis, so now my data are 98,295 SNPs from 54 individuals.

``` r
gen.matrix <- gen.matrix[, 3:56]
gen.matrix.raw <- t(gen.matrix)
dim(gen.matrix.raw)
```

    ## [1]    54 98295

Next, SpaceMix requires genotypes to be random with respect to the ancestral or derived allele. Because I gave the **-anc** flag when calling SNPs, calls of 1 or 2 are either heterozygous or homozygous, respectively, different from the ancestral 0 allele. So, let's make a function that can randomly choose an allele to use at the "ancestral" 0 state, and then change the other alleles accordingly. Thanks to Gideon Bradburg for help with this!

``` r
random.switcharoo <- function(x, sample.size = 2) {
    x <- ifelse(rep(runif(1) < 0.5, length(x)), x, sample.size - x)
    return(x)
}
switcharoo.data <- function(frequencies) {
    frequencies <- apply(frequencies, 2, random.switcharoo)
    return(frequencies)
}
```

Now run the function on our genotype matrix. We'll use the actual raw data to know which sites have missing data, so let's make a copy of the matrix, then run the function.

``` r
gen.matrix.copy <- gen.matrix.raw
gen.matrix.switcharoo <- switcharoo.data(gen.matrix.copy)
```

Right. So now we've got our switched data in gen.matrix.switcharoo and our raw data in gen.matrix.raw. SpaceMix uses a sample matrix to determine whether a site has missing data or not. So, for all the sites missing data in the genotype matrix, we need to replace the -1 with a 0.

``` r
gen.matrix.switcharoo[which(gen.matrix == -1)] <- 0
```

Making the sample matrix is now easy. It's just the raw data matrix, but with any site that was sucessfully genotyped as a 2, and any site that has missing data with a 0. (Technically, these are the number of alleles per site that were successfully genotyped!)

NB: because the switcharoo function didn't change whether or not a site was sequenced, just which allele was derived or ancestral, it doesn't matter if I use the raw data or the switcharoo data to make the samples matrix.

``` r
samples <- gen.matrix.raw
samples[samples >= 0] <- 2
samples[samples == -1] <- 0
```

As one last check, let's make sure that no individuals have lots of missing data. We can do this by summing the number of sampled alleles for each individual and making a histogram.

``` r
hist(rowSums(samples)/(2 * ncol(samples)))
```

![](SpaceMix_tutorial_files/figure-markdown_github/unnamed-chunk-7-1.png)

Well, that's not too good. Let's remove the individuals with more than 50% missing data from both the sample and genotype matrices so those individuals' missing data doesn't throw off the analysis.

``` r
newsamples <- samples[which(rowSums(samples)/(2 * ncol(samples)) > 0.5), ]
new.gen.matrix <- gen.matrix.switcharoo[which(rowSums(samples)/(2 * ncol(samples)) > 
    0.5), ]
hist(rowSums(newsamples)/(2 * ncol(newsamples)))
```

![](SpaceMix_tutorial_files/figure-markdown_github/unnamed-chunk-8-1.png)

Much better. Now we're down to 45 individuals.

``` r
dim(new.gen.matrix)
```

    ## [1]    45 98295

``` r
dim(newsamples)
```

    ## [1]    45 98295

Load the Geographic Locations
-----------------------------

Next, load in the latitude and longitude positions for each individual. Crucially, the order of this matrix must match the order of the individuals in the genotype matrix.

``` r
geo <- read.table("~/Documents/Berkeley/VoleProject/bestrad2016/Analysis/aligned_radtools_noclones/k90_outputs/spacemix/spatial_locations_allwildvoles.txt", 
    sep = "\t", header = T)
head(geo)
```

    ##        lat      long
    ## 1 35.88056 -116.2319
    ## 2 35.88037 -116.2318
    ## 3 35.88032 -116.2320
    ## 4 35.85505 -116.2324
    ## 5 35.85582 -116.2325
    ## 6 35.88500 -116.2342

I also have to remove the individuals that had too much missing data from these geographic positions.

``` r
newgeo <- geo[which(rowSums(samples)/(2 * ncol(samples)) > 0.5), ]
dim(newgeo)
```

    ## [1] 45  2

Run the SpaceMix Analysis
-------------------------

You now have all of the required input data for SpaceMix. Not so bad, was it? Now, let's run the analysis.

``` r
library(SpaceMix)
run.spacemix.analysis(n.fast.reps = 10, fast.MCMC.ngen = 10, fast.model.option = "source_and_target", 
    long.model.option = "source_and_target", data.type = "counts", counts = new.gen.matrix, 
    sample.sizes = newsamples, spatial.prior.X.coordinates = newgeo$long, spatial.prior.Y.coordinates = newgeo$lat, 
    round.earth = F, k = nrow(new.gen.matrix), loci = ncol(new.gen.matrix), ngen = 1.3e+07, 
    printfreq = 100, samplefreq = 1000, mixing.diagn.freq = 50, savefreq = 1e+05, 
    prefix = "admix_and_location_28Nov")
```

This code reduced my dataset to only 3,801 informative SNPs, and took a few days to run. I also had to tinker a bit with the parameters to ensure that I had adequate mixing and chain convergence. Definitely familiarize yourself with how to evaluate MCMC runs. [Gideon's excellent tutorial](https://github.com/gbradburd/SpaceMix/blob/master/vignettes/spacemix_vignette.Rmd) is a great place to start learning what to look for in the various plots.

I'll assume your model is working well, and continue instead with plotting the output.

Making a Pretty SpaceMix Plot
-----------------------------

SpaceMix requires a matrix of names for each of your individuals. Obviously, they should be in the same order as your genotype and lat/long matrices. SpaceMix is picky, so make sure everything is a matrix.

``` r
sample.names <- read.table("~/Documents/Berkeley/VoleProject/bestrad2016/Analysis/aligned_radtools_noclones/k90_outputs/spacemix/sample_names.txt", 
    sep = "\t")
sample.names <- sample.names[which(rowSums(samples)/(2 * ncol(samples)) > 0.5), ]
length(sample.names)
```

    ## [1] 45

``` r
sample.names <- as.matrix(sample.names)
newgeo <- as.matrix(newgeo[, c("long", "lat")])
```

Next we need a vector of coors for each individual. The tutorial gives a nifty way of coloring your individuals based on differing hues based on differences in admixture. However, I want to color mine by whether or not they belong to the northern clade (purple) or the southern clade (blue).

``` r
colors <- as.vector(c(rep("gray56", 26), rep("dodgerblue3", 5), rep("purple", 13), 
    rep("dodgerblue3", 1)))
color.vector <- as.vector(colors)
```

Finally, put it all together to make the map list with 95% confidence intervals around the geogenetic location. This command does take a few minutes to run on my machine.

``` r
spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "~/Documents/Berkeley/VoleProject/bestrad2016/Analysis/aligned_radtools_noclones/k90_outputs/spacemix/run_433172/admix_and_location_28Nov_LongRun/admix_and_location_28Nov_space_MCMC_output1.Robj", 
    geographic.locations = newgeo, name.vector = sample.names, color.vector = colors, 
    quantile = 0.95, burnin = 0)
```

For the actual plotting, instead of just throwing a whole box of code and adding comments to the side, I'll break the code down step by step without running it, then run it all at the end in a nice pretty plot.

First, make the actual plot.

``` r
make.spacemix.map(spacemix.map.list = spacemix.map.list, text = F, ellipses = T, 
    source.option = F, ylim = c(32.8, 37.5), xlim = c(-120.8, -114))
```

Next, to make it a bit nicer, we can take one individual from each population and plot their original location on the map. This will give us the geographic and geogenetic sampling locations. First, we'll make a list of latitude, longitude and colors to include for each population.

``` r
population.lats <- newgeo[c(1, 28, 33, 39, 43, 45), 2]
population.longs <- newgeo[c(1, 28, 33, 39, 43, 45), 1]
population.colors <- colors[c(1, 28, 33, 39, 43, 45)]
```

We can even add text labels next to those geographic sampling locations to label each population.

``` r
points(y = population.lats, x = population.longs, pch = 19, col = population.colors, 
    cex = 0.9)
text(x = population.longs + c(0, 0, -0.2, 0, -0.52, 0), y = population.lats + c(-0.1, 
    0.1, 0.1, -0.1, 0, 0.1), pch = 19, col = population.colors, labels = c("Tecopa", 
    "Victorville", "Deep Springs", "Lark Seep", "Lake Isabella", "Elizabeth Lake"), 
    cex = 0.9)
```

Finally, the outlines of California, Nevada and Arizona, to give the map context.

``` r
library(maps)
map("state", c("california", "nevada", "arizona"), add = T)
```

This last part is not very precise, but does work. We can add arrows going from the geographic sampling locations to the **approximate** geogenetic center of all of the individuals.

``` r
library(shape)
Arrows(x0 = population.longs, y0 = population.lats, x1 = c(-116.55, -117.25, -117.12, 
    -117.35, -117.3, -117.3), y1 = c(36.01, 34, 36.71, 36, 36.09, 34.4), col = population.colors, 
    lwd = 2, arr.length = 0.5)
```

![](SpaceMix_tutorial_files/figure-markdown_github/unnamed-chunk-21-1.png)
