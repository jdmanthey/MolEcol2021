## Estimating Genetic Differentiation from Genomic Data 30 September 2021

Today we will be working with the same genomic data as last week. Refer to the information about populations, etc. from last
week regarding the source of the data.

If you still have the directory with the data, set that to your current working directory. If not, repeat the download and 
set up from last week's exercise.

Again for today, you will need to coordinate with one or more people to make sure you can compare the patterns seen across
the two chromosomes we are looking at.

&nbsp;

Load the required package.
    
    library("PopGenome")

Read in the VCF file:

    x <- readVCF("chr1_1Mbp.recode.vcf.gz", frompos=1, topos=1000000, tid="Ca_0002_Tg_1", numcols=1000000)
    
    x <- readVCF("chr19_1Mbp.recode.vcf.gz", frompos=1, topos=1000000, tid="Ca_0029_Tg_19", numcols=1000000)

Set the populations (and make sure you have 4):

    x <- set.populations(x, list(c("7", "8", "9"), c("10", "11", "12"), c("19", "20", "21"), c("22", "23", "24")), diploid=T)
    
    x@populations
    
### Distribution

To give you a better understanding of where these samples are coming from in a taxonomic viewpoint, below is a map
showing the distribution of the Brown Creeper (_Certhia americana_). This species has two main lineages, denoted here by the
blue and red colors:

![distribution](https://github.com/jdmanthey/MolEcol2019/blob/master/05_genetic_differentiation/distribution.png)

### Whole region stats

Now that we have our populations defined, we can calculate differentiation within and between populations. This takes 
two commands to do:

    x <- F_ST.stats(x, mode="nucleotide")
    
    x <- diversity.stats(x)

Lets look again at nucleotide diversity between populations, which is another term for D<sub>XY</sub> (covered in lecture).
What we are calculating here is the mean nucleotide diversity between populations per bp, or the mean D<sub>XY</sub>. I would 
recommend you take notes of the values for D<sub>XY</sub> (and in the next step F<sub>ST</sub>) for comparison later on
with the sliding window analyses.

    x@nuc.diversity.between / x@n.sites

That was a measure of absolute genetic differentiation. Now let's look at our relative measure of genetic differentiation, 
F<sub>ST</sub>, between all subpopulations:

    x@nuc.F_ST.pairwise

We can plot the two measures of differentiation relative to one another:

    plot(x@nuc.F_ST.pairwise, x@nuc.diversity.between / x@n.sites, pch=19)

### Sliding windows

Now we will look at sliding windows, just like last week. First, we will  break up the entire segment we are 
investigating into two window objects:

    x_slide1 <- sliding.window.transform(x, width=25000, jump=25000, type=2)
    x_slide2 <- sliding.window.transform(x, width=100000, jump=100000, type=2)

Now, we'll calculate the same statistics as with the whole region for the windows:

    x_slide1 <- F_ST.stats(x_slide1, mode="nucleotide")
    x_slide1 <- diversity.stats(x_slide1)
    x_slide2 <- F_ST.stats(x_slide2, mode="nucleotide")
    x_slide2 <- diversity.stats(x_slide2)

Define the middle position of each of the sliding windows:
The numbers added to each item are half the window size to get the middle position.

    sliding_window_middle1 <- as.numeric(sapply(strsplit(x_slide1@region.names, " -"), "[[", 1)) + 12500
    sliding_window_middle2 <- as.numeric(sapply(strsplit(x_slide2@region.names, " -"), "[[", 1)) + 50000

Now, we'll do the plotting to explore the data. First, we will look at the patterns of D<sub>XY</sub>. The dotted lines are 
comparisons between the different lineages, and the solid lines are comparisons within lineages.

    # Set up plotting dimensions with 2 rows and 1 column
    par(mfrow=c(2,1))
    # Set up the plotting dimensions with axis and chart labels
    plot(c(1,1), col="white", ylim=c(0,0.01), xlim=c(0,1000000), xlab="Position(bp)", ylab="DXY", main="25 kbp windows")
    # Add lines for each of the four populations, each are colored by the col= command in their lines
    # the choice to use points or lines here is arbitrary
    lines(sliding_window_middle1, x_slide1@nuc.diversity.between[1,] / 25000, lwd=0.9, lty=1, col="black")
    lines(sliding_window_middle1, x_slide1@nuc.diversity.between[2,] / 25000, lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle1, x_slide1@nuc.diversity.between[3,] / 25000, lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle1, x_slide1@nuc.diversity.between[4,] / 25000, lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle1, x_slide1@nuc.diversity.between[5,] / 25000, lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle1, x_slide1@nuc.diversity.between[6,] / 25000, lwd=0.9, lty=1, col="black")
    # Set up the plotting dimensions and labels for the second set of windows
    plot(c(1,1), col="white", ylim=c(0,0.01), xlim=c(0,1000000), xlab="Position(bp)", ylab="DXY", main="100 kbp windows")
    # Add lines for each of the four populations, each are colored by the col= command in their lines
    lines(sliding_window_middle2, x_slide2@nuc.diversity.between[1,] / 100000, lwd=0.9, lty=1, col="black")
    lines(sliding_window_middle2, x_slide2@nuc.diversity.between[2,] / 100000, lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle2, x_slide2@nuc.diversity.between[3,] / 100000, lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle2, x_slide2@nuc.diversity.between[4,] / 100000, lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle2, x_slide2@nuc.diversity.between[5,] / 100000, lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle2, x_slide2@nuc.diversity.between[6,] / 100000, lwd=0.9, lty=1, col="black")

Next, we will look at patterns of F<sub>ST</sub>:

    # Set up plotting dimensions with 2 rows and 1 column
    par(mfrow=c(2,1))
    # Set up the plotting dimensions with axis and chart labels
    plot(c(1,1), col="white", ylim=c(-0.2,1), xlim=c(0,1000000), xlab="Position(bp)", ylab="FST", main="25 kbp windows")
    # Add lines for each of the four populations, each are colored by the col= command in their lines
    # the choice to use points or lines here is arbitrary
    lines(sliding_window_middle1, x_slide1@nuc.F_ST.pairwise[1,], lwd=0.9, lty=1, col="black")
    lines(sliding_window_middle1, x_slide1@nuc.F_ST.pairwise[2,], lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle1, x_slide1@nuc.F_ST.pairwise[3,], lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle1, x_slide1@nuc.F_ST.pairwise[4,], lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle1, x_slide1@nuc.F_ST.pairwise[5,], lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle1, x_slide1@nuc.F_ST.pairwise[6,], lwd=0.9, lty=1, col="black")
    # Set up the plotting dimensions and labels for the second set of windows
    plot(c(1,1), col="white", ylim=c(-0.2,1), xlim=c(0,1000000), xlab="Position(bp)", ylab="FST", main="100 kbp windows")
    # Add lines for each of the four populations, each are colored by the col= command in their lines
    lines(sliding_window_middle2, x_slide2@nuc.F_ST.pairwise[1,], lwd=0.9, lty=1, col="black")
    lines(sliding_window_middle2, x_slide2@nuc.F_ST.pairwise[2,], lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle2, x_slide2@nuc.F_ST.pairwise[3,], lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle2, x_slide2@nuc.F_ST.pairwise[4,], lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle2, x_slide2@nuc.F_ST.pairwise[5,], lwd=0.9, lty=2, col="black")
    lines(sliding_window_middle2, x_slide2@nuc.F_ST.pairwise[6,], lwd=0.9, lty=1, col="black")
    
Next, we will take a look at the matrices of the absolute and pairwise genetic differentiation for the larger windows. The 
output shows the differentiation measure for each pairwise comparison in rows and each of the windows in the 10 columns.
   
    x_slide2@nuc.F_ST.pairwise
    x_slide2@nuc.diversity.between / 100000

We can set up a bunch of plots to look at the relationship between the two measures by making a loop across the 6 comparisons.
Here, each point is a comparison of D<sub>XY</sub> (x-axis) and F<sub>ST</sub> (y-axis) for each window.
    
    #set up the plotting margins to be small to fit all plots
    par(mar=c(2,2,3,1))
    # set up number of plots
    par(mfrow=c(3,2))
    for(a in 1:6) {
      # Set up the plotting dimensions 
      plot(c(1,1), col="white", ylim=c(-0.2,1), xlim=c(0,0.01), main=rownames(x_slide2@nuc.diversity.between)[a])
      # plot the pairwise comparison for one of 'a' comparisons, fst on y axis, dxy on x axis
      points(x_slide2@nuc.diversity.between[a,]  / 100000, x_slide2@nuc.F_ST.pairwise[a,], pch=19)
    }

