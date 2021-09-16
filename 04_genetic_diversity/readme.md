## Estimating Genetic Diversity from Genomic Data 16 September 2021

Today we will be working with some genomic data in order to estimate population genetic diversity within and between populations.
First, do the following items:
1. Create a directory called 'diversity' somewhere (e.g., on the desktop).
2. Open up RStudio and set your working directory to the newly-created directory.
3. Download the 4 data files (they have .gz and .tbi suffixes) and put them in your new directory. 

Once you have your working directory set up and have your files there, you will need to install a package necessary for today:

    install.packages("PopGenome")

If the package is already installed it will tell you. Otherwise, you will see a bunch of information on screen as the package
is installed. After successful installation, you need to load the package:
    
    library("PopGenome")

One optional step (which I often prefer) is to remove scientific notation of large numbers. It is up to you if you want to do 
this, but I often find it helpful. This makes all the numbers always print out fully:

    options(scipen=999)

&nbsp;

For today's exercises, there are two potential datasets to use, for chromosomes 1 and 19 (the first 100 Mbp), of an organism
called the Brown Creeper (_Certhia americana_). For some of the questions on the exercises, you will need to communicate and 
compare results with one or more other people that have analyzed the opposite dataset. Check to make sure that your potential
partner(s) chooses the other dataset relative to you for comparison.

Once you have this figured out with your partner(s), read in your respective dataset with the function 'readVCF.' The function
uses several arguments:
1. the input vcf file (gzipped)
2. the starting position you want to read (here = 1)
3. the end position you want to read (here = 1,000,000)
4. the ID of the chromosome (full name)
5. the maximum number of SNPs you want to read (numcols)

&nbsp;

    # choose one or the other:
    x <- readVCF("chr1_1Mbp.recode.vcf.gz", frompos=1, topos=1000000, tid="Ca_0002_Tg_1", numcols=1000000)
    
    x <- readVCF("chr19_1Mbp.recode.vcf.gz", frompos=1, topos=1000000, tid="Ca_0029_Tg_19", numcols=1000000)

Normally in R, when you just type the name of the object (here we named it simply "x"), R displays the whole object. In 
contrast, objects created in PopGenome are often very large, and the code writers decided to make a summary of potential
commands, and then how you can look at the results, when you call an object. Try it here to see what that looks like:

    x

You can also look at the structure of an object, to see what components it has. In PopGenome, when you create a "genome-class"
object, it automatically makes many potential slots for statistics, but doesn't automatically calculate them. Check this out
with the following command:

    str(x)

We won't calculate all of these statistics, but it may be of interest to some of you that there are this many possibilities. 
All of the information for all of these stats can be found in the PopGenome package manual.

Next, we will define the populations sampled in the dataset. There are four populations each with 3 individuals: (1) Santa 
Rita Mountains in Arizona, (2) Morelos (state) in Mexico, (3) Pinal Mountains in Arizona, (4) Utah (state) in USA. The 
individuals are named with numbers in the vcf files. You can see the names of the individuals and set the populations with
the 2 following commands. Individuals marked with "-2" are the second allele for these individuals (they are diploid).

    get.individuals(x)
    
    x <- set.populations(x, list(c("7", "8", "9"), c("10", "11", "12"), c("19", "20", "21"), c("22", "23", "24")), diploid=T)

Now that we have our populations defined, we can calculate diversity within and between populations. This takes two commands
to do:

    x <- F_ST.stats(x, mode="nucleotide")
    
    x <- diversity.stats(x)

Now we can look at the mean diversity across the entire 1Mbp. We can look at the nucleotide diversity within populations
with the following command. Remember that the nucleotide diversity is the mean number of differences between individuals
sampled. We can look at the total nucleotide diversity (across 1,000,000 bp), and we can also look at the average per bp
nucleotide diversity:

    x@nuc.diversity.within
    
    x@nuc.diversity.within / x@n.sites

We can also look at nucleotide diversity between populations, or the average number of differences between one randomly 
sampled individual from each population. Again, in total and average per bp:

    x@nuc.diversity.between
    
    x@nuc.diversity.between / x@n.sites

Now that we've looked at the mean patterns across this whole length of sequence we are analyzing, we will look to see if there
are different patterns in sliding windows across the region. Sliding windows in the case we are setting up are looking at
non-overlapping segments of the chromosome. Here, we will set up two different-sized sliding windows: 25 kbp and 100 kbp. 
These are arbitrary numbers, but are representative of possible window sizes we may investigate across an entire genome. We
will first break up the entire segment we are investigating into two window objects:

    x_slide1 <- sliding.window.transform(x, width=25000, jump=25000, type=2)
    
    x_slide2 <- sliding.window.transform(x, width=100000, jump=100000, type=2)

We can look at the defined regions for both of our sliding window objects using the following commands:

    x_slide1@region.names
    
    x_slide2@region.names

Now, we'll calculate the same statistics as before, except for each of the sliding windows:

    x_slide1 <- F_ST.stats(x_slide1, mode="nucleotide")
    
    x_slide1 <- diversity.stats(x_slide1)
    
    x_slide2 <- F_ST.stats(x_slide2, mode="nucleotide")
    
    x_slide2 <- diversity.stats(x_slide2)

If we look at the within population diversity for each object, we can see the nucleotide diversity for each population for 
each of the sliding windows:
    
    x_slide1@nuc.diversity.within / 25000
    
    x_slide2@nuc.diversity.within / 100000

Now, we'll plot diversity across the sliding windows. First, we'll define the middle position of each of the sliding windows:
The numbers added to each item are half the window size to get the middle position.

    sliding_window_middle1 <- as.numeric(sapply(strsplit(x_slide1@region.names, " -"), "[[", 1)) + 12500
    
    sliding_window_middle2 <- as.numeric(sapply(strsplit(x_slide2@region.names, " -"), "[[", 1)) + 50000
    
Now, we'll do the plotting. I have added comments below to show each step. Remember, comments are preceded by the "#." If you 
are interested in any of the specific settings with plotting, remember the R tutorials and look up the plot, lines functions
to see what each of the commands do.

    # Set up plotting dimensions with 2 rows and 1 column
    par(mfrow=c(2,1))
    # Set up the plotting dimensions with axis and chart labels
    plot(c(1,1), col="white", ylim=c(0,0.01), xlim=c(0,1000000), xlab="Position(bp)", ylab="Nucleotide Diversity", main="25 kbp windows")
    # Add lines for each of the four populations, each are colored by the col= command in their lines
    # the choice to use points or lines here is arbitrary
    lines(sliding_window_middle1, x_slide1@nuc.diversity.within[,1] / 25000, lwd=0.9, col="darkblue")
    lines(sliding_window_middle1, x_slide1@nuc.diversity.within[,2] / 25000, lwd=0.9, col="darkorange2")
    lines(sliding_window_middle1, x_slide1@nuc.diversity.within[,3] / 25000, lwd=0.9, col="deeppink3")
    lines(sliding_window_middle1, x_slide1@nuc.diversity.within[,4] / 25000, lwd=0.9, col="black")
    # Set up the plotting dimensions and labels for the second set of windows
    plot(c(1,1), col="white", ylim=c(0,0.01), xlim=c(0,1000000), xlab="Position(bp)", ylab="Nucleotide Diversity", main="100 kbp windows")
    # Add lines for each of the four populations, each are colored by the col= command in their lines
    lines(sliding_window_middle2, x_slide2@nuc.diversity.within[,1] / 100000, lwd=0.9, col="darkblue")
    lines(sliding_window_middle2, x_slide2@nuc.diversity.within[,2] / 100000, lwd=0.9, col="darkorange2")
    lines(sliding_window_middle2, x_slide2@nuc.diversity.within[,3] / 100000, lwd=0.9, col="deeppink3")
    lines(sliding_window_middle2, x_slide2@nuc.diversity.within[,4] / 100000, lwd=0.9, col="black")







