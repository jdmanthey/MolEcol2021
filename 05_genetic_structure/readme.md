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




## Estimating Genetic Structure from SNPs

Next we will be working with a subset of genomic SNPs (n = 5000) from the Brown Creeper populations (the certhia_contact.stru file). Here is a map of the
subspecies again as well as a sampling map of all the populations used for the SNP dataset we are using today. In Panel B, 
green areas indicate places with dense vegetation, likely different types of forest.

![distribution](https://github.com/jdmanthey/MolEcol2019/blob/master/06_genetic_structure/sampling.png)

Open RStudio and set the working directory to the same one we have been using and download the .stru file from this 
directory on GitHub. 

We'll need to install a package for this exercises, do that here.

    install.packages("adegenet")
    
Load the libraries:
    
    library("adegenet")
    
### Using PCA

We will use a method called discriminant analysis of principal components (DAPC). This is an extension of the principle
components analysis (PCA) that we discussed in class. If you want more info about this method, a link to the paper is here:
https://bmcgenet.biomedcentral.com/articles/10.1186/1471-2156-11-94

To get started we load a structure-formatted file into R. If you are interested in what that looks like, you can open the 
file in a text editor and check it out.

    x <- read.structure("certhia_contact.stru",onerowperind=F,n.ind=24,n.loc=5000,ask=F,sep="\t")

Next, we'll use the DAPC program to identify the number of potential clusters in the genetic data. Here, we are setting the
maximum number of clusters to 6 (number of populations) and running for 1e5 iterations. Here, when the program asks you how 
many PCs to retain, choose 20. This is an amount that maintains most of the variation in the data and is less than the 
number of individuals we have sampled. Then, the program will show you a Bayesian Information Criterion (BIC) plot. The BIC
value will be _lowest_ where the number of genetic clusters is most supported. Choose that number of K (hopefully = 2). 

    grp <- find.clusters(x,max.n.clust=6,n.iter=1e5)
    
Next, we will choose the number of PCs again (choose 20 again) as well as the number of discriminant factors to include, which
has to be less than K. Choose the highest number you can based on the value of K you chose.
    
    dapc1 <- dapc(x,grp$grp)

Now, we'll plot the DAPC results in two ways. The first will be the principal components themselves, colored to the group that
each individual (here each point) belongs to:
    
    plot(dapc1$tab[grp$grp==1,1:2], xlim=c(min(dapc1$tab[,1]), max(dapc1$tab[,1])), ylim=c(min(dapc1$tab[,2]), max(dapc1$tab[,2])), col="blue", pch=19)
    points(dapc1$tab[grp$grp==2,1:2], col="orange2", pch=19)

We can also look at a STRUCTURE-like type of plot showing the assignment of each individual to the two genetic clusters:
    
    compoplot(dapc1, col=c("blue", "orange2"))

If you already looked at the text input file, you know we have three individuals per population from 8 localities:
    1. Utah (individuals 1-3)
    2. Pinal Mountains (4-6)
    3. Pinaleno Mountains (7-9)
    4. Santa Catalina Mountains (10-12)
    5. Chiricahua Mountains (13-15)
    6. Santa Rita Mountains (16-18)
    7. Huachuca Mountains (19-21)
    8. Central Mexico (22-24)

