#Population genetic stats for bluebonnet GBS SNP data
#1/5/2017
#tutorial at http://popgen.nescent.org/StartSNP.html#genetic-diversity-observed-and-expected-heterozygosity
#and http://popgen.nescent.org/DifferentiationSNP.html

library("adegenet")
library("hierfstat")
library("pegas")

bbsnp <- read.structure("dDoccat.FinalSNP.str")
# How many genotypes are there? 
# 187
# How many markers are there? 
# 11741
# Which column contains labels for genotypes ('0' if absent)? 
# 1
# Which column contains the population factor ('0' if absent)? 
# 2
#  Which other optional columns should be read (press 'return' when done)? 
# 1: 3
# 2: 4
# 3: 
#   Which row contains the marker names ('0' if absent)? ###file version on work laptop does not have this row!!!
# 1
# Are genotypes coded by a single row (y/n)? 
# n
# Converting data from a STRUCTURE .stru file to a genind object... 

bbsnp
#' /// GENIND OBJECT /////////
#'   
#'   // 187 individuals; 11,741 loci; 23,615 alleles; size: 22.4 Mb
#' 
#' // Basic content
#' @tab:  187 x 23615 matrix of allele counts
#' @loc.n.all: number of alleles per locus (range: 2-3)
#' @loc.fac: locus factor for the 23615 columns of @tab
#' @all.names: list of allele names for each locus
#' @ploidy: ploidy of each individual  (range: 2-2)
#' @type:  codom
#' @call: read.structure(file = "dDoccat.FinalSNP.str")
#' 
#' // Optional content
#' @pop: population of each individual (group size range: 7-8)
#' @other: a list containing: X 

bbsnp2 <- genind2hierfstat(bbsnp)

####obs and exp het####
div <- summary(bbsnp)
div[2] #group sizes
# $n.by.pop
# 1 10 11 12 13 14 15 16 17 18 19  2 20 21 22 23 24  3  4  5  6  7  8  9 
# 8  8  7  8  8  8  7  7  8  8  8  8  8  8  7  7  8  8  8  8  8  8  8  8

#div[6] #per site obs het
#div[7] #per site exp het

#plots of subset
plot(div$Hobs[1:12], xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")
plot(div$Hobs[1:12],div$Hexp[1:12], xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")

#test Hexp=Hobs
bartlett.test(list(div$Hexp, div$Hobs)) # a test : H0: Hexp = Hobs
# Bartlett test of homogeneity of variances
# 
# data:  list(div$Hexp, div$Hobs)
# Bartlett's K-squared = 2129.1, df = 1, p-value < 2.2e-16

basicstat <- basic.stats(bbsnp2, diploid = TRUE, digits = 6)
basicstat[7] #Estimates mean observed heterozygosities, 
# mean gene diversities within population Hs, 
# Gene diversities overall Ht and corrected Htp, and Dst, Dstp. 
# Finally, estimates Fst and Fstp as well as Fis following Nei (1987) 
# per locus and overall loci
# $overall
# Ho        Hs        Ht       Dst       Htp      Dstp       Fst      Fstp       Fis 
# 0.269910  0.257873  0.261543  0.003670  0.261702  0.003830  0.014032  0.014633 -0.046681 
# Dest 
# 0.005160 

boot.ppfis(bbsnp2) #Fis confidence interval for each pop;  FIS is the inbreeding coefficient of an individual (I) relative to the subpopulation (S)
# positive Fis values indicate that individuals in a population are more related than you would expect under a model of random mating. Inbreeding
# negative Fis values indicate that individuals in a population are less related than you would expect under a model of random mating. Mate choice?
# $call
# boot.ppfis(dat = bbsnp2)
# 
# $fis.ci
#       ll      hl
# 1  -0.0850 -0.0670 neg
# 2  -0.0140  0.0078 zero
# 3  -0.0522 -0.0291 neg
# 4  -0.0402 -0.0169 neg
# 5  -0.0164  0.0020 zero
# 6  -0.0148  0.0088 zero
# 7  -0.0269 -0.0062 neg
# 8  -0.0143  0.0096 zero
# 9  -0.0344 -0.0146 neg
# 10  0.0051  0.0256 pos - seeded
# 11 -0.0708 -0.0487 neg
# 12 -0.0392 -0.0171 neg
# 13 -0.0061  0.0165 zero
# 14 -0.0569 -0.0344 neg
# 15 -0.0984 -0.0770 neg
# 16 -0.0458 -0.0227 neg
# 17 -0.1131 -0.0943 neg
# 18 -0.0478 -0.0297 neg
# 19 -0.1436 -0.1245 neg
# 20 -0.1164 -0.0934 neg
# 21 -0.1012 -0.0771 neg
# 22 -0.1083 -0.0882 neg
# 23 -0.0581 -0.0371 neg
# 24 -0.0522 -0.0274 neg

boot.ppfst(bbsnp2) #conf int for Fst each pop; FST is the effect of subpopulations (S) compared to the total population (T)
#giant matrix/list thing. Takes forever!

boot.vc(bbsnp2)

x <- indpca(bbsnp2) 
plot(x, cex = 0.7)
