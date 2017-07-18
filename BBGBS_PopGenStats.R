#Population genetic stats for bluebonnet GBS SNP data
#1/5/2017
#tutorial at http://popgen.nescent.org/StartSNP.html#genetic-diversity-observed-and-expected-heterozygosity
#and http://popgen.nescent.org/DifferentiationSNP.html

library("adegenet")
library("hierfstat")
library("pegas")

bbsnp <- read.structure("dDoccat.FinalSNP.stru")
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
# 0
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
# boot.vc(gtrunchier[,1],gtrunchier[,-c(1:2)])$ci
boot.vc(bbsnp2[,1],bbsnp2[,-c(1:2)])$ci
#       H-Total(Ht) F-Pop/Total(Fst) F-Ind/Total  H-Pop(Hs)     F-Ind/Pop(Fis)   Hobs(Ho)
# 2.5%   0.2593      0.0142         -0.0399       0.2555   -0.0550      0.2662
# 50%    0.2618      0.0147         -0.0319       0.2580   -0.0473      0.2702
# 97.5%  0.2644      0.0152         -0.0237       0.2605   -0.0393      0.2739

x <- indpca(bbsnp2) 
plot(x, cex = 0.7)

####Calcuating Pi####
#per pop
pi1001 <- read.delim2("dDoccat.FinalSNP.1001Pi.windowed.pi", header = TRUE)
nrow(pi1001)
# [1] 3973
nrow(subset(pi1001, N_VARIANTS > 1))
# [1] 2475
nrow(subset(pi1001, N_VARIANTS > 1))/nrow(pi1001)
# [1] 0.6229549

pi1002 <- read.delim2("dDoccat.FinalSNP.1002Pi.windowed.pi", header = TRUE)
nrow(pi1002)
# [1] 3880
nrow(subset(pi1002, N_VARIANTS > 1))
# [1] 2388
nrow(subset(pi1002, N_VARIANTS > 1))/nrow(pi1002)
# [1] 0.6154639

pi1005 <- read.delim2("dDoccat.FinalSNP.1005Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1005))
# [1] 3908
(b <- nrow(subset(pi1005, N_VARIANTS > 1)))
# [1] 2385
b/a
# [1] 0.6102866

pi1006 <- read.delim2("dDoccat.FinalSNP.1006Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1006))
# [1] 3924
(b <- nrow(subset(pi1006, N_VARIANTS > 1)))
# [1] 2428
b/a
# [1] 0.6187564

pi1007 <- read.delim2("dDoccat.FinalSNP.1007Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1007))
# [1] 3802
(b <- nrow(subset(pi1007, N_VARIANTS > 1)))
# [1] 2254
b/a
# [1] 0.5928459

pi1008 <- read.delim2("dDoccat.FinalSNP.1008Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1008))
# [1] 3841
(b <- nrow(subset(pi1008, N_VARIANTS > 1)))
# [1] 2342
b/a
# [1] 0.609737

pi1101 <- read.delim2("dDoccat.FinalSNP.1101Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1101))
# [1] 3832
(b <- nrow(subset(pi1101, N_VARIANTS > 1)))
# [1] 2306
b/a
# [1] 0.6017745

pi1102 <- read.delim2("dDoccat.FinalSNP.1102Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1102))
# [1] 3845
(b <- nrow(subset(pi1102, N_VARIANTS > 1)))
# [1] 2361
b/a
# [1] 0.6140442

pi1103 <- read.delim2("dDoccat.FinalSNP.1103Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1103))
# [1] 3804
(b <- nrow(subset(pi1103, N_VARIANTS > 1)))
# [1] 2273
b/a
# [1] 0.5975289

pi1104 <- read.delim2("dDoccat.FinalSNP.1104Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1104))
# [1] 3911
(b <- nrow(subset(pi1104, N_VARIANTS > 1)))
# [1] 2390
b/a
# [1] 0.6110969

pi1105 <- read.delim2("dDoccat.FinalSNP.1105Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1105))
# [1] 3823
(b <- nrow(subset(pi1105, N_VARIANTS > 1)))
# [1] 2266
b/a
# [1] 0.5927282

pi1107 <- read.delim2("dDoccat.FinalSNP.1107Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1107))
# [1] 3859
(b <- nrow(subset(pi1107, N_VARIANTS > 1)))
# [1] 2335
b/a
# [1] 0.605079

pi1201 <- read.delim2("dDoccat.FinalSNP.1201Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1201))
# [1] 3893
(b <- nrow(subset(pi1201, N_VARIANTS > 1)))
# [1] 2355
b/a
# [1] 0.6049319

pi1202 <- read.delim2("dDoccat.FinalSNP.1202Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1202))
# [1] 3864
(b <- nrow(subset(pi1202, N_VARIANTS > 1)))
# [1] 2351
b/a
# [1] 0.6084369

pi1203 <- read.delim2("dDoccat.FinalSNP.1203Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1203))
# [1] 3835
(b <- nrow(subset(pi1203, N_VARIANTS > 1)))
# [1] 2293
b/a
# [1] 0.597914

pi1204 <- read.delim2("dDoccat.FinalSNP.1204Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1204))
# [1] 3817
(b <- nrow(subset(pi1204, N_VARIANTS > 1)))
# [1] 2280
b/a
# [1] 0.5973277

pi1205 <- read.delim2("dDoccat.FinalSNP.1205Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1205))
# [1] 3882
(b <- nrow(subset(pi1205, N_VARIANTS > 1)))
# [1] 2380
b/a
# [1] 0.613086

pi1301 <- read.delim2("dDoccat.FinalSNP.1301Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1301))
# [1] 3855
(b <- nrow(subset(pi1301, N_VARIANTS > 1)))
# [1] 2313
b/a
# [1] 0.6

pi1303 <- read.delim2("dDoccat.FinalSNP.1303Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1303))
# [1] 3826
(b <- nrow(subset(pi1303, N_VARIANTS > 1)))
# [1] 2294
b/a
# [1] 0.5995818

pi1304<- read.delim2("dDoccat.FinalSNP.1304Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1304))
# [1] 3890
(b <- nrow(subset(pi1304, N_VARIANTS > 1)))
# [1] 2359
b/a
# [1] 0.6064267

pi1601 <- read.delim2("dDoccat.FinalSNP.1601Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1601))
# [1] 3885
(b <- nrow(subset(pi1601, N_VARIANTS > 1)))
# [1] 2401
b/a
# [1] 0.618018

pi1602 <- read.delim2("dDoccat.FinalSNP.1602Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1602))
# [1] 3876
(b <- nrow(subset(pi1602, N_VARIANTS > 1)))
# [1] 2391
b/a
# [1] 0.6168731

pi1603 <- read.delim2("dDoccat.FinalSNP.1603Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1603))
# [1] 3760
(b <- nrow(subset(pi1603, N_VARIANTS > 1)))
# [1] 2206
b/a
# [1] 0.5867021

pi1604 <- read.delim2("dDoccat.FinalSNP.1604Pi.windowed.pi", header = TRUE)
(a <- nrow(pi1604))
# [1] 3911
(b <- nrow(subset(pi1604, N_VARIANTS > 1)))
# [1] 2402
b/a
# [1] 0.6141652

#P for all pops
(0.6141652+0.5867021+0.6168731+0.618018+0.6064267+0.5995818+0.6+0.613086+0.5973277+0.597914+0.6084369+0.6049319+0.605079+0.5927282+0.6110969+0.5975289+0.6140442+0.6017745+0.609737+ 0.5928459+0.6187564+0.6102866+0.6154639+0.6229549)/24
# 0.60649