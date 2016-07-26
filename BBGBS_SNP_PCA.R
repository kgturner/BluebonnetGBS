#PCA analysis of SNPs
#6/29/2016

source("https://bioconductor.org/biocLite.R")
biocLite("SNPRelate")
library(SNPRelate)
library(gdsfmt)
library(ggplot2)

#########dDocent snps, filtered, lane 7 only#########
FAM<-read.table(file="dDocL7.FinalSNPs.fam",sep=" ", header=FALSE,na="NA")
head(FAM)
dim(FAM)
unique(FAM$V1)

bim<-read.table(file="dDocL7.FinalSNPs.bim",sep="\t", header=FALSE,na="NA")
head(bim)
dim(bim) 
#only 10555 snps/rows in this file, should be 10719. WTF is plink doing?
#try straight from vcf?
# vcf.fn<-"~/xxx/tmp.vcf"
# snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")

POPINFO=read.table(file="bbpopmap.txt",header=F)
names(POPINFO) <- c("IndivID", "Population", "PopType")
POPINFO$Population <- as.factor(POPINFO$Population)
table(POPINFO$Population) #still includes individuals that were dropped because of poor coverage
POPINFO <- subset(POPINFO, IndivID%in%FAM$V2)
POPINFO$IndivID <- factor(POPINFO$IndivID) #droplevels() didn't work for some reason...

sum(POPINFO$IndivID!=FAM$V2)

bed.fn <- "C:/Users/Kat/Documents/GitHub/BluebonnetGBS/dDocL7.FinalSNPs.bed"
fam.fn <- "C:/Users/Kat/Documents/GitHub/BluebonnetGBS/dDocL7.FinalSNPs.fam"
bim.fn <- "C:/Users/Kat/Documents/GitHub/BluebonnetGBS/dDocL7.FinalSNPs.bim"

# convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "FinalSNPs.gds")
#summary
snpgdsSummary("FinalSNPs.gds")

## Open the GDS file
genofile <- snpgdsOpen("FinalSNPs.gds")
head(genofile)

head(read.gdsn(index.gdsn(genofile, "sample.id")))
head(read.gdsn(index.gdsn(genofile, "snp.id")))

## Now perform a PCA using a function from the SNPRelate package
pca <- snpgdsPCA(genofile, autosome.only=FALSE)


tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

##  by population
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
population=POPINFO$Population

tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(population)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1", main="PCA using all SNPs")
legend("bottomright", legend=levels(tab$pop), pch="o", col=1:(nlevels(tab$pop)))
#not enough colors for all of the populations

#seeded vs wild
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
PopType=POPINFO$PopType
population=POPINFO$Population

tab <- data.frame(sample.id = pca$sample.id,
                  population = factor(population)[match(pca$sample.id, sample.id)],
                  poptype = factor(PopType)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the third eigenvector
                  EV4 = pca$eigenvect[,4],    # the forth eigenvector
                  stringsAsFactors = FALSE)
head(tab)

plot(tab$EV2, tab$EV1, col=as.integer(tab$poptype), xlab="eigenvector 2", ylab="eigenvector 1", main="PCA using all SNPs")
legend("bottomright", legend=levels(tab$poptype), pch="o", col=1:(nlevels(tab$poptype)))

library(ggplot2)
png("FinalSNPsPCA_1v2.png")
ggplot(tab, aes(x=EV1, y=EV2, color=poptype))+
  geom_point(shape=16)+
  xlab("PC1")+ylab("PC2")+
  theme(legend.title=element_blank())
dev.off()

png("FinalSNPsPCA_3v4.png")
ggplot(tab, aes(x=EV3, y=EV4, color=poptype))+
  geom_point(shape=16)+
  xlab("PC3")+ylab("PC4")+
  theme(legend.title=element_blank())
dev.off()

## Now make scatterplots of the top 4 PCs with proportional variance explained included
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
png("FinalSNPsPCA_1thru4.png")
pairs(pca$eigenvect[,1:4], col=tab$poptype, labels=lbls)
dev.off()


## Do we need 10K SNPs for population structure infererence in this sample?
##Identify a subset of SNPs based on LD threshold of 0.2 
# snpset <- snpgdsLDpruning(genofile, ld.threshold=0.01, autosome.only=F)
# snpset.id <- unlist(snpset)
#not sure I can assess ld here...anyway ld.thresholds of 0.2, 0.5, 0.01 all return "0 SNPs are selected in total."

##Estimate proportional Ancestry from the PCA.
#hmmm, doesn't really work, since I only have two groups, not two groups + derived group
avgseed=mean(pca2$eigenvect[PopType=="seed",2])
avgwild=mean(pca2$eigenvect[PopType=="wild",2])

admix=(pca$eigenvect[,2])/(avgseed-avgwild)

tab3=cbind(admix,1-admix)
myorder=order(admix)
temp=t(as.matrix(tab3[myorder,]))
# temp=t(as.matrix(tab3))

png("ancestryaRt2.png")
barplot(temp, col=c("blue","green"),xlab="Individual ", ylab="Ancestry", border=NA,axisnames=FALSE,main="Ancestry of wild",ylim=c(0,1))
legend("bottomright", c("seed","wild"), lwd=4, col=c("blue","green"), bg="white",cex=0.85)
dev.off()

#example from SISG AssMap course notes
# ##   Estimate proportional Native American and European Ancestry for the MXL from the PCA.  ## ASSUME THAT MXL have negligible African Ancestry.
# avgCEU2=mean(pca2$eigenvect[population=="CEU",2])
# avgNAM2=mean(pca2$eigenvect[population=="NAM",2])
# MXLadmix=(pca2$eigenvect[population=="MXL",2]-avgNAM2)/(avgCEU2-avgNAM2)
# ### NOW MAKE A BARPLOT OF MXL  ESTIMATED ANCESTRY FROM THE PCA ###
# tab2=cbind(MXLadmix,1-MXLadmix)
# myorder=order(MXLadmix)
# temp=t(as.matrix(tab2[myorder,]))
# barplot(temp, col=c("blue","green"),xlab="Individual ", ylab="Ancestry", border=NA,axisnames=FALSE,main="Ancestry of MXL",ylim=c(0,1))
# legend("bottomright", c("European","Native American"), lwd=4, col=c("blue","green"), bg="white",cex=0.85)


####prelim dataset, pop3####
# PLINK BED files
pop3_filtered

bed.fn <- "C:/Users/Kat/Documents/GitHub/BluebonnetGBS/pop3_filtered.bed"
fam.fn <- "C:/Users/Kat/Documents/GitHub/BluebonnetGBS/pop3_filtered.fam"
bim.fn <- "C:/Users/Kat/Documents/GitHub/BluebonnetGBS/pop3_filtered.bim"
# convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "pop3filt.gds")
#summary
snpgdsSummary("pop3filt.gds")

# Principal Component Analysis
#
# open
genofile <- snpgdsOpen("pop3filt.gds")
RV <- snpgdsPCA(genofile, autosome.only=FALSE)
pc.percent <- RV$varprop*100
head(round(pc.percent, 2))


plot(RV$eigenvect[,2], RV$eigenvect[,1], xlab="PC 2", ylab="PC 1",
     col=rgb(0,0,150, 50, maxColorValue=255), pch=19)

pop3PC1 <- unlist(as.vector(RV$eigenvect[,1]))
pop3PC2 <- unlist(as.vector(RV$eigenvect[,2]))
sampID <- unlist(as.vector(RV$sample.id))
popType <- read.table("pop3order_origin.txt", header=FALSE, stringsAsFactor=FALSE)

pop3PC1_2 <- as.data.frame(cbind(sampID, popType,pop3PC1, pop3PC2,pop3PC3=unlist(as.vector(RV$eigenvect[,3])),pop3PC4=unlist(as.vector(RV$eigenvect[,4]))))
# pop3PC1_2 <- as.data.frame(cbind(sampID, popType,pop3PC1, pop3PC2))
ggplot(pop3PC1_2,aes(x=pop3PC1, y=pop3PC2, color=V2)) +geom_point(size=3)
ggplot(pop3PC1_2,aes(x=pop3PC1, y=pop3PC3, color=V2)) +geom_point(size=3)
ggplot(pop3PC1_2,aes(x=pop3PC1, y=pop3PC4, color=V2)) +geom_point(size=3)

ggplot(pop3PC1_2,aes(x=pop3PC1, y=pop3PC2, color=V1, shape=V2)) +geom_point(size=3,position=position_jitter(width=0.05,height=0.05))







# close the file
snpgdsClose(genofile)

#Fst
# open dataset
genofile <- snpgdsOpen("pop3filt.gds")
pop_code <- scan("pop3order_origin.txt", what=character())
table(pop_code)

popType <- read.table("pop3order_origin.txt", header=FALSE, stringsAsFactor=FALSE)

group <- as.factor(popType$V2)

# Fst estimation
snpgdsFst(genofile, autosome.only=FALSE, population=group, method="W&H02")
# Fst estimation on SNP genotypes:
#   Excluding 0 SNP (monomorphic: TRUE, < MAF: NaN, or > missing rate: NaN)
# Working space: 96 samples, 945 SNPs
# # of Populations: 2
# seed (48), wild (48)
# $Fst
# [1] 0.005359065
# 
# $Beta
# seed       wild
# seed -0.01054652 0.00000000
# wild  0.00000000 0.02126465

# or
snpgdsFst(genofile, autosome.only=FALSE,population=group, method="W&C84")
# Fst estimation on SNP genotypes:
#   Excluding 0 SNP (monomorphic: TRUE, < MAF: NaN, or > missing rate: NaN)
# Working space: 96 samples, 945 SNPs
# # of Populations: 2
# seed (48), wild (48)
# $Fst
# [1] 0.005304012


group <- as.factor(popType$V1)
snpgdsFst(genofile, autosome.only=FALSE, population=group, method="W&H02")
# Fst estimation on SNP genotypes:
#   Excluding 0 SNP (monomorphic: TRUE, < MAF: NaN, or > missing rate: NaN)
# Working space: 96 samples, 945 SNPs
# # of Populations: 12
# 1201wild (8), 1202seed (8), 1203wild (8), 1204seed (8), 1205wild (8), 1301seed (8), 1303wild (8), 1304seed (8), 1601wild (8), 1602seed (8), 1603wild (8), 1604seed (8)
# $Fst
# [1] 0.02961928
# 
# $Beta
# 1201wild     1202seed     1203wild      1204seed
# 1201wild  0.015680153 -0.010969768  0.005344490 -0.0112403963
# 1202seed -0.010969768  0.016084265  0.007696648 -0.0092018170
# 1203wild  0.005344490  0.007696648  0.029300396  0.0048178034
# 1204seed -0.011240396 -0.009201817  0.004817803  0.0216249859
# 1205wild -0.003330647  0.004240687  0.014152662  0.0004612191
# 1301seed -0.009094809 -0.004066844 -0.001290702 -0.0147017577
# 1303wild  0.005705748 -0.001806785  0.021315028  0.0076255213
# 1304seed -0.021800493 -0.018655850 -0.003336984 -0.0198390061
# 1601wild  0.007639725  0.009243477  0.025990180  0.0074650259
# 1602seed -0.011797972 -0.008526334  0.003943469 -0.0077961607
# 1603wild  0.007114369  0.010841822  0.022986747  0.0116770260
# 1604seed -0.026907319 -0.022910108 -0.017391522 -0.0167314089
# 1205wild     1301seed     1303wild     1304seed    1601wild
# 1201wild -0.0033306470 -0.009094809  0.005705748 -0.021800493 0.007639725
# 1202seed  0.0042406871 -0.004066844 -0.001806785 -0.018655850 0.009243477
# 1203wild  0.0141526625 -0.001290702  0.021315028 -0.003336984 0.025990180
# 1204seed  0.0004612191 -0.014701758  0.007625521 -0.019839006 0.007465026
# 1205wild  0.0396250949  0.008245670  0.023767066 -0.004323396 0.022139629
# 1301seed  0.0082456697  0.031708902 -0.002061151 -0.014427188 0.013963727
# 1303wild  0.0237670661 -0.002061151  0.054420054 -0.001777528 0.024742040
# 1304seed -0.0043233956 -0.014427188 -0.001777528 -0.007109178 0.007668669
# 1601wild  0.0221396289  0.013963727  0.024742040  0.007668669 0.059748417
# 1602seed  0.0046227722 -0.015151175  0.006977188 -0.019443543 0.011216057
# 1603wild  0.0265980552  0.019001545  0.017699787  0.006941429 0.036603576
# 1604seed -0.0152989450 -0.025403705 -0.020250325 -0.028468778 0.001052736
# 1602seed     1603wild     1604seed
# 1201wild -0.011797972  0.007114369 -0.026907319
# 1202seed -0.008526334  0.010841822 -0.022910108
# 1203wild  0.003943469  0.022986747 -0.017391522
# 1204seed -0.007796161  0.011677026 -0.016731409
# 1205wild  0.004622772  0.026598055 -0.015298945
# 1301seed -0.015151175  0.019001545 -0.025403705
# 1303wild  0.006977188  0.017699787 -0.020250325
# 1304seed -0.019443543  0.006941429 -0.028468778
# 1601wild  0.011216057  0.036603576  0.001052736
# 1602seed  0.008836328  0.002453111 -0.019225252
# 1603wild  0.002453111  0.075379547 -0.004727036
# 1604seed -0.019225252 -0.004727036  0.010132451

# #IBD using maximum likelihood estimation????
# sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
# s1202.id <- sample.id[popType$V1 == "1202seed"]
# ibd <- snpgdsIBDMLE(genofile, autosome.only=FALSE, sample.id=s1202.id,
#                     maf=0.05, missing.rate=0.05, num.thread=2)
# ibd.coeff <- snpgdsIBDSelection(ibd)
# 
# plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
#      xlab="k0", ylab="k1", main="1202seed samples (MLE)")
# lines(c(0,1), c(1,0), col="red", lty=2)

#IBS analysis
ibs <- snpgdsIBS(genofile, autosome.only=FALSE,num.thread=2)
pop.idx <- order(popType$V1)

image(ibs$ibs[pop.idx, pop.idx], col=terrain.colors(16))

#multi-dimensional scaling analysis (IBS)
loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[, 1]; y <- loc[, 2]
race <- as.factor(popType$V2)

plot(x, y, col=race, xlab = "", ylab = "",
     main = "Multidimensional Scaling Analysis (IBS)")
legend("bottomright", legend=levels(race), text.col=1:nlevels(race))

#cluster analysis(IBS)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, autosome.only=FALSE,num.thread=2))
rv <- snpgdsCutTree(ibs.hc)
plot(rv$dendrogram, leaflab="none", main="HapMap Phase II")

rv2 <- snpgdsCutTree(ibs.hc, samp.group=as.factor(popType$V2))
plot(rv2$dendrogram, leaflab="none", main="HapMap Phase II")
legend("bottomright", legend=levels(race), col=1:nlevels(race), pch=19, ncol=4)


# # close the genotype file
# snpgdsClose(genofile)


####example####
# Convert the PLINK BED file to the GDS file
#
# PLINK BED files
#or your own
bed.fn <- "C:/your_folder/your_plink_file.bed"
fam.fn <- "C:/your_folder/your_plink_file.fam"
bim.fn <- "C:/your_folder/your_plink_file.bim"
# convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "HapMap.gds")
####or convert vcf file####
# The VCF file
vcf.fn <- system.file("extdata", "sequence.vcf", package="SNPRelate")
cat(readLines(vcf.fn), sep="\n")
snpgdsVCF2GDS_R(vcf.fn, "test1.gds", method="biallelic.only")
snpgdsSummary("test1.gds")
snpgdsVCF2GDS_R(vcf.fn, "test2.gds", method="biallelic.only")
snpgdsSummary("test2.gds")
snpgdsVCF2GDS_R(vcf.fn, "test3.gds", method="copy.num.of.ref")
snpgdsSummary("test3.gds")
snpgdsVCF2GDS_R(vcf.fn, "test4.gds", method="copy.num.of.ref")
snpgdsSummary("test4.gds")

####################################################################
# Principal Component Analysis
#
# open
genofile <- snpgdsOpen("HapMap.gds")
RV <- snpgdsPCA(genofile)
plot(RV$eigenvect[,2], RV$eigenvect[,1], xlab="PC 2", ylab="PC 1",
     col=rgb(0,0,150, 50, maxColorValue=255), pch=19)
# close the file
snpgdsClose(genofile)
######################################################################
#Fst
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())
group <- as.factor(read.gdsn(index.gdsn(
  genofile, "sample.annot/pop.group")))
# Fst estimation
snpgdsFst(genofile, population=group, method="W&H02")
# or
snpgdsFst(genofile, population=group, method="W&C84")
# close the genotype file
snpgdsClose(genofile)