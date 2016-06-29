#PCA analysis of SNPs
#6/29/2016

source("https://bioconductor.org/biocLite.R")
biocLite("SNPRelate")
library(SNPRelate)
library(gdsfmt)
library(ggplot2)


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

# #Fst
# # open an example dataset (HapMap)
# genofile <- snpgdsOpen(snpgdsExampleFileName())
# group <- as.factor(read.gdsn(index.gdsn(
#   genofile, "sample.annot/pop.group")))
# # Fst estimation
# snpgdsFst(genofile, population=group, method="W&H02")
# # or
# snpgdsFst(genofile, population=group, method="W&C84")
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