#Standard errors for pop gen stats
#7/17/17

#se <- sd(x)/sqrt(length(x))

wvsFst <- read.delim("dDoccat.FinalSNP.wild_vs_seed.weir.fst")
(se <- sd(wvsFst$WEIR_AND_COCKERHAM_FST)/sqrt(length(wvsFst$WEIR_AND_COCKERHAM_FST)))
# [1] 8.582894e-05
mean(wvsFst$WEIR_AND_COCKERHAM_FST)
# [1] 0.0006407671 This is the same as unweighted estimate reported by vcftools

#confidence interval?
t.test(wvsFst$WEIR_AND_COCKERHAM_FST)

# One Sample t-test
# 
# data:  wvsFst$WEIR_AND_COCKERHAM_FST
# t = 7.4656, df = 10955, p-value = 8.919e-14
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   0.0004725269 0.0008090074
# sample estimates:
#   mean of x 
# 0.0006407671 This is the same as unweighted estimate reported by vcftools

asFst <- read.delim("dDoccat.FinalSNP.amongSeedFst.weir.fst")
(se <- sd(asFst$WEIR_AND_COCKERHAM_FST)/sqrt(length(asFst$WEIR_AND_COCKERHAM_FST)))
#[1] 0.0003516225
t.test(asFst$WEIR_AND_COCKERHAM_FST)
# One Sample t-test
# 
# data:  asFst$WEIR_AND_COCKERHAM_FST
# t = 37.069, df = 11244, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   0.01234508 0.01372356
# sample estimates:
#   mean of x 
# 0.01303432 

awFst <- read.delim("dDoccat.FinalSNP.amongWildFst.weir.fst")
(se <- sd(awFst$WEIR_AND_COCKERHAM_FST)/sqrt(length(awFst$WEIR_AND_COCKERHAM_FST)))
#[1] 0.0003945777
t.test(awFst$WEIR_AND_COCKERHAM_FST)
# One Sample t-test
# 
# data:  awFst$WEIR_AND_COCKERHAM_FST
# t = 45.323, df = 11062, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   0.01711006 0.01865695
# sample estimates:
#   mean of x 
# 0.01788351 

allFst <- read.delim("dDoccat.FinalSNP.allFst.weir.fst")
(se <- sd(allFst$WEIR_AND_COCKERHAM_FST)/sqrt(length(allFst$WEIR_AND_COCKERHAM_FST)))
#[1] 0.0002753403
t.test(allFst$WEIR_AND_COCKERHAM_FST)
# One Sample t-test
# 
# data:  allFst$WEIR_AND_COCKERHAM_FST
# t = 57.13, df = 10955, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   0.01519039 0.01626982
# sample estimates:
#   mean of x 
# 0.01573011 

#pairs
v1002Fst <- read.delim("dDoccat.FinalSNP.1002v1001Fst.weir.fst")
(se <- sd(v1002Fst$WEIR_AND_COCKERHAM_FST, na.rm=T)/sqrt(length(v1002Fst$WEIR_AND_COCKERHAM_FST)))
#[1]  0.000711596
t.test(v1002Fst$WEIR_AND_COCKERHAM_FST)
# One Sample t-test
# 
# data:  v1002Fst$WEIR_AND_COCKERHAM_FST
# t = 0.15615, df = 11065, p-value = 0.8759
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   -0.001315707  0.001543473
# sample estimates:
#   mean of x 
# 0.0001138832 

v1005Fst <- read.delim("dDoccat.FinalSNP.1005v1006Fst.weir.fst")
(se <- sd(v1005Fst$WEIR_AND_COCKERHAM_FST, na.rm=T)/sqrt(length(v1005Fst$WEIR_AND_COCKERHAM_FST)))
#[1] 0.000789743
t.test(v1005Fst$WEIR_AND_COCKERHAM_FST)
# One Sample t-test
# 
# data:  v1005Fst$WEIR_AND_COCKERHAM_FST
# t = 10.247, df = 11003, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   0.006725613 0.009907442
# sample estimates:
#   mean of x 
# 0.008316528 

v1007Fst <- read.delim("dDoccat.FinalSNP.1007v1008Fst.weir.fst")
(se <- sd(v1007Fst$WEIR_AND_COCKERHAM_FST, na.rm=T)/sqrt(length(v1007Fst$WEIR_AND_COCKERHAM_FST)))
#[1]  0.0008644499
t.test(v1007Fst$WEIR_AND_COCKERHAM_FST)
# One Sample t-test
# 
# data:  v1007Fst$WEIR_AND_COCKERHAM_FST
# t = 21.349, df = 10707, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   0.01742420 0.02094739
# sample estimates:
#   mean of x 
# 0.01918579 

v1101Fst <- read.delim("dDoccat.FinalSNP.1101v1102Fst.weir.fst")
(se <- sd(v1101Fst$WEIR_AND_COCKERHAM_FST, na.rm=T)/sqrt(length(v1101Fst$WEIR_AND_COCKERHAM_FST)))
#[1]  0.0008422158
t.test(v1101Fst$WEIR_AND_COCKERHAM_FST)
# One Sample t-test
# 
# data:  v1101Fst$WEIR_AND_COCKERHAM_FST
# t = 17.298, df = 10746, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   0.01337378 0.01679212
# sample estimates:
#   mean of x 
# 0.01508295 

v1103Fst <- read.delim("dDoccat.FinalSNP.1103v1104Fst.weir.fst")
(se <- sd(v1103Fst$WEIR_AND_COCKERHAM_FST, na.rm=T)/sqrt(length(v1103Fst$WEIR_AND_COCKERHAM_FST)))
#[1]  0.0008285444
t.test(v1103Fst$WEIR_AND_COCKERHAM_FST)
# One Sample t-test
# 
# data:  v1103Fst$WEIR_AND_COCKERHAM_FST
# t = 11.352, df = 10743, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   0.008051024 0.011411854
# sample estimates:
#   mean of x 
# 0.009731439 

v1105Fst <- read.delim("dDoccat.FinalSNP.1105v1107Fst.weir.fst")
(se <- sd(v1105Fst$WEIR_AND_COCKERHAM_FST, na.rm=T)/sqrt(length(v1105Fst$WEIR_AND_COCKERHAM_FST)))
#[1]  0.0008191593
t.test(v1105Fst$WEIR_AND_COCKERHAM_FST)
# One Sample t-test
# 
# data:  v1105Fst$WEIR_AND_COCKERHAM_FST
# t = 4.9691, df = 10614, p-value = 6.832e-07
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   0.002561663 0.005899322
# sample estimates:
#   mean of x 
# 0.004230492

v1201Fst <- read.delim("dDoccat.FinalSNP.1201v1202Fst.weir.fst")
(se <- sd(v1201Fst$WEIR_AND_COCKERHAM_FST, na.rm=T)/sqrt(length(v1201Fst$WEIR_AND_COCKERHAM_FST)))
#[1]  0.0007163044
t.test(v1201Fst$WEIR_AND_COCKERHAM_FST)
# One Sample t-test
# 
# data:  v1201Fst$WEIR_AND_COCKERHAM_FST
# t = -1.0038, df = 10816, p-value = 0.3155
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   -0.0021984202  0.0007093658
# sample estimates:
#   mean of x 
# -0.0007445272 

v1203Fst <- read.delim("dDoccat.FinalSNP.1203v1204Fst.weir.fst")
(se <- sd(v1203Fst$WEIR_AND_COCKERHAM_FST, na.rm=T)/sqrt(length(v1203Fst$WEIR_AND_COCKERHAM_FST)))
#[1]  0.000852277
t.test(v1203Fst$WEIR_AND_COCKERHAM_FST)
# One Sample t-test
# 
# data:  v1203Fst$WEIR_AND_COCKERHAM_FST
# t = 0.54599, df = 10590, p-value = 0.5851
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   -0.001253758  0.002221862
# sample estimates:
#   mean of x 
# 0.000484052 

v1205Fst <- read.delim("dDoccat.FinalSNP.1205v1301Fst.weir.fst")
(se <- sd(v1205Fst$WEIR_AND_COCKERHAM_FST, na.rm=T)/sqrt(length(v1205Fst$WEIR_AND_COCKERHAM_FST)))
#[1]  0.0007656696
t.test(v1205Fst$WEIR_AND_COCKERHAM_FST)
# One Sample t-test
# 
# data:  v1205Fst$WEIR_AND_COCKERHAM_FST
# t = 3.8186, df = 10782, p-value = 0.000135
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   0.001470536 0.004572723
# sample estimates:
#   mean of x 
# 0.003021629 

v1303Fst <- read.delim("dDoccat.FinalSNP.1303v1304Fst.weir.fst")
(se <- sd(v1303Fst$WEIR_AND_COCKERHAM_FST, na.rm=T)/sqrt(length(v1303Fst$WEIR_AND_COCKERHAM_FST)))
#[1]  0.0008501409
t.test(v1303Fst$WEIR_AND_COCKERHAM_FST)
# One Sample t-test
# 
# data:  v1303Fst$WEIR_AND_COCKERHAM_FST
# t = 14.211, df = 10837, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   0.01077698 0.01422561
# sample estimates:
#   mean of x 
# 0.01250129 

v1601Fst <- read.delim("dDoccat.FinalSNP.1601v1602Fst.weir.fst")
(se <- sd(v1601Fst$WEIR_AND_COCKERHAM_FST, na.rm=T)/sqrt(length(v1601Fst$WEIR_AND_COCKERHAM_FST)))
#[1]  0.0007434775
t.test(v1601Fst$WEIR_AND_COCKERHAM_FST)
# One Sample t-test
# 
# data:  v1601Fst$WEIR_AND_COCKERHAM_FST
# t = -0.61904, df = 10878, p-value = 0.5359
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   -0.001980375  0.001029752
# sample estimates:
#   mean of x 
# -0.0004753117

v1603Fst <- read.delim("dDoccat.FinalSNP.1603v1604Fst.weir.fst")
(se <- sd(v1603Fst$WEIR_AND_COCKERHAM_FST, na.rm=T)/sqrt(length(v1603Fst$WEIR_AND_COCKERHAM_FST)))
#[1]  [1] 0.0008573343
t.test(v1603Fst$WEIR_AND_COCKERHAM_FST)
# One Sample t-test
# 
# data:  v1603Fst$WEIR_AND_COCKERHAM_FST
# t = 14.143, df = 10786, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   0.01085099 0.01434272
# sample estimates:
#   mean of x 
# 0.01259685 
