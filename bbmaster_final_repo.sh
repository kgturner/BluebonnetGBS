#bbmaster_final, for github repository
#11/3/2016 - 10/17/2017

#for installation and prelim work, see bbmaster.sh

#gunzip raw files
#cat lane 6 and lane 7 R1 files together
cat HI.3413.006.BB_R1.fastq HI.3413.007.BB_R1.fastq >> BB_R1cat.fastq
#to save space, rezip/remove raw files
gzip HI.3413.006.BB_R1.fastq
gzip HI.3413.007.BB_R1.fastq
#demultiplex catted R1 file
# barcode table: bb.bars.names.txt
#in dir where you want files
mkdir demultiplexed_cat
cd demultiplexed_cat
perl ~/Bluebonnet/GBS_fastq_Demultiplexer_v9_2Enzyme2barcode.pl ~/Bluebonnet/bb.bars.names.txt ~/Bluebonnet/raw/BB_R1cat.fastq ~/Bluebonnet/raw/BB_R2cat.fastq BBcat

#gzip resulting file
cd demultiplexed_cat
gzip *

#check naming structure of previous run, 
#currently:
 BBcat_Pop1001_001_R1.fastq.gz
#should look like:
 1001_001.F.fq.gz
#demultiplex script makes 2 extra files (no bar)? Check and remove before renaming
ls -thor|less
mv BBcat_nobar* ~/Bluebonnet/temp_dDocent

#rename files for dDocent, wants Pop1_001.F.fq.gz, Pop1_001.R.fq.gz,Pop1_002.F.fq.gz,  Pop1_002.R.fq.gz
# -v for verbose
rename .fastq. .fq. *.gz
rename _R1. .F. *.fq.gz
rename _R2. .R. *.fq.gz
# remove a number of characters off the front, here "BBcat_Pop", because dDocent expects pop names to be 4 characters
for file in *gz; do mv "$file" "${file:9}"; done

#####running dDocent######
#need to update
#rm old dDocent dir
mkdir dDocent
#make sure dDocent in path

#copy install script into dDocent v2.2.6
chmod +x install_dDocent_requirements.sh
#run instal
sh install_dDocent_requirements.sh ~/dDocent

#update samtools to 1.3+ and add to path
mkdir samtool-1.3.1
#add samtools zipped file
#unzip
tar -jxvf samtools-1.3.1.tar.bz2

#build and install 
cd samtools-1.3.1   
make
make prefix=~/dDocent/samtools-1.3.1 install 

#FreeBayes, need version 1+ apparently, delete version installed by dDocent installer
rm -rf freebayes/
#get v1.1.0
git clone --recursive git://github.com/ekg/freebayes.git
cd freebayes/
make

#installed new version of cd-hit to fix seg fault
tar xvf cd-hit-v4.6.6-2016-0711.tar.gz -gunzip
cd  cd-hit-v4.6.6-2016-0711
make
cd cd-hit-auxtools
make
# in screen, add to path
#add stuff to path
export PATH=~/dDocent/samtools-1.3.1/bin:$PATH   # for sh or bash users
export PATH=~/dDocent/freebayes/bin:$PATH    # for sh or bash users
export PATH=~/Bluebonnet/bin/vcftools_0.1.11/bin:$PATH
export PATH=~/Bluebonnet/bin/cd-hit-v4.6.6-2016-0711:$PATH

#### check that it works by:
sh dDocent

#in screen, run dDocent in dir with demultiplexed reads
dDocent
#and answer prompts
#for 11/30/2016 run: 8 processors, 10G of memory (out of 23G free), trim, assemble, PE, don't change %similiary to cluster, map reads but dont' change parameters, call SNPs with freebayes, err on side of default values
#data cut-off 6x coverage with in indiv, data cut-off 19 individuals 

########using dDocent SNPs###########
#intense filtering tutorial here! https://github.com/jpuritz/dDocent/blob/master/tutorials/Filtering%20Tutorial.md
#filter in vcftools 0.1.15 (remove older vcftools from PATH or alternately open a new screen), vcflib, mawk
mkdir dDoc_cat_filtering
cp TotalRawSNPs.vcf ~/Bluebonnet/dDoc_cat_filtering/
cd ~/Bluebonnet/dDoc_cat_filtering/

#following tutorial
#keep variants that have been successfully genotyped in 50% of individuals, a minimum quality score of 30, and a minor allele count of 3.
vcftools --vcf TotalRawSNPs.vcf --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out dDoccat.g5mac3
# After filtering, kept 192 out of 192 Individuals
# Outputting VCF file...
# After filtering, kept 81662 out of a possible 180611 Sites

#minimum depth for a genotype call and a minimum mean depth; keep genotypes with 3+ reads
vcftools --vcf dDoccat.g5mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out dDoccat.g5mac3dp3
# After filtering, kept 192 out of 192 Individuals
# Outputting VCF file...
# After filtering, kept 81662 out of a possible 81662 Sites

# get rid of individuals that did not sequence well
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/filter_missing_ind.sh
chmod +x filter_missing_ind.sh
./filter_missing_ind.sh dDoccat.g5mac3dp3.recode.vcf DP3g95maf05
#85% cutoff/10% missing data leaves 166 individuals
#30% missing data leaves 187 individuals -> going with this for now...
# Excluding individuals in 'exclude' list
# After filtering, kept 187 out of 192 Individuals
# Outputting VCF file...
# After filtering, kept 81662 out of a possible 81662 Sites

#restrict the data to variants called in a high percentage of individuals and filter by mean depth of genotypes. This applied a genotype call rate (95%) across all individuals
vcftools --vcf DP3g95maf05.recode.vcf --max-missing 0.95 --maf 0.05 --recode --recode-INFO-all --out dDoccat.DP3g95maf05 --min-meanDP 20
# After filtering, kept 187 out of 187 Individuals
# Outputting VCF file...
# After filtering, kept 15041 out of a possible 81662 Sites
# vcftools --vcf dDoccat.g5mac3dp3.recode.vcf --max-missing 0.95 --maf 0.05 --recode --recode-INFO-all --out temp --min-meanDP 20

#This script will automatically filter a FreeBayes generated VCF file using criteria related to site depth,
# quality versus depth, strand representation, allelic balance at heterzygous individuals, and paired read representation.
# The script assumes that loci and individuals with low call rates (or depth) have already been removed.

curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/dDocent_filters
chmod +x dDocent_filters
./dDocent_filters
sh dDocent_filters.sh VCF_file Output_prefix
./dDocent_filters dDoccat.DP3g95maf05m.recode.vcf dDoccat.DP3g95maf05m.fils

# Number of sites filtered based on allele balance at heterozygous loci, locus quality, and mapping quality / Depth
 # 867 of 15041

# Number of additional sites filtered based on overlapping forward and reverse reads
 # 990 of 14174

# Is this from a mixture of SE and PE libraries? Enter yes or no.
# no
# Number of additional sites filtered based on properly paired status
 # 1021 of 13184

# Number of sites filtered based on high depth and lower than 2*DEPTH quality score
 # 422 of 12163


# VCFtools - 0.1.15
# (C) Adam Auton and Anthony Marcketta 2009

# Parameters as interpreted:
        # --vcf dDoccat.DP3g95maf05m.fils.fil.vcf
        # --exclude-positions dDoccat.DP3g95maf05m.recode.vcf.lowQDloci
        # --out dDoccat.DP3g95maf05m.recode.vcf
        # --site-depth
        # --remove-filtered NP

# After filtering, kept 187 out of 187 Individuals
# Outputting Depth for Each Site
# After filtering, kept 11929 out of a possible 12163 Sites
# Run Time = 1.00 seconds

                                             # Histogram of mean depth per site

  # 180 +-+-+---+---+---+---+---+---+--+---+---+---+---+---+---+---+---+---+---+---+---+---+--+---+---+---+---+---+-+-+
      # +   +   +   +   +  **   +   +  +   +   +   +   +   +   +   +   +   +   +   +   +   +  +   +   +   +   +   +   +
      # |              ** *****                             'meandepthpersite' using (bin($1,binwidth)):(1.0) ******* |
  # 160 +-+           *** *****                                                                                     +-+
      # |            ************                                                                                     |
  # 140 +-+          ************     *** ***            **                                                         +-+
      # |            **************   *** ***            **                                                           |
      # |       ************************* ***    **  **  **  *****                                                    |
  # 120 +-+     ***************************** ***** ***  ** ******                                                  +-+
      # |       ************************************************** **                                                 |
  # 100 +-+     ****************************************************** **                                           +-+
      # |       ****************************************************** **   *******                                   |
      # |       *********************************************************   *******       **                          |
   # 80 +-+     **********************************************************  *******   **  **                        +-+
      # |       *******************************************************************   **  **  **                      |
   # 60 +-+     ********************************************************************* *** **  **                    +-+
      # |       **************************************************************************** ****                     |
      # |       ********************************************************************************* **   **   **        |
   # 40 +-+     *************************************************************************************  **  ***      +-+
      # |       **********************************************************************************************        |
   # 20 +-+     ************************************************************************************************ *  +-+
      # |       *******************************************************************************************************
      # +   +   *******************************************************************************************************
    # 0 +-+-+---*******************************************************************************************************
      # 10  15  20  25  30  35  40  45 50  55  60  65  70  75  80  85  90  95 100 105 110 11 120 125 130 135 140 145 150
                                                        # Mean Depth

# If distrubtion looks normal, a 1.645 sigma cutoff (~90% of the data) would be 28230.2942
# The 95% cutoff would be 143
# Would you like to use a different maximum mean depth cutoff than 143, yes or no
# no

# VCFtools - 0.1.15
# (C) Adam Auton and Anthony Marcketta 2009

# Parameters as interpreted:
        # --vcf dDoccat.DP3g95maf05m.fils.fil.vcf
        # --recode-INFO-all
        # --max-meanDP 143
        # --out dDoccat.DP3g95maf05m.fils.FIL
        # --recode
        # --remove-filtered NP

# After filtering, kept 187 out of 187 Individuals
# Outputting VCF file...
# After filtering, kept 11340 out of a possible 12163 Sites
# Run Time = 16.00 seconds
# Number of sites filtered based on maximum mean depth
 # 823 of 12163


# VCFtools - 0.1.15
# (C) Adam Auton and Anthony Marcketta 2009

# Parameters as interpreted:
        # --vcf dDoccat.DP3g95maf05m.fils.fil.vcf
        # --exclude-positions dDoccat.DP3g95maf05m.recode.vcf.lowQDloci
        # --recode-INFO-all
        # --max-meanDP 143
        # --out dDoccat.DP3g95maf05m.fils.FIL
        # --recode
        # --remove-filtered NP

# After filtering, kept 187 out of 187 Individuals
# Outputting VCF file...
# After filtering, kept 11312 out of a possible 12163 Sites
# Run Time = 17.00 seconds
# Total number of sites filtered
 # 3729 of 15041

# Remaining sites
 # 11312

# Filtered VCF file is called Output_prefix.FIL.recode.vcf

# Filter stats stored in dDoccat.DP3g95maf05m.fils.filterstats

#filter departures from HWE within populations
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/filter_hwe_by_pop.pl
chmod +x filter_hwe_by_pop.pl
#filter our SNPs by population specific HWE
#decompose complex variant calls into phased SNP and INDEL genotypes and keep the INFO flags for loci and genotypes
vcfallelicprimitives dDoccat.DP3g95maf05m.fils.FIL.recode.vcf --keep-info --keep-geno > dDoccat.DP3g95maf05m.prim.vcf
#remove indels
vcftools --vcf dDoccat.DP3g95maf05m.prim.vcf --remove-indels --recode --recode-INFO-all --out SNP.dDoccat.DP3g95maf05m
mawk '!/#/' dDoccat.DP3g95maf05m.fils.FIL.recode.vcf | wc -l
#11312
mawk '!/#/' dDoccat.DP3g95maf05m.prim.vcf | wc -l
#12025
mawk '!/#/' SNP.dDoccat.DP3g95maf05m.recode.vcf | wc -l
#11741

#apply the HWE filter 
#default -h, Minimum cutoff for Hardy-Weinberg p-value (for test as implemented in vcftools) [Default: 0.001]
#default -c ,Proportion of all populations that a locus can be below HWE cutoff without being filtered. For example, choosing 0.5 willfilter SNPs that are below the p-value threshold in 50% or moreof the populations. [Default: 0.25]
./filter_hwe_by_pop.pl -v SNP.dDoccat.DP3g95maf05m.recode.vcf -p bbpopmap.txt -o SNP.dDoccat.DP3g95maf05m.HWE 
#at default h and c, and also at default h and c = 0.1, no loci filtered.
#at default c and h = 0.05, 615 loci removed. But this is probably too high an h

#estimate error
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/ErrorCount.sh
chmod +x ErrorCount.sh 
./ErrorCount.sh SNP.dDoccat.DP3g95maf05m.recode.vcf
# This script counts the number of potential genotyping errors due to low read depth
# It report a low range, based on a 50% binomial probability of observing the second allele in a heterozygote and a high range based on a 25% probability.
# Potential genotyping errors from genotypes from only 1 read range from 0 to 0.0
# Potential genotyping errors from genotypes from only 2 reads range from 0 to 0.0
# Potential genotyping errors from genotypes from only 3 reads range from 940 to 3161.34
# Potential genotyping errors from genotypes from only 4 reads range from 561 to 2840.208
# Potential genotyping errors from genotypes from only 5 reads range from 319 to 2420
# 187 number of individuals and 11741 equals 2195567 total genotypes
# Total genotypes not counting missing data 2185445
# Total potential error rate is between 0.000832782339524 and 0.00385347057464
# SCORCHED EARTH SCENARIO
# WHAT IF ALL LOW DEPTH HOMOZYGOTE GENOTYPES ARE ERRORS?????
# The total SCORCHED EARTH error rate is 0.0122309186459.


# ######Fst in vcftools
cd Bluebonnet
mkdir dDoc_Fst
cd dDoc_Fst/
cp ~/Bluebonnet/dDoc_cat_filtering/dDoccat.FinalSNP.vcf dDoccat.FinalSNP.vcf
 
# # #among all populations
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1201Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1203Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1205Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1303Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1601Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1603Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1202Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1204Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1301Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1304Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1602Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1604Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1001Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1002Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1005Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1006Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1007Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1008Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1101Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1102Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1103Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1104Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1105Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1107Indiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.allFst
# Weir and Cockerham mean Fst estimate: 0.01573
# Weir and Cockerham weighted Fst estimate: 0.014794

# # #wild vs. seed
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/wildIndiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/seedIndiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.wild_vs_seed
# Weir and Cockerham mean Fst estimate: 0.00064077
# Weir and Cockerham weighted Fst estimate: 0.00068971

# # #among wild
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1201Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1203Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1205Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1303Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1601Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1603Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1002Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1005Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1007Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1101Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1103Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1105Indiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.amongWildFst
# Weir and Cockerham mean Fst estimate: 0.017884
# Weir and Cockerham weighted Fst estimate: 0.017531

# # #among seed
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1202Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1204Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1301Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1304Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1602Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1604Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1001Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1006Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1008Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1102Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1104Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1107Indiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.amongSeedFst
# Weir and Cockerham mean Fst estimate: 0.013034
# Weir and Cockerham weighted Fst estimate: 0.012504

#between neighbors, wild vs seeded
# 1101 vs 1102
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1101Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1102Indiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.1101v1102Fst
# Weir and Cockerham mean Fst estimate: 0.015083
# Weir and Cockerham weighted Fst estimate: 0.023

# 1603 vs 1604
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1603Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1604Indiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.1603v1604Fst
# Weir and Cockerham mean Fst estimate: 0.012597
# Weir and Cockerham weighted Fst estimate: 0.019858

# 1601 vs 1602
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1601Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1602Indiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.1601v1602Fst
# Weir and Cockerham mean Fst estimate: -0.0004753
# Weir and Cockerham weighted Fst estimate: 0.0036834

# 1303 vs 1304
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1303Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1304Indiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.1303v1304Fst
# Weir and Cockerham mean Fst estimate: 0.012501
# Weir and Cockerham weighted Fst estimate: 0.020342

# 1205 vs 1301
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1205Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1301Indiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.1205v1301Fst
# Weir and Cockerham mean Fst estimate: 0.0030216
# Weir and Cockerham weighted Fst estimate: 0.0082165

# 1203 vs 1204
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1203Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1204Indiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.1203v1204Fst
# Weir and Cockerham mean Fst estimate: 0.00048405
# Weir and Cockerham weighted Fst estimate: 0.006448

# 1201 vs 1202
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1201Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1202Indiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.1201v1202Fst
# Weir and Cockerham mean Fst estimate: -0.00074452
# Weir and Cockerham weighted Fst estimate: 0.0012717

# 1105 vs 1107
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1105Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1107Indiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.1105v1107Fst
# Weir and Cockerham mean Fst estimate: 0.0042305
# Weir and Cockerham weighted Fst estimate: 0.010128

# 1103 vs 1104 (8 and 8)
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1103Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1104Indiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.1103v1104Fst
# Weir and Cockerham mean Fst estimate: 0.0097314
# Weir and Cockerham weighted Fst estimate: 0.016803

# 1007 vs 1008
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1007Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1008Indiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.1007v1008Fst
# Weir and Cockerham mean Fst estimate: 0.019186
# Weir and Cockerham weighted Fst estimate: 0.027048

# 1005 vs 1006
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1005Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1006Indiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.1005v1006Fst
# Weir and Cockerham mean Fst estimate: 0.0083165
# Weir and Cockerham weighted Fst estimate: 0.015101

# 1002 vs 1001
vcftools --vcf ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1002Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1001Indiv.txt --out ~/Bluebonnet/dDoc_Fst/dDoccat.FinalSNP.1002v1001Fst
# Weir and Cockerham mean Fst estimate: 0.00011389
# Weir and Cockerham weighted Fst estimate: 0.0037198


######use STRUCTURE with filtered SNPs
#use pgdspider on laptop to make structure input  file from filtered .vcf
#gave PGDSpider vcf file directly, columns worked out
#pgdspider asks if several columns are absent. Assuming there are 6 columns before the ATGCs start, the answer to all these questions should be 'no'. Group individuals acording to FamilyID (i.e. first column) Need to edit this column??? Specify data type in structure as SNP

#then edit structure file as tab/space delimted text file in excel(be sure to uncheck 'treat multiple delimters as one' option while importing to excel"
#add/edit columns as necessary. column order should be:
#1: individual id (2 rows for each individual)
#2: pop id (as integer 1-24)
#3: use pop flag (set all to 1)
#4: loc data for use in LOCPRIOR models, user set location or other info. Integer. Here wild = 1, seeded =2
#5 - n: SNP info. ATGC = 1 - 4, missing = -9

# #import into structure (on laptop) using "new project"
#187 indiv, 11741 loci, diploid, missing data value -9
# #build parameter sets, one set per group of models
# #Start job to run a batch of parameter sets, with several Ks
# #project: Bbcat


#####PCA#####
#in R
#

#################re-filtering following reviewer comments##########
# Lastly, I am a little concerned with the degree of apparent missing data... They also say that there was up to 30% missing data per individual. I think it 
# would be helpful if there were a statement to the effect of: Genotyping of SNPs ranged from X to Y % of individuals (mean = Z) and genotyping of individuals ranged from X to Y %
 # of loci (mean = Z). Given that there are so many SNPs, might having more stringent filtering affect the data?


########using dDocent SNPs###########
#intense filtering tutorial here! https://github.com/jpuritz/dDocent/blob/master/tutorials/Filtering%20Tutorial.md
#filter in vcftools 0.1.15 (remove older vcftools from PATH or alternately open a new screen), vcflib, mawk
mkdir dDoc_cat_Refiltering
cd  dDoc_cat_filtering
cp TotalRawSNPs.vcf ~/Bluebonnet/dDoc_cat_Refiltering/
cd ~/Bluebonnet/dDoc_cat_Refiltering/

# #following tutorial
# #keep variants that have been successfully genotyped in 50% of individuals, a minimum quality score of 30, and a minor allele count of 3.
# vcftools --vcf TotalRawSNPs.vcf --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out dDoccat.g5mac3
# # After filtering, kept 192 out of 192 Individuals
# # Outputting VCF file...
# # After filtering, kept 81662 out of a possible 180611 Sites

# #minimum depth for a genotype call and a minimum mean depth; keep genotypes with 3+ reads
# vcftools --vcf dDoccat.g5mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out dDoccat.g5mac3dp3
# # After filtering, kept 192 out of 192 Individuals
# # Outputting VCF file...
# # After filtering, kept 81662 out of a possible 81662 Sites

#####
#try this step again
cp dDoccat.g5mac3dp3.recode.vcf ~/Bluebonnet/dDoc_cat_Refiltering/
# get rid of individuals that did not sequence well
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/filter_missing_ind.sh
chmod +x filter_missing_ind.sh
./filter_missing_ind.sh dDoccat.g5mac3dp3.recode.vcf testMissing
#Previous: 85% cutoff/10% missing data leaves 166 individuals
#30% missing data leaves 187 individuals -> going with this for now...
# Excluding individuals in 'exclude' list
# After filtering, kept 187 out of 192 Individuals
# Outputting VCF file...
# After filtering, kept 81662 out of a possible 81662 Sites

# The 85% cutoff would be 0.102569
# Would you like to set a different cutoff, yes or no
# no
# After filtering, kept 166 out of 192 Individuals
# Outputting VCF file...
# After filtering, kept 81662 out of a possible 81662 Sites


./filter_missing_ind.sh dDoccat.g5mac3dp3.recode.vcf testMissing2
# The 85% cutoff would be 0.102569
# Would you like to set a different cutoff, yes or no
# yes
# Please enter new cutoff
# .3
# After filtering, kept 187 out of 192 Individuals
# Outputting VCF file...
# After filtering, kept 81662 out of a possible 81662 Sites

./filter_missing_ind.sh dDoccat.g5mac3dp3.recode.vcf testMissing3
# The 85% cutoff would be 0.102569
# Would you like to set a different cutoff, yes or no
# yes
# Please enter new cutoff
# .2
# After filtering, kept 183 out of 192 Individuals
# Outputting VCF file...
# After filtering, kept 81662 out of a possible 81662 Sites

###############move on with the default/more stringent cuttoff
#restrict the data to variants called in a high percentage of individuals and filter by mean depth of genotypes. This applied a genotype call rate (95%) across all individuals
vcftools --vcf testMissing.recode.vcf --max-missing 0.95 --maf 0.05 --recode --recode-INFO-all --out dDoccat.testMissing --min-meanDP 20
#Previous: After filtering, kept 15041 out of a possible 81662 Sites
#After filtering, kept 166 out of 166 Individuals
# Outputting VCF file...
# After filtering, kept 16157 out of a possible 81662 Sites

#when you have multiple localities being sampled You are also going to want to filter by a population specific call rate.
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/pop_missing_filter.sh
chmod +x pop_missing_filter.sh
./pop_missing_filter.sh
#Usage is pop_missing_filter vcffile popmap proportion_missing_per_pop number_of_pops_for_cutoff name_for_output
# ./pop_missing_filter.sh dDoccat.testMissing.recode.vcf bbpopmap.txt 0.3 1 dDoccat.testMissing.pop
# After filtering, kept 166 out of 166 Individuals
# Outputting VCF file...
# After filtering, kept 15875 out of a possible 16157 Sites
./pop_missing_filter.sh dDoccat.testMissing.recode.vcf bbpopmap.txt 0.1 2 dDoccat.testMissing.pop
# After filtering, kept 166 out of 166 Individuals
# Outputting VCF file...
# After filtering, kept 12135 out of a possible 16157 Sites

# ./pop_missing_filter.sh dDoccat.testMissing.recode.vcf bbpopmap.txt 0.3 1 dDoccat.testMissing.pop
# ./pop_missing_filter.sh dDoccat.testMissing.recode.vcf bbpopmap.txt 30 24 dDoccat.testMissing.pop4 #this called with bad arguments, does nothing
# mawk '!/#/' dDoccat.testMissing.fils.FIL.recode.vcf | wc -l

# # This script will automatically filter a FreeBayes generated VCF file using criteria related to site depth,
# # quality versus depth, strand representation, allelic balance at heterzygous individuals, and paired read representation.
# # The script assumes that loci and individuals with low call rates (or depth) have already been removed.

curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/dDocent_filters
chmod +x dDocent_filters
# ./dDocent_filters
# sh dDocent_filters.sh VCF_file Output_prefix
./dDocent_filters dDoccat.testMissing.pop.recode.vcf dDoccat.testMissing.fils

# # Number of sites filtered based on allele balance at heterozygous loci, locus quality, and mapping quality / Depth
 # #Previous: 867 of 15041
  772 of 12135

# # Number of additional sites filtered based on overlapping forward and reverse reads
 # #Previous: 990 of 14174
 851 of 11363

# # Is this from a mixture of SE and PE libraries? Enter yes or no.
# # no
# # Number of additional sites filtered based on properly paired status
 # # Previous:1021 of 13184
  881 of 10512

 
# # Number of sites filtered based on high depth and lower than 2*DEPTH quality score
 # # Previous:422 of 12163
 332 of 9631


                                             # # Previous Histogram of mean depth per site

  # # 180 +-+-+---+---+---+---+---+---+--+---+---+---+---+---+---+---+---+---+---+---+---+---+--+---+---+---+---+---+-+-+
      # # +   +   +   +   +  **   +   +  +   +   +   +   +   +   +   +   +   +   +   +   +   +  +   +   +   +   +   +   +
      # # |              ** *****                             'meandepthpersite' using (bin($1,binwidth)):(1.0) ******* |
  # # 160 +-+           *** *****                                                                                     +-+
      # # |            ************                                                                                     |
  # # 140 +-+          ************     *** ***            **                                                         +-+
      # # |            **************   *** ***            **                                                           |
      # # |       ************************* ***    **  **  **  *****                                                    |
  # # 120 +-+     ***************************** ***** ***  ** ******                                                  +-+
      # # |       ************************************************** **                                                 |
  # # 100 +-+     ****************************************************** **                                           +-+
      # # |       ****************************************************** **   *******                                   |
      # # |       *********************************************************   *******       **                          |
   # # 80 +-+     **********************************************************  *******   **  **                        +-+
      # # |       *******************************************************************   **  **  **                      |
   # # 60 +-+     ********************************************************************* *** **  **                    +-+
      # # |       **************************************************************************** ****                     |
      # # |       ********************************************************************************* **   **   **        |
   # # 40 +-+     *************************************************************************************  **  ***      +-+
      # # |       **********************************************************************************************        |
   # # 20 +-+     ************************************************************************************************ *  +-+
      # # |       *******************************************************************************************************
      # # +   +   *******************************************************************************************************
    # # 0 +-+-+---*******************************************************************************************************
      # # 10  15  20  25  30  35  40  45 50  55  60  65  70  75  80  85  90  95 100 105 110 11 120 125 130 135 140 145 150
                                                        # # Mean Depth

# # If distrubtion looks normal, a 1.645 sigma cutoff (~90% of the data) would be 28230.2942
# # The 95% cutoff would be 143
# # Would you like to use a different maximum mean depth cutoff than 143, yes or no
# # no
                                             Histogram of mean depth per site                                                                                                                
                                                                                                                                                                                             
  120 +-+-+---+---+---+*--+---+---+--+---+---+---+---+---+---+---+---+---+---+--                                                                                                             -+---+---+--+---+---+---+---+---+-+-+
      +   +   +   +   +*  +   +   +  +   +   +   +   +   +   +   +   +   +   +                                                                                                                +   +   +  +   +   +   +   +   +   +
      |               ***  *****  **                      'meandepthpersite' usi                                                                                                             ng (bin($1,binwidth)):(1.0) ******* |
      |            ** ***  *****  ****                             **                                                                                                                                                            |
  100 +-+          ** ********** ***** **    **         *          **                                                                                                                                                          +-+
      |            ** ********** ***** ***   **         *          **                                                                                                                                                            |
      |          *************** ***** ********         *          **                                                                                                                                                            |
      |        *********************** ********         *  **  ******                                                                                                                                                            |
   80 +-+      *********************** ******** **  **  *************  **  **                                                                                                                                                  +-+
      |        *********************** ***********  *********************  **                                                                                                                  ****                              |
      |        ********************************************************** *** *                                                                                                               *******                            |
   60 +-+     *************************************************************** *                                                                                                               ********   **       **           +-+
      |       *****************************************************************                                                                                                               ********   **       **             |
      |       ******************************************************************                                                                                                             *********   *****   ***             |
      |       ******************************************************************                                                                                                             ********** ******   ****            |
   40 +-+     ******************************************************************                                                                                                             *************************         +-+
      |       ******************************************************************                                                                                                             *************************        ***|
      |       ******************************************************************                                                                                                             ******************************   ****
      |       ******************************************************************                                                                                                             *************************************
   20 +-+     ******************************************************************                                                                                                             *************************************
      |       ******************************************************************                                                                                                             *************************************
      |       ******************************************************************                                                                                                             *************************************
      +   +   ******************************************************************                                                                                                             *************************************
    0 +-+-+---******************************************************************                                                                                                             *************************************
      10  15  20  25  30  35  40  45 50  55  60  65  70  75  80  85  90  95 100                                                                                                              105 110 11 120 125 130 135 140 145 150
                                                        Mean Depth                                                                                                                           
                                                                                                                                                                                             
If distrubtion looks normal, a 1.645 sigma cutoff (~90% of the data) would be 28                                                                                                             953.6902
The 95% cutoff would be 168
Would you like to use a different maximum mean depth cutoff than 168, yes or no
no

# # Number of sites filtered based on maximum mean depth
 # # Previous:823 of 12163
  615 of 9631

# # Total number of sites filtered
 # # Previous: 3729 of 15041
  3158 of 12135
# # Remaining sites
 # # Previous: 11312
 8977


#filter departures from HWE within populations
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/filter_hwe_by_pop.pl
chmod +x filter_hwe_by_pop.pl
#filter our SNPs by population specific HWE
#decompose complex variant calls into phased SNP and INDEL genotypes and keep the INFO flags for loci and genotypes
vcfallelicprimitives dDoccat.testMissing.fils.FIL.recode.vcf --keep-info --keep-geno > dDoccat.testMissing.prim.vcf
#remove indels
vcftools --vcf dDoccat.testMissing.prim.vcf --remove-indels --recode --recode-INFO-all --out SNP.dDoccat.testMissing
# After filtering, kept 166 out of 166 Individuals
# Outputting VCF file...
# After filtering, kept 12427 out of a possible 12742 Sites


mawk '!/#/' dDoccat.testMissing.fils.FIL.recode.vcf | wc -l
#Previous:11312
8977
mawk '!/#/' dDoccat.testMissing.prim.vcf | wc -l
#Previous:12025
9602
mawk '!/#/' SNP.dDoccat.testMissing.recode.vcf | wc -l
#Previous:11741
9374



# #apply the HWE filter 
# #default -h, Minimum cutoff for Hardy-Weinberg p-value (for test as implemented in vcftools) [Default: 0.001]
# #default -c ,Proportion of all populations that a locus can be below HWE cutoff without being filtered. For example, choosing 0.5 willfilter SNPs that are below the p-value threshold in 50% or moreof the populations. [Default: 0.25]
./filter_hwe_by_pop.pl -v SNP.dDoccat.testMissing.recode.vcf -p bbpopmap.txt -o SNP.dDoccat.testMissing.HWE 
# #at default h and c, no loci filtered.

#estimate error
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/ErrorCount.sh
chmod +x ErrorCount.sh 
./ErrorCount.sh SNP.dDoccat.testMissing.recode.vcf
# Previous: This script counts the number of potential genotyping errors due to low read depth
# It report a low range, based on a 50% binomial probability of observing the second allele in a heterozygote and a high range based on a 25% probability.
# Potential genotyping errors from genotypes from only 1 read range from 0 to 0.0
# Potential genotyping errors from genotypes from only 2 reads range from 0 to 0.0
# Potential genotyping errors from genotypes from only 3 reads range from 940 to 3161.34
# Potential genotyping errors from genotypes from only 4 reads range from 561 to 2840.208
# Potential genotyping errors from genotypes from only 5 reads range from 319 to 2420
# 187 number of individuals and 11741 equals 2195567 total genotypes
# Total genotypes not counting missing data 2185445
# Total potential error rate is between 0.000832782339524 and 0.00385347057464
# SCORCHED EARTH SCENARIO
# WHAT IF ALL LOW DEPTH HOMOZYGOTE GENOTYPES ARE ERRORS?????
# The total SCORCHED EARTH error rate is 0.0122309186459.
Potential genotyping errors from genotypes from only 1 read range from 0 to 0.0
Potential genotyping errors from genotypes from only 2 reads range from 0 to 0.0
Potential genotyping errors from genotypes from only 3 reads range from 136 to 458.22
Potential genotyping errors from genotypes from only 4 reads range from 84 to 425.336
Potential genotyping errors from genotypes from only 5 reads range from 57 to 432
166 number of individuals and 9374 equals 1556084 total genotypes
Total genotypes not counting missing data 1554985
Total potential error rate is between 0.000178136766593 and 0.000846024881269
SCORCHED EARTH SCENARIO
WHAT IF ALL LOW DEPTH HOMOZYGOTE GENOTYPES ARE ERRORS?????
The total SCORCHED EARTH error rate is 0.00274021935903.

#######Fst with more stringent filtering########
cd Bluebonnet/dDoc_cat_Refiltering
mkdir testMissingFst
cp SNP.dDoccat.testMissing.recode.vcf ~/Bluebonnet/dDoc_cat_Refiltering/testMissingFst/testMissing.FinalSNP.vcf
 
# # #among all populations
vcftools --vcf ~/Bluebonnet/dDoc_cat_Refiltering/testMissingFst/testMissing.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1201Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1203Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1205Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1303Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1601Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1603Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1202Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1204Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1301Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1304Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1602Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1604Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1001Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1002Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1005Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1006Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1007Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1008Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1101Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1102Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1103Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1104Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1105Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1107Indiv.txt --out ~/Bluebonnet/dDoc_cat_Refiltering/testMissingFst/testMissing.FinalSNP.allFst
# Previous:Weir and Cockerham mean Fst estimate: 0.01573
# Previous: Weir and Cockerham weighted Fst estimate: 0.014794
Weir and Cockerham mean Fst estimate: 0.016371
Weir and Cockerham weighted Fst estimate: 0.015353

# # # #wild vs. seed
vcftools --vcf ~/Bluebonnet/dDoc_cat_Refiltering/testMissingFst/testMissing.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/wildIndiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/seedIndiv.txt --out ~/Bluebonnet/dDoc_cat_Refiltering/testMissingFst/testMissing.FinalSNP.wild_vs_seed
# Previous:Weir and Cockerham mean Fst estimate: 0.00064077
# Previous:Weir and Cockerham weighted Fst estimate: 0.00068971
Weir and Cockerham mean Fst estimate: 0.00084713
Weir and Cockerham weighted Fst estimate: 0.00087882

# # #among wild
vcftools --vcf ~/Bluebonnet/dDoc_cat_Refiltering/testMissingFst/testMissing.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1201Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1203Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1205Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1303Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1601Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1603Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1002Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1005Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1007Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1101Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1103Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1105Indiv.txt --out ~/Bluebonnet/dDoc_cat_Refiltering/testMissingFst/testMissing.FinalSNP.amongWildFst
# Previous:Weir and Cockerham mean Fst estimate: 0.017884
# Previous:Weir and Cockerham weighted Fst estimate: 0.017531
Weir and Cockerham mean Fst estimate: 0.018926
Weir and Cockerham weighted Fst estimate: 0.018648


# # #among seed
vcftools --vcf ~/Bluebonnet/dDoc_cat_Refiltering/testMissingFst/testMissing.FinalSNP.vcf --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1202Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1204Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1301Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1304Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1602Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1604Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1001Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1006Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1008Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1102Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1104Indiv.txt --weir-fst-pop ~/Bluebonnet/dDoc_Fst/1107Indiv.txt --out ~/Bluebonnet/dDoc_cat_Refiltering/testMissingFst/testMissing.FinalSNP.amongSeedFst
# Weir and Cockerham mean Fst estimate: 0.013034
# Weir and Cockerham weighted Fst estimate: 0.012504
Weir and Cockerham mean Fst estimate: 0.013355
Weir and Cockerham weighted Fst estimate: 0.012534

####prepping files for Dryad data archive####
#datafile description text
Supplementary material for "Homogenization of populations in the wildflower Texas bluebonnet (Lupinus texensis)." Here we sample 187 individuals in 12 matched pairs of neighboring wild and seeded populations of the Texas bluebonnet (Lupinus texensis), a species popular in commercially available wildflower seed mixes used by both the Texas Department of Transportation and the public. We use genotyping by sequencing to identify 11,741 genome-wide single nucleotide polymorphisms, as well as a smaller number of SNPs from the chloroplast genome, to analyze population structure and genetic diversity within and between the populations.
