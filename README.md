# BluebonnetGBS

Population genetic analyses of seeded and wild populations of Texas bluebonnet (Lupinus texensis).

Code and results associated with:

Turner K, Huang D, Cronk Q, Rieseberg L Data from: Homogenization of populations in the wildflower Texas bluebonnet (Lupinus texensis). Journal of Heredity (in press).

Data for this paper and repo is deposited at Dryad:

(dryad citation)

bbmaster_final_repo.sh is the primary coding file and is not actually set up to run as a shell script (do not be fooled by the file extension). It includes code for demultiplexing raw reads, calling SNPs using dDocent pipeline, filtering SNPs, estimating Fst and other population genetic stats using vcftools, setting up a STRUCTURE run, and finally re-filtering SNPs extra stringently and reassessing Fst results. bb.bars.names.txt includes the barcodes for each individual in the raw sequencing reads, for use when demultiplexing. bbpopmap.txt indicates the each individuals population. Files named like "1201Indiv.txt" list all the individuals in that group (population, or seeded/wild). Files ending in .fst, .pi, .hwe, or .het are output files from vcftools analysis. Remaining files are R scripts and their outputs including PCA analysis of SNPs and population genetic stats using the R package hierfstat.

Manuscript abstract:

Wildflowers seeds are routinely spread along highways and thoroughfares throughout North America as part of federal beautification policy, but the genetic effect of the introduction of these cultivated populations on wild populations of the same species is unknown. Interbreeding may occur between these seeded and wild populations, resulting in several possible outcomes. Here we sample 187 individuals in 12 matched pairs of neighboring wild and seeded populations of the Texas bluebonnet (Lupinus texensis), a species popular in commercially available wildflower seed mixes used by both the Texas Department of Transportation and the public. We use genotyping by sequencing to identify 11,741 genome-wide single nucleotide polymorphisms, as well as a smaller number of SNPs from the chloroplast genome, to analyze population structure and genetic diversity within and between the populations. We find a striking lack of population structure both between wild and seeded populations and amongst wild populations. STRUCTURE analyses indicate that all populations are apparently panmictic. This pattern may be explained by extensive swamping of wild populations by seeded germplasm and increased dispersal of semi-domesticated seed across this species’ core native range by humans. We discuss the possible negative and positive ramifications of homogenization on the evolutionary future of this popular wildflower species.
