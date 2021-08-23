# LAVA
This tutorial will show you how to use Local Analysis of [co]Variant Association (LAVA, Werme et al.) Analysis to estimate local genetic correlation between coronary artery disease (CAD) and Alzheimer's disease (AD).

## Local genetic correlation between AD and CAD
Following the steps of the original [LAVA tutorial](https://github.com/josefin-werme/LAVA), I prepared the following files for input.
- ***Summary statistics*** for the traits of interest. For CAD, I downloaded the summary statistics from the [CARDIoGRAMplusC4D Consortium](http://www.cardiogramplusc4d.org/data-downloads/) 1000 Genomes-based GWAS. I used summary statistics for Alzheimer’s dementia from [Iris Jansen et al., 2019](https://www.nature.com/articles/s41588-018-0311-9), which are available for download from the [CNCR CTGlab](https://ctg.cncr.nl/software/summary_statistics). 
- ***Reference genome***. I used a reference genome of European ancestry, downloaded from the [CNCR CTGlab](https://ctg.cncr.nl/software/magma).  
- ***Sample overlap*** (optional). This file was left out of my analysis since there was no likely overlap between cohorts used for the two GWAS and the reference genome.
- ***Locus definition file***. 
- ***Input info file***
Note: the genomic build of all input files must match. Use [LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver/) for genomic build conversion. Here, both summary statistics and the reference genome are in build GRChr37/hg19.

## lava.R : analysis script for running LAVA across all genomic loci
This [R script](https://github.com/GaglianoTaliun-Lab/lava_protocol/blob/main/lava.R) runs a bivariate genetic correlation analysis across all genomic loci in the locus file which you've provided. It is best to use a cluster computer when analysing the local rg between a large amount of phenotypes. The script can be called from the command line:

```
lava.R "g1000_eur" "test.loci" "input.info.file" "sample.overlap.file" "cad;alz" "cad.alz"
```

In this example:
- ***g1000_eur*** : the reference genotype data (e.g. 1000 genomes). Pre-processed input data can be found [here](https://ctg.cncr.nl/software/magma).
- ***test.loci*** : the locus definition file
- ***input.info.file*** : used for easy processing of multiple phenotypes
- ***sample.overlap.file*** : (optional) provide if your GWAS studies have overlapping subjects
- ***"cad;alz"*** : in this example, we analyse the local correlation between these two phenotypes: coronary artery disease (cad) and Alzheimer's disease (alz)
- ***"cad.alz"*** : specifies the prefix of the output files



### **Complete R Script**

```
# command line arguments, specifying input/output file names and phenotype subset
arg = commandArgs(T); ref.prefix = arg[1]; loc.file = arg[2]; info.file = arg[3]; sample.overlap.file = arg[4]; phenos = unlist(strsplit(arg[5],";")); out.fname = arg[6]

### Load package
library(lava)

### Read in data
loci = read.loci(loc.file); n.loc = nrow(loci)
input = process.input(info.file, sample.overlap.file, ref.prefix, phenos)

print(paste("Starting LAVA analysis for",n.loc,"loci"))
progress = ceiling(quantile(1:n.loc, seq(.05,1,.05)))   # (if you want to print the progress)

### Analyse
u=b=list()
for (i in 1:n.loc) {
  if (i %in% progress) print(paste("..",names(progress[which(progress==i)])))     # (printing progress)
  locus = process.locus(loci[i,], input)                                          # process locus

  # It is possible that the locus cannot be defined for various reasons (e.g. too few SNPs), so the !is.null(locus) check is necessary before calling the analysis functions.
  if (!is.null(locus)) {
    # extract some general locus info for the output
    loc.info = data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)

    # run the univariate and bivariate tests
    loc.out = run.univ.bivar(locus, univ.thresh=1e-4)
    u[[i]] = cbind(loc.info, loc.out$univ)
    if(!is.null(loc.out$bivar)) b[[i]] = cbind(loc.info, loc.out$bivar)
  }
}

# save the output
write.table(do.call(rbind,u), paste0(out.fname,".univ.lava"), row.names=F,quote=F,col.names=T)
write.table(do.call(rbind,b), paste0(out.fname,".bivar.lava"), row.names=F,quote=F,col.names=T)

print(paste0("Done! Analysis output written to ",out.fname,".*.lava"))
```

# Sex Stratified Analysis

The .loci file for LAVA analysis can be updated by adding the locus descibed in the following articles:
- Deming et al. [article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6280657/) on sex specific predictors of Alzheimer’s disease.
- Dumitrescu et al. [article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6736148/) on sex differences in the genetic predictors of Alzheimer’s disease.
The genetic architecture of diseases can change based on sex. It would be more informative to perform genetic correlation analyses on men and women separately.

### Creating a Miami plot to visualize the differences in association signals between men and women, using EasyStrata and stratified GWAS results on the levels of tau found in the cerebospinal fluid of deceased Alzheimer's patients:

In R: 

```
library("EasyStrata")
EasyStrata("/lava/ecf/CSF_TAU_MIAMI.ecf")
```

The ecf file will take as input the summary statistics of the results to visualize. 
Example of the ecf file found [here](https://github.com/GaglianoTaliun-Lab/lava_protocol/blob/main/CSF_TAU_MIAMI.ecf).

##### Output: EasyStrata Miami Plot showing the differences in GWAS results between men and women
<img src="https://github.com/GaglianoTaliun-Lab/lava_protocol/blob/main/CSF_TAU_miami.png" width="600" height="300">

 
