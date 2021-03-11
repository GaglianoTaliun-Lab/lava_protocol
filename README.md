# lava_protocol
Local Analysis of [co]Variant Association (LAVA, Werme et al.) Analysis


# lava.R : analysis script for running LAVA across all genomic loci
This R script runs bivariate genetic correlation analysis across all genomic loci in the locus file which you've provided. If you are interested in analysing the local rg between a large amount of phenotypes, it is best to use a cluster computer. The script can be called from the command line:

```
lava.R "g1000_eur" "test.loci" "input.info.file" "sample.overlap.file" "cad;alz" "cad.alz"
```

Example input files:
- ***g1000_eur*** : the reference genotype data (e.g. 1000 genomes). Pre-processed input data can be found [here](https://ctg.cncr.nl/software/magma).
- ***test.loci*** : the locus definition file
- ***input.info.file*** : used for easy processing of multiple phenotypes
- ***sample.overlap.file*** : (optional) provide if your GWAS studies have overlapping subjects
- ***"cad;alz"*** : in this example, we analyse the local correlation between these two phenotypes: coronary artery disease (cad) and Alzheimer's disease (alz)
- ***"cad.alz"*** : specifies the prefix of the output files


**R Script**

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
