# lava_protocol
Local Analysis of [co]Variant Association (LAVA, Werme et al.) Analysis


# lava.R : analysis script for running LAVA across all genomic loci
This R script runs bivariate genetic correlation analysis across all genomic loci in the locus file which you've provided. If you are interested in analysing the local rg between a large amount of phenotypes, it is best to use a cluster computer. The script can be called from the command line:

```
lava.R "g1000_eur" "test.loci" "input.info.file" "sample.overlap.file" "cad;alz" "cad.alz"
```

your input files:
- ***g1000_eur*** : the reference genotype data
- ***test.loci*** : the locus definition file
- ***input.info.file*** : used for easy processing of multiple phenotypes
- ***sample.overlap.file*** : (optional) provide if your GWAS studies have overlapping subjects
- ***"cad;alz"*** : in this example, we analyse the local correlation between these two phenotypes: coronary artery disease (cad) and Alzheimer's disease (alz)
- ***"cad.alz"*** : specifies the prefix of the output files
