## cline estimates

## evaluate ability of bgchm to estimate cline parameters
## comparing hierarchical and non-hierarchical model to HIest
## with perfect knowledge of genotypes and fixed allele frequency differences

library(bgchm)
library(HIest)

## data and true parameter values
load("o_fixed_known.rda")

k<-1
## data set A

out_hiA<-est_genocl(Gx=o_SDA_fixed_known[[k]]$g,model="genotype",H=o_SDA_fixed_known[[k]]$h,
		p0=rep(0.001,100),p1=rep(.999,100),ploidy="diploid",hier=TRUE,sd0=2)
## data set B 

out_hiB<-est_genocl(Gx=o_SDB_fixed_known[[k]]$g,model="genotype",H=o_SDB_fixed_known[[k]]$h,
		p0=rep(0.001,100),p1=rep(.999,100),ploidy="diploid",hier=TRUE,sd0=2)


## data set C

out_hiC<-est_genocl(Gx=o_SDC_fixed_known[[k]]$g,model="genotype",H=o_SDC_fixed_known[[k]]$h,
		p0=rep(0.001,100),p1=rep(0.999,100),ploidy="diploid",hier=TRUE,sd0=2)

save(list=ls(),file="example1.rdat")
