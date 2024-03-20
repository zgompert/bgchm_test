## evaluate ability of bgchm to estimate h and Q
## example for illustration with large number of hybrids
library(bgchm)


## read in simulations
dat<-read.table("df.main",comment.char="#",header=TRUE,sep=",")

## focus on generation 200
sdat<-dat[dat$gen==200,]

## plot of h (q) vs q12 (het)
## looks like a nice distribution to work with

GenHybrids<-as.matrix(sdat[,-c(1:8)])
L<-100
inds<-floor(seq(1,500,length.out=100))
GenHybrids<-GenHybrids[inds,sample(1:dim(GenHybrids)[2],L,replace=FALSE)]

## first, assumes fixed differences, parental allele frequencies known
h_out<-est_hi(Gx=GenHybrids,p0=rep(0,L),p1=rep(1,L),model="genotype",ploidy="diploid")

## now same for Q
q_out<-est_Q(Gx=GenHybrids,p0=rep(0,L),p1=rep(1,L),model="genotype",ploidy="diploid")

save(list=ls(),file="bgc_hq_ex.rda")

