## cline estimates

## clines for stronger polygenic selection simulations

library(bgchm)

## read simulations
dat<-read.table("df_strongpoly.main",header=TRUE,sep=",",comment.char="#")

## get rep number
myargs<-commandArgs(trailingOnly=TRUE)
k<-as.numeric(myargs[1])

sdat<-dat[dat$gen==5000 & dat$rep==k,]

hyb<-sample(1:1500,100,replace=FALSE)

hi<-sdat$q[hyb]
ghyb<-as.matrix(dat[hyb,-c(1:8)])

out_est<-est_genocl(Gx=ghyb,H=hi,p0=rep(0.001,251),p1=rep(.999,251),ploidy="diploid",hier=TRUE,sd0=2)

out<-paste("out_dfcline_strpoly",k,".rdat",sep="")
save(list=ls(),file=out)
