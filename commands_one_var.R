## cline estimates

## evaluate ability of bgchm to estimate cline parameters
## comparing hierarchical and non-hierarchical model to HIest
## with perfect knowledge of genotypes and fixed allele frequency differences

library(bgchm)
library(HIest)

## data and true parameter values
load("o_fixed_known.rda")

## create objects for mean absolute error and cov
## non-hier, hier, HIest, and then each compared to genome-average 
mae_v<-matrix(NA,nrow=6,ncol=50)
mae_cent<-matrix(NA,nrow=6,ncol=50)

Lmae_v<-list(mae_v,mae_v,mae_v)
Lmae_cent<-list(mae_cent,mae_cent,mae_cent)

## non-hier, hier, coverage of true and then of genome-average
cov90_v<-matrix(NA,nrow=4,ncol=50)
cov90_cent<-matrix(NA,nrow=4,ncol=50)

Lcov90_v<-list(cov90_v,cov90_v,cov90_v)
Lcov90_cent<-list(cov90_cent,cov90_cent,cov90_cent)

## correlations with truth, non-hier, hier, HIest
cor_v<-matrix(NA,nrow=3,ncol=50)
cor_cent<-matrix(NA,nrow=3,ncol=50)

Lcor_v<-list(cor_v,cor_v,cor_v)
Lcor_cent<-list(cor_cent,cor_cent,cor_cent)

## hier params
est_SDC<-matrix(NA,nrow=50,ncol=3)
est_SDV<-matrix(NA,nrow=50,ncol=3)

L_SDC<-list(est_SDC,est_SDC,est_SDC)
L_SDV<-list(est_SDV,est_SDV,est_SDV)

## get rep number
myargs<-commandArgs(trailingOnly=TRUE)
k<-as.numeric(myargs[1])

#k<-8
#k<-3

## data set A
out_nh<-est_genocl(Gx=o_SDA_fixed_known[[k]]$g,model="genotype",H=o_SDA_fixed_known[[k]]$h,
		p0=rep(0.001,100),p1=rep(.999,100),ploidy="diploid",hier=FALSE,SDc=100,SDv=100)

out_hi<-est_genocl(Gx=o_SDA_fixed_known[[k]]$g,model="genotype",H=o_SDA_fixed_known[[k]]$h,
		p0=rep(0.001,100),p1=rep(.999,100),ploidy="diploid",hier=TRUE,sd0=2)
out_ll<-try(Cline.fit(Data=o_SDA_fixed_known[[k]]$g,S=o_SDA_fixed_known[[k]]$h,model="logit.logistic"),silent=TRUE)

X<-1

if(is(out_ll)[1]!="try-error"){
	Lmae_v[[X]][3,k]<-mean(abs(out_ll$logit.logistic$v-o_SDA_fixed_known[[k]]$v))
	Lmae_v[[X]][6,k]<-mean(abs(out_ll$logit.logistic$v-1))
	llc<-1/(1+exp(-1 * out_ll$logit.logistic$u/out_ll$logit.logistic$v))
	Lmae_cent[[X]][3,k]<-mean(abs(llc-o_SDA_fixed_known[[k]]$cent))
	Lmae_cent[[X]][6,k]<-mean(abs(llc-.5))
	Lcor_v[[X]][3,k]<-cor(out_ll$logit.logistic$v,o_SDA_fixed_known[[k]]$v)
	Lcor_cent[[X]][3,k]<-cor(llc,o_SDA_fixed_known[[k]]$cent)
}

Lmae_v[[X]][1,k]<-mean(abs(out_nh$gradient[,1]-o_SDA_fixed_known[[k]]$v))
Lmae_v[[X]][2,k]<-mean(abs(out_hi$gradient[,1]-o_SDA_fixed_known[[k]]$v))
Lmae_v[[X]][4,k]<-mean(abs(out_nh$gradient[,1]-1))
Lmae_v[[X]][5,k]<-mean(abs(out_hi$gradient[,1]-1))

Lmae_cent[[X]][1,k]<-mean(abs(out_nh$center[,1]-o_SDA_fixed_known[[k]]$cent))
Lmae_cent[[X]][2,k]<-mean(abs(out_hi$center[,1]-o_SDA_fixed_known[[k]]$cent))
Lmae_cent[[X]][4,k]<-mean(abs(out_nh$center[,1]-.5))
Lmae_cent[[X]][5,k]<-mean(abs(out_hi$center[,1]-.5))

Lcov90_v[[X]][1,k]<-mean(out_nh$gradient[,2] <= o_SDA_fixed_known[[k]]$v & out_nh$gradient[,3] >= o_SDA_fixed_known[[k]]$v)
Lcov90_v[[X]][2,k]<-mean(out_hi$gradient[,2] <= o_SDA_fixed_known[[k]]$v & out_hi$gradient[,3] >= o_SDA_fixed_known[[k]]$v)
Lcov90_v[[X]][3,k]<-mean(out_nh$gradient[,2] <= 1 & out_nh$gradient[,3] >= 1)
Lcov90_v[[X]][4,k]<-mean(out_hi$gradient[,2] <= 1 & out_hi$gradient[,3] >= 1)

Lcov90_cent[[X]][1,k]<-mean(out_nh$center[,2] <= o_SDA_fixed_known[[k]]$cent & out_nh$center[,3] >= o_SDA_fixed_known[[k]]$cent)
Lcov90_cent[[X]][2,k]<-mean(out_hi$center[,2] <= o_SDA_fixed_known[[k]]$cent & out_hi$center[,3] >= o_SDA_fixed_known[[k]]$cent)
Lcov90_cent[[X]][3,k]<-mean(out_nh$center[,2] <= .5 & out_nh$center[,3] >= .5)
Lcov90_cent[[X]][4,k]<-mean(out_hi$center[,2] <= .5 & out_hi$center[,3] >= .5)

Lcor_v[[X]][1,k]<-cor(out_nh$gradient[,1],o_SDA_fixed_known[[k]]$v)
Lcor_v[[X]][2,k]<-cor(out_hi$gradient[,1],o_SDA_fixed_known[[k]]$v)

Lcor_cent[[X]][1,k]<-cor(out_nh$center[,1],o_SDA_fixed_known[[k]]$cent)
Lcor_cent[[X]][2,k]<-cor(out_hi$center[,1],o_SDA_fixed_known[[k]]$cent)

L_SDC[[X]][k,]<-out_hi$SDc
L_SDV[[X]][k,]<-out_hi$SDv

## data set B 
out_nh<-est_genocl(Gx=o_SDB_fixed_known[[k]]$g,model="genotype",H=o_SDB_fixed_known[[k]]$h,
		p0=rep(0.001,100),p1=rep(.999,100),ploidy="diploid",hier=FALSE,SDc=100,SDv=100)

out_hi<-est_genocl(Gx=o_SDB_fixed_known[[k]]$g,model="genotype",H=o_SDB_fixed_known[[k]]$h,
		p0=rep(0.001,100),p1=rep(.999,100),ploidy="diploid",hier=TRUE,sd0=2)

out_ll<-try(Cline.fit(Data=o_SDB_fixed_known[[k]]$g,S=o_SDB_fixed_known[[k]]$h,model="logit.logistic"),silent=TRUE)

X<-2

if(is(out_ll)[1]!="try-error"){
	Lmae_v[[X]][3,k]<-mean(abs(out_ll$logit.logistic$v-o_SDB_fixed_known[[k]]$v))
	Lmae_v[[X]][6,k]<-mean(abs(out_ll$logit.logistic$v-1))
	llc<-1/(1+exp(-1 * out_ll$logit.logistic$u/out_ll$logit.logistic$v))
	Lmae_cent[[X]][3,k]<-mean(abs(llc-o_SDB_fixed_known[[k]]$cent))
	Lmae_cent[[X]][6,k]<-mean(abs(llc-.5))
	Lcor_v[[X]][3,k]<-cor(out_ll$logit.logistic$v,o_SDB_fixed_known[[k]]$v)
	Lcor_cent[[X]][3,k]<-cor(llc,o_SDB_fixed_known[[k]]$cent)
}

Lmae_v[[X]][1,k]<-mean(abs(out_nh$gradient[,1]-o_SDB_fixed_known[[k]]$v))
Lmae_v[[X]][2,k]<-mean(abs(out_hi$gradient[,1]-o_SDB_fixed_known[[k]]$v))
Lmae_v[[X]][4,k]<-mean(abs(out_nh$gradient[,1]-1))
Lmae_v[[X]][5,k]<-mean(abs(out_hi$gradient[,1]-1))

Lmae_cent[[X]][1,k]<-mean(abs(out_nh$center[,1]-o_SDB_fixed_known[[k]]$cent))
Lmae_cent[[X]][2,k]<-mean(abs(out_hi$center[,1]-o_SDB_fixed_known[[k]]$cent))
Lmae_cent[[X]][4,k]<-mean(abs(out_nh$center[,1]-.5))
Lmae_cent[[X]][5,k]<-mean(abs(out_hi$center[,1]-.5))

Lcov90_v[[X]][1,k]<-mean(out_nh$gradient[,2] <= o_SDB_fixed_known[[k]]$v & out_nh$gradient[,3] >= o_SDB_fixed_known[[k]]$v)
Lcov90_v[[X]][2,k]<-mean(out_hi$gradient[,2] <= o_SDB_fixed_known[[k]]$v & out_hi$gradient[,3] >= o_SDB_fixed_known[[k]]$v)
Lcov90_v[[X]][3,k]<-mean(out_nh$gradient[,2] <= 1 & out_nh$gradient[,3] >= 1)
Lcov90_v[[X]][4,k]<-mean(out_hi$gradient[,2] <= 1 & out_hi$gradient[,3] >= 1)

Lcov90_cent[[X]][1,k]<-mean(out_nh$center[,2] <= o_SDB_fixed_known[[k]]$cent & out_nh$center[,3] >= o_SDB_fixed_known[[k]]$cent)
Lcov90_cent[[X]][2,k]<-mean(out_hi$center[,2] <= o_SDB_fixed_known[[k]]$cent & out_hi$center[,3] >= o_SDB_fixed_known[[k]]$cent)
Lcov90_cent[[X]][3,k]<-mean(out_nh$center[,2] <= .5 & out_nh$center[,3] >= .5)
Lcov90_cent[[X]][4,k]<-mean(out_hi$center[,2] <= .5 & out_hi$center[,3] >= .5)

Lcor_v[[X]][1,k]<-cor(out_nh$gradient[,1],o_SDB_fixed_known[[k]]$v)
Lcor_v[[X]][2,k]<-cor(out_hi$gradient[,1],o_SDB_fixed_known[[k]]$v)

Lcor_cent[[X]][1,k]<-cor(out_nh$center[,1],o_SDB_fixed_known[[k]]$cent)
Lcor_cent[[X]][2,k]<-cor(out_hi$center[,1],o_SDB_fixed_known[[k]]$cent)

L_SDC[[X]][k,]<-out_hi$SDc
L_SDV[[X]][k,]<-out_hi$SDv


## data set C
out_nh<-est_genocl(Gx=o_SDC_fixed_known[[k]]$g,model="genotype",H=o_SDC_fixed_known[[k]]$h,
		p0=rep(0.001,100),p1=rep(0.999,100),ploidy="diploid",hier=FALSE,SDc=100,SDv=100)

out_hi<-est_genocl(Gx=o_SDC_fixed_known[[k]]$g,model="genotype",H=o_SDC_fixed_known[[k]]$h,
		p0=rep(0.001,100),p1=rep(0.999,100),ploidy="diploid",hier=TRUE,sd0=2)

out_ll<-try(Cline.fit(Data=o_SDC_fixed_known[[k]]$g,S=o_SDC_fixed_known[[k]]$h,model="logit.logistic"),silent=TRUE)



X<-3
if(is(out_ll)[1]!="try-error"){
	Lmae_v[[X]][3,k]<-mean(abs(out_ll$logit.logistic$v-o_SDC_fixed_known[[k]]$v))
	Lmae_v[[X]][6,k]<-mean(abs(out_ll$logit.logistic$v-1))
	llc<-1/(1+exp(-1 * out_ll$logit.logistic$u/out_ll$logit.logistic$v))
	Lmae_cent[[X]][3,k]<-mean(abs(llc-o_SDC_fixed_known[[k]]$cent))
	Lmae_cent[[X]][6,k]<-mean(abs(llc-.5))
	Lcor_v[[X]][3,k]<-cor(out_ll$logit.logistic$v,o_SDC_fixed_known[[k]]$v)
	Lcor_cent[[X]][3,k]<-cor(llc,o_SDC_fixed_known[[k]]$cent)
}

Lmae_v[[X]][1,k]<-mean(abs(out_nh$gradient[,1]-o_SDC_fixed_known[[k]]$v))
Lmae_v[[X]][2,k]<-mean(abs(out_hi$gradient[,1]-o_SDC_fixed_known[[k]]$v))
Lmae_v[[X]][4,k]<-mean(abs(out_nh$gradient[,1]-1))
Lmae_v[[X]][5,k]<-mean(abs(out_hi$gradient[,1]-1))

Lmae_cent[[X]][1,k]<-mean(abs(out_nh$center[,1]-o_SDC_fixed_known[[k]]$cent))
Lmae_cent[[X]][2,k]<-mean(abs(out_hi$center[,1]-o_SDC_fixed_known[[k]]$cent))
Lmae_cent[[X]][4,k]<-mean(abs(out_nh$center[,1]-.5))
Lmae_cent[[X]][5,k]<-mean(abs(out_hi$center[,1]-.5))

Lcov90_v[[X]][1,k]<-mean(out_nh$gradient[,2] <= o_SDC_fixed_known[[k]]$v & out_nh$gradient[,3] >= o_SDC_fixed_known[[k]]$v)
Lcov90_v[[X]][2,k]<-mean(out_hi$gradient[,2] <= o_SDC_fixed_known[[k]]$v & out_hi$gradient[,3] >= o_SDC_fixed_known[[k]]$v)
Lcov90_v[[X]][3,k]<-mean(out_nh$gradient[,2] <= 1 & out_nh$gradient[,3] >= 1)
Lcov90_v[[X]][4,k]<-mean(out_hi$gradient[,2] <= 1 & out_hi$gradient[,3] >= 1)

Lcov90_cent[[X]][1,k]<-mean(out_nh$center[,2] <= o_SDC_fixed_known[[k]]$cent & out_nh$center[,3] >= o_SDC_fixed_known[[k]]$cent)
Lcov90_cent[[X]][2,k]<-mean(out_hi$center[,2] <= o_SDC_fixed_known[[k]]$cent & out_hi$center[,3] >= o_SDC_fixed_known[[k]]$cent)
Lcov90_cent[[X]][3,k]<-mean(out_nh$center[,2] <= .5 & out_nh$center[,3] >= .5)
Lcov90_cent[[X]][4,k]<-mean(out_hi$center[,2] <= .5 & out_hi$center[,3] >= .5)

Lcor_v[[X]][1,k]<-cor(out_nh$gradient[,1],o_SDC_fixed_known[[k]]$v)
Lcor_v[[X]][2,k]<-cor(out_hi$gradient[,1],o_SDC_fixed_known[[k]]$v)

Lcor_cent[[X]][1,k]<-cor(out_nh$center[,1],o_SDC_fixed_known[[k]]$cent)
Lcor_cent[[X]][2,k]<-cor(out_hi$center[,1],o_SDC_fixed_known[[k]]$cent)

L_SDC[[X]][k,]<-out_hi$SDc
L_SDV[[X]][k,]<-out_hi$SDv


out<-paste("out_clinecomp",k,".rdat",sep="")
rm(out_nh)
rm(out_hi) ## just to shrink the data object
save(list=ls(),file=out)
