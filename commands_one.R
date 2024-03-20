## run hybrid index and q for one replicate
## do this with different amounts of ancestry information
## including with genotype likelihoods

## evaluate ability of bgchm to estimate h and Q
library(bgchm)

## read in simulations
dat<-read.table("repdf.main",comment.char="#",header=TRUE,sep=",")

## create objects for mean absolute error and cov
## rows are = fixed with genotypes, Fst = 0.5 or 0.1 with genotypes
## fixed with likelihoods, Fst = 0.5 or 0.1  with likelihoods
mae_h<-matrix(NA,nrow=6,ncol=50)
mae_q10<-matrix(NA,nrow=6,ncol=50)
cov90_h<-matrix(NA,nrow=6,ncol=50)
cov90_q10<-matrix(NA,nrow=6,ncol=50)

## get rep number
myargs<-commandArgs(trailingOnly=TRUE)
k<-as.numeric(myargs[1])

#k<-1

sdat<-dat[dat$rep==k,]

N<-dim(sdat)[1]
L<-dim(sdat)[2]-8

## sample 50 hybrids
hybs<-sample(1:N,50,replace=FALSE)

GenHybrids<-as.matrix(sdat[hybs,-c(1:8)])

## first, assumes fixed differences, parental allele frequencies known
h_out<-est_hi(Gx=GenHybrids,p0=rep(0,L),p1=rep(1,L),model="genotype",ploidy="diploid")
mae_h[1,k]<-mean(abs(h_out$hi[,1]-sdat$q[hybs]))
cov90_h[1,k]<-mean(round(h_out$hi[,3],2) <= sdat$q[hybs] & round(h_out$hi[,4],2) >= sdat$q[hybs])

q_out<-est_Q(Gx=GenHybrids,p0=rep(0,L),p1=rep(1,L),model="genotype",ploidy="diploid")
mae_q10[1,k]<-mean(abs(q_out$Q10[,1]-sdat$het[hybs]))
cov90_q10[1,k]<-mean(round(q_out$Q10[,2],2) <= sdat$het[hybs] & round(q_out$Q10[,3],2) >= sdat$het[hybs])


## now with afd = .5
p0<-rep(0.25,L)
p1<-rep(0.75,L)

GenHybridsFpt5<-GenHybrids
for(i in 1:L){for(j in 1:50){
	GenHybridsFpt5[j,i]<-sum(c(rbinom(n=2-GenHybrids[j,i],size=1,prob=p0[i]), rbinom(n=GenHybrids[j,i],size=1,prob=p1[i])))
}}

h_out<-est_hi(Gx=GenHybridsFpt5,p0=p0,p1=p1,model="genotype",ploidy="diploid")
mae_h[2,k]<-mean(abs(h_out$hi[,1]-sdat$q[hybs]))
cov90_h[2,k]<-mean(round(h_out$hi[,3],2) <= sdat$q[hybs] & round(h_out$hi[,4],2) >= sdat$q[hybs])

q_out<-est_Q(Gx=GenHybridsFpt5,p0=p0,p1=p1,model="genotype",ploidy="diploid")
mae_q10[2,k]<-mean(abs(q_out$Q10[,1]-sdat$het[hybs]))
cov90_q10[2,k]<-mean(round(q_out$Q10[,2],2) <= sdat$het[hybs] & round(q_out$Q10[,3],2) >= sdat$het[hybs])

## now with afd = .1
p0<-rep(0.45,L)
p1<-rep(0.55,L)

GenHybridsFpt1<-GenHybrids
for(i in 1:L){for(j in 1:50){
	GenHybridsFpt1[j,i]<-sum(c(rbinom(n=2-GenHybrids[j,i],size=1,prob=p0[i]), rbinom(n=GenHybrids[j,i],size=1,prob=p1[i])))
}}

h_out<-est_hi(Gx=GenHybridsFpt1,p0=p0,p1=p1,model="genotype",ploidy="diploid")
mae_h[3,k]<-mean(abs(h_out$hi[,1]-sdat$q[hybs]))
cov90_h[3,k]<-mean(round(h_out$hi[,3],2) <= sdat$q[hybs] & round(h_out$hi[,4],2) >= sdat$q[hybs])

q_out<-est_Q(Gx=GenHybridsFpt1,p0=p0,p1=p1,model="genotype",ploidy="diploid")
mae_q10[3,k]<-mean(abs(q_out$Q10[,1]-sdat$het[hybs]))
cov90_q10[3,k]<-mean(round(q_out$Q10[,2],2) <= sdat$het[hybs] & round(q_out$Q10[,3],2) >= sdat$het[hybs])

############# genotype likelihoods #######################
## assumes mean coverge of 7 (Poisson) and 1% error rate
cmat<-matrix(rpois(prod(dim(GenHybrids)),7),nrow=50,ncol=L)
Y<-GenHybrids/2
Y[Y==1]<-.99
Y[Y==0]<-.01
Gy<-rbinom(n=prod(dim(Y)),size=cmat,prob=Y)


GLik<-list(gl0=GenHybrids,gl1=GenHybrids,gl2=GenHybrids)
GLik$gl0<-matrix(dbinom(x=Gy,size=cmat,prob=0.01),nrow=50,ncol=L)
GLik$gl1<-matrix(dbinom(x=Gy,size=cmat,prob=0.5),nrow=50,ncol=L)
GLik$gl2<-matrix(dbinom(x=Gy,size=cmat,prob=0.99),nrow=50,ncol=L)

Gsum<-GLik$gl0+GLik$gl1+GLik$gl2
GLik$gl0<-GLik$gl0/Gsum
GLik$gl1<-GLik$gl1/Gsum
GLik$gl2<-GLik$gl2/Gsum

## first, assumes fixed differences, parental allele frequencies known
h_out<-est_hi(Gx=GLik,p0=rep(0,L),p1=rep(1,L),model="glik",ploidy="diploid")
mae_h[4,k]<-mean(abs(h_out$hi[,1]-sdat$q[hybs]))
cov90_h[4,k]<-mean(round(h_out$hi[,3],2) <= sdat$q[hybs] & round(h_out$hi[,4],2) >= sdat$q[hybs])

q_out<-est_Q(Gx=GLik,p0=rep(0,L),p1=rep(1,L),model="glik",ploidy="diploid")
mae_q10[4,k]<-mean(abs(q_out$Q10[,1]-sdat$het[hybs]))
cov90_q10[4,k]<-mean(round(q_out$Q10[,2],2) <= sdat$het[hybs] & round(q_out$Q10[,3],2) >= sdat$het[hybs])

## now with afd = .5
p0<-rep(0.25,L)
p1<-rep(0.75,L)
GenHybridsFpt5<-GenHybrids
for(i in 1:L){for(j in 1:50){
	GenHybridsFpt5[j,i]<-sum(c(rbinom(n=2-GenHybrids[j,i],size=1,prob=p0[i]), rbinom(n=GenHybrids[j,i],size=1,prob=p1[i])))
}}

## assumes mean coverge of 7 (Poisson) and 1% error rate
cmat<-matrix(rpois(prod(dim(GenHybridsFpt5)),7),nrow=50,ncol=L)
Y<-GenHybridsFpt5/2
Y[Y==1]<-.99
Y[Y==0]<-.01
Gy<-rbinom(n=prod(dim(Y)),size=cmat,prob=Y)


GLik<-list(gl0=GenHybrids,gl1=GenHybrids,gl2=GenHybrids)
GLik$gl0<-matrix(dbinom(x=Gy,size=cmat,prob=0.01),nrow=50,ncol=L)
GLik$gl1<-matrix(dbinom(x=Gy,size=cmat,prob=0.5),nrow=50,ncol=L)
GLik$gl2<-matrix(dbinom(x=Gy,size=cmat,prob=0.99),nrow=50,ncol=L)


h_out<-est_hi(Gx=GLik,p0=p0,p1=p1,model="glik",ploidy="diploid")
mae_h[5,k]<-mean(abs(h_out$hi[,1]-sdat$q[hybs]))
cov90_h[5,k]<-mean(round(h_out$hi[,3],2) <= sdat$q[hybs] & round(h_out$hi[,4],2) >= sdat$q[hybs])

q_out<-est_Q(Gx=GLik,p0=p0,p1=p1,model="glik",ploidy="diploid")
mae_q10[5,k]<-mean(abs(q_out$Q10[,1]-sdat$het[hybs]))
cov90_q10[5,k]<-mean(round(q_out$Q10[,2],2) <= sdat$het[hybs] & round(q_out$Q10[,3],2) >= sdat$het[hybs])

## now with afd = .1
p0<-rep(0.45,L)
p1<-rep(0.55,L)
GenHybridsFpt1<-GenHybrids
for(i in 1:L){for(j in 1:50){
	GenHybridsFpt1[j,i]<-sum(c(rbinom(n=2-GenHybrids[j,i],size=1,prob=p0[i]), rbinom(n=GenHybrids[j,i],size=1,prob=p1[i])))
}}

## assumes mean coverge of 7 (Poisson) and 1% error rate
cmat<-matrix(rpois(prod(dim(GenHybridsFpt1)),7),nrow=50,ncol=L)
Y<-GenHybridsFpt1/2
Y[Y==1]<-.99
Y[Y==0]<-.01
Gy<-rbinom(n=prod(dim(Y)),size=cmat,prob=Y)


GLik<-list(gl0=GenHybrids,gl1=GenHybrids,gl2=GenHybrids)
GLik$gl0<-matrix(dbinom(x=Gy,size=cmat,prob=0.01),nrow=50,ncol=L)
GLik$gl1<-matrix(dbinom(x=Gy,size=cmat,prob=0.5),nrow=50,ncol=L)
GLik$gl2<-matrix(dbinom(x=Gy,size=cmat,prob=0.99),nrow=50,ncol=L)

h_out<-est_hi(Gx=GLik,p0=p0,p1=p1,model="glik",ploidy="diploid")
mae_h[6,k]<-mean(abs(h_out$hi[,1]-sdat$q[hybs]))
cov90_h[6,k]<-mean(round(h_out$hi[,3],2) <= sdat$q[hybs] & round(h_out$hi[,4],2) >= sdat$q[hybs])

q_out<-est_Q(Gx=GLik,p0=p0,p1=p1,model="glik",ploidy="diploid")
mae_q10[6,k]<-mean(abs(q_out$Q10[,1]-sdat$het[hybs]))
cov90_q10[6,k]<-mean(round(q_out$Q10[,2],2) <= sdat$het[hybs] & round(q_out$Q10[,3],2) >= sdat$het[hybs])

out<-paste("out_hq",k,".rdat",sep="")
rm(h_out)
rm(q_out) ## just to shrink the data object
save(list=ls(),file=out)
