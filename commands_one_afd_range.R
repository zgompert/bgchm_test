## cline estimates

## evaluate ability of bgchm to estimate cline parameters
## comparing different degrees of minimum afd and uncertainty in genotypes

library(bgchm)

## data and true parameter values
load("o_range_genotypes.rda")

## create objects for mean absolute error and cov
## minimum  afd 1, 0.5, .1, genotypes then gl
mae_v<-matrix(NA,nrow=6,ncol=50)
mae_cent<-matrix(NA,nrow=6,ncol=50)

cov90_v<-matrix(NA,nrow=6,ncol=50)
cov90_cent<-matrix(NA,nrow=6,ncol=50)

cor_v<-matrix(NA,nrow=6,ncol=50)
cor_cent<-matrix(NA,nrow=6,ncol=50)

SDC<-matrix(NA,nrow=50,ncol=3)
SDV<-matrix(NA,nrow=50,ncol=3)

L_SDC<-list(SDC,SDC,SDC,SDC,SDC,SDC)
L_SDV<-list(SDV,SDV,SDV,SDV,SDV,SDV)



## get rep number
myargs<-commandArgs(trailingOnly=TRUE)
k<-as.numeric(myargs[1])

out<-vector("list",6)

## genotype models
out[[1]]<-est_genocl(Gx=o_afd_1.0[[k]]$g,model="genotype",H=o_afd_1.0[[k]]$h,
		p0=rep(0.001,100),p1=rep(.999,100),ploidy="diploid",hier=TRUE,sd0=2)
out[[2]]<-est_genocl(Gx=o_afd_0.5[[k]]$g,model="genotype",H=o_afd_0.5[[k]]$h,
		p0=p05,p1=p15,ploidy="diploid",hier=TRUE,sd0=2)
out[[3]]<-est_genocl(Gx=o_afd_0.1[[k]]$g,model="genotype",H=o_afd_0.1[[k]]$h,
		p0=p01,p1=p11,ploidy="diploid",hier=TRUE,sd0=2)


## gl models
N<-50;L<-100
## assumes mean coverge of 7 (Poisson) and 1% error rate
GenHybrids<-o_afd_1.0[[k]]$g
cmat<-matrix(rpois(prod(dim(GenHybrids)),7),nrow=N,ncol=L)
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

out[[4]]<-est_genocl(Gx=GLik,model="glik",H=o_afd_1.0[[k]]$h,
		p0=rep(0.001,100),p1=rep(.999,100),ploidy="diploid",hier=TRUE,sd0=2)

## assumes mean coverge of 7 (Poisson) and 1% error rate
GenHybrids<-o_afd_0.5[[k]]$g
cmat<-matrix(rpois(prod(dim(GenHybrids)),7),nrow=N,ncol=L)
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

out[[5]]<-est_genocl(Gx=GLik,model="glik",H=o_afd_0.5[[k]]$h,
		p0=p05,p1=p15,ploidy="diploid",hier=TRUE,sd0=2)

## assumes mean coverge of 7 (Poisson) and 1% error rate
GenHybrids<-o_afd_0.1[[k]]$g
cmat<-matrix(rpois(prod(dim(GenHybrids)),7),nrow=N,ncol=L)
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

out[[6]]<-est_genocl(Gx=GLik,model="glik",H=o_afd_0.1[[k]]$h,
		p0=p01,p1=p11,ploidy="diploid",hier=TRUE,sd0=2)

### save results
mae_v[1,k]<-mean(abs(out[[1]]$gradient[,1]-o_afd_1.0[[k]]$v))
mae_v[2,k]<-mean(abs(out[[2]]$gradient[,1]-o_afd_0.5[[k]]$v))
mae_v[3,k]<-mean(abs(out[[3]]$gradient[,1]-o_afd_0.1[[k]]$v))
mae_v[4,k]<-mean(abs(out[[4]]$gradient[,1]-o_afd_1.0[[k]]$v))
mae_v[5,k]<-mean(abs(out[[5]]$gradient[,1]-o_afd_0.5[[k]]$v))
mae_v[6,k]<-mean(abs(out[[6]]$gradient[,1]-o_afd_0.1[[k]]$v))

mae_cent[1,k]<-mean(abs(out[[1]]$center[,1]-o_afd_1.0[[k]]$cent))
mae_cent[2,k]<-mean(abs(out[[2]]$center[,1]-o_afd_0.5[[k]]$cent))
mae_cent[3,k]<-mean(abs(out[[3]]$center[,1]-o_afd_0.1[[k]]$cent))
mae_cent[4,k]<-mean(abs(out[[4]]$center[,1]-o_afd_1.0[[k]]$cent))
mae_cent[5,k]<-mean(abs(out[[5]]$center[,1]-o_afd_0.5[[k]]$cent))
mae_cent[6,k]<-mean(abs(out[[6]]$center[,1]-o_afd_0.1[[k]]$cent))

cov90_v[1,k]<-mean(out[[1]]$gradient[,2] <= o_afd_1.0[[k]]$v & out[[1]]$gradient[,3] >= o_afd_1.0[[k]]$v)
cov90_v[2,k]<-mean(out[[2]]$gradient[,2] <= o_afd_0.5[[k]]$v & out[[2]]$gradient[,3] >= o_afd_0.5[[k]]$v)
cov90_v[3,k]<-mean(out[[3]]$gradient[,2] <= o_afd_0.1[[k]]$v & out[[3]]$gradient[,3] >= o_afd_0.1[[k]]$v)
cov90_v[4,k]<-mean(out[[4]]$gradient[,2] <= o_afd_1.0[[k]]$v & out[[4]]$gradient[,3] >= o_afd_1.0[[k]]$v)
cov90_v[5,k]<-mean(out[[5]]$gradient[,2] <= o_afd_0.5[[k]]$v & out[[5]]$gradient[,3] >= o_afd_0.5[[k]]$v)
cov90_v[6,k]<-mean(out[[6]]$gradient[,2] <= o_afd_0.1[[k]]$v & out[[6]]$gradient[,3] >= o_afd_0.1[[k]]$v)

cov90_cent[1,k]<-mean(out[[1]]$center[,2] <= o_afd_1.0[[k]]$cent & out[[1]]$center[,3] >= o_afd_1.0[[k]]$cent)
cov90_cent[2,k]<-mean(out[[2]]$center[,2] <= o_afd_0.5[[k]]$cent & out[[2]]$center[,3] >= o_afd_0.5[[k]]$cent)
cov90_cent[3,k]<-mean(out[[3]]$center[,2] <= o_afd_0.1[[k]]$cent & out[[3]]$center[,3] >= o_afd_0.1[[k]]$cent)
cov90_cent[4,k]<-mean(out[[4]]$center[,2] <= o_afd_1.0[[k]]$cent & out[[4]]$center[,3] >= o_afd_1.0[[k]]$cent)
cov90_cent[5,k]<-mean(out[[5]]$center[,2] <= o_afd_0.5[[k]]$cent & out[[5]]$center[,3] >= o_afd_0.5[[k]]$cent)
cov90_cent[6,k]<-mean(out[[6]]$center[,2] <= o_afd_0.1[[k]]$cent & out[[6]]$center[,3] >= o_afd_0.1[[k]]$cent)

cor_v[1,k]<-cor(out[[1]]$gradient[,1],o_afd_1.0[[k]]$v)
cor_v[2,k]<-cor(out[[2]]$gradient[,1],o_afd_0.5[[k]]$v)
cor_v[3,k]<-cor(out[[3]]$gradient[,1],o_afd_0.1[[k]]$v)
cor_v[4,k]<-cor(out[[4]]$gradient[,1],o_afd_1.0[[k]]$v)
cor_v[5,k]<-cor(out[[5]]$gradient[,1],o_afd_0.5[[k]]$v)
cor_v[6,k]<-cor(out[[6]]$gradient[,1],o_afd_0.1[[k]]$v)

cor_cent[1,k]<-cor(out[[1]]$center[,1],o_afd_1.0[[k]]$cent)
cor_cent[2,k]<-cor(out[[2]]$center[,1],o_afd_0.5[[k]]$cent)
cor_cent[3,k]<-cor(out[[3]]$center[,1],o_afd_0.1[[k]]$cent)
cor_cent[4,k]<-cor(out[[4]]$center[,1],o_afd_1.0[[k]]$cent)
cor_cent[5,k]<-cor(out[[5]]$center[,1],o_afd_0.5[[k]]$cent)
cor_cent[6,k]<-cor(out[[6]]$center[,1],o_afd_0.1[[k]]$cent)


for(x in 1:6){
	L_SDC[[x]][k,]<-out[[x]]$SDc
	L_SDV[[x]][k,]<-out[[x]]$SDv
}

fout<-paste("out_range_gencomp",k,".rdat",sep="")
rm(out)
save(list=ls(),file=fout)
