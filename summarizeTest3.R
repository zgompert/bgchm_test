## combine and summarize the 50 cline models
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales)

ff<-list.files(pattern="out_range")

## afd 1, 0.5, 0.1, known then gl
Cmae_v<-matrix(NA,nrow=6,ncol=50)
Cmae_cent<-matrix(NA,nrow=6,ncol=50)

Ccov90_v<-matrix(NA,nrow=6,ncol=50)
Ccov90_cent<-matrix(NA,nrow=6,ncol=50)

Ccor_v<-matrix(NA,nrow=6,ncol=50)
Ccor_cent<-matrix(NA,nrow=6,ncol=50)

## hier params
est_SDC<-matrix(NA,nrow=50,ncol=3)
est_SDV<-matrix(NA,nrow=50,ncol=3)

C_SDC<-list(est_SDC,est_SDC,est_SDC,est_SDC,est_SDC,est_SDC)
C_SDV<-list(est_SDV,est_SDV,est_SDV,est_SDV,est_SDV,est_SDV)


for(fi in 1:length(ff)){
        cat(fi,"\n")
        load(ff[fi])
	Cmae_v[,k]<-mae_v[,k]
	Cmae_cent[,k]<-mae_cent[,k]
	Ccov90_v[,k]<-cov90_v[,k]
	Ccov90_cent[,k]<-cov90_cent[,k]
	Ccor_v[,k]<-cor_v[,k]
	Ccor_cent[,k]<-cor_cent[,k]
	for(j in 1:6){
		C_SDC[[j]][k,]<-L_SDC[[j]][k,]
		C_SDV[[j]][k,]<-L_SDV[[j]][k,]
	}
}
SDc<-.7
SDv<-.3

sdc<-c(C_SDC[[1]][,1],C_SDC[[2]][,1],C_SDC[[3]][,1],C_SDC[[4]][,1],C_SDC[[5]][,1],C_SDC[[6]][,1])
sdclb<-c(C_SDC[[1]][,2],C_SDC[[2]][,2],C_SDC[[3]][,2],C_SDC[[4]][,2],C_SDC[[5]][,2],C_SDC[[6]][,2])
sdcub<-c(C_SDC[[1]][,3],C_SDC[[2]][,3],C_SDC[[3]][,3],C_SDC[[4]][,3],C_SDC[[5]][,3],C_SDC[[6]][,3])

sdv<-c(C_SDV[[1]][,1],C_SDV[[2]][,1],C_SDV[[3]][,1],C_SDV[[4]][,1],C_SDV[[5]][,1],C_SDV[[6]][,1])
sdvlb<-c(C_SDV[[1]][,2],C_SDV[[2]][,2],C_SDV[[3]][,2],C_SDV[[4]][,2],C_SDV[[5]][,2],C_SDV[[6]][,2])
sdvub<-c(C_SDV[[1]][,3],C_SDV[[2]][,3],C_SDV[[3]][,3],C_SDV[[4]][,3],C_SDV[[5]][,3],C_SDV[[6]][,3])

xtable(cbind(round(apply(Cmae_v,1,mean),3),round(apply(Cmae_cent,1,mean),3),
	     round(apply(Ccor_v,1,mean),3),round(apply(Ccor_cent,1,mean),3),
	round(apply(Ccov90_v,1,mean),3),round(apply(Ccov90_cent,1,mean),3)))


cs<-rep(c("#f6e8c3","#d8b365","#8c510a","#c7eae5","#5ab4ac","#01665e"),each=50)

pdf("F_clineAFD.pdf",width=9,height=10.5)
par(mfrow=c(3,2))
par(mar=c(4.5,4.5,2.5,1))
cl<-1.62;ca<-1.1;cm<-1.575

plot(sdv,pch=19,ylim=c(min(sdvlb),max(sdvub)),xlab="Data set",ylab="Standard deviation for v",cex.lab=cl,cex.axis=ca,col=cs)
segments(1:300,sdvlb,1:300,sdvub,col=cs)
abline(h=SDv,lwd=1.2,lty=3)
title(main="(A) Cline v SD estimates",cex.main=cm)

plot(sdc,pch=19,ylim=c(min(sdclb),max(sdcub)),xlab="Data set",ylab="Standard deviation for center",cex.lab=cl,cex.axis=ca,col=cs)
segments(1:300,sdclb,1:300,sdcub,col=cs)
abline(h=SDc,lwd=1.2,lty=3)
title(main="(B) Cline center SD estiamtes",cex.main=cm)

cs<-unique(cs)

boxplot(t(Cmae_v),col="white",pch=NA,border=cs,axes=FALSE,cex.lab=cl,ylab="Absolute error",xlab="Minimum difference")
points(jitter(rep(1:6,50)),as.vector(Cmae_v),pch=19,col=rep(alpha(cs,.45),50))
axis(2,cex.axis=ca)
axis(1,at=1:6,c(1,.5,.1,1,.5,.1),cex.axis=ca)
title(main="(C) Mean absolute error for v",cex.main=cm)
box()

boxplot(t(Cmae_cent),col="white",pch=NA,border=cs,axes=FALSE,cex.lab=cl,ylab="Absolute error",xlab="Minimum difference")
points(jitter(rep(1:6,50)),as.vector(Cmae_cent),pch=19,col=rep(alpha(cs,.45),50))
axis(2,cex.axis=ca)
axis(1,at=1:6,c(1,.5,.1,1,.5,.1),cex.axis=ca)
title(main="(D) Mean absolute error for center",cex.main=cm)
box()

boxplot(t(Ccor_v),col="white",pch=NA,border=cs,axes=FALSE,ylim=c(0,1),cex.lab=cl,ylab="Correlation",xlab="Minimum difference")
points(jitter(rep(1:6,50)),as.vector(Ccor_v),pch=19,col=rep(alpha(cs,.45),50))
axis(2,cex.axis=ca)
axis(1,at=1:6,c(1,.5,.1,1,.5,.1),cex.axis=ca)
title(main="(E) Correlation for v",cex.main=cm)
box()
legend(0.4,0.25,c("1, known","0.5, known","0.1, known","1, uncertain","0.5, uncertain","0.1, uncertain"),col=cs,pch=19,bty='n',ncol=2,cex=ca*1.1)

boxplot(t(Ccor_cent),col="white",pch=NA,border=cs,axes=FALSE,ylim=c(0,1),cex.lab=cl,ylab="Absolute error",xlab="Minimum difference")
points(jitter(rep(1:6,50)),as.vector(Ccor_cent),pch=19,col=rep(alpha(cs,.45),50))
axis(2,cex.axis=ca)
axis(1,at=1:6,c(1,.5,.1,1,.5,.1),cex.axis=ca)
title(main="(F) Correlation for center",cex.main=cm)
box()

dev.off()

