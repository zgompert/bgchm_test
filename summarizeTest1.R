## combine and summarize the 50 cline models
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales)

ff<-list.files(path="o_clinecomp/",pattern="out_clinecomp",full.names=TRUE)

mae_v<-matrix(NA,nrow=6,ncol=50)
mae_cent<-matrix(NA,nrow=6,ncol=50)

Cmae_v<-list(mae_v,mae_v,mae_v)
Cmae_cent<-list(mae_cent,mae_cent,mae_cent)

## non-hier, hier, coverage of true and then of genome-average
cov90_v<-matrix(NA,nrow=4,ncol=50)
cov90_cent<-matrix(NA,nrow=4,ncol=50)

Ccov90_v<-list(cov90_v,cov90_v,cov90_v)
Ccov90_cent<-list(cov90_cent,cov90_cent,cov90_cent)

## correlations with truth, non-hier, hier, HIest
cor_v<-matrix(NA,nrow=3,ncol=50)
cor_cent<-matrix(NA,nrow=3,ncol=50)

Ccor_v<-list(cor_v,cor_v,cor_v)
Ccor_cent<-list(cor_cent,cor_cent,cor_cent)

## hier params
est_SDC<-matrix(NA,nrow=50,ncol=3)
est_SDV<-matrix(NA,nrow=50,ncol=3)

C_SDC<-list(est_SDC,est_SDC,est_SDC)
C_SDV<-list(est_SDV,est_SDV,est_SDV)


for(fi in 1:length(ff)){
        cat(fi,"\n")
        load(ff[fi])
	for(j in 1:3){
		Cmae_v[[j]][,k]<-Lmae_v[[j]][,k]
		Cmae_cent[[j]][,k]<-Lmae_cent[[j]][,k]
		Ccov90_v[[j]][,k]<-Lcov90_v[[j]][,k]
		Ccov90_cent[[j]][,k]<-Lcov90_cent[[j]][,k]
		Ccor_v[[j]][,k]<-Lcor_v[[j]][,k]
		Ccor_cent[[j]][,k]<-Lcor_cent[[j]][,k]
		C_SDC[[j]][k,]<-L_SDC[[j]][k,]
		C_SDV[[j]][k,]<-L_SDV[[j]][k,]
	}
}

SDc<-c(0.5,.8,1.2)
SDv<-c(.2,.4,.6)
save(list=ls(),file="example1.rdat")
#load("example1.rdat")
## make figure
cl<-1.6;ca<-1.15;cm<-1.4;ct<-1.25


pdf("F_clines.pdf",width=8.1,height=10.8)
layout(matrix(c(1,2,2,3,4,4,5,5,6,6,7,7,8,8,9,9),nrow=4,ncol=4,byrow=TRUE),widths=c(3,1.5,1.5,3),heights=c(3,3,3,3))
par(mar=c(4.5,5.5,2.5,1.5))

h <- seq(0, 1, 0.01)

center<-out_hiA$center[,1]
v<-out_hiA$gradient[,1]
u <- log(center/(1 - center)) * v
plot(h, h, type = "n", xlab = "Hybrid index", ylab = "Ancestry probability",xlim = c(0, 1), ylim = c(0, 1),cex.lab=cl,cex.axis=ca)
for (i in 1:length(v)) {
	phi <- (h^v[i])/(h^v[i] + (1 - h)^v[i] * exp(u[i]))
        lines(h, phi,col=alpha("darkgray",.5))
}
abline(a = 0, b = 1, lty = 2,lwd=1.5)
text(0.15,.87,expression(paste(sigma[c]," = 0.5"),sep=""),cex=ct,col="firebrick4")
text(0.15,.95,expression(paste(sigma[v]," = 0.2"),sep=""),cex=ct,col="firebrick4")

center<-out_hiB$center[,1]
v<-out_hiB$gradient[,1]
u <- log(center/(1 - center)) * v
plot(h, h, type = "n", xlab = "Hybrid index", ylab = "Ancestry probability",xlim = c(0, 1), ylim = c(0, 1),cex.lab=cl,cex.axis=ca)
for (i in 1:length(v)) {
	phi <- (h^v[i])/(h^v[i] + (1 - h)^v[i] * exp(u[i]))
        lines(h, phi,col=alpha("darkgray",.5))
}
abline(a = 0, b = 1, lty = 2,lwd=1.5)
text(0.15,.87,expression(paste(sigma[c]," = 0.8"),sep=""),cex=ct,col="firebrick4")
text(0.15,.95,expression(paste(sigma[v]," = 0.4"),sep=""),cex=ct,col="firebrick4")
title(main="(A) Example clines",cex.main=cm)

center<-out_hiC$center[,1]
v<-out_hiC$gradient[,1]
u <- log(center/(1 - center)) * v
plot(h, h, type = "n", xlab = "Hybrid index", ylab = "Ancestry probability",xlim = c(0, 1), ylim = c(0, 1),cex.lab=cl,cex.axis=ca)
for (i in 1:length(v)) {
	phi <- (h^v[i])/(h^v[i] + (1 - h)^v[i] * exp(u[i]))
        lines(h, phi,col=alpha("darkgray",.5))
}
abline(a = 0, b = 1, lty = 2,lwd=1.5)
text(0.15,.87,expression(paste(sigma[c]," = 1.2"),sep=""),cex=ct,col="firebrick4")
text(0.15,.95,expression(paste(sigma[v]," = 0.6"),sep=""),cex=ct,col="firebrick4")
#######################

lb<-min(c(C_SDV[[1]][,2],C_SDV[[2]][,2],C_SDV[[3]][,2])) - 0.025
ub<-max(c(C_SDV[[1]][,3],C_SDV[[2]][,3],C_SDV[[3]][,3])) + 0.025
cs<-rep(c("#bae4bc","#7bccc4","#2b8cbe"),each=50)
plot(c(C_SDV[[1]][,1],C_SDV[[2]][,1],C_SDV[[3]][,1]),pch=19,col=cs,ylim=c(lb,ub),cex.axis=ca,xlab="Data set",ylab="Standard deviation for v",cex.lab=cl)
segments(1:150,c(C_SDV[[1]][,2],C_SDV[[2]][,2],C_SDV[[3]][,2]),1:150,c(C_SDV[[1]][,3],C_SDV[[2]][,3],C_SDV[[3]][,3]),col=cs)
lines(c(1,50),c(SDv[1],SDv[1]),lty=3,lwd=1.5)
lines(c(51,100),c(SDv[2],SDv[2]),lty=3,lwd=1.5)
lines(c(101,150),c(SDv[3],SDv[3]),lty=3,lwd=1.5)
legend(1,.68,c(0.2,0.4,0.6),col=c("#bae4bc","#7bccc4","#2b8cbe"),pch=19,bty='n',cex=ca,title=expression(paste(sigma[v])))
title(main="(B) Cline v SD estimates",cex.main=cm)

lb<-min(c(C_SDC[[1]][,2],C_SDC[[2]][,2],C_SDC[[3]][,2])) - 0.05
ub<-max(c(C_SDC[[1]][,3],C_SDC[[2]][,3],C_SDC[[3]][,3])) + 0.05
plot(c(C_SDC[[1]][,1],C_SDC[[2]][,1],C_SDC[[3]][,1]),pch=19,col=cs,ylim=c(lb,ub),cex.axis=ca,xlab="Data set",ylab="Standard deviation for center",cex.lab=cl)
segments(1:150,c(C_SDC[[1]][,2],C_SDC[[2]][,2],C_SDC[[3]][,2]),1:150,c(C_SDC[[1]][,3],C_SDC[[2]][,3],C_SDC[[3]][,3]),col=cs)
lines(c(1,50),c(SDc[1],SDc[1]),lty=3,lwd=1.5)
lines(c(51,100),c(SDc[2],SDc[2]),lty=3,lwd=1.5)
lines(c(101,150),c(SDc[3],SDc[3]),lty=3,lwd=1.5)
legend(1,1.59,c(0.5,0.8,1.2),col=c("#bae4bc","#7bccc4","#2b8cbe"),pch=19,bty='n',cex=ca,title=expression(paste(sigma[c])))
title(main="(C) Cline center SD estimates",cex.main=cm)
####################
cs<-rep(c("#c2e699","#78c679","#238443"),3)
boxplot(t(rbind(Cmae_v[[1]][1:3,],Cmae_v[[2]][1:3,],Cmae_v[[3]][1:3,])),col="white",pch=NA,border=cs,axes=FALSE,cex.lab=cl,ylab="Absolute error",xlab="Cline standard deviation")
points(jitter(rep(1:9,50)),as.vector(rbind(Cmae_v[[1]][1:3,],Cmae_v[[2]][1:3,],Cmae_v[[3]][1:3,])),pch=19,col=rep(alpha(cs,.45),50))
axis(2,cex.axis=ca)
axis(1,at=c(2,5,8),c(0.2,0.4,0.6),cex.axis=ca)
legend(.45,27,c("bgchm simple","bgchm hierarchical","HIest"),col=alpha(cs[1:3],.6),pch=19,bty='n',cex=ca)
title(main="(D) Mean absolute error for v",cex.main=cm)
box()

boxplot(t(rbind(Cmae_cent[[1]][1:3,],Cmae_cent[[2]][1:3,],Cmae_cent[[3]][1:3,])),col="white",pch=NA,border=cs,axes=FALSE,cex.lab=cl,ylab="Absolute error",xlab="Cline standard deviation")
points(jitter(rep(1:9,50)),as.vector(rbind(Cmae_cent[[1]][1:3,],Cmae_cent[[2]][1:3,],Cmae_cent[[3]][1:3,])),pch=19,col=rep(alpha(cs,.45),50))
axis(2,cex.axis=ca)
axis(1,at=c(2,5,8),c(0.5,0.8,1.2),cex.axis=ca)
#legend(.55,.1,c("bgchm simple","bgchm hierarchical","HIest"),col=alpha(cs[1:3],.45),pch=19,bty='n')
title(main="(E) Mean absolute error for center",cex.main=cm)
box()

boxplot(t(rbind(Ccor_v[[1]][1:3,],Ccor_v[[2]][1:3,],Ccor_v[[3]][1:3,])),col="white",pch=NA,border=cs,axes=FALSE,cex.lab=cl,ylab="Correlation",xlab="Cline standard deviation",ylim=c(0,1))
points(jitter(rep(1:9,50)),as.vector(rbind(Ccor_v[[1]][1:3,],Ccor_v[[2]][1:3,],Ccor_v[[3]][1:3,])),pch=19,col=rep(alpha(cs,.45),50))
axis(2,cex.axis=ca)
axis(1,at=c(2,5,8),c(0.2,0.4,0.6),cex.axis=ca)
#legend(.55,1,c("bgchm simple","bgchm hierarchical","HIest"),col=alpha(cs[1:3],.6),pch=19,bty='n')
title(main="(F) Correlation with true v",cex.main=cm)
box()

boxplot(t(rbind(Ccor_cent[[1]][1:3,],Ccor_cent[[2]][1:3,],Ccor_cent[[3]][1:3,])),col="white",pch=NA,border=cs,axes=FALSE,cex.lab=cl,ylab="Correlation",xlab="Cline standard deviation",ylim=c(0,1))
points(jitter(rep(1:9,50)),as.vector(rbind(Ccor_cent[[1]][1:3,],Ccor_cent[[2]][1:3,],Ccor_cent[[3]][1:3,])),pch=19,col=rep(alpha(cs,.45),50))
axis(2,cex.axis=ca)
axis(1,at=c(2,5,8),c(0.5,0.8,1.2),cex.axis=ca)
#legend(.55,.1,c("bgchm simple","bgchm hierarchical","HIest"),col=alpha(cs[1:3],.45),pch=19,bty='n')
title(main="(F) Correlation with true center",cex.main=cm)
box()


dev.off()
