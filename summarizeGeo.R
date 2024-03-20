## combine and summarize geographic clines 
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales)
library(bgchm)

ff<-list.files(pattern="out_geo")

## create objects for mean absolute error and cov
Cmae_slope<-matrix(NA,nrow=3,ncol=100)
Cmae_cent<-matrix(NA,nrow=3,ncol=100)

Ccov90_slope<-matrix(NA,nrow=3,ncol=100)
Ccov90_cent<-matrix(NA,nrow=3,ncol=100)

## correlations with truth
Ccor_slope<-matrix(NA,nrow=3,ncol=100)
Ccor_cent<-matrix(NA,nrow=3,ncol=100)

## hier params
est_SDS<-matrix(NA,nrow=100,ncol=3)
C_SDS<-list(est_SDS,est_SDS,est_SDS)



for(fi in 1:length(ff)){
        cat(fi,"\n")
        load(ff[fi])
	Cmae_slope[,k]<-mae_slope[,k]
	Cmae_cent[,k]<-mae_cent[,k]
	Ccov90_slope[,k]<-cov90_slope[,k]
	Ccov90_cent[,k]<-cov90_cent[,k]
	Ccor_slope[,k]<-cor_slope[,k]
	Ccor_cent[,k]<-cor_cent[,k]
	for(j in 1:3){
		C_SDS[[j]][k,]<-L_SDS[[j]][k,]
	}
}

sds<-c(.1,.4,.8)

k<-1
out_low<-est_geocl(P=o_plow[[k]]$p,Geo=geo,prec=0.01,y_lb=-2,y_ub=2,hier=TRUE)
out_med<-est_geocl(P=o_pmed[[k]]$p,Geo=geo,prec=0.01,y_lb=-2,y_ub=2,hier=TRUE)
out_hi<-est_geocl(P=o_phi[[k]]$p,Geo=geo,prec=0.01,y_lb=-2,y_ub=2,hier=TRUE)


save(list=ls(),file="example_geo.rdat")
## make figure
cl<-1.6;ca<-1.15;cm<-1.4;ct<-1.25


pdf("SF_clines.pdf",width=8.1,height=10.8)
layout(matrix(c(1,2,2,3,4,4,4,4,5,5,6,6,7,7,8,8),nrow=4,ncol=4,byrow=TRUE),widths=c(3,1.5,1.5,3),heights=c(3,3,3,3))
par(mar=c(4.5,5.5,2.5,1.5))

cent<-out_low$cent[,1]
slope<-out_low$slope[,1]
plot(geo, rep(0,length(geo)), type = "n", xlab = "Location", ylab = "Logit allele frequency", ylim = c(-4,4),cex.lab=cl,cex.axis=ca)
for (i in 1:length(slope)) {
	y<-cent[i] + geo * slope[i]
	lines(geo, y,col=alpha("darkgray",.5))
}
text(-1.2,3.7,expression(paste(sigma," = 0.1"),sep=""),cex=ct,col="firebrick4")

cent<-out_med$cent[,1]
slope<-out_med$slope[,1]
plot(geo, rep(0,length(geo)), type = "n", xlab = "Location", ylab = "Logit allele frequency", ylim = c(-4,4),cex.lab=cl,cex.axis=ca)
for (i in 1:length(slope)) {
	y<-cent[i] + geo * slope[i]
	lines(geo, y,col=alpha("darkgray",.5))
}
text(-1.2,3.7,expression(paste(sigma," = 0.4"),sep=""),cex=ct,col="firebrick4")
title(main="(A) Example clines",cex.main=cm)

cent<-out_hi$cent[,1]
slope<-out_hi$slope[,1]
plot(geo, rep(0,length(geo)), type = "n", xlab = "Location", ylab = "Logit allele frequency", ylim = c(-4,4),cex.lab=cl,cex.axis=ca)
for (i in 1:length(slope)) {
	y<-cent[i] + geo * slope[i]
	lines(geo, y,col=alpha("darkgray",.5))
}
text(-1.2,3.7,expression(paste(sigma," = 0.8"),sep=""),cex=ct,col="firebrick4")

#######################

lb<-min(c(C_SDS[[1]][,2],C_SDS[[2]][,2],C_SDS[[3]][,2])) - 0.025
ub<-max(c(C_SDS[[1]][,3],C_SDS[[2]][,3],C_SDS[[3]][,3])) + 0.025
cs<-rep(c("#bae4bc","#7bccc4","#2b8cbe"),each=100)
plot(c(C_SDS[[1]][,1],C_SDS[[2]][,1],C_SDS[[3]][,1]),pch=19,col=cs,ylim=c(lb,ub),cex.axis=ca,xlab="Data set",ylab="Standard deviation",cex.lab=cl)
segments(1:300,c(C_SDS[[1]][,2],C_SDS[[2]][,2],C_SDS[[3]][,2]),1:300,c(C_SDS[[1]][,3],C_SDS[[2]][,3],C_SDS[[3]][,3]),col=cs)
lines(c(1,100),c(.1,.1),lty=3,lwd=1.5)
lines(c(101,200),c(.4,.4),lty=3,lwd=1.5)
lines(c(201,300),c(.8,.8),lty=3,lwd=1.5)
legend(1,1,c(0.1,0.4,0.8),col=c("#bae4bc","#7bccc4","#2b8cbe"),pch=19,bty='n',cex=ca,title=expression(paste(sigma)))
title(main="(B) Cline slope SD estimates",cex.main=cm)



####################
cs<-"#78c679"
boxplot(t(Cmae_slope),col="white",pch=NA,border=cs,axes=FALSE,cex.lab=cl,ylab="Absolute error",xlab="Cline standard deviation")
points(jitter(rep(1:3,100)),as.vector(Cmae_slope),pch=19,col=alpha(cs,.45))
axis(2,cex.axis=ca)
axis(1,at=1:3,sds,cex.axis=ca)
title(main="(C) Mean absolute error for slope",cex.main=cm)
box()

boxplot(t(Cmae_cent),col="white",pch=NA,border=cs,axes=FALSE,cex.lab=cl,ylab="Absolute error",xlab="Cline standard deviation")
points(jitter(rep(1:3,100)),as.vector(Cmae_cent),pch=19,col=alpha(cs,.45))
axis(2,cex.axis=ca)
axis(1,at=1:3,sds,cex.axis=ca)
title(main="(D) Mean absolute error for center",cex.main=cm)
box()

boxplot(t(Ccor_slope),col="white",pch=NA,border=cs,axes=FALSE,cex.lab=cl,ylab="Correlatin",xlab="Cline standard deviation",ylim=c(0,1))
points(jitter(rep(1:3,100)),as.vector(Ccor_slope),pch=19,col=alpha(cs,.45))
axis(2,cex.axis=ca)
axis(1,at=1:3,sds,cex.axis=ca)
title(main="(E) Correlation with true slope",cex.main=cm)
box()

boxplot(t(Ccor_cent),col="white",pch=NA,border=cs,axes=FALSE,cex.lab=cl,ylab="Correlatin",xlab="Cline standard deviation",ylim=c(0,1))
points(jitter(rep(1:3,100)),as.vector(Ccor_cent),pch=19,col=alpha(cs,.45))
axis(2,cex.axis=ca)
axis(1,at=1:3,sds,cex.axis=ca)
title(main="(F) Correlation with true center",cex.main=cm)
box()

dev.off()
