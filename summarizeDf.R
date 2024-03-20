library(bgchm)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales)

ff_neu<-list.files(pattern="out_dfcline_neu")
ff_oligo<-list.files(pattern="out_dfcline_oligo")
ff_poly<-list.files(pattern="out_dfcline_poly")
ff_strp<-list.files(pattern="out_dfcline_strpoly")

cent_neu<-vector("list",10)
cent_oligo<-vector("list",10)
cent_poly<-vector("list",10)
cent_spoly<-vector("list",10)
v_neu<-vector("list",10)
v_oligo<-vector("list",10)
v_poly<-vector("list",10)
v_spoly<-vector("list",10)

sdc_neu<-matrix(NA,nrow=10,ncol=3)
sdv_neu<-matrix(NA,nrow=10,ncol=3)

sdc_oligo<-matrix(NA,nrow=10,ncol=3)
sdv_oligo<-matrix(NA,nrow=10,ncol=3)

sdc_poly<-matrix(NA,nrow=10,ncol=3)
sdv_poly<-matrix(NA,nrow=10,ncol=3)

sdc_spoly<-matrix(NA,nrow=10,ncol=3)
sdv_spoly<-matrix(NA,nrow=10,ncol=3)

for(x in 1:10){
	cat(x,"\n")
	load(ff_neu[x])
	cent_neu[[x]]<-out_est$center
	v_neu[[x]]<-out_est$gradient
	sdc_neu[x,]<-out_est$SDc
	sdv_neu[x,]<-out_est$SDv
	
	load(ff_oligo[x])
	cent_oligo[[x]]<-out_est$center
	v_oligo[[x]]<-out_est$gradient
	sdc_oligo[x,]<-out_est$SDc
	sdv_oligo[x,]<-out_est$SDv
	
	load(ff_poly[x])
	cent_poly[[x]]<-out_est$center
	v_poly[[x]]<-out_est$gradient
	sdc_poly[x,]<-out_est$SDc
	sdv_poly[x,]<-out_est$SDv
	
	load(ff_strp[x])
	cent_spoly[[x]]<-out_est$center
	v_spoly[[x]]<-out_est$gradient
	sdc_spoly[x,]<-out_est$SDc
	sdv_spoly[x,]<-out_est$SDv
}

csn<-c(brewer.pal("Greens",n=9),"green")
cso<-c(brewer.pal("Oranges",n=9),"red")
csp<-c(brewer.pal("Blues",n=9),"cadetblue")
css<-c(brewer.pal("Purples",n=9),"violet")

dat_neu<-read.table("df_neu.main",header=TRUE,sep=",",comment.char="#")
dat_oligo<-read.table("df_oligo.main",header=TRUE,sep=",",comment.char="#")
dat_poly<-read.table("df_poly.main",header=TRUE,sep=",",comment.char="#")
dat_spoly<-read.table("df_strongpoly.main",header=TRUE,sep=",",comment.char="#")

cl<-1.4;ca<-1.1;cm<-1.4
pdf("F_dfclines.pdf",width=10,height=7.5)
par(mfrow=c(3,4))
par(mar=c(4,4.5,2,1.5))
## hybrid index boxplots
boxplot(dat_neu$q[dat_neu$gen==5000] ~ dat_neu$deme[dat_neu$gen==5000],pch=NA,xlab="Deme",ylab="Hybrid index",border=csn[5],col="white",cex.lab=cl,cex.axis=ca)
points(jitter(dat_neu$deme[dat_neu$gen==5000]),dat_neu$q[dat_neu$gen==5000],pch=19,col=alpha(csn[5],.05))
title(main="(A) HI neutral",cex.main=cm)

boxplot(dat_oligo$q[dat_oligo$gen==5000] ~ dat_oligo$deme[dat_oligo$gen==5000],pch=NA,xlab="Deme",ylab="Hybrid index",border=cso[5],col="white",cex.lab=cl,cex.axis=ca)
points(jitter(dat_oligo$deme[dat_oligo$gen==5000]),dat_oligo$q[dat_oligo$gen==5000],pch=19,col=alpha(cso[5],.05))
title(main="(B) HI oligogenic",cex.main=cm)

boxplot(dat_poly$q[dat_poly$gen==5000] ~ dat_poly$deme[dat_poly$gen==5000],pch=NA,xlab="Deme",ylab="Hybrid index",border=csp[5],col="white",cex.lab=cl,cex.axis=ca)
points(jitter(dat_poly$deme[dat_poly$gen==5000]),dat_poly$q[dat_poly$gen==5000],pch=19,col=alpha(csp[5],.05))
title(main="(C) HI weak polygenic",cex.main=cm)

boxplot(dat_spoly$q[dat_spoly$gen==5000] ~ dat_spoly$deme[dat_spoly$gen==5000],pch=NA,xlab="Deme",ylab="Hybrid index",border=css[5],col="white",cex.lab=cl,cex.axis=ca)
points(jitter(dat_spoly$deme[dat_spoly$gen==5000]),dat_spoly$q[dat_spoly$gen==5000],pch=19,col=alpha(css[5],.05))
title(main="(D) HI strong polygenic",cex.main=cm)

o<-sum2zero(center=cent_neu[[1]],v=v_neu[[1]],transform=TRUE)
sig<-(o$gradient[,2] > 1) +1
h <- seq(0, 1, 0.01)
center=o$cent[,1];v=o$gradient[,1]
u <- log(center/(1 - center)) * v
plot(h, h, type = "n", xlab = "Hybrid index", ylab = "Ancestry probability", xlim = c(0, 1), ylim = c(0, 1),cex.lab=cl)
for (i in 1:length(v)) {
	phi <- (h^v[i])/(h^v[i] + (1 - h)^v[i] * exp(u[i]))
        lines(h, phi,col=alpha(csn[5],sig[i]*.5),lwd=sig[i]*.4)
}
abline(a = 0, b = 1, lty = 2)
title(main="(E) Clines neutral",cex.main=cm)

o<-sum2zero(center=cent_oligo[[1]],v=v_oligo[[1]],transform=TRUE)
sig<-(o$gradient[,2] > 1) +1
h <- seq(0, 1, 0.01)
center=o$cent[,1];v=o$gradient[,1]
u <- log(center/(1 - center)) * v
plot(h, h, type = "n", xlab = "Hybrid index", ylab = "Ancestry probability", xlim = c(0, 1), ylim = c(0, 1),cex.lab=cl)
for (i in 1:length(v)) {
	phi <- (h^v[i])/(h^v[i] + (1 - h)^v[i] * exp(u[i]))
        lines(h, phi,col=alpha(cso[5],sig[i]*.5),lwd=sig[i]*.4)
}
abline(a = 0, b = 1, lty = 2)
title(main="(F) Clines oligoenic",cex.main=cm)

o<-sum2zero(center=cent_poly[[1]],v=v_poly[[1]],transform=TRUE)
sig<-(o$gradient[,2] > 1) +1
h <- seq(0, 1, 0.01)
center=o$cent[,1];v=o$gradient[,1]
u <- log(center/(1 - center)) * v
plot(h, h, type = "n", xlab = "Hybrid index", ylab = "Ancestry probability", xlim = c(0, 1), ylim = c(0, 1),cex.lab=cl)
for (i in 1:length(v)) {
	phi <- (h^v[i])/(h^v[i] + (1 - h)^v[i] * exp(u[i]))
        lines(h, phi,col=alpha(csp[5],sig[i]*.5),lwd=sig[i]*.4)
}
abline(a = 0, b = 1, lty = 2)
title(main="(G) Clines weak polygenic",cex.main=cm)

o<-sum2zero(center=cent_spoly[[1]],v=v_spoly[[1]],transform=TRUE)
sig<-(o$gradient[,2] > 1) +1
h <- seq(0, 1, 0.01)
center=o$cent[,1];v=o$gradient[,1]
u <- log(center/(1 - center)) * v
plot(h, h, type = "n", xlab = "Hybrid index", ylab = "Ancestry probability", xlim = c(0, 1), ylim = c(0, 1),cex.lab=cl)
for (i in 1:length(v)) {
	phi <- (h^v[i])/(h^v[i] + (1 - h)^v[i] * exp(u[i]))
        lines(h, phi,col=alpha(css[5],sig[i]*.5),lwd=sig[i]*.4)
}
abline(a = 0, b = 1, lty = 2)
title(main="(H) Clines strong polygenic",cex.main=cm)


plot(sdc_neu[,1],sdv_neu[,1],pch=19,xlim=c(.63,1.11),ylim=c(.125,.425),col=csn[5],xlab=expression(paste(SD[C])),ylab=expression(paste(SD[V])),cex.lab=cl,cex.axis=ca)
points(sdc_oligo[,1],sdv_oligo[,1],pch=19,col=cso[5])
points(sdc_poly[,1],sdv_poly[,1],pch=19,col=csp[5])
points(sdc_spoly[,1],sdv_spoly[,1],pch=19,col=css[5])
title(main="(I) Cline SDs",cex.main=cm)

ps<-c(20,19)

cors_o<-rep(NA,10)
sel<-251*c(.25,.75)
d<-(apply(cbind(abs(1:251 - sel[1]),abs(1:251 - sel[2])),1,min))
d<-100 * d/251 ## convert to cM
o<-sum2zero(center=cent_oligo[[1]],v=v_oligo[[1]],transform=TRUE)
sig<-(o$gradient[,2] > 1) +1
plot(d,log(o$gradient[,1]),ylim=c(-2,2),col=alpha(cso[1],.5),pch=ps[sig],xlab="Distance (cM)",ylab="Log gradient (v)",cex.lab=cl,cex.axis=ca)
cors_o[1]<-cor(d,log(o$gradient[,1]))
for(k in 2:10){
	o<-sum2zero(center=cent_oligo[[k]],v=v_oligo[[k]],transform=TRUE)
	sig<-(o$gradient[,2] > 1) +1
	points(d,log(o$gradient[,1]),col=alpha(cso[k],.5),pch=ps[sig])
	cors_o[k]<-cor(d,log(o$gradient[,1]))
}
mtext(paste("r = ",round(mean(cors_o),2),sep=""),1,line=-1.5,adj=.1)
title(main="(J) Oligogenic v",cex.main=cm)

cs<-alpha(c(brewer.pal(9,"Oranges"),"red"),.5)
ps<-c(20,19)

cors_p<-rep(NA,10)
sel<-251*seq(0,1,length.out=50)
d<-rep(NA,251)
for(i in 1:251){
	d[i]<-min(abs(sel-i))
}
d<-100 * d/251 ## convert to cM
o<-sum2zero(center=cent_poly[[1]],v=v_poly[[1]],transform=TRUE)
sig<-(o$gradient[,2] > 1) +1
plot(d,log(o$gradient[,1]),ylim=c(-2,2),col=alpha(csp[1],.5),pch=ps[sig],xlab="Distance (cM)",ylab="Log gradient (v)",cex.lab=cl,cex.axis=ca)
cors_p[1]<-cor(d,log(o$gradient[,1]))
for(k in 2:10){
	o<-sum2zero(center=cent_poly[[k]],v=v_poly[[k]],transform=TRUE)
	sig<-(o$gradient[,2] > 1) +1
	points(d,log(o$gradient[,1]),col=alpha(csp[k],.5),pch=ps[sig])
	cors_p[k]<-cor(d,log(o$gradient[,1]))
}
mtext(paste("r = ",round(mean(cors_p),2),sep=""),1,line=-1.5,adj=.1)
title(main="(K) Weak polygenic v",cex.main=cm)

cors_s<-rep(NA,10)
sel<-251*seq(0,1,length.out=50)
d<-rep(NA,251)
for(i in 1:251){
	d[i]<-min(abs(sel-i))
}
d<-100 * d/251 ## convert to cM
o<-sum2zero(center=cent_spoly[[1]],v=v_spoly[[1]],transform=TRUE)
sig<-(o$gradient[,2] > 1) +1
plot(d,log(o$gradient[,1]),ylim=c(-2,2),col=alpha(css[1],.5),pch=ps[sig],xlab="Distance (cM)",ylab="Log gradient (v)",cex.lab=cl,cex.axis=ca)
cors_s[1]<-cor(d,log(o$gradient[,1]))
for(k in 2:10){
	o<-sum2zero(center=cent_spoly[[k]],v=v_spoly[[k]],transform=TRUE)
	sig<-(o$gradient[,2] > 1) +1
	points(d,log(o$gradient[,1]),col=alpha(css[k],.5),pch=ps[sig])
	cors_s[k]<-cor(d,log(o$gradient[,1]))
}
mtext(paste("r = ",round(mean(cors_s),2),sep=""),1,line=-1.5,adj=.1)
title(main="(L) Strong polygenic v",cex.main=cm)



dev.off()

#####################
cors_n<-rep(NA,10)
sel<-251*seq(0,1,length.out=50)
d<-rep(NA,251)
for(i in 1:251){
	d[i]<-min(abs(sel-i))
}
o<-sum2zero(center=cent_neu[[1]],v=v_neu[[1]],transform=TRUE)
sig<-(o$gradient[,2] > 1) +1
plot(d,log(o$gradient[,1]),ylim=c(-2,2),col=cs[1],pch=ps[sig])
cors_n[1]<-cor(d,log(o$gradient[,1]))
for(k in 2:10){
	o<-sum2zero(center=cent_neu[[k]],v=v_neu[[k]],transform=TRUE)
	sig<-(o$gradient[,2] > 1) +1
	points(d,log(o$gradient[,1]),col=cs[k],pch=ps[sig])
	cors_n[k]<-cor(d,log(o$gradient[,1]))
}

library(xtable)
xtrable(round(cbind(cors_n,cors_o,cors_p,cors_s),3))

sigv_n<-rep(NA,10)
sigc_n<-rep(NA,10)
for(k in 1:10){
        o<-sum2zero(center=cent_neu[[k]],v=v_neu[[k]],transform=TRUE)
        sigv_n[k]<-sum(o$gradient[,2] > 1)
        sigc_n[k]<-sum(o$center[,2] > .5 | o$center[,3] < 0.5)
}


sigv_o<-rep(NA,10)
sigc_o<-rep(NA,10)
for(k in 1:10){
	o<-sum2zero(center=cent_oligo[[k]],v=v_oligo[[k]],transform=TRUE)
	sigv_o[k]<-sum(o$gradient[,2] > 1)
	sigc_o[k]<-sum(o$center[,2] > .5 | o$center[,3] < 0.5) 
}

sigv_p<-rep(NA,10)
sigc_p<-rep(NA,10)
for(k in 1:10){
	o<-sum2zero(center=cent_poly[[k]],v=v_poly[[k]],transform=TRUE)
	sigv_p[k]<-sum(o$gradient[,2] > 1)
	sigc_p[k]<-sum(o$center[,2] > .5 | o$center[,3] < 0.5) 
}

sigv_s<-rep(NA,10)
sigc_s<-rep(NA,10)
for(k in 1:10){
	o<-sum2zero(center=cent_spoly[[k]],v=v_spoly[[k]],transform=TRUE)
	sigv_s[k]<-sum(o$gradient[,2] > 1)
	sigc_s[k]<-sum(o$center[,2] > .5 | o$center[,3] < 0.5) 
}

####################
ps<-c(20,19)

pdf("F_ccenters.pdf",width=4,height=9)
par(mfrow=c(3,1))
par(mar=c(4.5,4.5,2.5,1))

cors_o<-rep(NA,10)
sel<-251*c(.25,.75)
d<-(apply(cbind(abs(1:251 - sel[1]),abs(1:251 - sel[2])),1,min))
d<-100 * d/251 ## convert to cM
o<-sum2zero(center=cent_oligo[[1]],v=v_oligo[[1]],transform=TRUE)
sig<-(o$center[,2] > .5 | o$center[,3] < .5) +1
plot(d,abs(log(o$center[,1]/(1-o$center[,1]))),ylim=c(0,4),col=alpha(cso[1],.5),pch=ps[sig],xlab="Distance (cM)",ylab="Logit center (c)",cex.lab=cl,cex.axis=ca)
cors_o[1]<-cor(d,abs(log(o$center[,1]/(1-o$center[,1]))))
for(k in 2:10){
	o<-sum2zero(center=cent_oligo[[k]],v=v_oligo[[k]],transform=TRUE)
	sig<-(o$center[,2] > .5 | o$center[,3] < .5) +1
	points(d,abs(log(o$center[,1]/(1-o$center[,1]))),col=alpha(cso[k],.5),pch=ps[sig])
	cors_o[k]<-cor(d,abs(log(o$center[,1]/(1-o$center[,1]))))
}
mtext(paste("r = ",round(mean(cors_o),2),sep=""),3,line=-1.5,adj=.1)
title(main="(A) Oligogenic center",cex.main=cm)

cors_p<-rep(NA,10)
sel<-251*seq(0,1,length.out=50)
d<-rep(NA,251)
for(i in 1:251){
	d[i]<-min(abs(sel-i))
}
o<-sum2zero(center=cent_poly[[1]],v=v_poly[[1]],transform=TRUE)
sig<-(o$center[,2] > .5 | o$center[,3] < .5) +1
plot(d,abs(log(o$center[,1]/(1-o$center[,1]))),ylim=c(0,4),col=alpha(csp[1],.5),pch=ps[sig],xlab="Distance (cM)",ylab="Logit center (c)",cex.lab=cl,cex.axis=ca)
cors_p[1]<-cor(d,abs(log(o$center[,1]/(1-o$center[,1]))))
for(k in 2:10){
	o<-sum2zero(center=cent_poly[[k]],v=v_poly[[k]],transform=TRUE)
	sig<-(o$center[,2] > .5 | o$center[,3] < .5) +1
	points(d,abs(log(o$center[,1]/(1-o$center[,1]))),col=alpha(csp[k],.5),pch=ps[sig])
	cors_p[k]<-cor(d,abs(log(o$center[,1]/(1-o$center[,1]))))
}
mtext(paste("r = ",round(mean(cors_p),2),sep=""),3,line=-1.5,adj=.1)
title(main="(B) Weak polygenic center",cex.main=cm)

cors_s<-rep(NA,10)
sel<-251*seq(0,1,length.out=50)
d<-rep(NA,251)
for(i in 1:251){
	d[i]<-min(abs(sel-i))
}
o<-sum2zero(center=cent_spoly[[1]],v=v_spoly[[1]],transform=TRUE)
sig<-(o$center[,2] > .5 | o$center[,3] < .5) +1
plot(d,abs(log(o$center[,1]/(1-o$center[,1]))),ylim=c(0,4),col=alpha(css[1],.5),pch=ps[sig],xlab="Distance (cM)",ylab="Logit center (c)",cex.lab=cl,cex.axis=ca)
cors_s[1]<-cor(d,abs(log(o$center[,1]/(1-o$center[,1]))))
for(k in 2:10){
	o<-sum2zero(center=cent_spoly[[k]],v=v_spoly[[k]],transform=TRUE)
	sig<-(o$center[,2] > .5 | o$center[,3] < .5) +1
	points(d,abs(log(o$center[,1]/(1-o$center[,1]))),col=alpha(css[k],.5),pch=ps[sig])
	cors_s[k]<-cor(d,abs(log(o$center[,1]/(1-o$center[,1]))))
}
mtext(paste("r = ",round(mean(cors_s),2),sep=""),3,line=-1.5,adj=.1)
title(main="(C) Strong polygenic center",cex.main=cm)

dev.off()
##############
cors_n<-rep(NA,10)
sel<-251*seq(0,1,length.out=50)
d<-rep(NA,251)
for(i in 1:251){
	d[i]<-min(abs(sel-i))
}
for(k in 1:10){
	o<-sum2zero(center=cent_neu[[k]],v=v_neu[[k]],transform=TRUE)
	sig<-(o$center[,2] > .5 | o$center[,3] < .5) +1
	cors_n[k]<-cor(d,abs(log(o$center[,1]/(1-o$center[,1]))))
}
library(xtable)
xtable(round(cbind(cors_n,cors_o,cors_p,cors_s),3))
