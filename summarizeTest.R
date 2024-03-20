## combine and summarize the 50 h/Q10 tests
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales)

ff<-list.files(pattern="out_hq")

cmae_h<-matrix(NA,nrow=6,ncol=50)
cmae_q10<-matrix(NA,nrow=6,ncol=50)
ccov90_h<-matrix(NA,nrow=6,ncol=50)
ccov90_q10<-matrix(NA,nrow=6,ncol=50)

for(fi in 1:length(ff)){
        cat(fi,"\n")
        load(ff[fi])
	cmae_h[,k]<-mae_h[,k]
	cmae_q10[,k]<-mae_q10[,k]
	ccov90_h[,k]<-cov90_h[,k]
	ccov90_q10[,k]<-cov90_q10[,k]

}

load("bgc_hq_ex.rda")



## make figure
cl<-1.6;ca<-1.15;cm<-1.35


pdf("F_hq10.pdf",width=7,height=10.5)
par(mfrow=c(3,2))
par(mar=c(4.5,5.5,2.5,1.5))

plot(h_out$hi[order(h_out$hi[,1]),1],xlab="Individual",ylab="Hybrid index",cex.lab=cl,cex.axis=ca,pch=19)
segments(1:100,h_out$hi[order(h_out$hi[,1]),3],1:100,h_out$hi[order(h_out$hi[,1]),4])
points(sdat$q[inds][order(h_out$hi[,1])],pch=19,col=alpha("firebrick",.3))
title(main="(A) Hybrid index estimates",cex.main=cm)


cx<-c("gray",brewer.pal(n=9,"BuPu")[-1],"black")
cs<-rep(cx[1],100)
for(i in 2:10){
	cs[q_out$Q10[,1] > (i-1)/10]<-cx[i]
}

plot(q_out$hi[,1],q_out$Q10[,1],xlab="Hybrid index",ylab="Interpopulation ancestry",ylim=c(0,1),cex.lab=cl,cex.axis=ca,pch=19,col=cs)
lines(c(0,.5),c(0,1))
lines(c(.5,1),c(1,0))
legend(0,1,(0:9)/10,pch=19,col=cx,bty='n')
title(main="(B) Ancestry triangle plot",cex.main=cm)

cs<-rep(c("#fdcc8a","#fc8d59","#d7301f"),2)
boxplot(t(cmae_h),col="white",pch=NA,border=cs,axes=FALSE,cex.lab=cl,ylab="Absolute error")
points(jitter(rep(1:6,50)),as.vector(cmae_h),pch=19,col=rep(alpha(cs,.45),50))
axis(2,cex.axis=ca)
axis(1,at=c(2,5),c("known genotypes","uncertain genotypes"),cex.axis=ca)
legend(.55,.155,c("1.0","0.5","0.1"),col=alpha(cs[1:3],.45),pch=19,bty='n',title="AFD")
title(main="(C) Mean absolute error for H",cex.main=cm)
box()

cs<-rep(c("#b3cde3","#8c96c6","#88419d"),2)
boxplot(t(cmae_q10),col="white",pch=NA,border=cs,axes=FALSE,cex.lab=cl,ylab="Absolute error")
points(jitter(rep(1:6,50)),as.vector(cmae_q10),pch=19,col=rep(alpha(cs,.45),50))
axis(2,cex.axis=ca)
axis(1,at=c(2,5),c("known genotypes","uncertain genotypes"),cex.axis=ca)
legend(.55,.205,c("1.0","0.5","0.1"),col=alpha(cs[1:3],.45),pch=19,bty='n',title="AFD")
title(main="(D) Mean absolute error for Q10",cex.main=cm)
box()

cs<-rep(c("#fdcc8a","#fc8d59","#d7301f"),2)
boxplot(t(ccov90_h),col="white",ylim=c(.5,1),pch=NA,border=cs,axes=FALSE,cex.lab=cl,ylab="CI coverage")
abline(h=.9,lty=2,col="darkgray")
points(jitter(rep(1:6,50)),as.vector(ccov90_h),pch=19,col=rep(alpha(cs,.45),50))
axis(2,cex.axis=ca)
axis(1,at=c(2,5),c("known genotypes","uncertain genotypes"),cex.axis=ca)
title(main="(E) 90% CI coverage for  H",cex.main=cm)
box()

cs<-rep(c("#b3cde3","#8c96c6","#88419d"),2)
boxplot(t(ccov90_q10),col="white",ylim=c(.5,1),pch=NA,border=cs,axes=FALSE,cex.lab=cl,ylab="CI coverage")
abline(h=.9,lty=2,col="darkgray")
points(jitter(rep(1:6,50)),as.vector(ccov90_q10),pch=19,col=rep(alpha(cs,.45),50))
axis(2,cex.axis=ca)
axis(1,at=c(2,5),c("known genotypes","uncertain genotypes"),cex.axis=ca)
title(main="(F) 90% CI coverage for Q10",cex.main=cm)
box()

dev.off()
