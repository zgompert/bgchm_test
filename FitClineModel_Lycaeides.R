## example analyses for Lycaeides
## load libraries
library(data.table)
library(bgchm)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting
library(ggplot2)
library(maps)

## first read in the genotype matrix, these are posteior point estimates, between 0 and 2, but not 
## constrained to integer values... I have estimates from four admixture models that I am combining here

gf<-list.files(pattern="geno_")
#[1] "geno_K2.txt" "geno_K3.txt" "geno_K4.txt" "geno_K5.txt"

gg<-vector("list",length(gf))
for(k in 1:4){
	gg[[k]]<-as.matrix(fread(gf[k],header=FALSE,sep=","))
}

## average
G<-(gg[[1]] + gg[[2]] + gg[[3]] + gg[[4]])/4


## read ids
ids<-read.table("Ids.txt",header=FALSE)

## PCA sanity check
pc<-prcomp(G,center=TRUE,scale=FALSE)
plot(pc$x[,1],pc$x[,2],type='n',xlab="PC1",ylab="PC2")
text(pc$x[,1],pc$x[,2],ids[,1])

## P1 = idas/JH, p2 = melissa, hybrids = dubois/DBS
## use BTB and BCR and BLD for P1
## use VIC and SIN for P2

## simple parental allele frequencies from genotypes
p1x<-which(ids[,1] %in% c("BTB","BCR","BLD"))
length(p1x) ## 166
P1<-apply(G[p1x,],2,mean)/2

p2x<-which(ids[,1] %in% c("VIC","SIN"))
length(p2x) ## 117
P2<-apply(G[p2x,],2,mean)/2

## ancestry informative
dp<-abs(P1-P2)
sum(dp > .3) ## 500 ancestry informative based on this definition, dif > .3 for these pops

anc<-which(dp > .3)
plot(P1[anc],P2[anc])

dbs<-which(ids[,1]=="DBS")
length(dbs)## 115 hybrids

######## plot with map #######
cl<-1.4;ca<-1.1;cm<-1.4

# Data for the six sites
sites <- data.frame(
  site = c("BCR", "BLD", "BTB", "SIN", "VIC", "DBS"),
  lat = c(43.30075, 43.5225, 43.63815, 41.8511, 43.65899, 43.5623),
  lon = c(-110.55297, -109.7161, -110.68197, -107.1131, -111.11138, -109.6991)
)

# Get map data for Wyoming
wyoming_map <- map_data("state", "wyoming")
idaho_map <- map_data("state", "idaho")


# Plot map of Wyoming and Idaho
pdf("map.pdf",width=4,height=4)
ggplot() +
  geom_polygon(data = wyoming_map, aes(x = long, y = lat, group = group), fill = "snow1", color = "black") +
  geom_polygon(data = idaho_map, aes(x = long, y = lat, group = group), fill = "snow1", color = "black") +
  geom_point(data = sites, aes(x = lon, y = lat, color = site), size = 3) +
  scale_color_manual(values = c("#FDBE85", "#FD8D3C", "#E6550D", "forestgreen", "#BDD7E7", "#3182BD")) +
  labs(title = "", x = "Longitude", y = "Latitude") +
  theme_minimal()
dev.off()


pdf("F_mapPca.pdf",width=8,height=4)
par(mfrow=c(1,2))
par(mar=c(4.5,5.5,2.5,1.5))

plot(0,0,xlab="",ylab="",axes=FALSE,pty='n')
title(main="(A) Map",cex.main=cm)

plot(pc$x[,1],pc$x[,2],type='n',xlab="PC1 (9.7%)",ylab="PC2 (1.7%)",cex.lab=cl,cex.axis=ca)
points(pc$x[-c(p1x,p2x,dbs),1],pc$x[-c(p1x,p2x,dbs),2],col=alpha("darkgray",.3),pch=19)
points(pc$x[p1x,1],pc$x[p1x,2],col=alpha("orange",.5),pch=19)
points(pc$x[p2x,1],pc$x[p2x,2],col=alpha("blue",.5),pch=19)
points(pc$x[dbs,1],pc$x[dbs,2],col=alpha("forestgreen",.5),pch=19)
legend(10,18,c("Parent 0","Parent 1","Hybrid zone"),pch=19,col=alpha(c("orange","blue","forestgreen"),.7),bty='n')
title(main="(B) PCA",cex.main=cm)
dev.off()
#############################


## genotype matrix for ancestry informative SNPs for hybrids
Ghyb<-G[dbs,anc]
GhybI<-round(Ghyb,0) ## using integer estimates, keeping it simple

## SNP/sex specific ploidy
ids<-read.table("Ids.txt",header=FALSE)
snps<-read.table("snps.txt",header=FALSE)
ZAnc<-which(snps[anc,1] == 1631)
Ploidy<-matrix(2,nrow=length(dbs),ncol=length(anc))
sex<-ids[dbs,2]
Ploidy[sex=="F",ZAnc]<-1
GhybI[Ploidy==1]<-round(GhybI[Ploidy==1]/2,0) 

GhybA<-GhybI[,-ZAnc]
Aanc<-anc[-ZAnc]

############# estiamte hybrid indexes##############
######### autosomal loci only ####################
out_h<-est_hi(Gx=GhybA,p0=P1[Aanc],p1=P2[Aanc],model="genotype",ploidy="diploid")
cor.test(out_h$hi[,1],pc$x[dbs,1])
#	Pearson's product-moment correlation
#
#data:  out_h$hi[, 1] and pc$x[dbs, 1]
#t = 77.34, df = 113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9865383 0.9935597
#sample estimates:
#      cor
#0.9906859

############# stanadrd cline analysis for autosomes ################
out_gcl_auto<-est_genocl(Gx=GhybA,p0=P1[Aanc],p1=P2[Aanc],H=out_h$hi[,1],model="genotype",ploidy="diploid",hier=TRUE,sd0=2)

############# autosomes with mean estimated, as a test ############
out_gcl_auto_mn<-est_genocl(Gx=GhybA,p0=P1[Aanc],p1=P2[Aanc],H=out_h$hi[,1],model="genotype",ploidy="diploid",hier=TRUE,sd0=2,estMu=TRUE)

############# Z linked, mean set to 0 ##############
GhybZ<-GhybI[,ZAnc]
Zploidy<-Ploidy[,ZAnc]
out_gcl_Z<-est_genocl(Gx=GhybZ,p0=P1[anc][ZAnc],p1=P2[anc][ZAnc],H=out_h$hi[,1],model="genotype",ploidy="mixed",pldat=Zploidy,hier=TRUE,sd0=2)

############## Z linked, with mean estiamted ###########
out_gcl_Z_mn<-est_genocl(Gx=GhybZ,p0=P1[anc][ZAnc],p1=P2[anc][ZAnc],H=out_h$hi[,1],model="genotype",ploidy="mixed",pldat=Zploidy,hier=TRUE,sd0=2,estMu=TRUE)

############## hybrid index, Z linked only ############
out_h_Z<-est_hi(Gx=GhybZ,p0=P1[anc][ZAnc],p1=P2[anc][ZAnc],model="genotype",ploidy="mixed",pldat=Zploidy)

out_gcl_Zhi<-est_genocl(Gx=GhybZ,p0=P1[anc][ZAnc],p1=P2[anc][ZAnc],H=out_h_Z$hi[,1],model="genotype",ploidy="mixed",pldat=Zploidy,hier=TRUE,sd0=2)

############# all together
out_h_all<-est_hi(Gx=GhybI,p0=P1[anc],p1=P2[anc],model="genotype",ploidy="mixed",pldat=Ploidy)
out_Q_all<-est_Q(Gx=GhybI,p0=P1[anc],p1=P2[anc],model="genotype",ploidy="mixed",pldat=Ploidy)
out_gcl_all<-est_genocl(Gx=GhybI,p0=P1[anc],p1=P2[anc],H=out_h_all$hi[,1],model="genotype",ploidy="mixed",pldat=Ploidy,hier=TRUE,sd0=2)


save(list=ls(),file="lycClines.rdat")

## read in lg sizes
lgsz<-read.table("LgSizes.txt",header=TRUE)

asnps<-snps[anc,1]
snpLgSz<-matrix(NA,nrow=length(anc),ncol=2)
snpLgSz<-as.data.frame(snpLgSz)
for(i in 1:length(anc)){
	a<-which(lgsz[,1]==asnps[i])
	snpLgSz[i,]<-lgsz[a,2:3]
}

lg<-which(is.na(snpLgSz[,1])==FALSE)

## snp annotation from Sam's paper
#cat /uufs/chpc.utah.edu/common/home/gompert-group4/data/lycaeides/dovetail_melissa_genome/Annotation/maker/melissa_round2.maker.output/snp_annotation/snp_annotation_outtable.txt | perl -p -i -e 's/[\t ]+/ /g' | cut -d " " -f 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16  > annot.txt

annot<-read.table("annot.txt",header=TRUE)
sannot<-matrix(NA,nrow=length(anc),ncol=16)
sannot<-as.data.frame(sannot)
for(i in 1:length(anc)){
	a<-which(annot[,1]==snps[i,1] & annot[,2]==snps[i,2])
	sannot[i,]<-annot[a,]
}

## within or near (within 1000 bp) of genes and te annotations
gene<-as.numeric(sannot[5]+sannot[6] > 0)
tapply(INDEX=gene,X=log(s2z_all$gradient[,1]),mean)
#          0           1
# 0.01845487 -0.04185517

null_gene<-rep(NA,1000)
for(i in 1:1000){
    a<-sample(1:500,1)
    gnull<-c(gene[a:500],gene[1:a])[-501]
    null_gene[i]<-tapply(INDEX=gnull,X=log(s2z_all$gradient[,1]),mean)[2]
}
mean(null_gene >= -0.04185517)
#[1] 0.792

te<-as.numeric(sannot[9]+sannot[16] > 0)
tapply(INDEX=te,X=log(s2z_all$gradient[,1]),mean)
#          0           1
#-0.01850953  0.20184204

null_te<-rep(NA,1000)
for(i in 1:1000){
    a<-sample(1:500,1)
    tenull<-c(te[a:500],te[1:a])[-501]
    null_te[i]<-tapply(INDEX=tenull,X=log(s2z_all$gradient[,1]),mean)[2]
}
# mean(null_te >= 0.20184204 )
#[1] 0.011

## checking need for sum2zoer
mean(log(out_gcl_all$gradient[,1]))
#[1] 0.2895367
mean(log(out_gcl_all$center[,1]/(1-out_gcl_all$center[,1])))
#[1] 0.1521315
s2z_all<-sum2zero(center=out_gcl_all$center,v=out_gcl_all$gradient)

sum(o$gradient[,2] > 1)
#[1] 40
sum(o$gradient[,3] < 1)
#[1] 60
sum(o$center[,2] > .5)
#[1] 48
sum(o$center[,3] < .5)
#[1] 70
sum(o$gradient[,2] > 1 | o$gradient[,3] < 1 | o$center[,2] > .5 | o$center[,3] < .5)
#[1] 182
## 218 total credible deviations, so some in two ways

## make plots ##
cl<-1.4;ca<-1.1;cm<-1.4
pdf("F_lycaeides.pdf",width=9,height=9)
par(mfrow=c(3,3))

## butterfly
plot(c(0,1),c(0,1),xlab="",ylab="",axes=FALSE,type='n')
title(main="(A) Lycaeides",cex.main=cm)

## hindex
N<-dim(out_h_all$hi)[1]
plot(out_h_all$hi[order(out_h_all$hi[,1]),1],xlab="Individual",ylab="Hybrid index",cex.lab=cl,cex.axis=ca,pch=20,col=alpha("black",.5))
segments(1:N,out_h_all$hi[order(out_h_all$hi[,1]),3],1:N,out_h_all$hi[order(out_h_all$hi[,1]),4],col=alpha("black",.5),lwd=.5)
title(main="(B) Hybrid indexes",cex.main=cm)

## Q plot
plot(out_Q_all$hi[,1],out_Q_all$Q10[,1],xlab="Hybrid index",ylab="Interpopulation ancestry",ylim=c(0,1),cex.lab=cl,cex.axis=ca,pch=20,col=alpha("black",.6),xlim=c(0,1))
segments(out_Q_all$hi[,1],out_Q_all$Q10[,2],out_Q_all$hi[,1],out_Q_all$Q10[,3],col=alpha("black",.5),lwd=.5)
segments(out_Q_all$hi[,2],out_Q_all$Q10[,1],out_Q_all$hi[,3],out_Q_all$Q10[,1],col=alpha("black",.5),lwd=.5)
lines(c(0,.5),c(0,1))
lines(c(.5,1),c(1,0))
title(main="(C) Ancestry triangle plot",cex.main=cm)

## cline plots
cs<-"cadetblue"
o<-s2z_all
sig<-(o$gradient[,2] > 1) +1
h <- seq(0, 1, 0.01)
center=o$cent[,1];v=o$gradient[,1]
u <- log(center/(1 - center)) * v
plot(h, h, type = "n", xlab = "Hybrid index", ylab = "Ancestry probability", xlim = c(0, 1), ylim = c(0, 1),cex.lab=cl)
for (i in seq(1,length(v),length.out=100)) {
        phi <- (h^v[i])/(h^v[i] + (1 - h)^v[i] * exp(u[i]))
        lines(h, phi,col=alpha(cs,sig[i]*.5),lwd=sig[i]*.4)
}
abline(a = 0, b = 1, lty = 2)
title(main="(D) Genomic clines",cex.main=cm)

sig<-1 + (o$gradient[order(snpLgSz[,1]),2] > 1)
ps<-c(21,20)
cs<-c(rep(alpha("gray30",.5),(500-length(ZAnc))),rep(alpha("firebrick",.5),length(ZAnc)))
plot(log(o$gradient[order(snpLgSz[,1]),1]),pch=ps[sig],col=cs,ylim=c(-2.2,2),cex.lab=cl,cex.axis=ca,xlab="SNP number",ylab="Log gradient (v)")
segments(1:500,log(o$gradient[order(snpLgSz[,1]),2]),1:500,log(o$gradient[order(snpLgSz[,1]),3]),col=cs,lwd=.3)
abline(h=0,lty=2)
#legend(0,-1.65,c("Autosome","Z chromosome"),fill=unique(cs),bty='n',cex=ca)
title(main="(E) Cline gradients (v)",cex.main=cm)

sig<-1 + (o$center[order(snpLgSz[,1]),2] > .5 | o$center[order(snpLgSz[,1]),3] < .5)
ps<-c(21,20)
cs<-c(rep(alpha("gray30",.5),(500-length(ZAnc))),rep(alpha("firebrick",.5),length(ZAnc)))
plot(log(o$center[order(snpLgSz[,1]),1]/(1-o$center[order(snpLgSz[,1]),1])),pch=ps[sig],col=cs,ylim=c(-3,5),cex.lab=cl,cex.axis=ca,xlab="SNP number",ylab="Logit center")
segments(1:500,log(o$center[order(snpLgSz[,1]),2]/(1-o$center[order(snpLgSz[,1]),2])),1:500,log(o$center[order(snpLgSz[,1]),3]/(1-o$center[order(snpLgSz[,1]),3])),col=cs,lwd=.3)
abline(h=0,lty=2)
legend(0,4.5,c("Autosome","Z chromosome"),fill=unique(cs),bty='n',cex=ca)
title(main="(F) Cline centers",cex.main=cm)

cs<-c(rep("gray30",18),"firebrick")
mns<-tapply(X=log(s2z_all$gradient[,1]),INDEX=snpLgSz[,1],mean)
sz<-tapply(X=snpLgSz[,2],INDEX=snpLgSz[,1],mean)
cnts<-table(INDEX=snpLgSz[,1])
plot(sz[cnts>5],mns[cnts>5],type='n',xlab="Chromosome length (bp)",ylab="Mean log gradient (v)",cex.lab=cl,cex.axis=ca)
text(sz[cnts>5],mns[cnts>5],names(mns)[cnts>5],col=cs)
oo<-lm(mns[cnts>5] ~ sz[cnts>5])
abline(oo$coefficients)
mtext(expression(paste(r^2," = 0.18",sep="")),line=-2.75,adj=.85,side=1,cex=.9)
mtext("P = 0.07",line=-1.5,adj=.85,side=1,cex=.9)
title(main="(G) Association with size",cex.main=cm)

#### difference in mean, auto vs Z
mnv_z<-extract(out_gcl_Z_mn$gencline_hmc,"muv")
mnv_auto<-extract(out_gcl_auto_mn$gencline_hmc,"muv")
dd<-density(mnv_z[[1]]-mnv_auto[[1]],adj=1.5)
plot(density(mnv_z[[1]]-mnv_auto[[1]],adj=1.5),type='n',xlim=c(-.4,.4),xlab="Difference in mean gradient",ylab="Density",cex.lab=cl,cex.axis=ca,main="")
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$x)),dd$y),border=alpha("green",.6),col=alpha("green",.3))
mean((mnv_z[[1]]-mnv_auto[[1]]) > 0)
abline(v=0,lty=2)
title(main="(H) Difference in means",cex.main=cm)
#[1] 0.99875

### SDs
sdv_z<-extract(out_gcl_Zhi$gencline_hmc,"sv")
sdv_auto<-extract(out_gcl_auto$gencline_hmc,"sv")
boxplot(sdv_auto[[1]],sdv_z[[1]],names=c("Autosomes","Z"),col="white",pch=NA,border=unique(cs),cex.lab=cl,cex.axis=ca,ylab=expression(paste(sigma[v])),xlab="Chromosome")
mean((sdv_z[[1]]-sdv_auto[[1]]) > 0)
points(jitter(rep(c(1,2),each=4000),2),c(sdv_auto[[1]],sdv_z[[1]]),pch=19,col=alpha(rep(unique(cs),each=4000),.05))
title(main="(I) Difference in SDs",cex.main=cm)
#[1] 0.74775
dev.off()

#### here jh and mel defined oppossie, deficient not excess 
sig_jh<-o$center[,2] > .5
sig_mel<-o$center[,3] < .5

zjh_obs<-sum(as.numeric(sig_jh)[ZAnc])
null_jh<-rep(NA,1000)
for(i in 1:1000){
    null_jh[i]<-sum(sample(as.numeric(sig_jh),length(ZAnc),replace=FALSE))zg
}

mean(null_jh >= zjh_obs)
#[1] 0.022
mean(null_jh)
#[1] 16.356
zjh_obs
#[1] 23


zmel_obs<-sum(as.numeric(sig_mel)[ZAnc])
null_mel<-rep(NA,1000)
for(i in 1:1000){
    null_mel[i]<-sum(sample(as.numeric(sig_mel),length(ZAnc),replace=FALSE))
}

mean(null_mel >= zmel_obs)
#[1] 0
mean(null_mel)
#[1] 23.838
zmel_obs
#[1] 37


sig_v<-o$gradient[,2] > 1

zv_obs<-sum(as.numeric(sig_v)[ZAnc])

null_v<-rep(NA,1000)
for(i in 1:1000){
    null_v[i]<-sum(sample(as.numeric(sig_v),length(ZAnc),replace=FALSE))
}
mean(null_v >= zv_obs)
#[1] 0
mean(null_v)
#[1] 13.685
zv_obs
#[1] 39 ## out of 40 = 97.5

sig_int<-o$gradient[,3] < 1

zint_obs<-sum(as.numeric(sig_int)[ZAnc])

null_int<-rep(NA,1000)
for(i in 1:1000){
    null_int[i]<-sum(sample(as.numeric(sig_int),length(ZAnc),replace=FALSE))
}
mean(null_int >= zint_obs)
#[1] 0.999
mean(null_int)
#[1] 40.471
zint_obs
#[1] 12

### difference in jh vs melissa, defintion here flipped to excess
sig_mel<-o$center[,2] > .5
sig_jh<-o$center[,3] < .5

dif_obs<-sum(as.numeric(sig_jh)[ZAnc])-sum(as.numeric(sig_mel)[ZAnc])
null_dif<-rep(NA,1000)
for(i in 1:1000){
    null_dif[i]<-sum(sample(as.numeric(sig_jh),length(ZAnc),replace=FALSE))-sum(sample(as.numeric(sig_mel),length(ZAnc),replace=FALSE))
}

mean(null_dif >= dif_obs)
#[1] 0.098
mean(null_dif)
#[1] 7.27
dif_obs
#[1] 14
