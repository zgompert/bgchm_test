## cline estimates

## evaluate ability of bgchm to geographic clines 
## considers different cline SDs

library(bgchm)

## data and true parameter values
load("geosims.rdat")

## create objects for mean absolute error and cov
mae_slope<-matrix(NA,nrow=3,ncol=100)
mae_cent<-matrix(NA,nrow=3,ncol=100)

cov90_slope<-matrix(NA,nrow=3,ncol=100)
cov90_cent<-matrix(NA,nrow=3,ncol=100)

## correlations with truth
cor_slope<-matrix(NA,nrow=3,ncol=100)
cor_cent<-matrix(NA,nrow=3,ncol=100)

## hier params
est_SDS<-matrix(NA,nrow=100,ncol=3)
L_SDS<-list(est_SDS,est_SDS,est_SDS)

## get rep number
myargs<-commandArgs(trailingOnly=TRUE)
k<-as.numeric(myargs[1])


out_low<-est_geocl(P=o_plow[[k]]$p,Geo=geo,prec=0.01,y_lb=-2,y_ub=2,hier=TRUE)
out_med<-est_geocl(P=o_pmed[[k]]$p,Geo=geo,prec=0.01,y_lb=-2,y_ub=2,hier=TRUE)
out_hi<-est_geocl(P=o_phi[[k]]$p,Geo=geo,prec=0.01,y_lb=-2,y_ub=2,hier=TRUE)

mae_slope[1,k]<-mean(abs(out_low$slope[,1]-o_plow[[k]]$slope))
mae_cent[1,k]<-mean(abs(out_low$cent[,1]-o_plow[[k]]$cent))
mae_slope[2,k]<-mean(abs(out_med$slope[,1]-o_pmed[[k]]$slope))
mae_cent[2,k]<-mean(abs(out_med$cent[,1]-o_pmed[[k]]$cent))
mae_slope[3,k]<-mean(abs(out_hi$slope[,1]-o_phi[[k]]$slope))
mae_cent[3,k]<-mean(abs(out_hi$cent[,1]-o_phi[[k]]$cent))

cov90_slope[1,k]<-mean(out_low$slope[,3] <= o_plow[[k]]$slope & out_low$slope[,4] >= o_plow[[k]]$slope)
cov90_cent[1,k]<-mean(out_low$cent[,3] <= o_plow[[k]]$cent & out_low$cent[,4] >= o_plow[[k]]$cent)
cov90_slope[2,k]<-mean(out_med$slope[,3] <= o_pmed[[k]]$slope & out_med$slope[,4] >= o_pmed[[k]]$slope)
cov90_cent[2,k]<-mean(out_med$cent[,3] <= o_pmed[[k]]$cent & out_med$cent[,4] >= o_pmed[[k]]$cent)
cov90_slope[3,k]<-mean(out_hi$slope[,3] <= o_phi[[k]]$slope & out_hi$slope[,4] >= o_phi[[k]]$slope)
cov90_cent[3,k]<-mean(out_hi$cent[,3] <= o_phi[[k]]$cent & out_hi$cent[,4] >= o_phi[[k]]$cent)

cor_slope[1,k]<-cor(out_low$slope[,1],o_plow[[k]]$slope)
cor_cent[1,k]<-cor(out_low$cent[,1],o_plow[[k]]$cent)
cor_slope[2,k]<-cor(out_med$slope[,1],o_pmed[[k]]$slope)
cor_cent[2,k]<-cor(out_med$cent[,1],o_pmed[[k]]$cent)
cor_slope[3,k]<-cor(out_hi$slope[,1],o_phi[[k]]$slope)
cor_cent[3,k]<-cor(out_hi$cent[,1],o_phi[[k]]$cent)


L_SDS[[1]][k,]<-out_low$sigma[c(1,3,4)]
L_SDS[[2]][k,]<-out_med$sigma[c(1,3,4)]
L_SDS[[3]][k,]<-out_hi$sigma[c(1,3,4)]

rm(out_low)
rm(out_med)
rm(out_hi)
out<-paste("out_geo",k,".rdat",sep="")
save(list=ls(),file=out)
