## simulate geographic cline data based on the actual assumed cline model
## this makes it possible to have "true" results that can be compared with
## cline estimates


simCline<-function(mu=1,sds=0.5,sdc=0.3,serr=0.1,geo=NA,L=100){

	N<-length(geo)

	## sample slopes and centers  
	## for individual clines
	slope<-rnorm(L,mu,sds)
	cent<-rnorm(L,0,sdc)
	
	## allele frequencye matrix
	p<-matrix(NA,nrow=N,ncol=L)

	## compute allele frequencies based on cline
	for(i in 1:L){
		lphat<-cent[i] + slope[i] * geo
		lp<-rnorm(N,lphat,serr)
		p[,i]<-1/(1+exp(-lp))
	}
	out<-list(mu=mu,sds=sds,sdc=sdc,serr=serr,slope=slope,cent=cent,p=p)
	return(out)
}

## 25 populations from a transect
geo<-sort(rnorm(25,0,1))
geo<-geo-mean(geo)


## low SD
o_plow<-vector("list",100)
for(k in 1:100){
	o_plow[[k]]<-simCline(mu=1.75,sds=.1,sdc=.3,serr=.1,geo=geo,L=100)	
}

## plot to check first one
plot(geo,o_plow[[1]]$p[,1],pch=19,ylim=c(0,1))
for(i in 2:25){
	points(geo,o_plow[[1]]$p[,i])
}


## med SD
o_pmed<-vector("list",100)
for(k in 1:100){
	o_pmed[[k]]<-simCline(mu=1.75,sds=.4,sdc=.3,serr=.1,geo=geo,L=100)	
}

## plot to check first one
plot(geo,o_pmed[[1]]$p[,1],pch=19,ylim=c(0,1))
for(i in 2:25){
	points(geo,o_pmed[[1]]$p[,i])
}

## hi SD
o_phi<-vector("list",100)
for(k in 1:100){
	o_phi[[k]]<-simCline(mu=1.75,sds=.8,sdc=.3,serr=.1,geo=geo,L=100)	
}

## plot to check first one
plot(geo,o_phi[[1]]$p[,1],pch=19,ylim=c(0,1))
for(i in 2:25){
	points(geo,o_phi[[1]]$p[,i])
}
save(list=ls(),file="geosims.rdat")
