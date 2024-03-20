## simulate genomic cline data based on the actual assumed cline model
## this makes it possible to have "true" results that can be compared with
## cline estimates


simCline<-function(sdc=0.1,sdv=0.1,h=NA,p0=rep(0,100),p1=rep(1,100)){
	N<-length(h)
	L<-length(p0)

	## sample values of v and center 
	## for individual clines
	l10v<-rnorm(L,0,sdv)
	lgcent<-rnorm(L,0,sdc)
	
	v<-10^l10v
	cent<-1/(1+exp(-lgcent))

	## compute u
	u<- log(cent/(1-cent)) * v

	## genotype matrix
	g<-matrix(NA,nrow=N,ncol=L)

	## compute genotypes based on cline
	for(i in 1:L){
		## ancestry prob
		phi<-(h^v[i])/(h^v[i] + (1-h)^v[i] * exp(u[i]))
		## sample ancestry
		z<-rbinom(n=N,size=2,prob=phi)
		## sample genotype
		for(j in 1:N){
			g[j,i]<-sum(c(rbinom(n=2-z[j],size=1,prob=p0[i]), rbinom(n=z[j],size=1,prob=p1[i])))
		}
	}
	out<-list(sdc=sdc,sdv=sdv,v=v,cent=cent,u=u,g=g,h=h)
	return(out)
}

## first test, low SDs, 50 hybrids, 100 loci, fixed differences, 50 data sets, all saved as one R workspace
o_SDA_fixed_known<-vector("list",50)
o_SDB_fixed_known<-vector("list",50)
o_SDC_fixed_known<-vector("list",50)

SDc<-c(0.5,.8,1.2)
SDv<-c(.2,.4,.6)

for(k in 1:50){
	hi<-runif(50,0,1)
	o_SDA_fixed_known[[k]]<-simCline(sdc=SDc[1],sdv=SDv[1],h=hi,p0=rep(0,100),p1=rep(1,100))
	o_SDB_fixed_known[[k]]<-simCline(sdc=SDc[2],sdv=SDv[2],h=hi,p0=rep(0,100),p1=rep(1,100))
	o_SDC_fixed_known[[k]]<-simCline(sdc=SDc[3],sdv=SDv[3],h=hi,p0=rep(0,100),p1=rep(1,100))
}

save(list=ls(),file="o_fixed_known.rda")

#ll_gcl<-function(h=NA,v=NA,cent=NA,p0=0.001,p1=0.999,g=NA){

#	u<-v*log(cent/(1-cent))
#	phi<-(h^v)/(h^v+(1-h)^v*exp(u))
#	x<-which(g==2)
#	pp2<-sum(log(phi[x]*p1 + (1-phi[x])*p0) + log(phi[x]*p1 + (1-phi[x])*p0))
#	x<-which(g==0)
#	pp0<-sum(log(phi[x]*(1-p1) + (1-phi[x])*(1-p0)) + log(phi[x]*(1-p1) + (1-phi[x])*(1-p0)))
#	x<-which(g==1)
#	pp1<-sum(log(phi[x]*p1 + (1-phi[x])*p0) + log(phi[x]*(1-p1) + (1-phi[x])*(1-p0)))
#	pp<-pp0+pp1+pp2
#	return(pp)
#	}

#k<-2

#ll<-rep(NA,100)
#v<-seq(.2,2,length.out=100)
#for(i in 1:100){
#ll[i]<-ll_gcl(h=o_SDpt1_fixed_known[[k]]$h,v=v[i],cent=0.5,g=o_SDpt1_fixed_known[[k]]$g[,4])
#}
