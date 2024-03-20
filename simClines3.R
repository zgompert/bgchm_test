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
o_afd_1.0<-vector("list",50)
o_afd_0.5<-vector("list",50)
o_afd_0.1<-vector("list",50)

SDc<-.7
SDv<-.3

dp5<-runif(100,.5,1)
dp1<-runif(100,.1,1)
p05<-0.5-dp5/2
p15<-0.5+dp5/2
p01<-0.5-dp1/2
p11<-0.5+dp1/2
for(k in 1:50){
	hi<-runif(50,0,1)
	o_afd_1.0[[k]]<-simCline(sdc=SDc,sdv=SDv,h=hi,p0=rep(0,100),p1=rep(1,100))
	o_afd_0.5[[k]]<-simCline(sdc=SDc,sdv=SDv,h=hi,p0=p05,p1=p15)
	o_afd_0.1[[k]]<-simCline(sdc=SDc,sdv=SDv,h=hi,p0=p01,p1=p11)
}



save(list=ls(),file="o_range_genotypes.rda")

