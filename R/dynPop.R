################################################
### Ecovirtual -  Population Dynamics Models ###
################################################
### Exponencial growth with Environmental Stochasticity
estExp <- function(N0,r,varr,tmax) 
{
## logical tests for initial conditions
#   N0 <- round(as.numeric(tclvalue(noVar)))
        if (is.na(N0) || N0 <= 0) 
        {
            stop("Number of individuos at the simulation start must be a positive integer")
#            return()
        }
#       tmax <- round(as.numeric(tclvalue(tmaxVar)))
        if (is.na(tmax) || tmax <= 0) 
        {
            stop("Number of simulations must be a positive integer")
#            return()
        }
#        varr <- as.numeric(tclvalue(varrVar))
        if (varr < 0)
        {
            stop(message = "r Variance must be zero (no stocatiscit) or positive value")
#            return()
        }
###############################################################################
resulta <- matrix(rep(NA,3*tmax),nrow=tmax)
resulta[,1] <- seq(0,tmax-1)
resulta[1,2:3] <- N0
for (t in 2:tmax) 
	{
	resulta[t,2] <- N0*(exp(r*(t-1)))
	re <- rnorm(1,r,sqrt(varr)) 
	resulta[t,3] <- resulta[t-1,3]*exp(re) 
		if (resulta[t,3]<1) 
		{
		resulta[t,3] <- 0
		}
	}
#x11()
plot(resulta[,1],resulta[,2],type="l",main="Exponential Population Growth", lty=2,xlab="Time(t)", ylab="Population size (N)",ylim=c(0,max(resulta[,2:3])))
lines(resulta[,1],resulta[,3], col="red", lty=2)
legend("topleft",c("deterministic","environment stocastic"),lty=2, col=c(1,2), bty="n")
text(x=tmax*0.4, y= resulta[(tmax/2),2], paste("r=", r), col="blue")
text(x=tmax*0.6, y= resulta[(tmax/2),2], paste("var=", varr), col="blue")
invisible(resulta)
}
#estExp(N0=1000,r=0.0488,varr=0.005,tmax=100) 
##################################################
### Simple Estocastic birth and death dynamics
estDem = function(tmax=10, n=0.2, m=0.2, N0=10, nsim=20, ciclo=5000)
{
ybord=matrix(runif(ciclo*nsim,0,1), ncol=ciclo, nrow=nsim,byrow=TRUE)
bdmat<-matrix(-1,ncol=ciclo, nrow=nsim)
bdmat[ybord <= n/(n+m)]<-1
bdmat<-cbind(N0, bdmat)
nindmat<-t(apply(bdmat, 1, cumsum))
##########
negat<-which(nindmat==0,arr.ind=TRUE)
uniq.row<-unique(negat[,1])
mat.uniq<-match(uniq.row,negat)
ext=0
if(length(negat)>0)
{
	for(i in 1:length(uniq.row))
	{
	nindmat[uniq.row[i],(negat[mat.uniq[i],2]):dim(nindmat)[2]]<-0	
	ext=ext+1
	}
}
ymat<-matrix(runif((ciclo)*nsim,0,1), ncol=ciclo, nrow=nsim,byrow=TRUE)
ymat<-cbind(NA,ymat)
smat<--(log(ymat))/(nindmat*(n+m))
smat[,1]<-0
smat[smat==Inf]<-NA
#ymat<-cbind(0,ymat)
cum.t<-t(apply(smat, 1, cumsum))
nciclo=1
	while(sum(cum.t[,dim(cum.t)[2]]<tmax, na.rm=TRUE)>=1)
	{
	ybord1=matrix(runif(ciclo*nsim,0,1), ncol=ciclo, nrow=nsim,byrow=TRUE)
	bdmat1<-matrix(-1,ncol=ciclo, nrow=nsim)
   bdmat1[ybord1 <= n/(n+m)]<-1
   bdmat1<-cbind(nindmat[,dim(nindmat)[2]], bdmat1)
	nindmat1<-t(apply(bdmat1, 1, cumsum))
	negat1<-which(nindmat1==0,arr.ind=TRUE)
	uniq.row1<-unique(negat1[,1])
	mat.uniq1<-match(uniq.row1,negat1)
		if(length(negat1)>0)
		{
			for(i in 1:length(uniq.row1))
			{
			nindmat1[uniq.row1[i],(negat1[mat.uniq1[i],2]):dim(nindmat1)[2]]<-0	
			ext=ext+1
			}
		}
	ymat1<-matrix(runif((ciclo)*nsim,0,1), ncol=ciclo, nrow=nsim,byrow=TRUE)
	ymat1<-cbind(NA,ymat1)
	smat1<--(log(ymat1))/(nindmat1*(n+m))
	smat1[,1]<-cum.t[,dim(cum.t)[2]]
	smat1[smat1==Inf]<-NA
#ymat<-cbind(0,ymat)
	cum.t1<-t(apply(smat1, 1, cumsum))
	cum.t<-cbind(cum.t, cum.t1[,-1])
	nindmat<-cbind(nindmat,nindmat1[,-1])
	nciclo=nciclo+1
	#cat(paste("iteration #", nciclo*ciclo, "\n"))
	}
#cum.t<-cbind(0,cum.t)
ind.tmax<-cum.t>(tmax)
nindmat[ind.tmax]<-NA
cum.t[ind.tmax]<-NA
p.maxt<-apply(cum.t, 1, which.max)
col.max<-max(p.maxt)
nindmat<-nindmat[,1:(col.max)]
cum.t<-cum.t[,1:col.max]
cores <- rainbow(nsim)
matplot(t(cum.t), t(nindmat),type="l", lty=2, main="Stocastich Simple Birth Death", xlab= "Time",col=cores, ylab="Population Size", cex.lab=1.2, cex.main=1.2, cex.axis=1.2, sub= paste("bird=",n, "dead =",m) )
curve(N0*exp((n-m)*x), add=TRUE, lwd=3)
######## numero de extinções
#p.maxt<-apply(cum.t, 1, which.max)
#n.end<-t(nindmat)[0:(dim(cum.t)[1]-1)*dim(cum.t)[2]+p.maxt]
#ext=sum(n.end==0, na.rm=TRUE)
n.min<-apply(nindmat,1,min, na.rm=TRUE)
n.ext<-sum(n.min==0)
legend("topleft",legend=paste("extinctions =", n.ext, "/", nsim), bty="n")
#points(temp.med,cresc.med, type="l", lwd=4)
invisible(list(time=cum.t, pop.size=nindmat))
}
#estDem(tmax=10, n=0.2, m=0.2, N0=100, nsim=20, ciclo=5000)
###################################
## Logistical Growth
crescLog=function(N0, r, K, tmax)
{
resulta=matrix(rep(NA,3*tmax),ncol=3)
colnames(resulta)=c("time", "Continuous Model ", "Discrete Model")
resulta[,1]=seq(0,(tmax-1))
resulta[1,2:3]=N0
#####################
	if (is.na(N0) || N0 <= 0) 
        {
        stop(message = "Number of individuals at the simulation start must be a positive integer")
        }
	if (is.na(tmax) || tmax <= 0) 
        {
            stop("Number of simulations must be a positive integer")
       }
	if (is.na(K) || K <= 0)
        {
            stop("Carrying Capacity (K) must be a positive integer")
        }
#####################
	for(t in 2:tmax)
	{
	nCont<-K/(1+((K-N0)/N0)*exp(-r*(t-1)))
	ifelse(nCont<0,resulta[t,2]<-0, resulta[t,2]<-nCont)
	lastN=resulta[t-1,3]
	nDisc<-lastN+r*lastN*(1-lastN/K)
	ifelse(nDisc<0,resulta[t,3]<-0, resulta[t,3]<-nDisc)
	}
resulta[,c(2,3)]
plot(resulta[,1],resulta[,2],type="l", lty=2,ylim=c(min(resulta[,c(2,3)])*0.8, max(c(K*1.2, resulta[,2], resulta[,3]))), xlab="Time (t)", main="Logistic Population Growth", ylab="Population size (N)", cex.lab=1.3, cex.axis=1.3, lwd=2, cex.main=1.5)	
lines(resulta[,1],resulta[,3], col="red", lwd=2)	
legend("bottomright", colnames(resulta)[2:3],lty=2,col=c(1,2),bty="n")
abline(h=K, lty=3, col="blue")
text(x=2, y=K, "Carrying capacity", col="blue",adj=c(0,0))
text(x=tmax*0.4, y= resulta[(tmax/2),2], paste("r=", r),pos=3)
invisible(resulta)
}
#crescLog(N0=10, r=0.05, K=80, tmax=100)
################################################
## Populational Model for structured populations
################################################
popStr=function(p.sj, p.jj, p.ja, p.aa, fec, ns, nj, na, rw, cl, tmax)
{
x11()
ncel=rw*cl
arena=matrix(0,nrow=rw,ncol=cl)
xy.sem=list()
pais=array(0,dim=c(rw, cl, tmax))
tab.fr=matrix(NA,ncol=4, nrow=tmax)
n0=rep(c(0,2,3), c((ncel-nj-na),nj, na))
arena[1:ncel]<-sample(n0)
image(0:rw, 0:cl, arena, main="Structured Population Dynamics", col=c("white", "green", "darkgreen") , breaks=c(-0.1,1.9,2.9,3.9), xlab="", ylab="")
grid(rw,cl)
xsem=sample(seq(0,cl,0.1), ns, replace=TRUE)
ysem=sample(seq(0,rw,0.1), ns, replace=TRUE)
ind.sem=floor(ysem)*cl + ceiling(xsem)
points(xsem,ysem, col="red", pch=16)
xy.sem[[1]]=cbind(x=xsem,y=ysem)
t.fr=table(arena)
tab.fr[1,as.numeric(names(t.fr))+1]<-t.fr[]
tab.fr[1,2]<-ns
pais[,,1]<-arena
	for (tc in 2:tmax)
	{
	j.vf=pais[,,(tc-1)]==2
		if(sum(j.vf)>0)
		{
		jovem=which(j.vf)
		pais[,,tc][jovem]<-sample(c(0,2,3),length(jovem),replace=TRUE, prob=c((1-(p.jj+p.ja)),p.jj,p.ja))
		}
	a.vf=pais[,,(tc-1)]==3
		if(sum(a.vf)>0)
		{
		adulto=which(a.vf)
		pais[,,tc][adulto]<-sample(c(0,3),length(adulto),replace=TRUE, prob=c((1-p.aa),p.aa))
		}
	n.fec=round(fec*sum(a.vf))
	vazio=which(pais[,,tc]==0)
	sv=vazio%in% ind.sem
		if(sum(sv)>0)
		{
		sem.vazio=vazio[sv]
		pais[,,tc][sem.vazio]<-sample(c(0,2),sum(sv),replace=TRUE, prob=c((1-p.sj),p.sj))
		}
	if(sum(pais[,,tc])==0 & n.fec==0)
	{
	image(0:rw,0:cl, matrix(0,nrow=rw,ncol=cl), col="white", xlab="", ylab="", add=TRUE)
	grid(rw,cl)
	text(rw/2, cl/2, "EXTINCTION", col="red", cex=4)
	break
	}
	image(0:rw,0:cl, matrix(0,nrow=rw,ncol=cl), col="white", xlab="", ylab="", add=TRUE)
	image(0:rw, 0:cl, pais[,,tc], col=c("white", "green", "darkgreen") ,breaks=c(0,1,2,3), xlab="", ylab="",  add=TRUE, sub=paste("simulation no. =",tc ))
	grid(rw,cl)
	xsem=sample(seq(0,cl,0.1), n.fec, replace=TRUE)
	ysem=sample(seq(0,rw,0.1), n.fec, replace=TRUE)
	xy.sem[[2]]=cbind(x=xsem,y=ysem)
	ind.sem=floor(ysem)*cl + ceiling(xsem)
	points(xsem,ysem, col="red", pch=16)
	Sys.sleep(.1)
	t.fr=table(pais[,,tc])
	tab.fr[tc,as.numeric(names(t.fr))+1]<-t.fr[]
	tab.fr[tc,2]<-n.fec
	}
tab.rel=tab.fr/apply(tab.fr,1,sum)
names(tab.rel)<-c("Empty", "Seed", "Juvenil", "Adult")
x11()
matplot(tab.rel, type="l",col=c("gray", "red", "green", "darkgreen"),lwd=2,main= "Stage Frequency", ylab="Frequency", xlab="Time (t)")
legend("topright",legend=c("Empty", "Seed", "Juvenil", "Adult") ,lty=1:4, col=c("gray", "red", "green", "darkgreen"), bty="n", cex=0.8 )
invisible(list(simula=pais, xy=xy.sem))
}
#popStr(p.sj=0.05, p.jj=0.99, p.ja=0, p.aa=1, fec=1.2, ns=100,nj=150,na=50, rw=20, cl=20, tmax=100)
#popStr(0.1,0.4,0.3,0.9,1.2,100,80,20, 20,20,100)
