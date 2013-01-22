######## Population Dynamcis for EcoVirtual Package #################### 

##########################################################################
######## crescimento exponencial
crescExp <- function(N0,lambda,r, tmax) 
{
#r=log(lambda)
resulta <- matrix(rep(NA,3*tmax),nrow=tmax)
colnames(resulta)<-c("time","Nd","Nc")
resulta[,1] <- seq(0,tmax-1)
resulta[1,2:3] <- N0
	for (t in 2:tmax) 
	{
	#print(t)
	#print(r)
	resulta[t,2] <- N0*(exp(r*(t-1))) 
	resulta[t,3] <- N0*(lambda^(t-1)) 
	}
x11()
plot(resulta[,1],resulta[,2],type="l",lty=2, main= "Exponential Growth", xlab="time(t)", ylab="population size (N)", col="red")
points(resulta[,1],resulta[,3])
legend("topleft",c("discrete growth","continuous"),lty=c(2,NA_integer_),pch=c(NA_integer_, 1), col=c(2,1), bty="n")
text(x=tmax*0.4, y= resulta[(tmax/2),2], paste("r=", r), col="blue")
invisible(resulta)
}

#crescExp(N0=100,lambda=1.05,r=log(1.05), tmax=20)


########################################################################
############# Exponencial com Estocasticidade Ambiental
estExp <- function(N0,r,varr,tmax) 
{
resulta <- matrix(rep(NA,3*tmax),nrow=tmax)
resulta[,1] <- seq(0,tmax-1)
resulta[1,2:3] <- N0
for (t in 2:tmax) 
	{
	resulta[t,2] <- N0*(exp(r*(t-1)))
	re <- rnorm(1,r,sqrt(varr)) 
	#cat("t =",t-1,"\n")
	#cat("r stocastic = ", re,"\n")
	resulta[t,3] <- resulta[t-1,3]*exp(re) 
		if (resulta[t,3]<1) 
		{
		resulta[t,3] <- 0
		}
	}
plot(resulta[,1],resulta[,2],type="l",main="Exponential Population Growth", lty=2,xlab="Time(t)", ylab="Population size (N)",ylim=c(0,max(resulta[,2:3])))
lines(resulta[,1],resulta[,3], col="red", lty=2)
legend("topleft",c("deterministic","environment stocastic"),lty=2, col=c(1,2), bty="n")
text(x=tmax*0.4, y= resulta[(tmax/2),2], paste("r=", r), col="blue")
text(x=tmax*0.6, y= resulta[(tmax/2),2], paste("var=", varr), col="blue")
invisible(resulta)
}

#estExp(N0=1000,r=0.0488,varr=0.005,tmax=100) 

##########################################################################
### crescimento populacional com taxas de nascimento e morte na pop

estDem=function(N0, b, d, tmax)
{
resulta=matrix(0,ncol=11, nrow=(tmax+1))
colnames(resulta)=c("time", "Nt.dt", "birth.dt", "death.dt", "b.dt", "d.dt", "Nt.st", "nasc.st",
"mort.st", "b.st", "d.st")
resulta[,1]=0:tmax
resulta[1,c(2,7)]=N0
pnasc=b/(b+d)
pmort=d/(b+d)
  for(t in 2:(tmax+1))
  {
  nt.dt=resulta[(t-1),2]
  nt.st=resulta[(t-1),7]
    if(nt.dt>0)
      {
      n.mort.dt=nt.dt*d
      n.nasc.dt=nt.dt*b
      nt1.dt=nt.dt+ n.nasc.dt -n.mort.dt
      resulta[(t-1),3 ]= n.nasc.dt
      resulta[(t-1),4 ]= n.mort.dt
      resulta[(t-1),5 ]=n.nasc.dt/nt.dt
      resulta[(t-1),6 ]=n.mort.dt/nt.dt
        if(nt1.dt>0)
        { 
        resulta[t,2]=nt1.dt
        }
      }
    if(nt.st>0)
    {
    n.nasc.st=rbinom(1,nt.st,pnasc)
    n.mort.st=rbinom(1,nt.st,pmort)
    nt1.st=nt.st + n.nasc.st - n.mort.st
    resulta[(t-1),8]=n.nasc.st
    resulta[(t-1),9]=n.mort.st
    resulta[(t-1),10]=n.nasc.st/nt.st
    resulta[(t-1),11]=n.mort.st/nt.st
        if(nt1.st>0)
        {
        resulta[t,7]=nt1.st
        }
        else
        {
        resulta[t,7]=0
        }
    }
  }
y.max=max(resulta[,c(5,6,10,11)])* 1.2
y.min=min(resulta[,c(5,6,10,11)])* 0.8
nt=dim(resulta)[1]-1

x11()
old=par(mfrow=c(2,1), mar=c(5,4,2,2), cex.lab=0.7, cex.axis=0.7,cex.main=0.8)
plot(resulta[1:nt,"time"] ,resulta[1:nt,"b.dt"],type="l", main="Birth and Death Rates",cex.main=0.7,ylim=c(y.min,y.max),xlab="Time", ylab="Rate", lty=2)
lines(resulta[1:nt,"time"], resulta[1:nt,"d.dt"], col="red",lty=2) 
lines(resulta[1:nt,"time"], resulta[1:nt,"d.st"], col="red")
lines(resulta[1:nt,"time"], resulta[1:nt,"b.st"])
legend("bottomright", legend=c("basal birth", "basal death","stochastic bird", "stochastic dead"),lty=c(2,2,1,1),col=c(1,2,1,2), bty="n", cex=0.7)

plot(resulta[1:nt,"time"] , resulta[1:nt,"Nt.dt"],main="Population Growth", type="l",xlab="Time", ylab="Population size (N)",lty=2)
lines(resulta[1:nt,"time"], resulta[1:nt,"Nt.st"], col="red") 
legend("bottomright", legend=c("Deterministic Model", "Stochastic Model"),lty=c(1,2),col=c(1,2), bty="n", cex=0.7)
par(old)
invisible(resulta)
}

#estDem(N0=100, b=0.55, d=0.5, tmax=50)

#######################################################################
#### crescimento logístico
crescLog=function(N0, r, K, tmax)
{
resulta=matrix(rep(NA,3*tmax),ncol=3)
colnames(resulta)=c("time", "Continuous Model ", "Discrete Model")
resulta[,1]=seq(0,(tmax-1))
resulta[1,2:3]=N0
	for(t in 2:tmax)
	{
	#cat("t= ", t-1,"\n")
	resulta[t,2]=K/(1+((K-N0)/N0)*exp(-r*(t-1)))
	lastN=resulta[t-1,3]
	resulta[t,3]=lastN+r*lastN*(1-lastN/K)
	}
plot(resulta[,1],resulta[,2],type="l", lty=2,ylim=c(min(resulta[,c(2,3)])*0.8, K*1.2), xlab="Time (t)", main="Logistic Population Growth", ylab="Population size (N)", sub=paste("Intrinsic rate (r) =", r, sep=""))	
lines(resulta[,1],resulta[,3], col="red")	
legend("bottomright", colnames(resulta)[2:3],lty=2,col=c(1,2),bty="n")
abline(h=K, lty=3, col="blue")
text(x=2, y=K, "Carying capacity", col="blue",adj=c(0,0))
text(x=tmax*0.5, y= resulta[(tmax/2),2], paste("r=", r),pos=3)
invisible(resulta)
}

#crescLog(N0=10, r=0.05, K=80, tmax=100)

################################################################

discrLog<-function(N0, rd, K, tmax)
  {
  Nt=c(N0,numeric(tmax))
    for(t in 2: (tmax+1))
    {
    Nt[t]=Nt[t-1] + (rd * Nt[t-1]* (1-(Nt[t-1]/K))) 
    }
return(Nt)
}

#discrLog(N0=10, rd=0.05, K=80, tmax=100)

#############################################################
#####
atrBif=function(N0, K, tmax, nrd,maxrd=3)
{
rd.s=seq(1,maxrd,length=nrd)
r1=sapply(rd.s, function(x){discrLog(N0=N0, rd=x, K=100,tmax=200)})
r2=stack(as.data.frame(r1))
names(r2)=c("N", "old.col")
r2$rd=rep(rd.s,each=tmax+1)
r2$time=rep(0:tmax, nrd)
res.bif=subset(r2, time>0.5*tmax)
plot(N~rd, data=res.bif, pch=".", cex=2)
}

#atrBif(N0=50,K=100,tmax=200,nrd=500, maxrd=3)

########################################################################
crescAtr<-function( N0, lambda,varl,rd,K, tmax)
  {
  resulta=matrix(0, ncol=3, nrow=tmax)
  resulta[,1]=1:tmax
  colnames(resulta)=c("time", "exp.estocastic", "logist.discrete")
  resulta[1,c(2,3)]= N0
    for(t in 2: tmax)
    {
    nt.exp=resulta[(t-1),2]
    nt.log=resulta[(t-1),3]
    resulta[t,3]=nt.log + (rd * nt.log* (1-(nt.log/K)))
    lamb.est=rnorm(1,lambda,sd=sqrt(varl))
    resulta[t,2]=nt.exp * lamb.est
      if(resulta[t,2]<1)
      {
      resulta[t,2]=0
      }
    }
  x11()
op <- par(mfrow = c(2, 2)) # 2 x 2 pictures on one plot
plot(resulta[,1],resulta[,2],main="Exponential growth", sub=paste("lamb= ", lambda, "  var= ",
varl), type="l",lty=2,xlab="Time (t)", ylab="Population size (N)")
#tmax <- dim(resulta)[1] 
plot(resulta[1:(tmax-1),2],resulta[2:tmax,2],type="l",col="red",xlab="N[t]",ylab="N[t+1]")
points(resulta[1:(tmax-1),2],resulta[2:tmax,2],pch=20)
plot(resulta[,1],resulta[,3],type="l", main="Logistic Growth", sub=paste(" rd= ",  rd), xlab="Time (t)", ylab="Population size (N)")
plot(resulta[1:(tmax-1),3],resulta[2:tmax,3],type="l",col="red",xlab="N[t]",ylab="N[t+1]")
points(resulta[1:(tmax-1),3],resulta[2:tmax,3],pch=20)
par(op)
invisible(resulta)
}

#crescAtr(N0=610, lambda=1.1,varl=0.05,rd=2.99,K=600, tmax=100)

###################################################
#############################################################
popStr=function(p.sj, p.jj, p.ja, p.aa, fec, ns,nj,na, ln, cl, tmax)
{
ncel=ln*cl
arena=matrix(0,nrow=ln,ncol=cl)
xy.sem=list()
#sem[1]=ns
pais=array(0,dim=c(ln, cl, tmax))
tab.fr=matrix(NA,ncol=4, nrow=tmax)
#image(0:ln,0:cl, matrix(0,nrow=ln,ncol=cl), col="white", xlab="", ylab="")
#grid(ln,cl)
n0=rep(c(0,2,3), c((ncel-nj-na),nj, na))
arena[1:ncel]<-sample(n0)
image(0:ln, 0:cl, arena, main="Structured Population Dynamics", col=c("white", "green", "darkgreen") , breaks=c(-0.1,1.9,2.9,3.9), xlab="", ylab="")
grid(ln,cl)
xsem=sample(seq(0,cl,0.1), ns, replace=TRUE)
ysem=sample(seq(0,ln,0.1), ns, replace=TRUE)
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
		#match(vazio,ind.sem)
	sv=vazio%in% ind.sem
		if(sum(sv)>0)
		{
		sem.vazio=vazio[sv]
		pais[,,tc][sem.vazio]<-sample(c(0,2),sum(sv),replace=TRUE, prob=c((1-p.sj),p.sj))
		}
	if(sum(pais[,,tc])==0 & n.fec==0)
	{
	image(0:ln,0:cl, matrix(0,nrow=ln,ncol=cl), col="white", xlab="", ylab="", add=TRUE)
	grid(ln,cl)
	text(ln/2, cl/2, "EXTINCTION", col="red", cex=4)
	break
	}
#	if(sum[pais[,,tc]==0 & n.fec>0)
#	{
#	image(0:ln,0:cl, matrix(0,nrow=ln,ncol=cl), col="white", xlab="", ylab="", add=TRUE)
#	#text(ln/2, cl/2, "EXTINTION", col="red", cex=4)
#	grid(ln,cl)
#	stop()
#	}
	image(0:ln,0:cl, matrix(0,nrow=ln,ncol=cl), col="white", xlab="", ylab="", add=TRUE)
	image(0:ln, 0:cl, pais[,,tc], col=c("white", "green", "darkgreen") ,breaks=c(0,1,2,3), xlab="", ylab="",  add=TRUE, sub=paste("simulation no. =",tc ))
	grid(ln,cl)
	xsem=sample(seq(0,cl,0.1), n.fec, replace=TRUE)
	ysem=sample(seq(0,ln,0.1), n.fec, replace=TRUE)
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
#t.sim=apply(pais,3, table)
invisible(list(simula=pais, xy=xy.sem))
}

#popStr(p.sj=0.05, p.jj=0.99, p.ja=0, p.aa=1, fec=1.2, ns=100,nj=150,na=50, ln=20, cl=20, tmax=100)
#popStr(0.1,0.4,0.3,0.9,1.2,100,80,20, 20,20,100)


###############################################################
sobrevive=function(p.mort,N0)
{
res=rep(0,N0)
  for(i in 1: N0)
  {
  conta=0
  while(sample(c("m","v"),prob=c(p.mort,1-p.mort),size=1)=="v")
    {
    conta=conta+1
    }
  res[i]=conta  
}
return(res)
}

#sobrevive(0.5,200)

##########################################	
