##########################################################
## Metapopulations models functions of EcoVirtual R Package
## Alexandre Adalardo de Oliveira - 17 february 2011
###########################################################
metapop <-function(tf,cl,ln,fi,pc,pe)
{
	paisag=array(0,dim=c(ln,cl,tf))
   nmanchas=cl*ln
	paisag[,,1]=matrix(sample(c(1,0),nmanchas,prob=c(fi,1-fi), replace=TRUE),ln,cl)
	resultado=numeric()
		for(tc in 2:tf)
		{
	       paisag[,,tc][paisag[,,(tc-1)]==1]<-sample(c(0,1),sum(paisag[,,(tc-1)]), replace=TRUE, prob=c(pe,1-pe))
	       paisag[,,tc][paisag[,,(tc-1)]==0]<-sample(c(0,1),cl*ln-sum(paisag[,,(tc-1)]), replace=TRUE, prob=c(1-pc,pc))
	       resultado[tc-1]=sum(paisag[,,tc])/(cl*ln)	
	   }
	meta.anima2(paisag)
	#x11()
	graf.fim(paisag)
	x11()
	F=pc/(pc+pe)
	plot(1:tf,c(fi,resultado),type="l",xlab="Time",ylab="Proportion of occupation",
	ylim=c(0,1),main=paste("Propagulus rain","\n c=",cl," l=",ln," fi=",fi," pc=",pc," pe=",pe),font.lab=2,lwd=2)
	abline(h=F,col=2,lwd=2,lty=2)
	legend("topright", legend=("expected equilibrium"), lty=2, col="red", bty="n")
   invisible(paisag)
}
######################################

meta.spac <-function(tf,cl,ln,fi,pe, pc, canto=FALSE)
{
paisag=array(0,dim=c(ln,cl,tf))
pais.ind=which(paisag[,,1]==0, arr.ind=TRUE)
nmanchas=cl*ln
image(0:ln, 0:cl, paisag[,,1], col=c("white", "green") , breaks=c(0,0.99,5),main="Metapopulation Dynamics", sub=paste("time = 0/", tf, sep=""), xlab="", ylab="")	
	grid(ln,cl)
dist.max=sqrt(cl^2+ln^2)
dif.x=outer(pais.ind[,1],pais.ind[,1],"-")
dif.y=outer(pais.ind[,2],pais.ind[,2],"-")
dist.all=sqrt(dif.x^2 + dif.y^2)
n0ocup=nmanchas*fi
	if(canto==TRUE)
	{
	ocup.l=sqrt(n0ocup)
	rest.fr=(ocup.l-floor(ocup.l))
	rest=(2*rest.fr*ocup.l)-(rest.fr*rest.fr)
	paisag[1:round(ocup.l),1:round(ocup.l),1]=1
		if(rest>0)
		{
		paisag[1:rest,round(ocup.l+1),1]=1
		}
	}
	else
	{
	paisag[,,1]=matrix(sample(c(rep(1,fi*nmanchas), rep(0,round((1-fi)*nmanchas)))), ncol=cl)
	}
resultado=numeric()
	for(tc in 2:tf)
	{
	vazio=which(paisag[,,(tc-1)]==0)
	nvazio=length(vazio)
	ocupa=which(paisag[,,(tc-1)]==1)
	nocupa=length(ocupa)
	d_vz.oc=dist.all[vazio,ocupa]
	dist.media=apply(d_vz.oc,1,mean)
	dist.rel=(dist.max-dist.media)/dist.max
	pc.aj=pc*dist.rel
	paisag[,,tc][paisag[,,(tc-1)]==1]<-sample(c(0,1),nocupa,replace=TRUE,prob=c(pe,1-pe))
   paisag[,,tc][paisag[,,(tc-1)]==0]<-rbinom(nvazio,1,prob=pc.aj)
	resultado[tc-1]=sum(paisag[,,tc])/nmanchas
	}
#F=1-(pe/i)
meta.anima2(paisag)
	#x11()
graf.fim(paisag)
x11()
plot(1:tf,c(fi,resultado),type="l",xlab="Time",ylab="Patch Proportion",
ylim=c(0,1),main="Patch Dynamics",font.lab=2,lwd=2)
#abline(h=F,col=2,lwd=2,lty=2)
return(paisag)
}
################################

meta.inter <-function(tf,cl,ln,fi,i,pe)
{
paisag=array(0,dim=c(ln,cl,tf))
nmanchas=cl*ln
paisag[,,1]=matrix(sample(c(rep(1,fi*nmanchas), rep(0,round((1-fi)*nmanchas)))), ncol=cl)
resultado=numeric()
	for(tc in 2:tf)
	{
	pc=i*sum(paisag[,,tc-1])/(cl*ln)
	paisag[,,tc][paisag[,,(tc-1)]==1]<-sample(c(0,1),sum(paisag[,,(tc-1)]), replace=TRUE,prob=c(pe,1-pe))
	paisag[,,tc][paisag[,,(tc-1)]==0]<-sample(c(0,1),cl*ln-sum(paisag[,,(tc-1)]), replace=TRUE,prob=c(1-pc,pc))
   resultado[tc-1]=sum(paisag[,,tc])/nmanchas
   }
meta.anima2(paisag)
	#x11()
graf.fim(paisag)
x11()
F=1-(pe/i)
plot(1:tf,c(fi,resultado),type="l",xlab="Time",ylab="Proportion of ocupation",
ylim=c(0,1),main=paste("Internal Colonization","\n cols=",cl," rowns=",ln," fi=",fi," i=",i," pe=",pe),font.lab=2,lwd=2)
abline(h=F,col=2,lwd=2,lty=2)
legend("topright", legend=("expected equilibrium"), lty=2, col="red", bty="n")
return(paisag)
}
##########################################
meta.er <-function(tf,cl,ln,fi,pc,e)
{
nmanchas=cl*ln
paisag=array(0,dim=c(ln,cl,tf))
paisag[,,1]=matrix(sample(c(1,0),nmanchas,prob=c(fi,1-fi), replace=TRUE),ln,cl)
resultado=numeric()
res=numeric()
	for(tc in 2:tf)
	{
	pe=e*(1-sum(paisag[,,tc-1])/nmanchas)
	paisag[,,tc][paisag[,,(tc-1)]==1]<-sample(c(0,1),sum(paisag[,,(tc-1)]), replace=TRUE, prob=c(pe,1-pe))
	paisag[,,tc][paisag[,,(tc-1)]==0]<-sample(c(0,1),cl*ln-sum(paisag[,,(tc-1)]), replace=TRUE, prob=c(1-pc,pc))
	resultado[tc-1]=sum(paisag[,,tc])/nmanchas
	res[tc-1]=pe
	}
meta.anima2(paisag)
	#x11()
graf.fim(paisag)
x11()
F=pc/e
plot(1:tf,c(fi,resultado),type="l",xlab="Time",ylab="Proportion of occupancy", ylim=c	(0,1),main=paste("Propagulus Rain and Rescue Effect","\n cols=",cl," rowns=",ln," fi=",fi," pc=",pc," e=",e),font.lab=2,lwd=2) 
abline(h=F,col=2,lwd=2,lty=2) # equilibrio F
points(1:tf,c(e*(1-fi),res),type='l',lwd=2,col="blue") # pe observado
abline(h=e-pc,col="green",lwd=2,lty=2) # pe equilibrio
legend("topright", legend=c("proportion of occupancy", "equilibrium F", "extintion probability(pe)", "pe equilibrium"), lty=c(1,2,1,2), col=c("black","red","blue", "green"), bty="n")
return(paisag)
}
########################################################
meta.cier <-function(tf,cl,ln,fi,i,e)
{
nmanchas=cl*ln
paisag=array(0,dim=c(ln,cl,tf))
paisag[,,1]=sample(c(rep(0,round(nmanchas-fi*nmanchas)),rep(1,round(fi*nmanchas))))
resultado=numeric()
rese=numeric()
resi=numeric()
	for(tc in 2:tf)
	{
	
	pe=e*(1-(sum(paisag[,,tc-1])/nmanchas))
	pc=i*sum(paisag[,,tc-1])/nmanchas
	paisag[,,tc][paisag[,,(tc-1)]==1]<-sample(c(0,1),sum(paisag[,,tc-1]),replace=TRUE,prob=c(pe,1-pe))
	paisag[,,tc][paisag[,,(tc-1)]==0]<-sample(c(0,1),nmanchas-sum(paisag[,,tc-1]),replace=TRUE,prob=c(1-pc,pc))
	resultado[tc-1]=sum(paisag[,,tc])/nmanchas
	rese[tc-1]=pe
	resi[tc-1]=pc
	}
meta.anima2(paisag)
	#x11()
graf.fim(paisag)
x11()
plot(1:tf,c(fi,resultado),type="l",xlab="Time",ylab="Occupancy proportion", ylim=c(0,1),main=paste("Internal colonization","\n cols=",cl," rowns=",ln," fi=",fi," i=",i, "e=",e),font.lab=2,lwd=2)
abline(h=0,lty=2)
points(1:tf,c(e*(1-fi),rese),type='l',lwd=2,col=4,lty=3)
points(1:tf,c(i*fi,resi),type='l',lwd=2,col=6,lty=3)
legend("topright", legend=c("patchs occupance", "colonization", "extintion"), lty=c(1,3,3), col=c(1,6,4), bty="n")
return(paisag)
}
#######################################################
compete=function(n01,n02,tmax,r1,r2,k1,k2,alfa,beta)
{
resulta=matrix(0, ncol=3, nrow=tmax)
resulta[,1]=0:(tmax-1)
resulta[1,c(2,3)]=c(n01,n02)
  for(t in 2:tmax)
  {
   nsp1=resulta[(t-1),2]
   nsp2=resulta[(t-1),3]
   resulta[t,2]=nsp1 + r1*nsp1*((k1-nsp1-alfa*nsp2)/k1)
   resulta[t,3]=nsp2 + r2*nsp2*((k2-nsp2-beta*nsp1)/k2)
     if (resulta[t,2]<1)  
     {
     resulta[t,2]=0
     }
     if (resulta[t,3]<1)  
     {
     resulta[t,3]=0
     }
  }
par(mfrow=c(1,2))
plot(resulta[,1],resulta[,2],ylim=c(0,max(na.omit(resulta[,2:3]))),type="l",lty=1,xlab="time (t)",
ylab="population size (N1,N2)")
lines(resulta[,1],resulta[,3],lty=2)
plot(resulta[,2],resulta[,3],type="l",col="red",xlab="N1",ylab="N2",ylim=c(0,max(c(na.omit(resulta[,3
]),k1/alfa,k2))),xlim=c(0,max(c(na.omit(resulta[,2]),k2/beta,k1))))
segments(0,k1/alfa,k1,0,lty=1)
segments(0,k2,k2/beta,0,lty=2)
return(resulta)
}
########################################################

