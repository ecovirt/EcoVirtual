##########################################################
## Metapopulations models functions of EcoVirtual R Package
## Alexandre Adalardo de Oliveira - 17 february 2011
###########################################################


### chuva de propagulos seed rain 
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
          x11()
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

##ex
metapop(tf=25,cl=20,ln=20,fi=.01,pc=0.2,pe=0.5)


######################################

#################################
### colonização interna ##

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
x11()
meta.anima2(paisag)
	#x11()
graf.fim(paisag)
x11()
F=1-(pe/i)
plot(1:tf,c(fi,resultado),type="l",xlab="Time",ylab="Proportion of ocupation",
ylim=c(0,1),main=paste("Internal Colonization","\n cols=",cl," rows=",ln," fi=",fi," i=",i," pe=",pe),font.lab=2,lwd=2)
abline(h=F,col=2,lwd=2,lty=2)
legend("topright", legend=("expected equilibrium"), lty=2, col="red", bty="n")
return(paisag)
}

##ex
meta.inter(tf=100,cl=10,ln=10,fi=.1,i=1,pe=0.5)

##########################################

## efeito resgate
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
x11()
meta.anima2(paisag)
	#x11()
graf.fim(paisag)
x11()
F=pc/e
plot(1:tf,c(fi,resultado),type="l",xlab="Time",ylab="Proportion of occupancy", ylim=c	(0,1),main=paste("Propagulus Rain and Rescue Effect","\n cols=",cl," rows=",ln," fi=",fi," pc=",pc," e=",e),font.lab=2,lwd=2) 
abline(h=F,col=2,lwd=2,lty=2) # equilibrio F
points(1:tf,c(e*(1-fi),res),type='l',lwd=2,col="blue") # pe observado
abline(h=e-pc,col="green",lwd=2,lty=2) # pe equilibrio
legend("topright", legend=c("proportion of occupancy", "equilibrium F", "extintion probability(pe)", "pe equilibrium"), lty=c(1,2,1,2), col=c("black","red","blue", "green"), bty="n")
return(paisag)
}
########################################################


#efeito resgate com colonização interna

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
x11()
meta.anima2(paisag)
	#x11()
graf.fim(paisag)
x11()
plot(1:tf,c(fi,resultado),type="l",xlab="Time",ylab="Occupancy proportion", ylim=c(0,1),main=paste("Internal colonization","\n cols=",cl," rows=",ln," fi=",fi," i=",i, "e=",e),font.lab=2,lwd=2)
abline(h=0,lty=2)
points(1:tf,c(e*(1-fi),rese),type='l',lwd=2,col=4,lty=3)
points(1:tf,c(i*fi,resi),type='l',lwd=2,col=6,lty=3)
legend("topright", legend=c("patchs occupance", "colonization", "extintion"), lty=c(1,3,3), col=c(1,6,4), bty="n")
return(paisag)
}

#ex
meta.cier(100,10,10,0.5,0.5,0.5)

#######################################################

