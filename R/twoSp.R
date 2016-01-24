###################################################
### Ecovirtual - Two Species Competition Models ###
###################################################
## Lotka-Volterra competition, populational growth and isoclines
compLV=function(n01,n02,tmax,r1,r2,k1,k2,alfa,beta)
{
resulta=matrix(0, ncol=3, nrow=tmax, dimnames=list(NULL, c("time", "Nsp1","Nsp2")))
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
dev.new()
old=par(mfrow=c(1,2), mar=c(4,4,2,1))
plot(resulta[,1],resulta[,2],ylim=c(0,max(na.omit(resulta[,2:3]))),type="l",lty=4,xlab="time (t)",ylab="Population size", main="Population Growth", col="blue", lwd=1.5 )
legend("topleft", legend=c("Sp. 1", "Sp. 2"), lty=4, col=c("blue", "green"), bty="n", cex=0.8)
lines(resulta[,1],resulta[,3], col="green", lty=4, lwd=1.5)
plot(resulta[,2],resulta[,3],type="l",col="red",xlab="N1",ylab="N2",ylim=c(0,max(c(na.omit(resulta[,3]),k1/alfa,k2))),xlim=c(0,max(c(na.omit(resulta[,2]),k2/beta,k1))), main="Isoclines")
segments(0,k1/alfa,k1,0,lty=4, lwd=1.5, col="blue")
segments(0,k2,k2/beta,0,lty=4,lwd=1.5, col="green" )
legend("topleft", title="Equilibrium without habitat destruction",legend=c("isocline sp.1 ", "Isocline sp. 2", "Populations trajectory"), lty=c(4,4,1), col=c("blue", "green", "red"), bty="n", cex=0.8)
invisible(resulta)
}

#compLV(n01=10, n02=10,r1=0.05, r2=0.03, k1=80, k2=50, alfa=1.2, beta=0.5, tmax=200)


## Metapopulation competition - patch occupancy between superior and inferior competitors
metaComp<-function(tmax,rw,cl,f01,f02,i1,i2,pe,D=0, anima=TRUE)
{
pais<-array(0, dim=c(rw,cl,tmax))
F1 <- 1-(pe/i1)
F2 <- pe/i1-i1/i2
if(F1<=0) 
    {
        F1=0
        F2 <- 1-(pe/i2)
    }
Nt <- rw*cl
N <- floor(Nt*(1-D))
resultado=matrix(nrow=tmax,ncol=3)
n1 <- floor(f01*N)
n2 <- floor(f02*N)
n0 <- N-(n1+n2)
antes <- sample(rep(c(2,1,0),c(n2,n1,n0)))
nD=rep(-0.5,Nt-N)
pais[,,1]<-c(antes,nD)
resultado[,1] <- 1:tmax
resultado[1,2:3] <- c(sum(antes==1),sum(antes==2))/N
for(tc in 2:tmax)
    {
        depois <- rep(0,N) 
        pi1=i1*sum(antes==1)/Nt
        pi2=i2*sum(antes==2)/Nt
        if(pi1>1)
            {
                pi1= 1
            }
    	if(pi2>1)
            {
                pi2= 1
            }
        depois[antes==1]<-sample(c(0,1),sum(antes==1),replace=TRUE,prob=c(pe,1-pe))
        depois[antes==2]<-sample(c(0,2),sum(antes==2),replace=TRUE,prob=c(pe,1-pe))
        depois[antes==0] <- sample(c(0,2),sum(antes==0),replace=TRUE,prob=c(1-pi2,pi2))
        d1<-sample(c(0,1),sum(antes!=1),replace=TRUE,prob=c(1-pi1,pi1))
        depois[antes!=1][d1==1] <- 1
        resultado[tc,2:3]=c(sum(depois==1),sum(depois==2))/Nt
        pais[,,tc]<-c(depois,nD)
        antes <- depois
    }
dev.new()
if(anima==TRUE)
    {
        animaMetaComp(pais)
    }
dev.new()  
op = par(mar=c(5.5,5,4,2), las=1,mgp= c(3.2,1,0))
plot(1:tmax,resultado[,2],type="l",xlab="Time",ylab="Patch occupancy", ylim=c(0,max(resultado[,c(2,3)]*1.1)),main="Competition and Internal Colonization", sub=paste("cl=",cl,"; rw=",rw,";  f01=",f01,";  f02=", f02,";  i1=",i1,";  i2=",i2,";  pe=",pe,";  D=",D, sep=""),cex.lab=1.3, cex.main=1.4, cex.axis=1.2,lwd=2, cex.sub=1.1, col="blue")
lines(1:tmax,resultado[,3],col="green", lwd=2)
abline(h=F1,col="blue",lwd=1.5,lty=2)
if(F2>0)abline(h=F2,col="green",lwd=1.5,lty=2)
if(F2<0)abline(h=0, col="green",lwd=1.5,lty=2)
if(D>0)abline(h=1-D,lty=2)
legend("topright",legend= c("Superior competitor", "Inferior competitor"),col=c("blue","green"),lty=2, bty="n", title="Equilibrium without Habitat destruction")
invisible(pais)
par(op)
}
#metaComp(tmax=100,cl=20,rw=20,f01=0.4,f02=0.4,i1=0.1,i2=0.1,pe=0.05, D=0,anima=TRUE)
#metaComp(tmax=100, cl=100, rw=100, f01=0.1, f02=0.4, i1=0.4, i2=0.5, pe=0.25, D=0)
###################END###############################
