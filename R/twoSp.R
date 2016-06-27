###################################################
### Ecovirtual - Two Species Competition Models ###
###################################################
## Lotka-Volterra competition, populational growth and isoclines


##' Lotka-Volterra Competition Model
##' 
##' Simulate the Lotka-Volterra competition model for two populations.
##' 
##' The Lotka-Volterra competition model follows the equations: \itemize{ \item
##' SP1: %%%% NOTE: deqn requires two arguments: \deqn{LaTeX}{ascii}
##' \deqn{\frac{dN_1}{dt}=r_1N_1\left(\frac{K_1-N_1-\alpha N_2}{K_1}\right)}{%
##' dN_1/dt = r_1*N_1*((K_1-N_1-alpha*N_2)/K_1)} \item SP2:
##' \deqn{\frac{dN_2}{dt}=r_2N_2\left(\frac{K_2-N_2-\beta N_1}{K_2}\right)}{%
##' dN_2/dt = r_2*N_2*((K_2-N_2-beta*N_1)/K_2)} }
##' 
##' @param n01 initial population for the superior competitor species.
##' @param n02 initial population for the inferior competitor species.
##' @param tmax maximum simulation time.
##' @param r1 intrinsic growth rate for the superior competitor species.
##' @param r2 intrinsic growth rate for the inferior competitor species.
##' @param k1 carrying capacity for the superior competitor species.
##' @param k2 carrying capacity for the inferior competitor species.
##' @param alfa alfa coefficient.
##' @param beta beta coefficient
##' @return 'compLV' returns a graph of the population size in time, and a
##' graph with the isoclines of the equilibrium for both species. 'compLV' also
##' returns an invisible matrix with the population size of each species in
##' time.
##' @author Alexandre Adalardo de Oliveira \email{ecovirtualpackage@@gmail.com}
##' @seealso \url{http://ecovirtual.ib.usp.br}
##' @references Gotelli, N.J. 2008. A primer of Ecology. 4th ed. Sinauer
##' Associates, 291pp. Hastings, A. 1980. Disturbance, coexistence, history and
##' competition for space. Theoretical Population Biology, 18:363-373. Stevens,
##' M.H.H. 2009. A primer in ecology with R. New York, Springer.
##' @keywords simulation competition
##' @examples
##' 
##' \dontrun{
##' 	compLV(n01=10, n02=10,r1=0.05, r2=0.03, k1=80, k2=50, alfa=1.2, beta=0.5, tmax=200)
##' }
##' 
##' @export compLV
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


##' Metapopulation Competition Model
##' 
##' Simulate a metapopulation dynamics with two competing species, a superior
##' and an inferior competitor. Includes the possibility of habitat destruction
##' in the model.
##' 
##' This function uses the metapopulationa model with internal colonization
##' (see function metaCi in \code{\link{metapopulation}}) for the superior
##' competitor. The inferior competitor can only occupy empty patches and is
##' displaced by the superior competitor if it occupies the same patch.
##' 
##' The argument 'D' inserts the influences of habitat destruction in the
##' model.
##' 
##' The number of patches in the simulated landscape is defined by rw*cl.
##' 
##' @param tmax maximum simulation time.
##' @param cl number of columns for the simulated landscape.
##' @param rw number of rowns for the simulated landscape.
##' @param f01 initial fraction of patches occupied by the superior competitor.
##' @param f02 initial fraction of patches occupied by the inferior competitor.
##' @param i1 colonization coefficient for the superior competitor.
##' @param i2 colonization coefficient for the inferior competitor.
##' @param pe probability of extinction (equal for both species).
##' @param D proportion of habitat destroyed.
##' @param anima logical; if TRUE, show simulation frames.
##' @return 'metaComp' returns a graphic with the simulated landscapes and the
##' results of the proportion of patch occupied by both species.
##' 
##' This function also return an invisible array with the simulation results.
##' @author Alexandre Adalardo de Oliveira and Paulo Inacio Prado
##' \email{ecovirtualpackage@@gmail.com}
##' @seealso \code{\link{comCompete}}, \url{http://ecovirtual.ib.usp.br}
##' @references Stevens, M.H.H. 2009. A primer in ecology with R. New York,
##' Springer.
##' 
##' Gotelli, N.J. 1991. Metapopulation models: the rescue effect, the propagule
##' rain, and the core-satellite hypothesis. The American Naturalist
##' 138:768-776.
##' @keywords metapopulation simulation
##' @examples
##' 
##' \dontrun{
##' metaComp(tmax=100,cl=20,rw=20,f01=0.1,f02=0.4,i1=0.4,i2=0.5,pe=0.25)
##' metaComp(tmax=100,cl=20,rw=20,f01=0.1,f02=0.4,i1=0.4,i2=0.5,pe=0.25, D=0.1)
##' }
##' 
##' @export metaComp
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
if(F2 > 0)abline(h=F2,col="green",lwd=1.5,lty=2)
if(F2 <= 0)abline(h=0, col="green",lwd=1.5,lty=2)
if(D > 0)abline(h=1-D,lty=2)
legend("topright",legend= c("Superior competitor", "Inferior competitor"),col=c("blue","green"),lty=2, bty="n", title="Equilibrium without Habitat destruction")
invisible(pais)
par(op)
}
#metaComp(tmax=100,cl=20,rw=20,f01=0.4,f02=0.4,i1=0.1,i2=0.1,pe=0.05, D=0,anima=TRUE)
#metaComp(tmax=100, cl=100, rw=100, f01=0.1, f02=0.4, i1=0.4, i2=0.5, pe=0.25, D=0)
###################END###############################
