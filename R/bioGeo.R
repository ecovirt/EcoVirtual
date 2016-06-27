##################################################################
### Ecovirtual - Island Biogeography and Neutral Theory Models ###
##################################################################
###########################
### Island Biogeography ###
###########################
## Relationship between extinction and colonization rates for the species richness
##############################
##' Colonization and Extinction balance in the Island Biogeography Equilibrium
##' model
##' 
##' Simulate the balance between extinction and colonization rates given the
##' equilibrium number of species in a island, based on the Island Biogeography
##' Equilibrium model.
##' 
##' The number of species is the balance between extinction and colonization
##' rates at the equilibrium.
##' 
##' @param min between 0-1. The minimum value of the extinction and
##' colonization rates.
##' @param max between 0-1. The maximum value of the extinction and
##' colonization rates.
##' @param cycles number of cycles in the simulation.
##' @param Ext a string representing the extinction rate. This can be 'crs' for
##' an increasing extinction rate, 'fix' for a fixed extinction rate in 0.5, or
##' 'dcr' for a decreasing extinction rate.
##' @param Col a string representing the colonization rate. This can be 'crs'
##' for an increasing colonization rate, 'fix' for a fixed colonization rate in
##' 0.5, or 'dcr' for a decreasing colonization rate.
##' @return 'animaColExt' returns a graph of the extinction and colonization
##' rates varying one or both rates in relation with the number of species of
##' an island.
##' @author Alexandre Adalardo de Oliveira \email{ecovirtualpackage@@gmail.com}
##' @seealso \code{\link{archip}} \code{\link{bioGeoIsl}}
##' \url{http://ecovirtual.ib.com.br}
##' @references Gotelli, N.J. 2008. A primer of Ecology. 4th ed. Sinauer
##' Associates, 291pp.
##' @keywords simulation Island Biogeography
##' @examples
##' 
##' \dontrun{
##' animaColExt(Ext='fix', Col="fix")
##' }
##' 
##' @export animaColExt
animaColExt=function(min=0.01, max=1, cycles=100, Ext="crs", Col="dcr")
{
          a=seq(from=min,to=max,length.out=cycles)
          b=seq(from=max, to=min, length.out=cycles)
          nt=length(a)
          if(Ext=="fix"){ext=rep(0.5,nt)}
          if(Ext=="crs"){ext=a}
          if(Ext=="dcr"){ext=b}
          if(Col=="fix"){col=rep(0.5,nt)}
          if(Col=="crs"){col=a}
          if(Col=="dcr"){col=b}
          dev.new()
          grColExt(E=ext,I=col,P=100, area=1)
}

#animaColExt(Ext='crs', Col="dcr")



## Species colonization and species-area relationship in archipelagoes


##' Species Colonization and Species-Area Relationship in Archipelagos
##' 
##' Simulate species colonization from mainland to islands with different
##' sizes.
##' 
##' The mainland has richness (S) and the evenness can be controled argument
##' abund. The 'abund' argument can be one of these 3 options: \enumerate{
##' \item a vector with the same length of the species richness, meaning the
##' proportion of each species population; \item a single value more than 1
##' representing equal abundance of each species (maximum evenness); \item a
##' single value between 0 and 1, meaning the model of geometric species
##' rank-abundance distribution. The model is: abund*(1-abund)*((1:S)-1), where
##' S is the number of species. }
##' 
##' @param n.isl numeric, number of islands.
##' @param ar.min numeric, area of the smallest island.
##' @param ar.max numeric, area of the biggest island.
##' @param S numeric, number of species (species richness from mainland).
##' @param seed.rain numeric, seed rain. Number of seeds colonizing islands on
##' each time.
##' @param abund numeric, abundance of each species in the seed rain.
##' @param tmax numeric, maximum time for the simulations.
##' @param anima logical; if TRUE, show simulation frames.
##' @return 'archip' returns 3 graphics: \itemize{ \item The species-area
##' relationship: number of species x island area at the end of the simulation.
##' It also returns the coefficients c and z from species-area relationship
##' \eqn{S=cA^z}.
##' 
##' \item Colonization rate curves: colonization (number of species per cycle)
##' x number of species for each island.
##' 
##' \item Passive colonization: number of species x time for each island.
##' 
##' 'archip' also returns an invisible array with the simulation results. }
##' @author Alexandre Adalardo de Oliveira \email{ecovirtualpackage@@gmail.com}
##' @seealso \code{\link{animaColExt}}, \code{\link{bioGeoIsl}},
##' \url{http://ecovirtual.ib.usp.br}
##' @references Gotelli, N.J. 2008. A primer of Ecology. 4th ed. Sinauer
##' Associates, 291pp.
##' @keywords simulation island biogeography Species-area relationship
##' @examples
##' 
##' 
##' \dontrun{
##' archip(n.isl=10,ar.min=10, ar.max=100, S=1000, seed.rain=100, abund=10, tmax=100, anima=TRUE)
##' 
##' archip(n.isl=10,ar.min=10, ar.max=100, S=1000, seed.rain=100, abund=0.5, tmax=100, anima=TRUE)
##' }
##' 
##' 
##' @export archip
archip=function(n.isl,ar.min, ar.max, S, seed.rain, abund, tmax=100, anima=TRUE)
{
	ar.ampl=ar.max -ar.min
	ar.isl= seq(ar.min, ar.max, length.out=n.isl)
	spp=1:S
	cena=array(0, dim=c(S,n.isl, tmax)) 
	local=seq(0, ar.max , len=n.isl*10)
	local[c(1,n.isl*10)]=local[c(2, n.isl*10-1)] 
	locxy<-list()
	sprain<-list()
		if(length(abund)==S) {abund=abund/sum(abund)}else{
		cat("\n abundance vector length is different from the number of species, only the first value considered\n")
			abund=abund[1]
			if(abund==0 | abund>=1){abund=rep(1/S, S);cat("\n maximum eveness\n")}else{
				if(abund<=1 & abund>0){abund = abund*(1-abund)^((1:S)-1); cat("\n geometric species rank-abundance distribution\n")}
				}
		}## modelo tilman geometrico ## todas especies igualmente contribuem para a chuva
	for(i in 2:tmax)
		{
		cena[,,i]<-cena[,,(i-1)]
		chuva=sample(spp, seed.rain, prob=abund, replace=TRUE)
		loc.x=sample(local, seed.rain, replace=TRUE)
		loc.y=sample(local, seed.rain, replace=TRUE)

		locxy[[i]]<-cbind(loc.x, loc.y)
		sprain[[i]]<-chuva
		for(l in 1:n.isl)
			{
			v_x=loc.x<=ar.isl[l]
			v_y=loc.y<=ar.isl[l]
			v_spp=unique(chuva[v_x & v_y])
			cena[v_spp,l,i]<-1
			}
	          }
riq.tempo=t(apply(cena, c(2,3), sum))
	if(i>1 & anima==TRUE)
	{
	animaIsl(riq.tempo, ar.isl, locxy, sprain, S=S)
	}
dev.new()
layout(matrix(data=c(1,2), nrow=2, ncol=1), widths=c(1,1), heights=c(5,1))
old<-par(mar=c(5,4,3,3))#, oma=c(0,0,0,0))			
matplot(riq.tempo, type="l", col=rainbow(n.isl), bty="l", cex.lab=1.2, xlab="Time", ylab="Number of Species", cex.axis=1.2, main="Passive Colonization", cex.main=1.2 )
par(mar=c(2,2,1,2))
image(x=1:n.isl, y=1, matrix(data=1:n.isl, nrow=n.isl,ncol=1),col=rainbow(n.isl), ylab="",xlab="", xaxt="n", yaxt="n", main="Island Size (Area)", cex.main=0.8)
pos.x=1:(n.isl)
area.isl=round(ar.isl^2,0)
axis(1,at=pos.x, area.isl, cex.axis=0.8)
dev.new()
par(mfrow=c(2,1), mar=c(5,5,4,2))
riq.final<-riq.tempo[tmax,]
mod1<-lm(log10(riq.final)~log10(area.isl))
plot(area.isl,riq.final,log="xy",pch=16,col=rainbow(n.isl),bty="l",main=paste("N Islands=",n.isl,"; N spp=",S,"; Time=",tmax), sub=paste("c=",round(10^coef(mod1)[1],2),"; z=",round(coef(mod1)[2],2)),xlab="Island Area",ylab="Number of species",ylim=c(1,max(riq.final)))
abline(mod1, lty=2)
rqz<-apply(cena, c(2,3), sum)
clz<-diff(riq.tempo)
matplot(riq.tempo[2:100,],clz, type="l", col=rainbow(n.isl), bty="l", xlab="Number of species", ylab="Colonization\n (Number of species/cycle)", main="Colonization Rate Curves", cex.main=1.2 )
invisible(cena)
}
#archip(n.isl=10,ar.min=10, ar.max=100, S=1000, seed.rain=100, abund=rep(10,1000), tmax=100, anima=TRUE) #abund 'NORMAL'
#archip(n.isl=10,ar.min=10, ar.max=100, S=1000, seed.rain=100, abund=0.5, tmax=100, anima=TRUE) #abund igual RCMDR, em eveness
##################################################
## Island Biogeography, rates of colonization and extinctions for islands of different sizes and distances to continent.
##################################################
##' Island Biogeographical Model
##' 
##' Simulates island biogeographical models, with rates of colonization and
##' extinction for islands of different sizes and distances to the mainland.
##' 
##' 
##' @param area a vector with the sizes of the island areas. It must have the
##' same length as 'dist'
##' @param dist a vector with the distances of the islands to the mainland. It
##' must have the same length as 'areas'
##' @param P the number of species in the mainland (species richness of the
##' pool).
##' @param weight.A ratio between the area and distance effects. Should be a
##' number between 0 to 1. When the ration is 1 the extinction is only affected
##' by size and colonization only by distance. The default ratio is 0.5,
##' meaning that distance and size equally influence colonization and
##' extinction.
##' @param a.e basal extiction coefficient for area.
##' @param b.e extinction/area coefficient.
##' @param c.i basal colonization coefficient for distance.
##' @param d.i numeric, colonization/distance coefficient.
##' @param e.i basal colonization coefficient for area.
##' @param f.i colonization/area coefficient.
##' @param g.e basal extinction coefficient for distance.
##' @param h.e extinction/distance coefficient.
##' @return 'bioGeoIsl' returns a graph with the rates of colonization and
##' extinction in relation with the species richness for each island.
##' 
##' 'bioGeoIsl' also returns a invisible data frame with the values for area,
##' distance and species richness (S) for each island.
##' @author Alexandre Adalardo de Oliveira \email{ecovirtualpackage@@gmail.com}
##' @seealso \code{\link{animaColExt}} \code{\link{archip}},
##' \url{http://ecovirtual.ib.usp.br}
##' @references Gotelli, N.J. 2008. A primer of Ecology. 4th ed. Sinauer
##' Associates, 291pp.
##' @keywords simulation island biogeography
##' @importFrom stats coef lm 
##' @examples
##' 
##' \dontrun{
##' bioGeoIsl(area=c(5,10,50,80), dist=c(10,100,100,10), P=100, weight.A=.5, a=1, 
##' b=-0.01, c=1, d=-0.01, e=0, f=.01, g=0, h=.01)
##' }
##' 
##' @export bioGeoIsl
bioGeoIsl=function(area, dist , P , weight.A=.5 , a.e=1, b.e=-.01, c.i=1, d.i=-.01,e.i=0, f.i=.01,g.e=0, h.e=.01)
{
dev.new()
nf <- layout(matrix(c(1,2), 2, 1),widths=c(1), heights=c(4,1))
def_par<-par(mar=c(4,7,3,7))
  E=((a.e+b.e*area)*weight.A+(g.e+h.e*dist)*(1-weight.A))
  I=((c.i+d.i*dist)*weight.A+(e.i+f.i*area)*(1-weight.A))
I[I<=0]<-0.001
E[E<=0]<-0.001
S=I*P/(I+E)
T=I*E/(I+E)
nIsl=length(area)
grColExt(E=E , I=I , P=P, area=area)
ex=data.frame(area=area,dist=dist, S=S)
par(mar=c(0,0,0,0))
plot(1:10, 1:10, type="n", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")
points(rep(4,nIsl), 2:(nIsl+1), col=rainbow(nIsl))
text(c (5, 6),c(nIsl+3,nIsl+3), c("Size","Distance"))
text(rep(5,nIsl),2:(nIsl+1), area)
text(rep(6,nIsl),2:(nIsl+1), dist)
segments(4.5,nIsl+2, 6.5, nIsl+2)
segments(4.5, nIsl+3, 4.5, 1)
par(def_par)

invisible(ex)

}

#bioGeoIsl(area=c(5,10,50,80) , dist=c(10,100,100,10), P=100 , weight.A=.5 , a=1, b=-.01, c=1, d=-.01, e=0, f=.01, g=0, h=.01)
######################
### Neutral Theory ###
######################
## Null models - random walk simuation
######################################
##' Random Walk Simulations
##' 
##' Simulates random walk models.
##' 
##' Random walk is a stochastic process of a succession of random steps.
##' 
##' Zero is the absorption surface. When an individual simulation reaches zero,
##' it means that the individual is dead.
##' 
##' See \url{ http://en.wikipedia.org/wiki/Random_walk}.
##' 
##' @param S number of individuals.
##' @param step step size (number of steps on each time)
##' @param tmax maximum simulation time.
##' @param x1max maximum initial distance from absorption surface.
##' @param alleq logical; if TRUE, all initial distance are equal. if FALSE,
##' initial distances for each individual is a sample between 1 and maximum
##' initial distance(x1max).
##' @return 'randWalk' returns a graphic with the simulated trajectories of
##' each individual.
##' 
##' 'randWalk' also returns an invisible matrix with the distance from de edge
##' for each individual on each time.
##' @author Alexandre Adalardo de Oliveira and Paulo Inacio Prado
##' \email{ecovirtualpackage@@gmail.com}
##' @seealso \code{\link{extGame}}, \code{\link{simHub}},
##' \url{http://ecovirtual.ib.usp.br}
##' @references \url{ http://en.wikipedia.org/wiki/Random_walk}
##' @keywords simulation Neutral Theory
##' @examples
##' 
##' \dontrun{
##' randWalk(S=100,step=2,tmax=2e5)
##' randWalk(S=10,step=1,tmax=1e4, x1max=300, alleq=TRUE)
##' }
##' 
##' @export randWalk
randWalk <- function(S=1,step=1,tmax=1e5,x1max=200, alleq=FALSE){
  cont=round(tmax/100)
  sleep=0.01
  if(cont>5e4){sleep=0}
      if(alleq){
                x1=rep(x1max,S)  
               }else{
                    x1 <- sample(1:x1max,S,replace=TRUE)
                    }
  results <- matrix(NA,nrow=1+tmax/cont,ncol=S) 
  results[1,] <- x1
  X <- x1
  for(i in 2:(1+tmax/cont)){
    for(j in 1:cont){
      X[X<=0] <- NA
      X <- X +sample(c(step,-1*step),S,replace=TRUE)
    }
    results[i,] <- X
  }
  results[is.na(results)] <- 0
  time <- seq(0,tmax,by=cont)
  dev.new()
  animaRandWalk(rwData=results, time= time, sleep=sleep)
  invisible(results)
}
## Zero Sum Game


##' Zero-sum game
##' 
##' Simulates a zero-sum game between two competitors with a fixed amount of
##' resource.
##' 
##' A zero-sum game is a mathematical representation of a situation in which a
##' participant's gain (or loss) of resource is exactly balanced by the losses
##' (or gains) of the resource of the other participant(s). If the total gains
##' of the participants are added up, and the total losses are subtracted, they
##' will sum to zero.
##' 
##' @param bet bet size of each competitor on each time.
##' @param total total amount of resource.
##' @param tmax maximum game time.
##' @return 'extGame' returns a graphic with the amount of resource of each
##' competitor on each time.
##' 
##' 'extGame' also returns an invisible vector with the results of the loser on
##' each time.
##' @author Alexandre Adalardo de Oliveira and Paulo Inacio Prado
##' \email{ecovirtualpackage@@gmail.com}
##' @seealso \code{\link{simHub}}, \code{\link{randWalk}},
##' \url{http://ecovirtual.ib.usp.br}
##' @references \url{http://en.wikipedia.org/wiki/Zero-sum_game}
##' @keywords simulation neutral theory
##' @examples
##' 
##' \dontrun{
##' extGame(bet=1,total=20)
##' extGame(bet=1,total=100)
##' }
##' 
##' @export extGame
extGame <- function(bet=1,total=100, tmax=2){
  X <- total/2
  results <- X
  t0=Sys.time()
  while(X>0&X<total){
    X <- X+sample(c(bet,-1*bet),1)
    results <- c(results,X)
    ti=Sys.time()
    time.sim=round(difftime(ti, t0, units="min"), 1)
    if(time.sim>tmax)
    {
    cat("\ntimeout, no losers!\n")
    break()
    }
  }
  dev.new()
  animaGame(results, total)
  invisible(results)
}
#####################################
##' Neutral Theory of Biogeography
##' 
##' Simulates Community Dynamics as in the Neutral Theory of Biogeography
##' 
##' 'simHub1' is the model without immigration.
##' 
##' 'simHub2' incorporates immigration rate from the metacommunity
##' 
##' 'simHub3' incorporates immigration and speciation rates in the
##' metacommunity.
##' @name simHub
##' @aliases simHub simHub1 simHub2 simHub3
##' @param S number of species in the community.
##' @param Sm number of species in the metacommunity.
##' @param j individuals per species in the metacommunity.
##' @param jm individuals per species in the metacommunity.
##' @param D number of deaths per cycle.
##' @param cycles number of cycles in the simulation.
##' @param m colonization/immigration rate.
##' @param nu speciation rate.
##' @param m.weights Mortality weights for each species. Mortality rates of
##' individuals of each species is proportional to species' abundances
##' multiplied by these weights as in Yu et al. (1998). In neutral dynamics all
##' weigths are equal. If length(m.weights)<S then species are divided in
##' groups of (approximately) S/length(m.weights) and species of each group
##' have a value in m.weights. This allows to create groups of species with
##' different mortality probabilities and compare to the neutral dynamics.
##' @param anima logical; if TRUE, the simulation frames of the metacommunity
##' are shown.
##' @return These functions returns a graph with the number of species in time
##' (cycles) in the metacommunity.
##' 
##' They also return an invisible matrix with the results of species richness
##' on each community per time.
##' @author Alexandre Adalardo de Oliveira and Paulo Inacio Prado
##' \email{ecovirtualpackage@@gmail.com}
##' @seealso \code{\link{extGame}}, \code{\link{randWalk}},
##' \url{http://ecovirtual.ib.usp.br}
##' @references Hubbell, S.P. 2001. The Unified Neutral Theory of Biodiversity
##' and Biogeography. Princeton University Pres, 448p.
##' 
##' Yu, D. W., Terborgh, J. W., and Potts, M. D. 1998. Can high tree species
##' richness be explained by Hubbell's null model?. Ecology Letters, 1(3):
##' 193--199.
##' @keywords simulation Neutral Theory
##' @examples
##' 
##' \dontrun{
##' simHub1(S=10,j=10, D=1, cycles=5e3)
##' simHub2(j=2,cycles=2e4,m=0.1)
##' simHub3(Sm=200, jm=20, S= 10, j=100, D=1, cycles=1e4, m=0.01, nu=0.001, anima=TRUE)
##' }
##' 
##################################################
## Hubbell Neutral Model without imigration
##################################################
##' @rdname simHub
simHub1=function(S= 100, j=10, D=1, cycles=1e4, m.weights=1, anima=TRUE){
    if(length(m.weights)>S){
        m.weigths <- m.weights[1:S]
        warning("length of m.weights truncated at end to match number of species (S)")
    }
    if(length(m.weights)<S){
        w1 <- rep(NA,S)
        w1 <- rep(m.weights, each=floor(S/length(m.weights)))
        w1[is.na(w1)] <- m.weights[length(m.weights)]
    }
    if(length(m.weights)==S) w1 <- m.weights
    if(cycles<200){cycles=200; cat("\n Minimum number of cycles: 200\n")}
    stepseq=round(seq(101, cycles+1, len=100))
    step=stepseq[2]- stepseq[1]
    ## Size of the community
    J <- S*j
    ## Matrices to save results
    ind.mat=matrix(nrow=J,ncol=100+length(stepseq)) 
    ## Initial conditions##
    ## All species start with the same number of individuals
    ind.mat[,1] <- rep(1:S,each=j)
    cod.sp <- ind.mat[,1]
#################################################################
###########      incluindo 100 primeiros ciclos          ########
#################################################################
    for(k in 2:100)
        {
            ##Indice dos individuos que morrem
            mortek <- sample(1:J,D, prob=w1[cod.sp])
            ##Indice dos individuos que produzem filhotes para substituir os mortos
            novosk <- sample(1:J,D,replace=TRUE)
            ##Substituindo
            cod.sp[mortek]<-cod.sp[novosk]
            ind.mat[,k] <- cod.sp
        }
###########################
    cont=100
    tempo=0:99
    ##Aqui comecam as simulacoes
    if(!is.null(stepseq))
        {
            for(i in 1:length(stepseq)){
                cont=cont+1
                for(j in 1:step){
                    ##Indice dos individuos que morrem
                    morte <- sample(1:J,D, prob=w1[cod.sp])
                    ##Indice dos individuos que produzem filhotes para substituir os mortos
                    novos <- sample(1:J,D,replace=TRUE)
                    ##Substituindo
                    cod.sp[morte]<-cod.sp[novos]
                }
                ## A cada step ciclos os resultados sao gravados
                ind.mat[,cont] <- cod.sp
            }
            tempo <- c(tempo,stepseq)
        }
    colnames(ind.mat) <- tempo
    if(anima==TRUE)
        {
            dev.new()
            animaHub(dadoHub=ind.mat)
        }
    dev.new()
    op = par(mar=c(6,5,5,2), las=1)
    plot(as.numeric(colnames(ind.mat)),apply(ind.mat,2,rich), xlab="Time (cycles)", ylab="Number of species",ylim=c(0,S), cex.lab=1.2, type="l", col="red", lty=2,  main=paste("Neutral Model Without Colonization", "\n S=",S," J=",J),  sub=paste("Mean extintion=",(S-rich(ind.mat[,ncol(ind.mat)]))/cycles,"sp/cycle"), cex.sub=1, cex.axis=1.2, cex.lab=1.2, lwd=2) 
    invisible(ind.mat)
}
##################################################
## Hubbell Neutral Model with immigration from a Metacommunity
##' @rdname simHub
simHub2=function(S= 100, j=10, D=1, cycles=1e4, m=0.01, anima=TRUE)
{ 
if(cycles<200){cycles=200; cat("\n Minimum number of cycles: 200\n")}
  stepseq=round(seq(101, cycles+1, len=100))
  step=stepseq[2]- stepseq[1]
  ## Tamanho da comunidade
  J <- S*j
  ##Matrizes para guardar os resultados
  ## matriz da especie de cada individuo por ciclo
  ind.mat=matrix(nrow=J,ncol=100+length(stepseq))
  ##CONDICOES INICIAIS##
  ## Todas as especies comecam com o meamo numero de individuos (j=J/S)
  ## Rotulo de especies para cada um dos individuos
  ind.mat[,1] <- rep(1:S,each=j)
  ## Repetindo este rotulo no vetor que sofrera modificacoes
  cod.sp <- ind.mat[,1]
####################################
#### primeiras 100 simulacoes  #####
####################################
    for(k in 2:100)
    {
      ##Indice dos individuos que morrem
      morte <- sample(1:J,D)
      ## Indice dos individuos mortos que serao repostos por migrantes
      defora <- sample(c(TRUE,FALSE),size=D,replace=TRUE,prob=c(m,1-m))
      ##Indice dos individuos que produzem filhotes para substituir os mortos
      novosd <- sample(1:J,D-sum(defora),replace=TRUE)
      novosf <- sample(1:J,sum(defora),replace=TRUE)
      ##Substituindo
      ## Mortos por propagulos de dentro
      if(length(novosd)>0){
        cod.sp[morte[!defora]]<-cod.sp[novosd]
      }
      ## Mortos por propagulos de fora
      if(length(novosf)>0){
        cod.sp[morte[defora]]<-ind.mat[,1][novosf]
      }
    ## A cada step ciclos os resultados sao gravados
    ind.mat[,k] <- cod.sp
  }
#####################################################
	cont=100
  ##Aqui comecam as simulacoes
  for(i in 1:length(stepseq)){
  cont=cont+1
    for(j in 1:step)
    {
      ##Indice dos individuos que morrem
      morte <- sample(1:J,D)
      ## Indice dos individuos mortos que serao repostos por migrantes
      defora <- sample(c(TRUE,FALSE),size=D,replace=TRUE,prob=c(m,1-m))
      ##Indice dos individuos que produzem filhotes para substituir os mortos
      novosd <- sample(1:J,D-sum(defora),replace=TRUE)
      novosf <- sample(1:J,sum(defora),replace=TRUE)
      ##Substituindo
      ## Mortos por propagulos de dentro
      if(length(novosd)>0){
        cod.sp[morte[!defora]]<-cod.sp[novosd]
      }
      ## Mortos por propagulos de fora
      if(length(novosf)>0){
        cod.sp[morte[defora]]<-ind.mat[,1][novosf]
      }
    }
    ## A cada step ciclos os resultados sao gravados
    ind.mat[,cont] <- cod.sp
  }
  tempo <- c(0:99,stepseq)
  colnames(ind.mat) <- tempo
  if(anima==TRUE)
      {
          dev.new()
          animaHub(dadoHub=ind.mat)
      }
  ########### grafico interno ###############
dev.new()
op = par(las=1, mar=c(6,5,4,2))
  plot(tempo,apply(ind.mat,2,rich), xlab="Time (cycles)", ylab="Number of species", type="l",
       main="Neutral Dynamics - Original Community Colonization",sub=paste( "S=",S," J=",J," m=",m,"Mean Extintion rate =",(S-rich(ind.mat[,ncol(ind.mat)]))/cycles,"sp/cycle"),ylim=c(0,S), cex.sub=1, lwd=2, col="blue",cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
par(op)
  invisible(ind.mat)
}
##################################################
## Hubbel Neutral Model with Immigration and speciation from a metacommunity
##################################################
##' @rdname simHub
simHub3=function(Sm=200, jm=20, S= 100, j=10, D=1, cycles=1e4, m=0.01, nu=0.001, anima=TRUE)
{
if(cycles<200){cycles=200; cat("\n Minimum number of cycles: 200\n")}
  stepseq=round(seq(101, cycles+1, len=100))
  step=stepseq[2]- stepseq[1]
  ## Tamanho da metacomunidade
  Jm <- Sm*jm
  ## Tamanho da comunidade
  J <- S*j
  ## Matrizes para guardar os resultados
  ## matriz da especie de cada individuo por ciclo
  ## Na metacomunidade
  meta.mat=matrix(nrow=Jm,ncol=100+length(stepseq)) 
  ## Na comunidade
  ind.mat=matrix(nrow=J,ncol=100+length(stepseq))
  ##CONDICOES INICIAIS##
  ## Todas as especies comecam com o mesmo numero de individuos (j=J/S)
  ## METACOMUNIDADE
  meta.mat[,1] <- rep(1:Sm,each=jm)
  ## Repetindo este rotulo no vetor que sofrera modificacoes
  meta.sp <- meta.mat[,1]
  ## COMUNIDADE
  ## Rotulo de especies para cada um dos individuos
  ind.mat[,1] <- rep(1:S,each=j)
  ## Repetindo este rotulo no vetor que sofrera modificacoes
  cod.sp <- ind.mat[,1]
###################################
#### primeiras 100 simulacoes #####
###################################
    for(k in 2:100)
    {
      ##Indice dos individuos que morrem
      ## Na comunidade
      morte <- sample(1:J,D)
      ## Na metacomunidade
      meta.morte <- sample(1:Jm,D)
      ## Indice dos individuos mortos da comunidade que serao repostos por migrantes 
      defora <- sample(c(TRUE,FALSE),size=D,replace=TRUE,prob=c(m,1-m))
      ## Indice dos individuos mortos da metacomunidade que serao repostos por novas especies 
      meta.defora <- sample(c(TRUE,FALSE),size=D,replace=TRUE,prob=c(nu,1-nu))
      ##Indice dos individuos que produzem filhotes para substituir os mortos da comunidade
      novosd <- sample(1:J,D-sum(defora),replace=TRUE)
      novosf <- sample(1:Jm,sum(defora),replace=TRUE)
      ##Indice dos individuos que produzem filhotes para substituir os mortos da metacomunidade
      meta.novosd <- sample(1:Jm,D-sum(meta.defora),replace=TRUE)
      meta.novosf <- sample(1:Jm,sum(meta.defora),replace=TRUE)
      ##Substituindo
      ## Na metacomunidade ##
      ## Mortos por propagulos de dentro
      if(length(meta.novosd)>0){
        meta.sp[meta.morte[!meta.defora]]<-meta.sp[meta.novosd]
      }
      ## Mortos por novas especies
      if(length(meta.novosf)>0){
        meta.sp[meta.morte[meta.defora]]<-max(meta.sp)+1
      }
      ## Na comunidade ##
      ## Mortos por propagulos de dentro
      if(length(novosd)>0){
        cod.sp[morte[!defora]]<-cod.sp[novosd]
      }
      ## Mortos por propagulos de fora
      if(length(novosf)>0){
        cod.sp[morte[defora]]<-meta.sp[novosf]
      }
      ind.mat[,k] <- cod.sp
      meta.mat[,k] <- meta.sp
    }
#####################################################
	cont=100
  ##Aqui comecam as simulacoes
  for(i in 1:length(stepseq)){
  cont=cont+1
    for(j in 1:step)
    {
      ##Indice dos individuos que morrem
      ## Na comunidade
      morte <- sample(1:J,D)
      ## Na metacomunidade
      meta.morte <- sample(1:Jm,D)
      ## Indice dos individuos mortos da comunidade que serao repostos por migrantes 
      defora <- sample(c(TRUE,FALSE),size=D,replace=TRUE,prob=c(m,1-m))
      ## Indice dos individuos mortos da metacomunidade que serao repostos por novas especies 
      meta.defora <- sample(c(TRUE,FALSE),size=D,replace=TRUE,prob=c(nu,1-nu))
      ##Indice dos individuos que produzem filhotes para substituir os mortos da comunidade
      novosd <- sample(1:J,D-sum(defora),replace=TRUE)
      novosf <- sample(1:Jm,sum(defora),replace=TRUE)
      ##Indice dos individuos que produzem filhotes para substituir os mortos da metacomunidade
      meta.novosd <- sample(1:Jm,D-sum(meta.defora),replace=TRUE)
      meta.novosf <- sample(1:Jm,sum(meta.defora),replace=TRUE)
      ##Substituindo
      ## N metacomunidade ##
      ## Mortos por propagulos de dentro
      if(length(meta.novosd)>0){
        meta.sp[meta.morte[!meta.defora]]<-meta.sp[meta.novosd]
      }
      ## Mortos por novas especies
      if(length(meta.novosf)>0){
        meta.sp[meta.morte[meta.defora]]<-max(meta.sp)+1
      }
      ## Na comunidade ##
      ## Mortos por propagulos de dentro
      if(length(novosd)>0){
        cod.sp[morte[!defora]]<-cod.sp[novosd]
      }
      ## Mortos por propagulos de fora
      if(length(novosf)>0){
        cod.sp[morte[defora]]<-meta.sp[novosf]
      }
    }
    ## A cada step ciclos os resultados sao gravados
    ind.mat[,cont] <- cod.sp
    meta.mat[,cont] <- meta.sp
  }
  tempo <- c(0:99,stepseq)
  colnames(ind.mat) <- tempo
  colnames(meta.mat) <- tempo
  resultados <- list(metacomunidade=meta.mat,comunidade=ind.mat)
if(anima==TRUE)
    {
        dev.new()
        animaHub(dadoHub=resultados$comunidade)
    }
  ########### grafico interno ###############
  dev.new()
  mrich<-apply(meta.mat,2,rich)
  crich<-apply(ind.mat,2,rich)
  ymax<-max(c(mrich,crich))
ymax=ymax*1.1
op <- par(las=1, mar=c(6,5,4,2))
  plot(tempo,apply(meta.mat,2,rich), xlab="Time (cycles)", ylab="Number of species", type="l",  main="Neutral Dynamics - Metacommunity Colonization" ,sub=paste( "Jm=",Jm,"; nu=",nu,"; Theta=",2*Jm*nu, "; S=",S,"; J=",J," m=",m, "; Mean Extintion Rate=",(S-rich(ind.mat[,ncol(ind.mat)]))/cycles,"sp/cycle", sep=""), col="blue",  ylim=c(0,ymax), cex.sub=.75, cex.lab=1.2, cex.main=1.2, cex.axis=1.2, lwd=2)
  lines(tempo,apply(ind.mat,2,rich),col="red", lwd=2)
  text(tempo[length(tempo)*.6] ,crich[length(tempo)*.6]*1.1, "Community", col="red")
text(tempo[length(tempo)*.9] ,mrich[length(tempo)*.9]*1.1, "Metacommunity", col="blue")
par(op)
  invisible(resultados)
}
#####################END#############################
