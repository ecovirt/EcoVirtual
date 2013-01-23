##################################################################
### Ecovirtual - Island Biogeography and Neutral Theory Models ###
##################################################################

###########################
### Island Biogeography ###
###########################

## Species colonization and species-area relationship in arquipelagoes
arquip=function(nIsl,ar.min, ar.max, Nspp, seed.rain, abund, tmax, anima=TRUE)
{
	ar.ampl=ar.max -ar.min
	ar.isl= seq(ar.min, ar.max, length.out=nIsl)
	#lado.isl=sqrt(ar.isl)
	spp=1:Nspp
	cena=array(0, dim=c(Nspp,nIsl, tmax)) 
	local=seq(0, ar.max , len=nIsl*10)
	local[c(1,nIsl*10)]=local[c(2, nIsl*10-1)] 
	locxy<-list()
	sprain<-list()
		if(length(abund)==Nspp) {abund=abund/sum(abund)}else{
			cat("\n valores de abundância não correspondem ao número de espécie, apenas o primeiro valor foi considerado\n")
			abund=abund[1]
			if(abund==0 | abund>=1){abund=rep(1/Nspp, Nspp);cat("\n distribuição equitativa de abundância\n")}else{
				if(abund<=1 & abund>0){abund = abund*(1-abund)^((1:Nspp)-1); cat("\n modelo de distribuição geometrica de abundância\n")}
				}
		}## modelo tilman geometrico ## todas especies igualmente contribuem para a chuva
	for(i in 2:tmax)
		{
		cena[,,i]<-cena[,,(i-1)]
		chuva=sample(spp, seed.rain, prob=abund, replace=TRUE)
		loc.x=sample(local, seed.rain, replace=TRUE)
		loc.y=sample(local, seed.rain, replace=TRUE)
#		nsemIlh=function(x,y,...){comp.isl<=x & comp.ils}
#		outer(loc.x, loc.y, sum  )
		#v.x=loc.x<ar.isl[l]
		#v.y=loc.y<ar.isl[l]
		#v.spp=unique(chuva[v.x & v.y])
		#cena[v.spp,l,1]<-1
		locxy[[i]]<-cbind(loc.x, loc.y)
		sprain[[i]]<-chuva
		for(l in 1:nIsl)
			{
			v_x=loc.x<=ar.isl[l]
			v_y=loc.y<=ar.isl[l]
			v_spp=unique(chuva[v_x & v_y])
			cena[v_spp,l,i]<-1
			}
#		if(i>1 & anima==TRUE)	
#		animaIsl(cena[,,i],ar.isl, Nspp, loc.x, loc.y, chuva,i)
		}
riq.tempo=t(apply(cena, c(2,3), sum))
          x11()
	if(i>1 & anima==TRUE)
	{
	animaIsl(riq.tempo, ar.isl, locxy, sprain)
	}
x11()
layout(matrix(data=c(1,2), nrow=2, ncol=1), widths=c(1,1), heights=c(5,1))
old<-par(mar=c(5,4,3,3))#, oma=c(0,0,0,0))			
matplot(riq.tempo, type="l", col=rainbow(nIsl), bty="l", cex.lab=1.2, xlab="Time", ylab="Number of Species", cex.axis=1.2, main="Passive Colonization", cex.main=1.2 )
par(mar=c(2,2,1,2))
image(x=1:nIsl, y=1, matrix(data=1:nIsl, nrow=nIsl,ncol=1),col=rainbow(nIsl), ylab="",xlab="", xaxt="n", yaxt="n", main="Island Size (Area)", cex.main=0.8)
pos.x=1:(nIsl)
area.isl=round(ar.isl^2,0)
axis(1,at=pos.x, area.isl, cex.axis=0.8)

x11()
par(mfrow=c(2,1))
riq.final<-riq.tempo[tmax,]
mod1<-lm(log10(riq.final)~log10(area.isl))
plot(area.isl,riq.final,log="xy",pch=16,col=rainbow(nIsl),bty="l",main=paste("Nº Islands=",nIsl,"; Nº spp=",Nspp,"; Time=",tmax), sub=paste("c=",round(10^coef(mod1)[1],2),"; z=",round(coef(mod1)[2],2)),xlab="Island Area",ylab="Number of species",ylim=c(1,max(riq.final)))
abline(mod1, lty=2)
rqz<-apply(cena, c(2,3), sum)
clz<-diff(riq.tempo)
matplot(riq.tempo[2:100,],clz, type="l", col=rainbow(nIsl), bty="l", cex.lab=1.2, xlab="Species Number", ylab="Colonization (species/cicle)", cex.axis=1.2, main="Colonization Rate Curves", cex.main=1.2 )

invisible(cena)
}

#arquip(nIsl=10,ar.min=10, ar.max=100, Nspp=1000, seed.rain=100, abund=10, tmax=100, anima=TRUE)


## Relationship between extinction and colonization rates for the species richness
animaColExt=function(minimo=0.01, maximo=1, ciclos=100, Ext="crs", Col="dcr")
{
a=seq(from=minimo,to=maximo,length.out=ciclos)
b=seq(from=maximo, to=minimo, length.out=ciclos)
nt=length(a)
if(Ext=="fix"){ext=rep(0.5,nt)}
if(Ext=="crs"){ext=a}
if(Ext=="dcr"){ext=b}
if(Col=="fix"){col=rep(0.5,nt)}
if(Col=="crs"){col=a}
if(Col=="dcr"){col=b}
x11()
grColExt(E=ext,I=col,P=100, areas=1)
}

#animaColExt(Ext='crs', Col="dcr")


## Island Biogeography, rates of colonization and extinctions for islands of different sizes and distances to continent.
bioGeoIsl=function(areas , dist , P , peso.A=.5 , a=1, b=-.01, c=1, d=-.01,e=0, f=.01,g=0, h=.01)
{
x11()
nf <- layout(matrix(c(1,2), 2, 1),widths=c(1), heights=c(4,1))
#layout.show(nf)
def_par<-par(mar=c(4,7,3,7))
  E=((a+b*areas)*peso.A+(g+h*dist)*(1-peso.A))
  I=((c+d*dist)*peso.A+(e+f*areas)*(1-peso.A))
#E=((b*areas)*peso.A + (h*dist)*(1-peso.A))
#I=((d*dist)*peso.A+(f*areas)*(1-peso.A))
I[I<=0]<-0.001
E[E<=0]<-0.001
#E=((b*r.areas) * peso.A) + ((h*r.dist)*(1-peso.A))
#I= ((d*r.dist) * (1-peso.A)) + (f*r.areas*peso.A)
S=I*P/(I+E)
T=I*E/(I+E)
nIsl=length(areas)
grColExt(E=E , I=I , P=P, areas=areas)
ex=data.frame(areas=areas,spp=S,dist=dist)
par(mar=c(0,0,0,0))
plot(1:10, 1:10, type="n", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")
points(rep(4,nIsl), 2:(nIsl+1), col=rainbow(nIsl))
text(c (5, 6),c(nIsl+3,nIsl+3), c("Size","Distance"))
text(rep(5,nIsl),2:(nIsl+1), areas)
text(rep(6,nIsl),2:(nIsl+1), dist)
segments(4.5,nIsl+2, 6.5, nIsl+2)
segments(4.5, nIsl+3, 4.5, 1)
par(def_par)
invisible(ex)
}

#bioGeoIsl(areas=c(5,10,50,80) , dist=c(10,100,100,10), P=100 , peso.A=.5 , a=1, b=-.01, c=1, d=-.01, e=0, f=.01, g=0, h=.01)



######################
### Neutral Theory ###
######################

## Null models - random walk simuation
randWalk <- function(n=1,step=1,ciclo=1e5,x1max=200, alleq=FALSE){
  cont=round(ciclo/100)
  sleep=0.01
  if(cont>5e4){sleep=0}
      if(alleq){
                x1=rep(x1max,n)  
               }else{
                    x1 <- sample(1:x1max,n,replace=TRUE)
                    }
  results <- matrix(NA,nrow=1+ciclo/cont,ncol=n) 
  results[1,] <- x1
  X <- x1
  for(i in 2:(1+ciclo/cont)){
    for(j in 1:cont){
      X[X<=0] <- NA
      X <- X +sample(c(step,-1*step),n,replace=TRUE)
    }
    results[i,] <- X
  }
  results[is.na(results)] <- 0
  time <- seq(0,ciclo,by=cont)
  x11()
  animaRandWalk(rwData=results, time= time, sleep=sleep)
  invisible(results)
#  matplot(time,results,type="l", col=rainbow(n),lwd=2, xlab="Steps",  main="Randon Walk",ylab="Distance from the edge")
#  abline(h=0,lwd=4)
}

#randWalk(n=10,step=10,ciclo=1e4)
#randWalk(n=10,step=1,ciclo=1e4)
#randWalk(n=10,step=1,ciclo=1e4, x1max=300, alleq=TRUE)
#rand.walk(n=100,step=2,ciclo=2e5)


## Zero Sum Game
extGame <- function(aposta=1,total=100, tmax=2){
  X <- total/2
  results <- X
  t0=Sys.time()
  while(X>0&X<total){
    X <- X+sample(c(aposta,-1*aposta),1)
    results <- c(results,X)
    ti=Sys.time()
    time.sim=round(difftime(ti, t0, units="min"), 1)
    if(time.sim>tmax)
    {
    cat("\ntimeout, no losers!\n")
    break()
    }
  }
  x11()
  animaGame(results, total)
  invisible(results)
}

#old<-par(mfrow=c(2,2))
#extGame(aposta=1,total=20)
#extGame(aposta=1,total=50)
#extGame(aposta=1,total=100)
#extGame(aposta=1,total=200)
#par(old)


## Hubbell Neutral Model without imigration
simHub1=function(S= 100, j=10, D=1, ciclo=1e4, anima=TRUE)
{
if(ciclo<200){ciclo=200; cat("\n Minimum number of ciclos: 200\n")}
  stepseq=round(seq(101, ciclo+1, len=100))
  step=stepseq[2]- stepseq[1]
  ## Tamanho da comunidade
  J <- S*j
  ##Matrizes para guardar os resultados
  #simHub1(S=10,j=10, D=1, ciclo=2e4, anima=TRUE)# matriz da especie de cada individuo por ciclo
  ind.mat=matrix(nrow=J,ncol=100+length(stepseq)) 
  ##CONDICOES INICIAIS##
  ##Deduzidas de acordo com o modelo de Hubbell:
  ## Todas as especies comecam com o mesmo numero de individuos (j=J/S)
  ind.mat[,1] <- rep(1:S,each=j)
  cod.sp <- ind.mat[,1]
#################################################################
###########      incluindo 100 primeiros ciclos          ########
#################################################################
      for(k in 2:100)
      {
      ##Indice dos individuos que morrem
      mortek <- sample(1:J,D)
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
      morte <- sample(1:J,D)
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
x11()
if(anima==TRUE)
  {
  animaHub(dadoHub=ind.mat)
  }
  x11()
    plot(as.numeric(colnames(ind.mat)),apply(ind.mat,2,rich), xlab="Time (cicles)", ylab="Number of species",ylim=c(0,S), cex.lab=1.2, type="l", col="red", lty=2,  main=paste("Neutral Model Without Colonization", "\n S=",S," J=",J), sub=paste("Mean extintion=",(S-rich(ind.mat[,ncol(ind.mat)]))/ciclo,"sp/cicle"), cex.sub=0.8) 
  invisible(ind.mat)
}

#par(mfrow=c(2,2))
#simHub1(S=10,j=10, D=1, ciclo=5e3, anima=FALSE)
#simHub1(j=5,ciclo=2e4)
#simHub1(j=10,ciclo=2e4)
#simHub1(j=20,ciclo=2e4)
#par(mfrow=c(1,1))


## Hubbell Neutral Model with immigration from a Metacommunity
simHub2=function(S= 100, j=10, D=1, ciclo=1e4, step=1000, m=0.01, anima=TRUE)
{ 
if(ciclo<200){ciclo=200; cat("\n Minimum number of ciclos: 200\n")}
  stepseq=round(seq(101, ciclo+1, len=100))
  step=stepseq[2]- stepseq[1]
  ## Tamanho da comunidade
  J <- S*j
  ##Matrizes para guardar os resultados
  ## matriz da especie de cada individuo por ciclo
  ind.mat=matrix(nrow=J,ncol=100+length(stepseq))
  ##CONDICOES INICIAIS##
  ## Todas as especies comecam com o meamo numero de individuos (j=J/S)
  ## Rotulo de especies para cada um dos inividuos
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
x11()
  if(anima==TRUE)
  {
  animaHub(dadoHub=ind.mat)
  }
  ########### grafico interno ###############
  x11()
  plot(tempo,apply(ind.mat,2,rich), xlab="Time (cicles)", ylab="Number of species", type="l",
       main="Neutral Dynamics - Original Community Colonization",sub=paste( "S=",S," J=",J," m=",m,"Mean Extintion rate =",(S-rich(ind.mat[,ncol(ind.mat)]))/ciclo,"sp/cicle"),ylim=c(0,S), cex.sub=0.7)
  invisible(ind.mat)
}

#simHub2(j=2,ciclo=2e4,step=1e2,m=0.1)


## Hubbel Neutral Model with Immigration and speciation from a metacommunity
simHub3=function(Sm=200, jm=20, S= 100, j=10, D=1, ciclo=1e4, m=0.01, nu=0.001, anima=TRUE)
{
if(ciclo<200){ciclo=200; cat("\n Minimum number of ciclos: 200\n")}
  stepseq=round(seq(101, ciclo+1, len=100))
  step=stepseq[2]- stepseq[1]
  ## Tamanho da metacomunidade
  Jm <- Sm*jm
  #cores<-c("#FFFFFF", topo.colors(Sm)) # preto = #000000
  ## Tamanho da comunidade
  J <- S*j
  ##Matrizes para guardar os resultados
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
  ##COMUNIDADE
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
x11()
if(anima==TRUE)
  {
  animaHub(dadoHub=resultados$comunidade)
  }
  ########### grafico interno ###############
  x11()
  mrich<-apply(meta.mat,2,rich)
  crich<-apply(ind.mat,2,rich)
  ymax<-max(c(mrich,crich))
  ymax=ymax*1.1
  plot(tempo,apply(meta.mat,2,rich), xlab="Time (cicles)", ylab="Number of species", type="l",
       main="Neutra Dynamics - Metacomunity Colonization" ,sub=paste( "Jm=",Jm," nu=",nu," Theta=",2*Jm*nu, "S=",S," J=",J," m=",m, " Mean Extintion Rate=",(S-rich(ind.mat[,ncol(ind.mat)]))/ciclo,"sp/cicle"), col="blue",  ylim=c(0,ymax), cex.sub=0.7)
  lines(tempo,apply(ind.mat,2,rich),col="red")
  text(tempo[length(tempo)*.6] ,crich[length(tempo)*.6]*1.1, "Community", col="red")
  text(tempo[length(tempo)*.9] ,mrich[length(tempo)*.9]*1.1, "Metacommunity", col="blue")
  invisible(resultados)
}

#simHub3(Sm=200, jm=20, S= 10, j=100, D=1, ciclo=1e4, m=0.01, nu=0.001, anima=TRUE)
#simHub3(j=10, ciclo=2e4,m=0.1, anima=FALSE)
#simHub3(j=2, ciclo=2e3,nu=0.00001,m=0.1)
