###### Biogeography functions - ECOVIRTUAL PACKAGE 
############# Specie Area relationship ###############
######################################################################
################# Arquipelago: ilhas de diferentes tamanhos ###########
### Thu 17 Nov 2011 05:32:27 PM BRST Alexandre Adalardo

arquip=function(nIsl,ar.min, ar.max, Nspp, chuva.total, abund, tmax, anima=TRUE)
{
	#n.ilhas=nIsl
	ar.ampl=ar.max -ar.min
	ar.isl= seq(ar.min, ar.max, length.out=nIsl)
	spp=1:Nspp
	cena=array(0, dim=c(Nspp,nIsl, tmax)) 
	local=1:100
	
	for(i in 1:tmax)
		{
		if(i>1)
			{
			cena[,,i]<-cena[,,(i-1)]
			}
		if(length(abund)==1 | sum(abund)==0){abund=rep(chuva.total/Nspp,Nspp)}
			
		chuva=sample(spp, chuva.total, prob=abund, replace=TRUE)
		loc.x=sample(local, chuva.total, replace=TRUE)
		loc.y=sample(local, chuva.total, replace=TRUE)
		#v.x=loc.x<ar.isl[l]
		#v.y=loc.y<ar.isl[l]
		#v.spp=unique(chuva[v.x & v.y])
		#cena[v.spp,l,1]<-1
		for(l in 1:nIsl)
			{
			v.x=loc.x<ar.isl[l]
			v.y=loc.y<ar.isl[l]
			v.spp=unique(chuva[v.x & v.y])
			cena[v.spp,l,i]<-1
			}
		if(i>1 & anima==TRUE)	
		anima.isl(cena[,,i],nIsl,ar.max,ar.isl, Nspp, loc.x, loc.y, chuva,i)
		}

riq.tempo=t(apply(cena, c(2,3), sum))	
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
plot(area.isl,riq.final,log="xy",pch=16,col=2,bty="l",main=paste("Nº Islands=",nIsl,"; Nº spp=",Nspp,"; Time=",tmax),xlab="Island Area",ylab="Number of species",ylim=c(1,max(riq.final)))
abline(mod1, lty=2)
rqz<-apply(cena, c(2,3), sum)
clz<-diff(riq.tempo)
matplot(riq.tempo[2:100,],clz, type="l", col=rainbow(nIsl), bty="l", cex.lab=1.2, xlab="Species Number", ylab="Colonization (species)", cex.axis=1.2, main="Colonization Rate Curves", cex.main=1.2 )

invisible(cena)
}
#######################################################
#cena<-arquip(nIsl=10,ar.min=10, ar.max=100, Nspp=1000, chuva.total=100, abund=10, tmax=100, anima=FALSE)
########################################################
anima.isl=function(dados, nIsl,ar.max, ar.isl, Nspp, loc.x, loc.y, chuva,i)
{
nspp=apply(dados,2, sum)
layout(matrix(data=c(1,2), nrow=2, ncol=1), widths=c(1,1), heights=c(5,1))
old<-par(mar=c(0,2,3,2), oma=c(0,0,0,0))
plot(0:ar.max, 0:ar.max, usr=c(0,ar.max,0,ar.max), type="n", yaxt="n", xaxt="n", xlab="", ylab="", bty="n", main="Passive Sampling and Area ")
#grid(ar.max,ar.max)   
segments(x0=c(0,0,ar.max,0), y0=c(0,0,0,ar.max), x1=c(0,rep(ar.max,3)), y1=c(ar.max,0,ar.max,ar.max))
segments(x0=c(rep(0,nIsl), ar.isl), y0=c(ar.isl,rep(0,nIsl)), x1=c(ar.isl,ar.isl), y1=c(ar.isl,ar.isl))
col.spp=rainbow(Nspp)
col.func=colorRamp(c("white", "green3"))
col.riq=rgb(col.func(seq(0,1, length.out=Nspp)), max=255)
#par(new=TRUE)
		for(f in nIsl:1)
			{
			vert=ar.isl[f]
			polygon(x=c(0,vert, vert,0),y=c(0,0,vert,vert), col=col.riq[nspp[f]] )
			}
points(loc.x,loc.y, col=col.spp[chuva], pch=16)			
par(mar=c(2,2,1,2))#, oma=c(0,0,0,0))#, mgp=c(1,0,0), omd=c(0,0,0,0))
image(x=1:Nspp, y=1, matrix(data=1:Nspp, nrow=Nspp,ncol=1),col=col.riq, ylab="",xlab=paste("time", i), xaxt="n", yaxt="n", main="Richness")
axis(3, at=c(1.5,Nspp),tick=FALSE, labels=c("0", Nspp))
polygon(x=c(0,Nspp,Nspp,0), y=c(0,0,Nspp,Nspp))
par(old)
Sys.sleep(.1)
}

##########
graColExt=function(E , I , P, areas=areas)
{
	S = I*P/(I+E) ; T = I*E/(I+E)
	nIsl=length(E)
	corIsl=rainbow(nIsl)
	curve(I[1]-I[1]*x/P[1],0,P[1],bty="n",xlab="Number of Species", ylab="Rate",xaxt="n",yaxt="n", font.lab=2,lwd=2,ylim=c(0,1),  main="Island Biogeography", col=corIsl[1])
	curve((E[1]/P[1])*x,0,P,lwd=2,add=TRUE, col=corIsl[1], lty=2) #xlim=c(0,1),
	legend("top", legend=c("Colonization", "Extinction"),  bty="n",lty=c(1,2))
	abline(v=0)
	abline(h=0)
	mtext("St",side=1,at=P,font=2, line=1)
	linhas=seq(0,1.5, length.out=nIsl)
	for(i in 1:nIsl)	
	{
	curve(I[i]-I[i]*x/P,0,P,lwd=2,add=TRUE, col=corIsl[i], lty=1)
	curve((E[i]/P)*x,0,P,lwd=2,add=TRUE, col=corIsl[i], lty=2)
	mtext(paste("S", i, sep=""),side=1,at=S[i], cex=0.8,font=2,col=corIsl[i], line=linhas[i])
	mtext(paste("T", i, sep=""),side=2,at=T[i],cex=0.8,font=2,las=1,col=corIsl[i], line=linhas[i])
	points(S[i],T[i],col=corIsl[i],pch=16,cex=1)
	if(length(unique(areas))>1)
		{
		siz.ar=2 +(areas/max(areas))
		points(S[i],T[i],col=corIsl[i],cex=siz.ar[i])
		}
	segments(S[i],T[i],S[i],0,lty=3,col=corIsl[i])
	segments(S[i],T[i],0,T[i],lty=3,col=corIsl[i])
	Sys.sleep(0.5)
	}	
#	mtext("I",side=2,at=I,font=2,las=1, line=2)
#	mtext("E",side=4,at=E,font=2,las=1)
}
# testando... graColExt(E = .5 , I = .5 , P = 100)
####################################
#######################################
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
	for(i in 1:nt)
	{
		graColExt(E=ext[i],I=col[i],P=100, areas=1)
		Sys.sleep(.01)
	}
}

#animaColExt(Ext='crs', Col="dcr")
####################################
###### Mon 21 Nov 2011 12:58:57 PM BRST Alexandre Adalardo
bioGeoIsl=function(areas , dist , P , peso.A=.5 , a=1, b=-.01, c=1, d=-.01,e=0, f=.01,g=0, h=.01)
{
x11()
nf <- layout(matrix(c(1,2), 2, 1),widths=c(1), heights=c(4,1))
#layout.show(nf)
def.par<-par(mar=c(4,7,3,7))
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
graColExt(E=E , I=I , P=P, areas=areas)
ex=data.frame(areas=areas,spp=S,dist=dist)
par(mar=c(0,0,0,0))
plot(1:10, 1:10, type="n", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")
points(rep(4,nIsl), 2:(nIsl+1), col=rainbow(nIsl))
text(c (5, 6),c(nIsl+3,nIsl+3), c("Size","Distance"))
text(rep(5,nIsl),2:(nIsl+1), areas)
text(rep(6,nIsl),2:(nIsl+1), dist)
segments(4.5,nIsl+2, 6.5, nIsl+2)
segments(4.5, nIsl+3, 4.5, 1)
par(def.par)
return(ex)
}

######################################
###teste
#bioGeoIsl(areas=c(5,10,50,80) , dist=c(10,100,100,10), P=100 , peso.A=.5 , a=1, b=-.01, c=1, d=-.01, e=0, f=.01, g=0, h=.01)
###########################

spp.area=function(c , z){
	curve(expr = c*x^z , from=1, to=10^10, xlab="Area", ylim=c(1,500),
	 ylab="Species number", font.lab=2, lwd=2, col=2, main=paste("c = ",c,"; z = ",z))
	curve(expr = c*x^z , from=1, to=10^10, xlab="Area", ylim=c(1,500),
	 ylab="Species number", font.lab=2, lwd=2, col=2, main=paste("c = ",c,"; z = ",z), log="xy")
	}
#par(mfrow=c(2,2))
#spp.area(c = 1.5 , z = .25)
#spp.area(c = 2.1 , z = .25)

iRain=function(Nspp , chuva , abund , tempo){
	spp=paste("sp.",1:Nspp)
	ilha=numeric()
	riq=numeric()
	for(i in 1:tempo){
		ilha=union(unique(sample(spp,chuva,replace=TRUE,prob=abund/sum(abund))),ilha)
		riq[i]=length(ilha)
		}
	plot(0:tempo,c(0,riq),type="l",lwd=2,bty='n',xlab="time",ylab="number of species",
	 font.lab=2,col=2,ylim=c(0,Nspp),main=c("Island richness",paste("(Nspp=",Nspp," ; rain=",chuva,")")))
	abline(h=Nspp,lty=3)
	return(riq)
	}

# faca um teste:
#iRain(Nspp=10, chuva=10, abund=c(10,10,10,10,10,10,10,10,10,10), tempo=10)

###### tamanho da ilha x spp
iCol=function(areas, Nspp, chuva.total, abund, tempo){
	n.ilhas=length(areas)
	spp=paste("sp.",1:Nspp)
	ilhas=paste("Island",1:n.ilhas)
	chuva=round(chuva.total*areas/sum(areas))
	riq=numeric()
	par(mfrow=c(ceiling(sqrt(n.ilhas)),ceiling(sqrt(n.ilhas))))
	for(i in 1:n.ilhas){
		riq[i]=iRain(Nspp,chuva[i],abund,tempo)[tempo]
		}
	names(riq)=ilhas
	mod=lm(log10(riq)~log10(areas))
	x11()
	plot(areas,riq,log="xy",font.lab=2,pch=16,col=2,bty="l",
		main=paste("N Islands=",n.ilhas,"; N spp=",Nspp,"; time=",tempo),
		xlab="Island area",ylab="Number of species",ylim=c(1,Nspp))
	abline(mod,lty=2)
	cat("c=",mod[[1]][1],"z=",mod[[1]][2],"\n")
	return(riq)
	}

#arquip(areas=c(10,20,40,80),Nspp=1000,chuva.total=100,abund=rep(10,1000),tempo=10)

iColExt=function(Nspp, chuva, abund, tempo, tx.ext){
	spp=paste("sp.",1:Nspp,sep="")
	ilha=numeric()
	riq=numeric()
	ilha=unique(sample(spp,chuva,replace=TRUE,prob=abund/sum(abund)))
	riq[1]=length(ilha)

	for(i in 2:tempo){
		rr=sample(c("V","M"),length(ilha),replace=TRUE,prob=c((1-tx.ext),tx.ext))
		roleta=data.frame(ilha,rr)
		vivos=roleta[roleta$rr=="V",1]
		ilha=union(sample(spp,chuva,replace=TRUE,prob=abund/sum(abund)),vivos)
		riq[i]=length(unique(ilha))
		}
	plot(0:tempo,c(0,riq),type="l",lwd=2,bty='n',xlab="time",ylab="number of species",
	 font.lab=2,col=2,ylim=c(0,Nspp),
	 main=c("Island Richness",paste("Nspp =",Nspp," ; rain =",chuva,"; tx.ext = ",tx.ext)))
	abline(h=Nspp,lty=3)
	return(riq)
	}
#iColExt(Nspp=100, chuva=5, abund=rep(100,100), tempo=100, tx.ext=.1)

#### biogeografia de Ilhas
MW=function(areas , dist , P , a=1, b=-.01, c=1, d=-.01){
  par(mfrow=c(1,2))
  E=a+b*areas
  I=c+d*dist
  S=numeric()
  for(i in 1:length(areas)){S[i]=I[i]*P/(I[i]+E[i])}
  Tn=numeric()
	for(i in 1:length(areas)){Tn[i]=I[i]*E[i]/(I[i]+E[i])}
  
  curve(I[1]-I[1]*x/P,0,P,bty="n",xaxt="n",yaxt="n",xlab="Species number",
        ylab="Rates",font.lab=2,lwd=2,ylim=c(0,1))
  curve((E[1]/P)*x,0,P,lwd=2,add=TRUE)
  abline(v=0)
  abline(h=0)
  mtext("P",side=1,at=P,font=2)
  mtext("I1",side=2,at=I[1],font=2,las=1)
  mtext("E1",side=4,at=E[1],font=2,las=1)
  mtext("S1",side=1,at=S[1],font=2,col=2)
  mtext("T1",side=2,at=Tn[1],font=2,las=1,col=2)
  points(S[1],Tn[1],col=2,pch=16,cex=1.3)
  segments(S[1],Tn[1],S[1],0,lty=3,col=2)
  segments(S[1],Tn[1],0,Tn[1],lty=3,col=2)
  
  for(i in 2:length(areas)){
    curve(I[i]-I[i]*x/P,0,P,lwd=2,ylim=c(0,1),add=TRUE,lty=i)
    curve((E[i]/P)*x,0,P,lwd=2,add=TRUE,lty=i)
    mtext(paste("I",i,sep=""),side=2,at=I[i],font=2,las=1)
    mtext(paste("E",i,sep=""),side=4,at=E[i],font=2,las=1)
    mtext(paste("S",i,sep=""),side=1,at=S[i],font=2,col=i+1)
    mtext(paste("T",i,sep=""),side=2,at=Tn[i],font=2,las=1,col=i+1)
    points(S[i],Tn[i],col=i+1,pch=16,cex=1.3)
    segments(S[i],Tn[i],S[i],0,lty=3,col=i+1)
    segments(S[i],Tn[i],0,Tn[i],lty=3,col=i+1)
  }
  
  ex=data.frame(areas=areas,spp=S,dist=dist)
  y=lm(S~areas)[[1]][1]
  z=lm(S~areas)[[1]][2]	
  plot(spp~areas,data=ex,log="xy",font.lab=2,pch=as.character(1:length(areas)),col=2,bty="l",
       main=c("Equilibrium",paste("c = ",round(y,2),"; z = ",round(z,2))),
       xlab="Island area",ylab="Species number",ylim=c(1,P))
  abline(lm(log10(spp)~log10(areas),data=ex),lty=3)
  
  return(ex)
  par(mfrow=c(1,2))
}

#################################################
### Island Biog. Plus Rescue Effect and Internal Colonization  
MW.2.0=function(areas , dist , P , peso.A=.5 , a=1, b=-.01, c=1, d=-.01, e=0, f=.01, g=0, h=.01){
	E=((a+b*areas)*peso.A+(g+h*dist)*(1-peso.A))/(peso.A+(1-peso.A))
	I=((c+d*dist)*peso.A+(e+f*areas)*(1-peso.A))/(peso.A+(1-peso.A))
	S=numeric()
	for(i in 1:length(areas)){S[i]=I[i]*P/(I[i]+E[i])}
	Tn=numeric()
	for(i in 1:length(areas)){Tn[i]=I[i]*E[i]/(I[i]+E[i])}

	curve(I[1]-I[1]*x/P,0,P,bty="n",xaxt="n",yaxt="n",xlab="Species number",
	 ylab="Rates",font.lab=2,lwd=2,ylim=c(0,1),main="Equilibrium")
	curve((E[1]/P)*x,0,P,lwd=2,add=TRUE)
	abline(v=0)
	abline(h=0)
	mtext("P",side=1,at=P,font=2)
	mtext("I1",side=2,at=I[1],font=2,las=1)
	mtext("E1",side=4,at=E[1],font=2,las=1)
	mtext("S1",side=1,at=S[1],font=2,col=2)
	mtext("T1",side=2,at=Tn[1],font=2,las=1,col=2)
	points(S[1],Tn[1],col=2,pch=16,cex=1.3)
	segments(S[1],Tn[1],S[1],0,lty=3,col=2)
	segments(S[1],Tn[1],0,Tn[1],lty=3,col=2)

	for(i in 2:length(areas)){
		curve(I[i]-I[i]*x/P,0,P,lwd=2,ylim=c(0,1),add=TRUE,lty=i)
		curve((E[i]/P)*x,0,P,lwd=2,add=TRUE,lty=i)
		mtext(paste("I",i,sep=""),side=2,at=I[i],font=2,las=1)
		mtext(paste("E",i,sep=""),side=4,at=E[i],font=2,las=1)
		mtext(paste("S",i,sep=""),side=1,at=S[i],font=2,col=i+1)
		mtext(paste("T",i,sep=""),side=2,at=Tn[i],font=2,las=1,col=i+1)
		points(S[i],Tn[i],col=i+1,pch=16,cex=1.3)
		segments(S[i],Tn[i],S[i],0,lty=3,col=i+1)
		segments(S[i],Tn[i],0,Tn[i],lty=3,col=i+1)
		}

	ex=data.frame(areas=areas,spp=S,dist=dist)
	y=lm(S~areas)[[1]][1]
	z=lm(S~areas)[[1]][2]

	plot(spp~areas,data=ex,log="xy",font.lab=2,pch=16,col=2,bty="l",
	 main=c("Species-area relationship",paste("c = ",round(y,2),"; z = ",round(z,2))),
	 xlab="Island area",ylab="Species number",ylim=c(1,P))
	abline(lm(log10(spp)~log10(areas),data=ex),lty=3)

	return(ex)
	}
########### MODELOS NULO ###############
rand.walk <- function(n=1,step=1,ciclo=1e5,cont=1e3,x1=NULL){
  if(is.null(x1)){
    x1 <- sample(1:200,n,replace=TRUE)
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
  matplot(time,results,type="l", col=rainbow(n),lwd=2, xlab="Steps",  main="Randon Walk",ylab="Distance from the edge")
  abline(h=0,lwd=4)
}
#rand.walk(n=10,step=10,ciclo=1e4,cont=1e3)
#rand.walk(n=10,step=2,ciclo=1e4,cont=1e3)
#rand.walk(n=10,step=2,ciclo=5e4,cont=1e3)

##### Game
ext.game <- function(aposta=1,total=100){
  X <- total/2
  results <- X
  while(X>0&X<total){
    X <- X+sample(c(aposta,-1*aposta),1)
    results <- c(results,X)
  }
  plot(1:length(results),results, type="l", col="blue",ylim=c(0,total), xlab="Cicle", ylab="Number of Individuals")
  lines(1:length(results),total-results, col="red")
  abline(h=c(0,total),lty=2)
  legend(total*0.9, legend=c("sp 1", "sp2"), lty=1, col=c("red", "blue")  )
}

#MELINA MODIFICOU A LEGENDA DA FUNCAO, LEGENDA ANTIGA:
#legend(total*0.5,dim(results)*0.8, legend=c("sp 1", "sp2"), lty=1, col=c("red", "blue")  )

#old<-par(mfrow=c(2,2))
#ext.game(aposta=1,total=20)
#ext.game(aposta=1,total=50)
#ext.game(aposta=1,total=100)
#ext.game(aposta=1,total=200)
#par(old)


##Modelo neutro sem imigracao 
sim.hub1=function(S= 100, j=10, D=1, ciclo=1e4, step=1000){ 
  ## Tamanho da comunidade
  rich <- function(x)length(unique(x))
  J <- S*j
  ##Matrizes para guardar os resultados
  ## matriz da especie de cada individuo por ciclo
  ind.mat=matrix(nrow=J,ncol=1+ciclo/step) 
  ##CONDICOES INICIAIS##
  ##Deduzidas de acordo com o modelo de Hubbell:
  ## Todas as especies comecam com o mesmo numero de individuos (j=J/S)
  ind.mat[,1] <- rep(1:S,each=j)
  cod.sp <- ind.mat[,1]
  ##Aqui comecam as simulacoes
  for(i in 2:(1+ciclo/step)){
    for(j in 1:step){
      ##Indice dos individuos que morrem
      morte <- sample(1:J,D)
      ##Indice dos individuos que produzem filhotes para substituir os mortos
      novos <- sample(1:J,D,replace=TRUE)
      ##Substituindo
      cod.sp[morte]<-cod.sp[novos]
    }
    ## A cada step ciclos os resultados sao gravados
    ind.mat[,i] <- cod.sp
  }
  tempo <- seq(0,ciclo,by=step)
  colnames(ind.mat) <- tempo
  invisible(ind.mat)
  plot(tempo,apply(ind.mat,2,rich), xlab="Time (cicles)", ylab="Number of species",ylim=c(0,S), type="l", main=paste("Neutral Model Without Colonization", "\n S=",S," J=",J), sub=paste("Mean extintion=",(S-rich(ind.mat[,ncol(ind.mat)]))/ciclo,"sp/ciclo"), cex.sub=0.7) 
}

#par(mfrow=c(2,2))
#sim.hub1(j=2,ciclo=2e4,step=1e2)
#sim.hub1(j=5,ciclo=2e4,step=1e2)
#sim.hub1(j=10,ciclo=2e4,step=1e2)
#sim.hub1(j=20,ciclo=2e4,step=1e2)
#par(mfrow=c(1,1))
##


## Com migracao de uma metacomunidade com a composicao inicial
sim.hub2=function(S= 100, j=10, D=1, ciclo=1e4, step=1000, m=0.01){ 
  ## Tamanho da comunidade
  rich <- function(x)length(unique(x))
  J <- S*j
  ##Matrizes para guardar os resultados
  ## matriz da especie de cada individuo por ciclo
  ind.mat=matrix(nrow=J,ncol=1+ciclo/step) 
  ##CONDICOES INICIAIS##
  ## Todas as especies comecam com o mesmo numero de individuos (j=J/S)
  ## Rotulo de especies para cada um dos inividuos
  ind.mat[,1] <- rep(1:S,each=j)
  ## Repetindo este rotulo no vetor que sofrera modificacoes
  cod.sp <- ind.mat[,1]
  ##Aqui comecam as simulacoes
  for(i in 2:(1+ciclo/step)){
    for(j in 1:step){
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
    ind.mat[,i] <- cod.sp
  }
  tempo <- seq(0,ciclo,by=step)
  colnames(ind.mat) <- tempo
  invisible(ind.mat)
  plot(tempo,apply(ind.mat,2,rich), xlab="Time (cicles)", ylab="Number of species", type="l",
       main="Neutral Dynamics - Original Community Colonization",sub=paste( "S=",S," J=",J," m=",m,"Mean Extintion rate =",(S-rich(ind.mat[,ncol(ind.mat)]))/ciclo,"sp/ciclo"),ylim=c(0,S), cex.sub=0.7)
  }
#teste2 <- sim.hub2(j=2,ciclo=2e4,step=1e2,m=0.1)

## Com migracao de uma metacomunidade com especiacao
### funcao elaborada por Paulo Inacio Prado e modificada por Alexandre Adalardo Seg 21 Nov 2011 20:56:53 BRST 
sim.hub3=function(Sm=200, jm=20, S= 100, j=10, D=1, ciclo=1e4, step=1000, m=0.01, nu=0.001)
{ 
  rich <- function(x)length(unique(x))
  ## Tamanho da metacomunidade
  Jm <- Sm*jm
  cores<-c("#FFFFFF", topo.colors(Sm))
  ## Tamanho da comunidade
  J <- S*j
  ##Matrizes para guardar os resultados
  ## matriz da especie de cada individuo por ciclo
  ## Na metacomunidade
  meta.mat=matrix(nrow=Jm,ncol=1+ciclo/step) 
  ## Na comunidade
  ind.mat=matrix(nrow=J,ncol=1+ciclo/step)
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
  ##Aqui comecam as simulacoes
  for(i in 2:(1+ciclo/step))
  {
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
    ind.mat[,i] <- cod.sp
    meta.mat[,i] <- meta.sp
    quad<-sqrt(J)
    if(is.integer(quad))
    {
    image(matrix(cod.sp,ncol=quad, nrow=quad), col=topo.colors(Sm),xaxt="n", yaxt="n", main=" Null Community Dynamics", sub=paste("time =",(i-1) * ciclo/step ))
    }
    else
    {
    lad1=ceiling(quad)
    lad2=floor(quad)
    are=lad1^2
    vaz=are-J
    cod.aj<-c(rep(1,vaz),cod.sp+1)
    
    image(matrix(cod.aj, ncol=lad1, nrow=lad2), col=cores,xaxt="n", yaxt="n", main=" Null Community Dynamics", sub=paste("time =",(i-1) * ciclo/step ))
	}
  }
  tempo <- seq(0,ciclo,by=step)
  colnames(ind.mat) <- tempo
  colnames(meta.mat) <- tempo
 
  resultados <- list(metacomunidade=meta.mat,comunidade=ind.mat)
  ## Graficos
  x11()
  plot(tempo,apply(meta.mat,2,rich), xlab="Time (cicles)", ylab="Number of species", type="l",
       main="Neutra Dynamics - Metacomunity Colonization" ,sub=paste( "Jm=",Jm," nu=",nu," Theta=",2*Jm*nu, "S=",S," J=",J," m=",m, " Mean Extintion Rate=",(S-rich(ind.mat[,ncol(ind.mat)]))/ciclo,"sp/ciclo"), col="blue",  ylim=c(0,max(apply(meta.mat,2,rich))), cex.sub=0.7)
  lines(tempo,apply(ind.mat,2,rich),col="red")
  text(ciclo/2 ,length(unique(ind.mat[,round(dim(ind.mat)[2]/2)]))*1.3, "Community", col="red")
  text(ciclo/2 ,length(unique(meta.mat[,round(dim(ind.mat)[2]/2)]))*0.95, "Metacommunity", col="blue")
  invisible(resultados)
  }

#teste3 <- sim.hub3(j=10, ciclo=2e4,step=1e2,m=0.1)
#teste3 <- sim.hub3(j=2, ciclo=2e5,step=1e3,nu=0.00001,m=0.1)
