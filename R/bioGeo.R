######################################################
###### Biogeography functions - ECOVIRTUAL PACKAGE ###
############# Specie Area relationship ###############
######################################################################
################# Arquipelago: ilhas de diferentes tamanhos ###########
### Thu 17 Nov 2011 05:32:27 PM BRST Alexandre Adalardo
## fuction Rich
rich <- function(x)length(unique(x))
####################################
arquip=function(nIsl,ar.min, ar.max, Nspp, chuva.total, abund, tmax, anima=TRUE)
{
	#n.ilhas=nIsl
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
		chuva=sample(spp, chuva.total, prob=abund, replace=TRUE)
		loc.x=sample(local, chuva.total, replace=TRUE)
		loc.y=sample(local, chuva.total, replace=TRUE)
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
plot(area.isl,riq.final,log="xy",pch=16,col=2,bty="l",main=paste("Nº Islands=",nIsl,"; Nº spp=",Nspp,"; Time=",tmax),xlab="Island Area",ylab="Number of species",ylim=c(1,max(riq.final)))
abline(mod1, lty=2)
rqz<-apply(cena, c(2,3), sum)
clz<-diff(riq.tempo)
matplot(riq.tempo[2:100,],clz, type="l", col=rainbow(nIsl), bty="l", cex.lab=1.2, xlab="Species Number", ylab="Colonization (species)", cex.axis=1.2, main="Colonization Rate Curves", cex.main=1.2 )

invisible(cena)
}
#######################################################
#cena<-arquip(nIsl=10,ar.min=10, ar.max=100, Nspp=1000, chuva.total=100, abund=10, tmax=100, anima=TRUE)
########################################################
animaIsl=function(riq.tempo, ar.isl, locxy, sprain)
{
Nspp=max(riq.tempo)
maxt=dim(riq.tempo)[1]
nIsl<-length(ar.isl)
comp.max<-max(ar.isl)
tempo=length(riq.tempo)
col_spp=rainbow(max(riq.tempo))
col_func=colorRamp(c("white", "green3"))
col_riq=rgb(col_func(seq(0,1, length.out=Nspp)), max=255)
## aqui inicia o grafico
layout(matrix(data=c(2,1), nrow=2, ncol=1), widths=c(1,1), heights=c(5,1))
old<-par(mar=c(2,2,1,3))
image(x=1:Nspp, y=1, matrix(data=1:Nspp, nrow=Nspp,ncol=1),col=col_riq, ylab="",xlab=paste("cicle", 1:length(maxt)), xaxt="n", yaxt="n", main="Richness")
axis(3, at=c(1.5,Nspp),tick=FALSE, labels=c("0", Nspp), mgp=c(0,0,0))
polygon(x=c(1.5,1.5,Nspp,Nspp), y=c(0.6,1.4,1.4,0.6), lwd=2)
plot(0:comp.max, 0:comp.max, usr=c(0,comp.max,0,comp.max), type="n", yaxt="n", xaxt="n", xlab="", ylab="", bty="n", main="Passive Sampling and Area ",mar=c(0,2,3,2), oma=c(0,0,0,0))
segments(x0=c(0,0,comp.max,0), y0=c(0,0,0,comp.max), x1=c(0,rep(comp.max,3)), y1=c(comp.max,0,comp.max,comp.max))
segments(x0=c(rep(0,nIsl), ar.isl), y0=c(ar.isl,rep(0,nIsl)), x1=c(ar.isl,ar.isl), y1=c(ar.isl,ar.isl))
	for (i in 2:maxt)
	{
	lxy=locxy[[i]]
	nspp=riq.tempo[i,]
		for(f in nIsl:1)
			{
			vert=ar.isl[f]
			polygon(x=c(0,vert, vert,0),y=c(0,0,vert,vert), col=col_riq[nspp[f]] )
			}
	points(lxy[,1],lxy[,2], col=col_spp[sprain[[i]]], pch=16)
	Sys.sleep(.1)
	}
par(old)
}
#########################################
grColExt=function(E , I , P, areas)
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
		siz_ar=2 +(areas/max(areas))
		points(S[i],T[i],col=corIsl[i],cex=siz_ar[i])
		}
	segments(S[i],T[i],S[i],0,lty=3,col=corIsl[i])
	segments(S[i],T[i],0,T[i],lty=3,col=corIsl[i])
	Sys.sleep(0.5)
	}	
#	mtext("I",side=2,at=I,font=2,las=1, line=2)
#	mtext("E",side=4,at=E,font=2,las=1)
}
# testando... grColExt(E = .5 , I = .5 , P = 100, areas=1:10)
####################################
#######################################
animaColExt=function(minimo=0.01, maximo=1, ciclos=100, Ext="crs", Col="dcr")
{
a=seq(from=minimo,to=maximo,length.out=ciclos)
b=seq(from=maximo, to=minimo, length.out=ciclos)
#nt=length(a)
if(Ext=="fix"){ext=rep(0.5,nt)}
if(Ext=="crs"){ext=a}
if(Ext=="dcr"){ext=b}
if(Col=="fix"){col=rep(0.5,nt)}
if(Col=="crs"){col=a}
if(Col=="dcr"){col=b}
grColExt(E=ext,I=col,P=100, areas=1)
}
#animaColExt(Ext='crs', Col="dcr")
####################################
###### Mon 21 Nov 2011 12:58:57 PM BRST Alexandre Adalardo
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

######################################
###teste
#bioGeoIsl(areas=c(5,10,50,80) , dist=c(10,100,100,10), P=100 , peso.A=.5 , a=1, b=-.01, c=1, d=-.01, e=0, f=.01, g=0, h=.01)
###########################
sppArea=function(c , z){
	curve(expr = c*x^z , from=1, to=10^10, xlab="Area", ylim=c(1,500),
	 ylab="Species number", font.lab=2, lwd=2, col=2, main=paste("c = ",c,"; z = ",z))
	curve(expr = c*x^z , from=1, to=10^10, xlab="Area", ylim=c(1,500),
	 ylab="Species number", font.lab=2, lwd=2, col=2, main=paste("c = ",c,"; z = ",z), log="xy")
	}
#par(mfrow=c(2,2))
#spp_area(c = 1.5 , z = .25)
#spp_area(c = 2.1 , z = .25)

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
	invisible(riq)
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
	invisible(riq)
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
	invisible(riq)
	}
#iColExt(Nspp=100, chuva=5, abund=rep(100,100), tempo=100, tx.ext=.1)

#### biogeografia de ilhas
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
  invisible(ex)
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
	curve(I[1]-I[1]*x/P,0,P,bty="n", xaxt="n",yaxt="n",xlab="Species number",	 ylab="Rates",font.lab=2,lwd=2,ylim=c(0,1),main="Equilibrium")
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
	invisible(ex)
	}
########################################
########### MODELOS NULO ###############
########################################
randWalk <- function(n=1,step=1,ciclo=1e5,x1max=200, alleq=FALSE){
  cont=round(ciclo/100)
  sleep=1/cont 
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
#rd1<-randWalk(n=10,step=10,ciclo=1e4)
#randWalk(n=10,step=1,ciclo=1e4)
#randWalk(n=10,step=1,ciclo=1e4, x1max=300, alleq=TRUE)
#rand.walk(n=100,step=2,ciclo=2e5)
###################
## animaRandWalk
################### 
animaRandWalk = function(rwData, time=2, sleep=0.1)
{
#par( )
xplus=max(time)*0.1
ymax=max(apply(rwData, 2, max))[1]
plot(time, rwData[,which.max(apply(rwData, 2, max))[1]], xlab="Steps", ylab="Distance from the edge",cex.axis=1.2, cex.lab=1.2,ylim=c(-.1* ymax,ymax), main="Randon Walk", cex.main=1.5, type="n", xlim=c(0,max(time)))

polygon(x=c(-xplus, -xplus, max(time)+xplus, max(time)+xplus), y=c(ymax*-0.15,0,0,ymax*-0.15), col="gray")
text(max(time)/2, -0.05* ymax, labels="Absortion Surface", col="red", cex=1.5)
n=dim(rwData)[2]
#ncolors= terrain.colors(n)
ncolors= rainbow(n)
	for(i in 2:length(time))
	{
		for(j in 1:n)
		{
		lines(time[1:i], rwData[1:i,j], col=ncolors[j], lty=j )
		}
 	 Sys.sleep(sleep)
	}
}
###################
## Zero Sum Game ##
###################
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
animaGame = function(xGame, total, sleep=0.01)
{
xmax=length(xGame)
xseq=1:xmax
	if(xmax>1e3){sleep=0}
	if(xmax>1e4)
	{
	indx=ceiling(seq(1,xmax, len=1000)) 
	xGame=xGame[indx]
	xseq=xseq[indx]
	}
plot(0:xmax, seq(0,total, len=xmax+1), xlab="Cicle", ylab="Number of Individuals",cex.axis=1.2, cex.lab=1.2, ylim=c(-.1* total,total+total*0.1), main="Zero Sum Game", cex.main=1.5, type="n", sub=paste("Maximum number of individuals = ", total), cex.sub=0.9)
abline(h=total/2, lty=2, col="red")
cores= c("blue","black")
#n=dim(rwData)[2]
	for(i in 2:xmax)
	{
		lines(xseq[1:i], xGame[1:i], col=cores[1], lty=2)
		lines(xseq[1:i], total - xGame[1:i], col=cores[2], lty=3)
 	 Sys.sleep(sleep)
	}
	polygon(x=c(-.2* xmax, -.2* xmax, xmax+ 0.1*xmax, xmax+ 0.1*xmax), y=c(-.2*total,0,0,-.2* total), col="gray")
	polygon(x=c(-.2*xmax, -.2*xmax, xmax+ 0.1*xmax, xmax+ 0.1*xmax), y=c(total,total+total*.5,total +total*.5,total), col="gray")
	text(xmax/2, - 0.05* total, labels="Loser", col="red", cex=1.5)
	text(xmax/2, total + 0.05* total, labels="Winner", col="green", cex=1.5)
}

################################################################
############## Neutral Model without Imigration ################
################################################################
simHub1=function(S= 100, j=10, D=1, ciclo=1e4, anima=TRUE)
{
if(ciclo<200){ciclo=200; cat("\n Minimum number of ciclos: 200\n")}
  stepseq=round(seq(101, ciclo+1, len=100))
  step=stepseq[2]- stepseq[1]
  ## Tamanho da comunidade
  rich <- function(x)length(unique(x))
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
if(anima==TRUE)
  {
  animaHub1(dadoHub=ind.mat)
  }
  x11()
    plot(as.numeric(colnames(ind.mat)),apply(ind.mat,2,rich), xlab="Time (cicles)", ylab="Number of species",ylim=c(0,S), cex.lab=1.2, type="l", col="red", lty=2,  main=paste("Neutral Model Without Colonization", "\n S=",S," J=",J), sub=paste("Mean extintion=",(S-rich(ind.mat[,ncol(ind.mat)]))/ciclo,"sp/cicle"), cex.sub=0.8) 
  invisible(ind.mat)
}
#par(mfrow=c(2,2))

#simHub1(S=10,j=10, D=1, ciclo=2e4, anima=TRUE)
#simHub1(j=5,ciclo=2e4)
#simHub1(j=10,ciclo=2e4)
#simHub1(j=20,ciclo=2e4)
#par(mfrow=c(1,1))
#################################################################
## Null Model Simulation with imigration from a Metacommunity  ##
#################################################################
simHub2=function(S= 100, j=10, D=1, ciclo=1e4, step=1000, m=0.01, anima=TRUE)
{ 
if(ciclo<200){ciclo=200; cat("\n Minimum number of ciclos: 200\n")}
  stepseq=round(seq(101, ciclo+1, len=100))
  step=stepseq[2]- stepseq[1]
  ## Tamanho da comunidade
  rich <- function(x)length(unique(x))
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
  if(anima==TRUE)
  {
  animaHub1(dadoHub=ind.mat)
  }
  ########### grafico interno ###############
  x11()
  plot(tempo,apply(ind.mat,2,rich), xlab="Time (cicles)", ylab="Number of species", type="l",
       main="Neutral Dynamics - Original Community Colonization",sub=paste( "S=",S," J=",J," m=",m,"Mean Extintion rate =",(S-rich(ind.mat[,ncol(ind.mat)]))/ciclo,"sp/cicle"),ylim=c(0,S), cex.sub=0.7)
  invisible(ind.mat)
}
#teste2 <- simHub2(j=2,ciclo=2e4,step=1e2,m=0.1)
###########################################################
########## Null Model Simulation 3 ########################
## Imigracao and speciation from a metacommunity ##########
###########################################################
### funcao elaborada por Paulo Inacio Prado e modificada por Alexandre Adalardo Seg 21 Nov 2011 20:56:53 BRST 
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
  if(anima==TRUE)
  {
  animaHub1(dadoHub=resultados$comunidade)
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
# hub3<-simHub3(Sm=200, jm=20, S= 10, j=100, D=1, ciclo=1e4, m=0.01, nu=0.001, anima=TRUE)
#dadoHub <- simHub3(j=10, ciclo=2e4,m=0.1, anima=FALSE)
#teste3 <- simHub3(j=2, ciclo=2e3,nu=0.00001,m=0.1)
#############################################################
############### animation for neutral models ################
#############################################################
animaHub1=function(dadoHub, sleep=0.1)
{
library(tcltk)
#nsp=length(unique(dadoHub[,1]))
maxsp=max(dadoHub)[1]
uniqsp=unique(as.numeric(dadoHub))
nind=dim(dadoHub)[1]
#nindsp=table(dadoHub[,1])[[1]]
nsim=dim(dadoHub)[2]
ciclo=as.numeric(colnames(dadoHub))
pb = tkProgressBar(title = "Simulation Progress", max = nsim)
riq=apply(dadoHub, 2, rich)
## definindo o tamanho do retangulo
    lado<-round(sqrt(nind))
    lado2<-ceiling(nind/lado)
    lastLine=lado*lado2 - nind
cormix=sample(rainbow(maxsp+10))
#if(lastLine !=0){cor=c("#000000", cor)}
#ffffff
cor=c("#FFFFFF", cormix)
mcor<-c("#FFFFFF00","#000000")
spcol<-c(rep(0, lastLine),dadoHub[,1])
############ escala das especies da metapopulacao ########
layout(matrix(data=c(2,1), nrow=2, ncol=1), widths=c(1,1), heights=c(5,1))
old<-par(mar=c(2,2,1,2))
image(x=1:maxsp, y=1, matrix(data=1:maxsp, nrow=maxsp,ncol=1),col=rainbow(maxsp), ylab="",xlab="", xaxt="n", yaxt="n", main="Metacommunity Species colors", cex.main=0.8)
axis(3, at = c(1,maxsp), labels = c(1, maxsp), tick = FALSE, mgp=c(1,0,0), cex.axis=0.8)
hmat=matrix(spcol,ncol=lado, nrow=lado2)
#cormat=matrix(cor[factor(spcol, levels=0:maxsp)], ncol=lado, nrow=lado2)
par(mar=c(2,2,2,2))
image(hmat, col=cor[sort(unique(as.numeric(hmat)))], xaxt="n", yaxt="n")
#mtext(text="simulation ", side=1, adj=0)
grid(nx=lado2, ny=lado)
#mtext(text="                   1", side=1, col="white",adj=0)
	for (i in 2:nsim)
	{
	#if(riq[i]==1){cor=cor[unique(dadoHub[,i])[1]]}
#	mtext(text=paste("                    ", ciclo[i-1]), side=1, col="white", adj=0)
	mvf=dadoHub[,i-1]!=dadoHub[,i]
	matm<-matrix(c(rep(TRUE, lastLine),mvf ),ncol=lado, nrow=lado2)
	image(matm,col=mcor, add=TRUE)
	Sys.sleep(sleep)
	spcol<-c(rep(0, lastLine),dadoHub[,i] )
	cores=cor[sort(unique(spcol)+1)]
	scol<-sort(unique(spcol))
	lcol<-length(scol)
	mcol<-match(spcol, scol)
	hmat=(matrix(mcol,ncol=lado, nrow=lado2))
	image(hmat, col=cores, add=TRUE)
	grid(nx=lado2, ny=lado)
#	mtext(text=paste("                    ", ciclo[i]), side=1, adj=0)
	setTkProgressBar(pb, value = i, label = paste("Simulation #", ciclo[i], sep="")) 
	}
	close(pb)
}
#animaHub1(dadoHub)


