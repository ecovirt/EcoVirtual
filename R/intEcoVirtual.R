#######################################
### EcoVirtual - multispecies internal functions
###############
############################
### graficos metacomunidade
#####################
metacomp.anima=function(dados)
{
nsim=dim(dados)[3]
ln=dim(dados)[1]
cl=dim(dados)[2]
op=par(mar=c(1,2,2,2))
layout(matrix(c(2,1), ncol=1, nrow=2), heights=c(5,1),widths=c(1,1))
plot(1:10,1:10,xaxt="n", yaxt="n", xlab="", ylab="", cex=0.8,type="n", , bty="n")
legend(2,10,ncol=4, legend=c("not available", "empty", "species 1", "species 2"), pch=c(15,22,15,15), title="Patches legend", col=c("red","black", "blue", "green"),bty="n")
#layout.show()
#par(op)
image(0:ln, 0:cl, dados[,,1], col=c("red", "white","blue" ,"green") , breaks=c(-0.9,-0.001,0.1,1.5,2.9),main="Metapopulations Competition",  xlab="", ylab="")
grid(ln,cl)
Sys.sleep(.5)
	for(i in 2:nsim)
	{
	par(new=TRUE)
image(0:ln, 0:cl, dados[,,i], col=c("red", "white","blue" ,"green") , breaks=c(-0.1,-0.001,0.1,1.9,2.9), xlab="", ylab="")
#	,main="Metapopulations Competition", main=paste("Simulation no.", i,"; total", nsim,   sep="")
grid(ln,cl)
	Sys.sleep(.1)
	}
}
###############################################


###################
anima <-function(dados)
{
x11()
	for(i in 1:dim(dados)[3])
        {
	image(dados[,,i], main=("Metapopulation"),sub=paste("simulation no.= ", i), col=c("white","red"), bty="n",xaxt='n',yaxt='n')
	grid(dim(dados)[1],dim(dados)[2])
	Sys.sleep(.2)
	}
}
########################################################
meta.anima2=function(dados)
{
nsim=dim(dados)[3]
ln=dim(dados)[1]
cl=dim(dados)[2]
image(0:ln, 0:cl, dados[,,1], col=c("white", "green") , breaks=c(0,0.99,5),main="Metapopulation Dynamics", sub=paste("Initial configuration from", nsim," simulations",  sep=""), xlab="", ylab="")	
grid(ln,cl)
Sys.sleep(.5)
conta12=dados[,,1]+ (2*dados[,,2])
image(0:ln, 0:cl, conta12, col=c("white","red","lightgreen", "darkgreen") , breaks=c(0,0.9,1.9,2.9,3.9),main="Metapopulation Dynamics", sub=paste("red= extintion; light green= colonization; dark green = permanence \n maximum time = ", nsim, sep=""), xlab="", ylab="")
	for(i in 3:nsim)
	{
	conta12=dados[,,(i-1)]+ (2*dados[,,i])
	image(0:ln, 0:cl, conta12, col=c("white","red","lightgreen", "darkgreen") , breaks=c(0,0.9,1.9,2.9,3.9), xlab="", ylab="", add=TRUE)
	Sys.sleep(.1)
	}
}
###############################################
graf.fim=function(dados)
{
op=par(mfrow=c(2,2))
nsim=dim(dados)[3]
ln=dim(dados)[1]
cl=dim(dados)[2]
image(0:ln, 0:cl, dados[,,1], col=c("white", "green") , breaks=c(0,0.99,5),main="Metapopulation Dynamics", sub=paste("time = 1/", nsim, sep=""), xlab="", ylab="")	
grid(ln,cl)
	for(ts in c(4,2,1))
	{
	sim=round(nsim/ts)
	conta12=dados[,,(sim-1)]+ (2*dados[,,sim])
	image(0:ln, 0:cl, conta12, col=c("white","red","lightgreen", "darkgreen") , breaks=c(0,0.9,1.9,2.9,3.9),main="Metapopulation Dynamics", sub=paste("red= extintion; light green= colonization;\n dark green = permanence \t time = ", sim, "/", nsim, sep=""), xlab="", ylab="")
#	Sys.sleep(.5)
	}
par(op)
}
###############################
#Trade-off Multispecies Graphic
### 
gr.toff=function(rq, fsp1,pe,add=FALSE,...)
{
#	rq <- as.numeric(tclvalue(rqVar))
#	fsp1 <- as.numeric(tclvalue(fsp1Var))
#	pe <- as.numeric(tclvalue(peVar))
	rank=1:rq
	ci= pe/(1-fsp1)^(2*rank-1)
	px= fsp1*(1-fsp1)^(rank-1)
		if(add==FALSE)
		{
		toff<-x11( width=5, height=5)
		}
	old<-par(mar=c(3,3,3,3))
	plot(ci~rank, col="red",ylim=c(0,max(ci)*1.2), type="b", ann=FALSE, axes=FALSE,)
	axis(4, cex.axis=0.8, col.axis='red', col='red')#, yaxp=c(0,3,3))
	par(new=TRUE)
	plot(px~rank, ylim=c(0,fsp1),type="b", bty="n",  ann=FALSE, cex.axis=0.8)#yaxt="n", xaxp=c(0,10,5))
	#axis(2, cex.axis=0.8)#, yaxp=c(0,0.2,4))
	mtext("Specie competitive rank", 1, 2, cex=0.9)
	mtext("Abundance", 2, 2, cex=0.9)
	mtext("Colonization rate", 4, 2, cex=0.9, col='red')
	mtext("Trade-off Species Rank ", 3, 0, cex=1.2)
	par(old)
}

#ex:
gr.toff(rq =  10 , fsp1 =  0.2 , pe =  0.1 ,add=FALSE)

############################
### Sucessional Niche Graphic
############################
anima.cena=function(dados)
{
nt=dim(dados)[3]
x11()
op=par(mfrow=c(5,5),  mar=c(0.1,0.1,0.1,0.1))
	for(i in 1:nt)
	{
	image(dados[,,i], main="",  bty="n",xaxt='n',yaxt='n', col=c("white", "yellow", "orange", "blue", "green"))
	grid(dim(dados)[2],dim(dados)[1])
	}
x11()
par(mfrow=c(2,2))
image(dados[,,1], main= paste("Patches occupancy\n \ttime=", 1 ),  bty="n",xaxt='n',yaxt='n',col=c("white", "yellow", "orange", "blue", "green"))
grid(dim(dados)[2],dim(dados)[1])
image(dados[,,round(nt/3)], main= paste("Patches occupancy\n \t time=", round(nt/3) ),  bty="n",xaxt='n',yaxt='n',col=c("white", "yellow", "orange", "blue", "green"))
grid(dim(dados)[2],dim(dados)[1])
image(dados[,,round(2*nt/3)], main= paste("Patches occupancy\n \ttime=", round(2*nt/3) ),  bty="n",xaxt='n',yaxt='n',col=c("white", "yellow", "orange", "blue", "green"))
grid(dim(dados)[2],dim(dados)[1])
image(dados[,,nt], main= paste("Patches occupancy\n \t time=", nt ),  bty="n",xaxt='n',yaxt='n',col=c("white", "yellow", "orange", "blue", "green"))
grid(dim(dados)[2],dim(dados)[1])
par(op)
}

##################################
graf.com=function(medias, desvios, minimo, maximo)
{
  nsp=length(medias)
  cor=rainbow(nsp)
  curve(dnorm.trunc(x, medias[1], desvios[1], maximo=maximo, minimo=minimo),from=minimo, to=maximo, ylim=c(0,1), ylab="Population Density", xlab="Gradient Value", main="Species Distribution", col=cor[1])
	for (i in 2:nsp)
 	{
 	curve(dnorm.trunc(x, medias[i], desvios[i], maximo=maximo, minimo=minimo),from=minimo, to=maximo,add=TRUE, col=cor[i], lty=2)
 	} 
text(medias+1, dnorm.trunc(medias, medias, desvios,maximo=maximo,minimo=minimo)+0.5, labels=(paste("sp",(1:(nsp)),sep="_")), col=cor, cex=0.8)
}  
###################################
dnorm.trunc=function(x, minimo=-Inf, maximo=Inf, media=0, desvio=1)
{
res=numeric(length(x))
x.prov=dnorm(x,mean=media, sd=desvio)
ampl.norm=pnorm(maximo,mean=media, sd=desvio)-pnorm(minimo,mean=media, sd=desvio)
x.prov/ampl.norm
}
###################################
####################################
### f para truncar a amostra ###
###################################
pnorm.trunc=function(x,minimo=-Inf, maximo=Inf, media=0, desvio=1)
{
denom <- pnorm(maximo, mean=media, sd=desvio) - pnorm(minimo, mean=media, sd=desvio)
qtmp <- pnorm(x, mean=media, sd=desvio) - pnorm(minimo, mean=media, sd=desvio)
qtmp/denom
}
##### Proportion of species at each sample
prob.ssp=function(medias, desvios, amostra, minimo, maximo)
{
nsp=length(medias)
namostra=length(amostra)
resulta=matrix(NA, nrow=nsp, ncol=namostra)
rownames(resulta)=paste("sp", 1:nsp, sep="_")
colnames(resulta)=paste("plot", 1:namostra, sep="_")
  for(k in 1:namostra)
  {
	for(i in 1:nsp)
	{
	resulta[i,k]= pnorm.trunc(amostra[k]+1,minimo=minimo, maximo=maximo, media=medias[i], desvio=desvios[i])- pnorm.trunc(amostra[k],minimo=minimo, maximo=maximo, media=medias[i], desvio=desvios[i] )
	}
  }
invisible(resulta)
}
#############################
######### Func Distancia Bray-Curtis
####################################
ddis.bc<-function(dados)
	{
	nplot=dim(dados)[2]
	similar=matrix(NA,ncol=nplot,nrow=nplot)
	rownames(similar)<-paste("plot", c(1:nplot))
	colnames(similar)<-paste("plot", c(1:nplot))
		for(i in 1:(nplot-1))
		{
		m=i+1
			for(m in m:nplot)
			{
			bc.dist=sum(abs(dados[,i]-dados[,m]))/(sum (dados[,c(i,m)]))
			similar[m,i]=bc.dist
			}
		}
	invisible(round(similar,3))
	}
###########################
####Polar-Ordination Func
########################
ordena.polar=function(dist)
{
somadist1.cont=apply(dist, 1, sum, na.rm=TRUE) + apply(dist,2,sum, na.rm=TRUE)
nomes.parc=names(somadist1.cont)
parc.ax=nomes.parc[somadist1.cont==max(somadist1.cont)][1]
dist.ax=rbind(dist[,parc.ax], dist[parc.ax,])
dist.ax=apply(dist.ax,2,sum, na.rm=TRUE)
max.ax=max(dist.ax)
parc.bx=nomes.parc[dist.ax==max.ax]
	if(length(parc.bx)>1)
	{
	  somamax.bx=max(somadist1.cont[parc.bx])
	  parc.bx=nomes.parc[somadist1.cont==somamax.bx][1]
	  parc.bx
	}
dist.bx=rbind(dist[,parc.bx], dist[parc.bx,])
dist.bx=apply(dist.bx,2,sum, na.rm=TRUE)
xi= (max.ax^2 + dist.ax^2 - dist.bx^2)/(2*max.ax)
yi=sqrt((dist.ax)^2-xi^2)
yi[parc.bx]=max(dist.ax)
op.xy=data.frame(xi,yi)
plot(op.xy, pch=19, col=rainbow(length(xi)), xlim=c(-0.1, 1), ylim=c(-0.1,1), main="Polar Ordination", sub="Bray-Curtis distance")
text(op.xy-0.05, labels=rownames(op.xy))
invisible(op.xy)
}
############################
### Matrix Similarity
############################
sim<-function(dados, indice="bc")
	{
	nplot=dim(dados)[2]
	similar=matrix(1,ncol=nplot,nrow=nplot)
	rownames(similar)<-paste("plot", c(1:nplot))
	colnames(similar)<-paste("plot", c(1:nplot))
		for(i in 1:(nplot-1))
		{
		m=i+1
		for(m in m:nplot)
		{
		if(indice=="jacc")
			{
			dados[dados>0]=1
			co.oc=sum(dados[,i]>0 & dados[,m]>0)
			total.sp=sum(dados[,i])+sum(dados[,m])-co.oc
			similar[i,m]=co.oc/total.sp 
			similar[m,i]=co.oc/total.sp
			}
		if(indice=="bc") 
			{
			bc.sim=sum(apply(dados[,c(i,m)], 1, min))/(sum (dados[,c(i,m)]))
			similar[i,m]=bc.sim
			similar[m,i]=bc.sim
			}
		}
		}
	invisible(round(similar,3))
	}
#################################
############################
### hcluster
############################
#  clas.cont1=hclust(as.dist(1-sim.cont1), method="average")
#  dend.cont1=as.dendrogram(clas.cont1, hang=-1)
#  plot(dend.cont1)


#############################
rich <- function(x)length(unique(x))
#####################
###### Grafico biog ilha
grafeq=function(E , I , P){
	S = I*P/(I+E) ; T = I*E/(I+E)
	curve(I-I*x/P,0,P,bty="n",xaxt="n",yaxt="n",xlab="Species number",
	 ylab="Taxas",font.lab=2,lwd=2,ylim=c(0,1))
	curve((E/P)*x,0,P,lwd=2,add=T)
	abline(v=0)
	abline(h=0)
	mtext("P",side=1,at=P,font=2)
	mtext("I",side=2,at=I,font=2,las=1)
	mtext("E",side=4,at=E,font=2,las=1)
	mtext("S",side=1,at=S,font=2,col=2)
	mtext("T",side=2,at=T,font=2,las=1,col=2)
	points(S,T,col=2,pch=16,cex=1.3)
	segments(S,T,S,0,lty=3,col=2)
	segments(S,T,0,T,lty=3,col=2)
	}
############################################
