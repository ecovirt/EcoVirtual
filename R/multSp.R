### EcoVirtual - multispecies functions
###############
############################
### Sucessional Niche
############################
regNicho=function(tmax, ln, cl, c1,c2, ec, dst,  Er, Sc, Mx, Rs) 
{
N=cl*ln
V=1-Er-Sc-Mx-Rs
cena=array(NA,dim=c(ln,cl,tmax))
cena[,,1]<-sample(c(0:4), N, prob=c(V,Er,Sc,Mx, Rs), replace=TRUE)
resulta=matrix(0, ncol=5, nrow = tmax)
conta=table(cena[,,1])/N
resulta[1,(as.numeric(names(conta))+1)]<-conta
	for (t in 2:tmax)
	{
	Vvf<-cena[,,t-1]==0
	nV=sum(Vvf)
	Ervf<-cena[,,t-1]==1
	nEr=sum(Ervf)
	Scvf<-cena[,,t-1]==2
	nSc=sum(Scvf)
	Mxvf<-cena[,,t-1]==3
	nMx=sum(Mxvf)
	Rsvf<-cena[,,t-1]==4
	nRs=sum(Rsvf)
	p_col1=c1*(nSc+nMx+nRs)/N
	p_col2=c2*(nEr+nMx)/N
	p_ncol=1-p_col1- p_col2
	p_permEr = 1- (dst + p_col1)
	p_permSc = 1- (dst + p_col2 + ec)
	p_permMx = 1 - (dst + ec)
	if(p_ncol<0){p_ncol=0}
	if(p_permEr<0){p_permEr=0}
	if(p_permSc<0){p_permSc=0}
	if(p_permMx<0){p_permMx=0}
	cena[,,t][Vvf]<-sample(c(0,1,2),nV, replace=TRUE, prob=c(p_ncol,p_col2,p_col1)) 
	cena[,,t][Ervf]<-sample(c(0,1,3), nEr, replace=TRUE, prob=c(dst,p_permEr , p_col1))
	cena[,,t][Scvf]<-sample(c(0,2,3,4), nSc, replace=TRUE, prob=c(dst,p_permSc, p_col2, ec))
	cena[,,t][Mxvf]<-sample(c(0,3,4), nMx, replace=TRUE, prob=c(dst,p_permMx, ec))
	cena[,,t][Rsvf]<-sample(c(0,4), nRs, replace=TRUE, prob=c(dst,1 - dst))
	conta=table(cena[,,t])/N
	resulta[t,(as.numeric(names(conta))+1)]<-conta
	}
animaCena(cena)
x11()
matplot( 1:tmax,resulta[,2:5], type="l", main="Niche Regeneration Model" , xlab="time", ylab="State proportion", lty=2:5, col=2:5)
legend("topright", c("Early", "Susceptible", "Mixed", "Resistant"), bty="n", lty=2:5, col=2:5, cex=0.7)
invisible(cena)
}

#regNicho(tmax=50, ln=100, cl=100, c1=0.2, c2=0.8, ec=0.5, dst=0.04,  Er=0.08, Sc=0.02, Mx=0, Rs=0)

##########################
# Trade-off
#########################
comCompete = function(tmax,ln,cl, rq, fi, fsp1, pe,fr=0,int=0)
{
rank=1:rq
vetor_dist=rep("n", tmax)
  if(fr>0 & int>0)
  {
#  n_dist=round(tmax*fr,0)
#  max_dist=n_dist*fr
  vetor_dist[round(seq(0,tmax, length.out=fr*tmax),0)[-1]]="d"
  }
vetor_dist=vetor_dist[-1]
ci= pe/(1-fsp1)^(2*rank-1)
N <- ln*cl
resulta=matrix(nrow=rq,ncol=tmax)
attributes(resulta)=c(attributes(resulta),list(tempo=tmax, riqueza=rq, n_manchas=N, inicial= fi, comp=fsp1, coloniza=ci, freq_dist=fr, int=int))
n_ocup= sum(fi)*N
temp=1
	if(length(fi)==rq)
	{ 
	resulta[,1]=fi
	antes=sample(c(rep(0,N-n_ocup),rep(rank,each=fi*N)))
	}
	if(length(fi)==1)
	{
	antes <- sample(c(1:rq, sample(1:rq, (n_ocup-rq),replace=TRUE), rep(0, N-n_ocup)))
	t_antes=table(antes)
	n_sp=rep(0,rq)
	names_antes=match(as.numeric(names(t_antes)), rank)
	rank_antes=match(rank,as.numeric(names(t_antes)))
	n_sp[na.omit(names_antes)]=t_antes[na.omit(rank_antes)]
	resulta[,1]<-n_sp/N
	}
	for (f in vetor_dist)
		{
		temp=temp+1
		depois <- rep(0,N)
		#ncum=cumsum(resulta[,1])-resulta[1,1] ## o que é isso?
		pi=ci*resulta[,temp-1]
		pi[pi>1]=0.999
			for(rs in rq:1)
			{
			depois[antes==rs]<-sample(c(0,rs),sum(antes==rs),replace=TRUE,prob=c(pe,1-pe))
			d1<-sample(c(0,rs),sum(antes>rs | antes==0),replace=TRUE,prob=c(1-pi[rs],pi[rs]))
			depois[antes>rs | antes==0][d1==rs] <- rs
			}
		n_sp=rep(0,rq)
		if(f=="d")
		{
		depois[sample(1:N,N*int)]=0
		}
		t_depois=table(depois)
		names_match=match(as.numeric(names(t_depois)), rank)
		rank_match=match(rank,as.numeric(names(t_depois)))
		n_sp[na.omit(names_match)]=t_depois[na.omit(rank_match)]
		resulta[,temp]<-n_sp/N
		antes<-depois
		}
#grafico
x11()
layout(matrix(data=c(1,2), nrow=2, ncol=1), widths=c(1,1), heights=c(4,1))
matplot(1:tmax,t(resulta),type="l", lty=3, col=rainbow(rq),bty="n", lwd=2,xlab="Time", ylab="Patch occupance", main="Competition/Colonization Trade-off", sub=paste("\n fsp1 =",fsp1,"; pe=",pe, "; fr=",fr, "int=", int), cex.sub=0.7) 
old<-par(mar=c(3,5,2,4))
image(x=1:rq, y=1, matrix(data=1:rq, nrow=rq,ncol=1),col=rainbow(rq), xlab="competition/colonizatio scale", ylab="",xaxt="n", yaxt="n")
axis(1, at=c(1.5,9.5),tick=FALSE, labels=c("better competitor", "better colonizator"))
par(old)
##resultado
invisible(resulta)
}


#comCompete(tmax=1000,ln=100,cl=100, rq=10, fi=1, fsp1=0.20, pe=0.01,fr=0,int=0)


#########################
##### Sucession
####################
sucMatrix=function(mat_trans, init_prop, ln, cl, tmax)
{
mat_trans=as.matrix(mat_trans)
porc1=apply(mat_trans,2,sum)
if(sum(porc1!=1)>0)
{
stop("the transition for each fase should sum 1: there is no extintion of area under the model")
}
if(sum(init_prop)!=1 | length(init_prop) != dim(mat_trans)[2])
{
stop("the initial proportion of ocupance should sum 1 and the number of stages should be iqual to transition matrix")
}
nfase=dim(mat_trans)[1]
ncel=ln*cl
fase_n=round(init_prop*ncel)
cl_fase=colorRampPalette(c("gray","yellow", "orange","green"))
#cores=c("#ffffff",colors(nfase-1))
#cores=terrain.colors(nfase)
arena=matrix(NA,nrow=ln,ncol=cl)
#resulta=matrix(0,ncol=nfase, nrow=tmax)
pais=array(0,dim=c(ln, cl, tmax))
#image(0:ln,0:cl, matrix(0,nrow=ln,ncol=cl), col="white", xlab="", ylab="")
#grid(ln,cl)
n0=sample(rep(0:(nfase-1), fase_n))
arena[1:ncel]<-n0
#t.n0=table(n0)
#resulta[1,as.numeric(names(t.n0))+1]<-t.n0
pais[,,1]<-arena
image(0:ln, 0:cl, arena, col=cl_fase(nfase) , breaks=c(-0.1,(1:nfase)-0.1), xlab="", ylab="", main="Sucessional Model")
grid(ln,cl)
	for (tc in 2:tmax)
	{
		for(nf in 0:(nfase-1))
		{
		nf_vf=pais[,,(tc-1)]==nf
		contn=sum(nf_vf)
		pais[,,tc][nf_vf]<-sample(0:(nfase-1),contn,replace=TRUE, prob=as.numeric(mat_trans[,(nf+1)]))
		}
#	resulta[tc,as.numeric(names(t_n0))+1]<-t_n0
	image(0:ln, 0:cl, pais[,,tc], col=cl_fase(nfase) , breaks=c(-0.1,(1:nfase)-0.1), add=TRUE)
	Sys.sleep(.1)
	}
x11()
op=par(mfrow=c(2,2))
image(0:ln, 0:cl, arena, col=cl_fase(nfase) , breaks=c(-0.1,(1:nfase)-0.1), xlab="", ylab="", main="Ernitial Conditions")
for(ts in c(4,2,1))
	{
	image(0:ln, 0:cl, pais[,,round(tc/ts)], col=cl_fase(nfase) , breaks=c(-0.1,(1:nfase)-0.1), main=paste("Cicle", round(tc/ts)))
	}
par(op)
x11()
resulta=t(apply(pais,3, table))
matplot(resulta, type="l", ylim=c(min(resulta)*0.8, max(resulta)*1.1), main="Stage Distribution",xlab="Number of patchs", col=cl_fase(nfase), lty=2, lwd=2)
legend("topright", legend=paste("Stage", 1:nfase), lty=2, lwd=2, col=cl_fase(nfase), bty="n", cex=0.8)
#resulta=as.data.frame(resulta)
eigs_st=eigen(mat_trans)
dom_pos=which.max(Re(eigs_st$values))
stage_v<- Re(eigs_st[["vectors"]][, dom_pos])
stage_stable=(stage_v/sum(stage_v))*ncel
abline(h=stage_stable, col=cl_fase(nfase), lwd=0.8)
legend("topleft", legend=paste("Stage Stable", 1:nfase), lty=1, lwd=0.9, col=cl_fase(nfase), bty="n", cex=0.8)
invisible(pais)
}

