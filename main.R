
# Set environments --------------------------------------------------------
rm(list=ls()); graphics.off()
library(latex2exp)


# Functions ---------------------------------------------------------------
trapezoidal_fn = function(x=NULL,lb=0,ub=1,m1=0.5,m2=0.5){
  #It generalizes: (i) Triangular case (m1=m2), (ii) Rectangular case (lb=m1,ub=m2)
  if(is.null(x)){x=seq(from=lb,to=ub,length.out=101)}
  y=rep(1e-99,length(x))
  y[x>=lb&x<m1] = (x[x>=lb&x<m1]-lb)/(m1-lb)
  y[x>=m1&x<=m2] = 1
  y[x>m2&x<=ub] = (ub-x[x>m2&x<=ub])/(ub-m2)
  return(y)
}


logLikel_fun = function(Alpha,sigma_eta,mu_eta=0,Y,L,R,H=21,Tm,Dm){
  W=statmod::gauss.quad(n = H,kind = "hermite")
  FY = matrix(NA,I*J,H)
  
  for(h in 1:H){
    Delta_h = matrix(sqrt(2)*sigma_eta*W$nodes[h]+mu_eta,I,N)
    
    PY = matrix(NA,I*J,M) #JIxM (J nested within I)
    #ETA = kronecker(Eta,matrix(1,J,1)) #JIxN (rep each Eta J times)
    ETA = kronecker(Delta_h,matrix(1,J,1)) #JIxN (rep each Eta J times)
    ALPHA = kronecker(matrix(1,I,1),Alpha) #JIxN (rep each Alpha I times)
    for(m in 1:M){
      TM = matrix(1,I*J,1)%*%Tm[m,] 
      DM = matrix(1,I*J,1)%*%Dm[m,]
      PY[,m] = apply( (exp((ETA+ALPHA)*TM)/(matrix(1,J*I,N)+exp((ETA+ALPHA))))^DM, 1, prod)
    }
    
    Omega_t = 1-diag(PY[,Y])
    JJD=t(apply(Y,1,function(x){as.numeric(1:M<x)}))
    p_s = apply(PY*JJD,1,sum) / Omega_t
    
    vs = apply(PY,1,function(x){sum(x*(1:M-sum(x*1:M))^2)})
    #csi = vs/(1+vs)
    csi = apply(PY,1,function(x){-sum(x*log(x))/log(M)})
    
    K = ((factorial(Y-1)*factorial(M-Y)) / (factorial(L)*factorial(Y-1-L) * factorial(R)*factorial(M-Y-R)))
    FY[,h] = diag(PY[,Y]) * ( csi*K*( p_s^(L+M-Y-R) * (1-p_s)^(R+Y-1-L) - 0^(L+R) * 1^(M-R-1-L)) + K*0^(L+R) * 1^(M-R-1-L)) #over i=1..I
    
  }
  
  likx = apply(FY*matrix(W$weights,I,H,byrow = TRUE),1,sum)*(1/sqrt(pi))
  likx[likx==0 | likx<0] = 1
  loglikx = sum(log(likx))
  #loglikx = sum(log(apply(FY*matrix(W$weights,I,H,byrow = TRUE),1,sum)*(1/sqrt(pi))))
  return(loglikx)
}

add_legend = function(...) {
  #From: https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}



# Case study --------------------------------------------------------------
load("datard.Rdata")
J=1 #number of items
I=NROW(datard) #number of subjects
M=4 #number of response categories
N=M-1 #number of nodes
y = datard$RDB_c; l = datard$RDB_l; r = datard$RDB_r

#Mapping matrix: response categories (rows) x nodes (cols)
#It represents a linear response tree
Tm = matrix(c(0,NA,NA,1,0,NA,1,1,0,1,1,1),M,N,byrow=TRUE); rownames(Tm) = paste0(c(1:M));colnames(Tm) = paste0("node",1:N); print(Tm)
Dm=matrix(1,M,N); for(n in 1:N){Dm[is.na(Tm[,n]),n]=0} #Dm matrix (delta_mn values)

## Model 1
res1 = optim(par = c(0.1,1.5),fn = function(x){logLikel_fun(matrix(c(x[1],x[1],x[1]),J,N),exp(x[2]),0,matrix(y,I,1),matrix(l,I,1),matrix(r,I,1),H=21,Tm,Dm)},hessian = FALSE,
             method = "BFGS",control = list(fnscale=-1,trace=3))
res1$par[1] #alpha
exp(res1$par[3]) #sigma
bic1 = -2*res1$value + length(res1$par)*log(I)

## Model 2
D = model.matrix(~datard$sex)
res2 = optim(par = c(1,0,1.5),fn = function(x){logLikel_fun(matrix(c(x[1],x[1],x[1]),J,N),exp(x[3]),D%*%c(0,x[2]),matrix(y,I,1),matrix(l,I,1),matrix(r,I,1),H=21,Tm,Dm)},hessian = FALSE,
             method = "BFGS",control = list(fnscale=-1,trace=3))
res2$par[1:2] #alpha, beta
exp(res2$par[3]) #sigma
bic2 = -2*res2$value + length(res2$par)*log(I)

## Model 3
D = model.matrix(~datard$sex+datard$DAS_slowd_m)
res3 = optim(par = c(1,0,1,1.5),fn = function(x){logLikel_fun(matrix(c(x[1],x[1],x[1]),J,N),exp(x[4]),D%*%c(0,x[2],x[3]),matrix(y,I,1),matrix(l,I,1),matrix(r,I,1),H=21,Tm,Dm)},hessian = TRUE,
             method = "BFGS",control = list(fnscale=-1,trace=3))
res3$par[1:3] #alpha, beta
exp(res3$par[4]) #sigma
bic3 = -2*res3$value + length(res3$par)*log(I)
se_pars = sqrt(1/diag(-res3$hessian))

## Model 4
D = model.matrix(~datard$sex*datard$DAS_slowd_m)
res4 = optim(par = c(1,0,1,1,1.5),fn = function(x){logLikel_fun(matrix(c(x[1],x[1],x[1]),J,N),exp(x[5]),D%*%c(0,x[2],x[3],x[4]),matrix(y,I,1),matrix(l,I,1),matrix(r,I,1),H=21,Tm,Dm)},hessian = FALSE,
             method = "BFGS",control = list(fnscale=-1,trace=3))
res4$par[1:4] #alpha, beta
exp(res4$par[5]) #sigma
bic4 = -2*res4$value + length(res4$par)*log(I)

## Model 5 (model 3 with a different tree)
Tm = matrix(c(0,0,NA,0,1,NA,1,NA,0,1,NA,1),M,N,byrow=TRUE) #nested response tree
rownames(Tm) = paste0(c(1:M));colnames(Tm) = paste0("node",1:N); print(Tm)
Dm=matrix(1,M,N); for(n in 1:N){Dm[is.na(Tm[,n]),n]=0}

D = model.matrix(~datard$sex+datard$DAS_slowd_m)
res5 = optim(par = c(1,0,0,1,1.5),fn = function(x){logLikel_fun(matrix(c(x[1],x[1],x[1]),J,N),exp(x[4]),D%*%c(0,x[2],x[3]),matrix(y,I,1),matrix(l,I,1),matrix(r,I,1),H=21,Tm,Dm)},hessian = FALSE,
              method = "BFGS",control = list(fnscale=-1,trace=3))
res5$par[1:4] #alpha, beta
exp(res5$par[5])
bic5 = -2*res5$value + length(res5$par)*log(I)

## Table 1
X=rbind(
  c("M1: linear tree","-",length(res1$par),round(res1$value,3),round(bic1,3)),
  c("M2: linear tree","\\texttt{sex}",length(res2$par),round(res2$value,3),round(bic2,3)),
  c("M3: linear tree","\\texttt{sex}, \\texttt{DAS}",length(res3$par),round(res3$value,3),round(bic3,3)),
  c("M4: linear tree","\\texttt{sex}, \\texttt{DAS}, \\texttt{sex:DAS}",length(res4$par),round(res4$value,3),round(bic4,3)),
  c("M5: nested tree","\\texttt{sex}, \\texttt{DAS}",length(res5$par),round(res5$value,3),round(bic5,3))
)
colnames(X) = c("Model","Covariates","$p$","L","BIC")
Xtab_tex = xtable::xtable(X,digits = rep(4,NCOL(X)+1))
attributes(Xtab_tex)$caption = "Application: "
attributes(Xtab_tex)$label = "tab1"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x},include.rownames = FALSE)

## Table 2
X=cbind(c("$\\alpha$","$\\beta_{\\texttt{sex}}$","$\\beta_{\\texttt{das}}$","$\\sigma_\\eta$"),
        c(round(res3$par[1:3],3),round(exp(res3$par[4]),3)),round(se_pars,3))
colnames(X) = c("Parameter","Estimate","Std. Error")
Xtab_tex = xtable::xtable(X,digits = rep(3,NCOL(X)+1))
attributes(Xtab_tex)$caption = "Application: "
attributes(Xtab_tex)$label = "tab2"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x},include.rownames = FALSE)


## Marginal effects
pi_y.F = mapply(function(m){prod((exp((res3$par[1]+(res3$par[3]*0))*Tm[m,])/(1+exp(res3$par[1]+(res3$par[3]*0))))^Dm[m,])},1:M)
pi_y.M = mapply(function(m){prod((exp((res3$par[1]+(res3$par[2]+res3$par[3]*0))*Tm[m,])/(1+exp(res3$par[1]+(res3$par[2]+res3$par[3]*0))))^Dm[m,])},1:M)
Out0.F=Out0.M=matrix(NA,M,3)
Out0.F[,1]=pi_y.F; Out0.M[,1]=pi_y.M
for(y0 in 1:M){
  omega_t = 1-pi_y.F[y0]
  jjd = as.numeric(1:M<y0)
  Out0.F[y0,2] = sum(pi_y.F*jjd) / omega_t
  Out0.F[y0,3] = -sum(pi_y.F*log(pi_y.F))/log(M)
  
  omega_t = 1-pi_y.M[y0]
  jjd = as.numeric(1:M<y0)
  Out0.M[y0,2] = sum(pi_y.M*jjd) / omega_t
  Out0.M[y0,3] = -sum(pi_y.M*log(pi_y.M))/log(M)
}

pi_y.F = mapply(function(m){prod((exp((res3$par[1]+(res3$par[3]*min(D[,3])))*Tm[m,])/(1+exp(res3$par[1]+(res3$par[3]*min(D[,3])))))^Dm[m,])},1:M)
pi_y.M = mapply(function(m){prod((exp((res3$par[1]+(res3$par[2]+res3$par[3]*min(D[,3])))*Tm[m,])/(1+exp(res3$par[1]+(res3$par[2]+res3$par[3]*min(D[,3])))))^Dm[m,])},1:M)
Out_min.F=Out_min.M=matrix(NA,M,3)
Out_min.F[,1]=pi_y.F; Out_min.M[,1]=pi_y.M
for(y0 in 1:M){
  omega_t = 1-pi_y.F[y0]
  jjd = as.numeric(1:M<y0)
  Out_min.F[y0,2] = sum(pi_y.F*jjd) / omega_t
  Out_min.F[y0,3] = -sum(pi_y.F*log(pi_y.F))/log(M)
  
  omega_t = 1-pi_y.M[y0]
  jjd = as.numeric(1:M<y0)
  Out_min.M[y0,2] = sum(pi_y.M*jjd) / omega_t
  Out_min.M[y0,3] = -sum(pi_y.M*log(pi_y.M))/log(M)
}

pi_y.F = mapply(function(m){prod((exp((res3$par[1]+(res3$par[3]*mean(D[,3])))*Tm[m,])/(1+exp(res3$par[1]+(res3$par[3]*mean(D[,3])))))^Dm[m,])},1:M)
pi_y.M = mapply(function(m){prod((exp((res3$par[1]+(res3$par[2]+res3$par[3]*mean(D[,3])))*Tm[m,])/(1+exp(res3$par[1]+(res3$par[2]+res3$par[3]*mean(D[,3])))))^Dm[m,])},1:M)
Out_mean.F=Out_mean.M=matrix(NA,M,3)
Out_mean.F[,1]=pi_y.F; Out_mean.M[,1]=pi_y.M
for(y0 in 1:M){
  omega_t = 1-pi_y.F[y0]
  jjd = as.numeric(1:M<y0)
  Out_mean.F[y0,2] = sum(pi_y.F*jjd) / omega_t
  Out_mean.F[y0,3] = -sum(pi_y.F*log(pi_y.F))/log(M)
  
  omega_t = 1-pi_y.M[y0]
  jjd = as.numeric(1:M<y0)
  Out_mean.M[y0,2] = sum(pi_y.M*jjd) / omega_t
  Out_mean.M[y0,3] = -sum(pi_y.M*log(pi_y.M))/log(M)
}

pi_y.F = mapply(function(m){prod((exp((res3$par[1]+(res3$par[3]*max(D[,3])))*Tm[m,])/(1+exp(res3$par[1]+(res3$par[3]*max(D[,3])))))^Dm[m,])},1:M)
pi_y.M = mapply(function(m){prod((exp((res3$par[1]+(res3$par[2]+res3$par[3]*max(D[,3])))*Tm[m,])/(1+exp(res3$par[1]+(res3$par[2]+res3$par[3]*max(D[,3])))))^Dm[m,])},1:M)
Out_max.F=Out_max.M=matrix(NA,M,3)
Out_max.F[,1]=pi_y.F; Out_max.M[,1]=pi_y.M
for(y0 in 1:M){
  omega_t = 1-pi_y.F[y0]
  jjd = as.numeric(1:M<y0)
  Out_max.F[y0,2] = sum(pi_y.F*jjd) / omega_t
  Out_max.F[y0,3] = -sum(pi_y.F*log(pi_y.F))/log(M)
  
  omega_t = 1-pi_y.M[y0]
  jjd = as.numeric(1:M<y0)
  Out_max.M[y0,2] = sum(pi_y.M*jjd) / omega_t
  Out_max.M[y0,3] = -sum(pi_y.M*log(pi_y.M))/log(M)
}

## Figure 3
#x11()
tikzDevice::tikz(file='fig3.tex',width=6.5,height=5.5)
par(mfrow=c(2,2),mai=c(1.15, 0.85, 0.45, 0.15))
thr=0.25; cols = c("orangered3","orchid4","palegreen4","peru")
# row=1,column=1
plot(1:M,Out0.M[,1],bty="n",xlab="",ylab=TeX("$\\pi^y$"),ylim=c(0,0.9),xlim=c(1,M+0.75),axes = FALSE,pch=20,lwd=3,col=cols[1],cex.axis=2,cex.lab=1.25,main="sex = M")
segments(x0 = 1:M,x1 = 1:M,y0 = 0,y1 = Out0.M[,1],lty = 1,lwd = 2,col=cols[1])
#axis(side = 1,at = seq(1,M+1,length=5)+((thr*3)/2),labels = c(1,2,3,4,"")); 
axis(side = 2,at = seq(0,1,length=5),labels = seq(0,1,length=5))
points(1:M+thr,Out_min.M[,1],bty="n",xlab="",ylab=TeX("$\\pi^y$"),ylim=c(0,1),pch=20,lwd=3,col=cols[2])
segments(x0 = 1:M+thr,x1 = 1:M+thr,y0 = 0,y1 = Out_min.M[,1],lty = 1,lwd = 2,col=cols[2])
points(1:M+thr*2,Out_mean.M[,1],bty="n",xlab="",ylab=TeX("$\\pi^y$"),ylim=c(0,1),pch=20,lwd=3,col=cols[3])
segments(x0 = 1:M+thr*2,x1 = 1:M+thr*2,y0 = 0,y1 = Out_mean.M[,1],lty = 1,lwd = 2,col=cols[3])
points(1:M+thr*3,Out_max.M[,1],bty="n",xlab="",ylab=TeX("$\\pi^y$"),ylim=c(0,1),pch=20,lwd=3,col=cols[4])
segments(x0 = 1:M+thr*3,x1 = 1:M+thr*3,y0 = 0,y1 = Out_max.M[,1],lty = 1,lwd = 2,col=cols[4])
abline(v = c(1.9,2.9,3.9,4.9),lty=2); text(x=1:M+((thr*3)/2),y=0.85,c("m=1","m=2","m=3","m=4"),cex = 1.25)
# row=1,column=2
plot(1:M,Out0.F[,1],bty="n",xlab="",ylab="",ylim=c(0,0.9),xlim=c(1,M+0.75),axes = FALSE,pch=20,lwd=3,col=cols[1],cex.axis=2,cex.lab=1.25,main="sex = F")
segments(x0 = 1:M,x1 = 1:M,y0 = 0,y1 = Out0.F[,1],lty = 1,lwd = 2,col=cols[1])
#axis(side = 1,at = seq(1,M+1,length=5)+((thr*3)/2),labels = c(1,2,3,4,"")); 
axis(side = 2,at = seq(0,1,length=5),labels = seq(0,1,length=5))
points(1:M+thr,Out_min.F[,1],bty="n",xlab="",ylab=TeX("$\\pi^y$"),ylim=c(0,1),pch=20,lwd=3,col=cols[2])
segments(x0 = 1:M+thr,x1 = 1:M+thr,y0 = 0,y1 = Out_min.F[,1],lty = 1,lwd = 2,col=cols[2])
points(1:M+thr*2,Out_mean.F[,1],bty="n",xlab="",ylab=TeX("$\\pi^y$"),ylim=c(0,1),pch=20,lwd=3,col=cols[3])
segments(x0 = 1:M+thr*2,x1 = 1:M+thr*2,y0 = 0,y1 = Out_mean.F[,1],lty = 1,lwd = 2,col=cols[3])
points(1:M+thr*3,Out_max.F[,1],bty="n",xlab="",ylab=TeX("$\\pi^y$"),ylim=c(0,1),pch=20,lwd=3,col=cols[4])
segments(x0 = 1:M+thr*3,x1 = 1:M+thr*3,y0 = 0,y1 = Out_max.F[,1],lty = 1,lwd = 2,col=cols[4])
abline(v = c(1.9,2.9,3.9,4.9),lty=2); text(x=1:M+((thr*3)/2),y=0.85,c("m=1","m=2","m=3","m=4"),cex = 1.25)
# row=2,column=1
plot(1:M,Out0.M[,2],bty="n",xlab="",ylab=TeX("$\\pi^s$"),ylim=c(0,1.20),xlim=c(1,M+0.75),axes = FALSE,pch=20,lwd=3,col=cols[1],cex.axis=2,cex.lab=1.25)
segments(x0 = 1:M,x1 = 1:M,y0 = 0,y1 = Out0.M[,2],lty = 1,lwd = 2,col=cols[1])
#axis(side = 1,at = seq(1,M+1,length=5)+((thr*3)/2),labels = c(1,2,3,4,"")); 
axis(side = 2,at = seq(0,1,length=5),labels = seq(0,1,length=5))
points(1:M+thr,Out_min.M[,2],bty="n",xlab="",ylab=TeX("$\\pi^y$"),ylim=c(0,1),pch=20,lwd=3,col=cols[2])
segments(x0 = 1:M+thr,x1 = 1:M+thr,y0 = 0,y1 = Out_min.M[,2],lty = 1,lwd = 2,col=cols[2])
points(1:M+thr*2,Out_mean.M[,2],bty="n",xlab="",ylab=TeX("$\\pi^y$"),ylim=c(0,1),pch=20,lwd=3,col=cols[3])
segments(x0 = 1:M+thr*2,x1 = 1:M+thr*2,y0 = 0,y1 = Out_mean.M[,2],lty = 1,lwd = 2,col=cols[3])
points(1:M+thr*3,Out_max.M[,2],bty="n",xlab="",ylab=TeX("$\\pi^y$"),ylim=c(0,1),pch=20,lwd=3,col=cols[4])
segments(x0 = 1:M+thr*3,x1 = 1:M+thr*3,y0 = 0,y1 = Out_max.M[,2],lty = 1,lwd = 2,col=cols[4])
abline(v = c(1.9,2.9,3.9,4.9),lty=2); text(x=1:M+((thr*3)/2),y=1.20,c("m=1","m=2","m=3","m=4"),cex = 1.25)
# row=2,column=2
plot(1:M,Out0.F[,2],bty="n",xlab="",ylab="",ylim=c(0,1.20),xlim=c(1,M+0.75),axes = FALSE,pch=20,lwd=3,col=cols[1],cex.axis=2,cex.lab=1.25)
segments(x0 = 1:M,x1 = 1:M,y0 = 0,y1 = Out0.F[,2],lty = 1,lwd = 2,col=cols[1])
#axis(side = 1,at = seq(1,M+1,length=5)+((thr*3)/2),labels = c(1,2,3,4,"")); 
axis(side = 2,at = seq(0,1,length=5),labels = seq(0,1,length=5))
points(1:M+thr,Out_min.F[,2],bty="n",xlab="",ylab=TeX("$\\pi^y$"),ylim=c(0,1),pch=20,lwd=3,col=cols[2])
segments(x0 = 1:M+thr,x1 = 1:M+thr,y0 = 0,y1 = Out_min.F[,2],lty = 1,lwd = 2,col=cols[2])
points(1:M+thr*2,Out_mean.F[,2],bty="n",xlab="",ylab=TeX("$\\pi^y$"),ylim=c(0,1),pch=20,lwd=3,col=cols[3])
segments(x0 = 1:M+thr*2,x1 = 1:M+thr*2,y0 = 0,y1 = Out_mean.F[,2],lty = 1,lwd = 2,col=cols[3])
points(1:M+thr*3,Out_max.F[,2],bty="n",xlab="",ylab=TeX("$\\pi^y$"),ylim=c(0,1),pch=20,lwd=3,col=cols[4])
segments(x0 = 1:M+thr*3,x1 = 1:M+thr*3,y0 = 0,y1 = Out_max.F[,2],lty = 1,lwd = 2,col=cols[4])
abline(v = c(1.9,2.9,3.9,4.9),lty=2); text(x=1:M+((thr*3)/2),y=1.20,c("m=1","m=2","m=3","m=4"),cex = 1.25)

add_legend("bottom",fill = cols,legend = c("das=0","min(das)","mean(das)","max(das)"),border = FALSE,bty = "n",ncol = 4,cex=1.5)
dev.off()


## Figure 4
cols = c("skyblue4","tomato4")
#x11()
tikzDevice::tikz(file='fig4.tex',width=4.5,height=3.5)
par(mfrow=c(1,1),mai=c(1.15, 0.85, 0.45, 0.15))
plot(x = c(0,min(D[,3]),mean(D[,3]),max(D[,3])),
     y = c(Out0.M[1,3],Out_min.M[1,3],Out_mean.M[1,3],Out_max.M[1,3]),
     bty="n",type="b",xlab="das",ylab=TeX("$\\xi$"),ylim=c(0.35,0.99),lwd=2,col=cols[1])
points(x = c(0,min(D[,3]),mean(D[,3]),max(D[,3])),
       y = c(Out0.F[1,3],Out_min.F[1,3],Out_mean.F[1,3],Out_max.F[1,3]),type="b",lwd=2,col=cols[2])

add_legend("bottom",fill = cols,legend = c("sex=M","sex=F"),border = FALSE,bty = "n",ncol = 2,cex=1.25)
dev.off()






