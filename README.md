[<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/banner.png" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **InternerZins** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)


```yaml
Name of Quantlet: InternerZins 

Description: 'Approximiert den internen Zins in Aufgabe 7.' 

Author:   Awdesch Melzer 

Example: 'Interner Zins' 

```
![Picture1](08.jpeg)

```R
# WERMssvdeFASTEC

# clear history
rm(list = ls(all = TRUE))
graphics.off()

setwd(Sys.glob("~/desktop/WERM/WERM_MERRA")) 

# install and load packages
libraries = c("latticeExtra","gplots","splines", "maps","mapdata","maptools","scales","marmap","mapproj", "matrixStats","RNetCDF","zoo","lattice","viridis","ncdf4","akima","foreach","parallel","snow","locpol","nlmrt","doParallel","ssvd","np","iterators","softImpute")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

## Function: FISTA algorithm with ALS-SVD
mer = function(Y, X, tau, lambda, epsilon = 10^(-6), itt = 2000) {
  ## Initialize
  m       = ncol(Y)
  n       = nrow(Y)
  p       = ncol(X)
  L       = 2 * (m * n)^(-1) * max(tau, (1 - tau)) * norm(X, type = "F")^2
  Omega   = matrix(0, nrow = p, ncol = m)
  delta   = 1  # step size
  error   = 1e+07
  L.error = 1e+10
  it      = 1
  ## Output
  A       = matrix(0, nrow = p, ncol = m)
  A_      = matrix(0, nrow = p, ncol = m)
  ## Main iteration
  while(it<itt & error > epsilon){
    sS         = svd.als(Omega - L^(-1) * G.er.new(Omega, Y, X, tau), rank.max=5)
 #   sS         = svd.als(Omega - L^(-1) * G.qr(Omega, Y, X, tau,kappa=1e-04,m,n), rank.max=5)
    temp.sv   = sS$d - (lambda/L)
    temp.sv[temp.sv < 0] = 0
    A         = sS$u %*% diag(temp.sv) %*% t(sS$v)
    delta.for = (1 + sqrt(1 + 4 * delta^2))/2
    Omega     = A + (delta - 1)/(delta.for) * (A - A_)
    error     = L.error - (sum((tau - matrix(as.numeric(Y - X %*% A < 0), n, m)) * (Y - X %*% A)^2) + lambda * sum(temp.sv))
    L.error   = sum((tau - matrix(as.numeric(Y - X %*% A < 0), n, m)) * (Y - X %*%A)^2) + lambda * sum(temp.sv)
    A_        = A
    delta     = delta.for
    it        = it + 1
    print(it)
    print(c(error, delta, sum((tau - matrix(as.numeric(Y - X %*% A < 0), n, m)) * (Y - X %*% A)^2), sum(temp.sv)))
    # if(it < 10){error=1000000}
  }
  list(Gamma = A, d = sS$d, U = sS$u, V = sS$v, error = error, loss = sum((tau - matrix(as.numeric(Y - X %*% A < 0), n, m)) * (Y - X %*% A)^2), norm = sum(temp.sv), lambda = lambda, 
       iteration = it)
}
## Function: Computing the gradient of the expectile loss function
  G.er.new = function(A, Y, X, tau) {
  m = ncol(Y)
  n = nrow(Y)
  p = ncol(X)
#  w = 0
  G = matrix(0, p, m)
  u   = Y - X%*%A
  w   = matrix(0,nrow(u),ncol(u))
  w[which(u > 0,arr.ind=T)] = 2 * tau*  u[which(u > 0,arr.ind=T)]
  w[which(u <= 0,arr.ind=T)] = 2 * (1 - tau)*  u[which(u <= 0,arr.ind=T)]
  G = (m * n)^(-1) * -t(X)%*%w
  return(G)
  }
  ## Function: Computing the gradient of the expectile loss function
G.er.old = function(A, Y, X, tau) {
    m = ncol(Y)
    n = nrow(Y)
    p = ncol(X)
    w = 0
    G = matrix(0, p, m)
    for (r in 1:p) {
        for (s in 1:m) {
            for (i in 1:n) {
                u   = Y[i, s] - X[i, ] %*% A[, s]
                w   = w + (if (u > 0) {
                  2 * tau * u
                } else {
                  2 * (1 - tau) * u
                }) * X[i, r]
            }
            G[r, s] = (m * n)^(-1) * w
            w       = 0
        }
    }
    G
}

########

## Function: FISTA algorithm with sSVD
mer.ssvd = function(Y, X, tau, lambda, epsilon = 10^(-6), itt = 2000,r=10) {
  ## Initialize
  m       = ncol(Y)
  n       = nrow(Y)
  p       = ncol(X)
  L       = 2 * (m * n)^(-1) * max(tau, (1 - tau)) * norm(X, type = "F")^2
  Omega   = matrix(0, nrow = p, ncol = m)
  delta   = 1  # step size
  error   = 1e+07
  L.error = 1e+10
  it      = 1
  ## Output
  A       = matrix(0, nrow = p, ncol = m)
  A_      = matrix(0, nrow = p, ncol = m)
  ## Main iteration
  while(it<itt & error > epsilon){
    sS         = ssvd(Omega - L^(-1) * G.er.new(Omega, Y, X, tau), r=r)
    temp.sv   = sS$d - (lambda/L)
    temp.sv[temp.sv < 0] = 0
    A         = sS$u %*% diag(temp.sv) %*% t(sS$v)
    delta.for = (1 + sqrt(1 + 4 * delta^2))/2
    Omega     = A + (delta - 1)/(delta.for) * (A - A_)
    error     = L.error - (sum((tau - matrix(as.numeric(Y - X %*% A < 0), n, m)) * 
                                 (Y - X %*% A)^2) + lambda * sum(temp.sv))
    L.error   = sum((tau - matrix(as.numeric(Y - X %*% A < 0), n, m)) * (Y - X %*% 
                                                                           A)^2) + lambda * sum(temp.sv)
    A_        = A
    delta     = delta.for
    it        = it + 1
    print(it)
    print(c(error, delta, sum((tau - matrix(as.numeric(Y - X %*% A < 0), n, m)) * (Y - X %*% A)^2), sum(temp.sv)))
    # if(it < 10){error=1000000}
  }
  list(Gamma = A, d = sS$d, U = sS$u, V = sS$v, error = error, loss = sum((tau - matrix(as.numeric(Y - X %*% A < 0), n, m)) * (Y - X %*% A)^2), norm = sum(temp.sv), lambda = lambda, 
       iteration = it)
}

## Function: FISTA algorithm SVD
mer.svd = function(Y, X, tau, lambda, epsilon = 10^(-6), itt = 2000,r=10) {
  ## Initialize
  m       = ncol(Y)
  n       = nrow(Y)
  p       = ncol(X)
  L       = 2 * (m * n)^(-1) * max(tau, (1 - tau)) * norm(X, type = "F")^2
  Omega   = matrix(0, nrow = p, ncol = m)
  delta   = 1  # step size
  error   = 1e+07
  L.error = 1e+10
  it      = 1
  ## Output
  A       = matrix(0, nrow = p, ncol = m)
  A_      = matrix(0, nrow = p, ncol = m)
  ## Main iteration
  while(it<itt & error > epsilon){
    sS         = svd(Omega - L^(-1) * G.er.new(Omega, Y, X, tau), nu = p, nv = m)
    temp.sv   = sS$d - (lambda/L)
    temp.sv[temp.sv < 0] = 0
    A         = sS$u %*% diag(temp.sv,nrow=p,ncol=m) %*% t(sS$v)
    delta.for = (1 + sqrt(1 + 4 * delta^2))/2
    Omega     = A + (delta - 1)/(delta.for) * (A - A_)
    error     = L.error - (sum((tau - matrix(as.numeric(Y - X %*% A < 0), n, m)) * 
                                 (Y - X %*% A)^2) + lambda * sum(temp.sv))
    L.error   = sum((tau - matrix(as.numeric(Y - X %*% A < 0), n, m)) * (Y - X %*% 
                                                                           A)^2) + lambda * sum(temp.sv)
    A_        = A
    delta     = delta.for
    it        = it + 1
    print(it)
    print(c(error, delta, sum((tau - matrix(as.numeric(Y - X %*% A < 0), n, m)) * (Y - X %*% A)^2), sum(temp.sv)))
    # if(it < 10){error=1000000}
  }
  list(Gamma = A, d = sS$d, U = sS$u, V = sS$v, error = error, loss = sum((tau - matrix(as.numeric(Y - X %*% A < 0), n, m)) * (Y - X %*% A)^2), norm = sum(temp.sv), lambda = lambda, 
       iteration = it)
}

#######  

## Function: FISTA algorithm
mer.old = function(Y, X, tau, lambda, epsilon = 10^(-6), itt = 2000) {
    ## Initialize
    m = ncol(Y)
    n = nrow(Y)
    p = ncol(X)
    w = 0
    D = matrix(0, p, m)
    for (r in 1:p) {
        for (s in 1:m) {
            for (i in 1:n) {
                w = w + max(2 * tau, 2 * (1 - tau)) * X[i, r]^2
            }
            D[r, ] = rep((m * n)^(-1) * w, 1, m)
            w = 0
        }
    }
    L       = norm(D, type = "E")
    Omega   = matrix(0, nrow = p, ncol = m)
    delta   = 1  # step size
    error   = 1e+07
    L.error = 1e+10
    it      = 1
    ## Output
    A       = matrix(0, nrow = p, ncol = m)
    A_      = matrix(0, nrow = p, ncol = m)
    ## Main iteration
    while (it < itt & error > epsilon) {
        S         = svd(Omega - L^(-1) * G.er.old(Omega, Y, X, tau), nu = p, nv = m)
        temp.sv   = S$d - (lambda/L)
        temp.sv[temp.sv < 0] = 0
        A         = S$u %*% diag(temp.sv, nrow = p, ncol = m) %*% t(S$v)
        delta.for = (1 + sqrt(1 + 4 * delta^2))/2
        Omega     = A + (delta - 1)/(delta.for) * (A - A_)
        error     = L.error - (sum((tau - matrix(as.numeric(Y - X %*% A < 0), n, m)) * 
            (Y - X %*% A)^2) + lambda * sum(temp.sv))
        L.error   = sum((tau - matrix(as.numeric(Y - X %*% A < 0), n, m)) * (Y - X %*% 
            A)^2) + lambda * sum(temp.sv)
        A_        = A
        delta     = delta.for
        it        = it + 1
        print(c(error, delta, sum((tau - matrix(as.numeric(Y - X %*% A < 0), n, m)) * 
            (Y - X %*% A)^2), sum(temp.sv)))
            print(it)
        # if(it < 10){error=1000000}
    }
    list(Gamma = A, d = S$d, U = S$u, V = S$v, error = error, loss = sum((tau - matrix(as.numeric(Y - 
        X %*% A < 0), n, m)) * (Y - X %*% A)^2), norm = sum(temp.sv), lambda = lambda, 
        iteration = it)
}

###########

setwd(Sys.glob("~/desktop/WERM/WERM_MERRA")) 
load("DElat.RData")
load("DElon.RData")

matrix.lat = rep(DE.lat,length(DE.lon))
matrix.lon = rep(DE.lon,each=length(DE.lat))

# load("Yseas.RData")
# load("Xseas.RData")
load("resid.RData")
#load("resx.RData")
Y = resid;rm(resid)
#X = res.x;rm(res.x)

######

set.seed(1)
sam=1:4650
#sam = sample(1:4650,1000,replace=F)#1:4650
YY      = Y[,sam];rm(Y)  # transform from data frame to matrix
#XX      = X[,sam];rm(X)
load("W100M.RData")
W100 = array(0,c(dim(W100M)[1],dim(W100M)[2],dim(W100M)[3]*dim(W100M)[4]))
for(i in 1:dim(W100M)[1]){
	for(j in 1:dim(W100M)[2]){
		W100[i,j,] = c(W100M[i,j,,])
	}
}
W100M = matrix(0,365*24,dim(W100)[3])
 for(i in 1:dim(W100)[3]){
	 W100M[,i] = matrix(W100[,,i],nrow=365*24,ncol=1,byrow=T)
 }
W100M = W100M[,sam]
tau1 = sum(W100M<4)/length(c(W100M))
tau3 = 1-sum(W100M>26)/length(c(W100M))
rm(W100M)

setwd(Sys.glob("~/desktop/WERM/WERM_ERA"))
# load("interpolt2m.RData")
# X.i = interpol.t2m
# t2m = array(0,c(dim(X.i)[1],dim(X.i)[2],dim(X.i)[3]*dim(X.i)[4]))
# for(i in 1:dim(X.i)[1]){
	# for(j in 1:dim(X.i)[2]){
		# t2m[i,j,] = c(X.i[i,j,,])
	# }
# }
# T2M = matrix(0,365*24,dim(t2m)[3])
 # for(i in 1:dim(t2m)[3]){
	 # T2M[,i] = matrix(t2m[,,i],nrow=365*24,ncol=1,byrow=T)
 # }
# T2 = T2M[,sam]
	  # bwreg.f     = list() 
	  # t2m.res = matrix(0,8760,length(sam))
# for(i in 1:length(sam)){
# # bw[[i]]        = npregbw(xdat=1:8760, ydat=c(T2[,i]), regtype="ll", bwmethod="cv.aic")
 # bwreg.f[[i]]     = npreg(txdat=1:8760, tydat=c(T2[,i]),bws=1.423915,ckertype="epanechnikov")
 # t2m.res[,i] = T2[,i] - bwreg.f[[i]]$mean
 # print(i)
# }
# save(t2m.res, file="t2mres.RData")
# save(T2,file="T2.RData")
n      = nrow(YY)
p      = ceiling(365^0.4)
xx     = seq(0, 1, length = 365)
X.fac  = bs(xx, df = p, intercept = TRUE)
X      = matrix(0, nrow = n, ncol = 0)
for (i in 1:p) {
    X  = cbind(X, X.fac[, i])
}

#load("t2mres.RData")
load("T2.RData")

#XX = cbind(X, t2m.res[,sam]); rm(t2m.res)
XX = cbind(X, T2[,sam]); rm(T2)
n = nrow(YY)
p = ncol(XX)
m = dim(YY)[2]
sig_x  = sqrt(norm(XX, type = "F")/n)



alpha            = 0.1
sim.lambda.tau1  = numeric(0)
sim.lambda.tau3  = numeric(0)
B                = 1000
set.seed(1001)

# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)
registerDoParallel(cl)
tau = c(tau1,tau3)

sim.lambda.tau1 <- foreach (i = 1:B,.combine="c")%dopar%{
  W.temp.tau1         = matrix(as.numeric(runif(n * m) < tau1) - tau1, nrow = n, ncol = m)
  temp.score.tau1     = norm(t(XX) %*% W.temp.tau1/n, type = "E")
  temp.score.tau1/m
}
lamb1            = 2 * quantile(sim.lambda.tau1, p = 0.9)
sim.lambda.tau3 <- foreach (i = 1:B,.combine="c")%dopar%{
  W.temp.tau3         = matrix(as.numeric(runif(n * m) < tau3) - tau3, nrow = n, ncol = m)
  temp.score.tau3     = norm(t(XX) %*% W.temp.tau3/n, type = "E")
  temp.score.tau3/m
}
lamb3            = 2 * quantile(sim.lambda.tau3, p = 0.9)

stopCluster(cl)
###########
lambdas = c(lamb1, lamb3)
#save(lambdas,file="lambdas.RData")

# ALS-SVD
system.time(fit11.1000 <- mer(YY,XX,tau=tau[1],lambda=lamb1,itt=1000))
#    user  system elapsed 
#   1.069   0.077   1.141
#     user   system  elapsed 
# 2827.745   22.365 2865.702 
save(fit11.1000, file="fit111000.RData")#tau = 0.069
system.time(fit12.1000 <- mer(YY,XX,tau=tau[2],lambda=lamb3,itt=1000))
#    user  system elapsed 
#  41.765   4.602  46.096 
#      user    system   elapsed 
# 11997.364   114.758 25025.501 
save(fit12.1000, file="fit121000.RData") # tau = 0.867
### ssvd
system.time(fit21.1000 <- mer.ssvd(YY,XX,tau=tau[1],lambda=lamb1,itt=1000,r=10))
#    user  system elapsed 
#   1.081   0.078   1.152
#     user   system  elapsed 
# 2781.837   14.415 2798.319 
save(fit21.1000, file="fit211000.RData")  
system.time(fit22.1000 <- mer.ssvd(YY,XX,tau=tau[2],lambda=lamb3,itt=1000,r=10))
#    user  system elapsed 
#  32.089   3.568  35.439
#      user    system   elapsed 
# 11342.515    74.216 11377.154  
save(fit22.1000, file="fit221000.RData")
#### new, standard svd
 system.time(fit31.1000 <- mer.svd(YY,XX,tau=tau[1],lambda=lamb1,itt=1000))
#    user  system elapsed 
#   0.637   0.077   0.708
#     user   system  elapsed 
# 3432.875   18.930 3442.779 
save(fit31.1000, file="fit311000.RData")
 system.time(fit32.1000 <- mer.svd(YY,XX,tau=tau[2],lambda=lamb3,itt=1000))
#    user  system elapsed 
#   2.620   0.393   2.994
#      user    system   elapsed 
# 14140.207    75.971 14178.131 
save(fit32.1000, file="fit321000.RData")

# load("fit111000.RData")
# load("fit121000.RData")

pc1            = XX %*% fit11.1000$U %*% diag(fit11.1000$d)
pc3            = XX %*% fit12.1000$U %*% diag(fit12.1000$d)
varipc1        = colVars(pc1)
varipc3        = colVars(pc3)

score1 = fit11.1000$U %*% diag(fit11.1000$d)
score1.ord = score1[,order(varipc1,decreasing=T)]
score1.ord = score1.ord[order(score1.ord[,1],score1.ord[,2],decreasing=T),]
score3 = fit12.1000$U %*% diag(fit12.1000$d)
score3.ord = score3[,order(varipc3,decreasing=T)]
score3.ord = score3.ord[order(score3.ord[,1],score3.ord[,2],decreasing=T),]

load1     = (fit11.1000$V[, order(varipc1, decreasing = TRUE)])
load3     = (fit12.1000$V[, order(varipc3, decreasing = TRUE)])

xx = 1:8760
days = 365
samp = 1:(24*days)
col = c("blue3","red3")#viridis(2)
par(mfrow      = c(2, 2),         # 2x2 layout
    oma        = c(0, 1, 0, 0),   # two rows of text at the outer left and bottom margin
    mar        = c(2, 1.8, 1, 0), # space for one row of text at ticks and to separate plots
    mgp        = c(1.5, 0.5, 0),  # axis label at 2 rows distance, tick labels at 1 row
    xpd        = NA )             # allow content to protrude into outer margin (and beyond)
plot(xx[samp], pc1[samp, 1], type = "l", lwd = 1.2, main = "1st factor", xlab = "", ylab = "", 
  ylim = c(min(pc1[samp, 1],pc3[samp, 1]), max(pc1[samp, 1],pc3[samp, 1])),col=col[2] )
lines(xx[samp], pc3[samp, 1], lty = 1, lwd = 1, col = col[1])
lines(xx[samp], pc1[samp, 1], lty = 1, lwd = 1, col = col[2])

plot(xx[samp], pc1[samp, 2], type = "l", lwd = 1.2, main = "2nd factor", xlab = "", ylab = "", 
  ylim = c(min(pc1[samp, 2],pc3[samp, 2]), max(pc1[samp, 2],pc3[samp, 2])),col=col[2] )
lines(xx[samp], pc3[samp, 2], lty = 1, lwd = 1,col=col[1] )
lines(xx[samp], pc1[samp, 2], lty = 1, lwd = 1, col = col[2])

plot(xx[samp], pc1[samp, 3], type = "l", lwd = 1.2, main = "3rd factor", xlab = "", ylab = "", 
  ylim = c(min(pc1[samp, 3],pc3[samp, 3]), max(pc1[samp, 3],pc3[samp, 3])),col=col[2] )
lines(xx[samp], pc3[samp, 3], lty = 1, lwd = 1,col=col[1] )
lines(xx[samp], pc1[samp, 3], lty = 1, lwd = 1, col = col[2])

plot(xx[samp], pc1[samp, 4], type = "l", lwd = 1.2, main = "4th factor", xlab = "", ylab = "", 
  ylim = c(min(pc1[, 4],pc3[samp, 4]), max(pc1[samp, 4],pc3[samp, 4])),col=col[2] )
lines(xx[samp], pc3[samp, 4], lty = 1, lwd = 1.2, col=col[1] )
lines(xx[samp], pc1[samp, 4], lty = 1, lwd = 1, col = col[2])

pdf("WERMalsWPDpec.pdf",width=8,height=5)
par(mfrow = c(2, 1),           # 2x2 layout
rm(list=ls(all=T))
graphics.off()


epsilon = 0.05                # Größe des erlaubten Fehlers 
                             # (entscheidet über Anzahl der Iterationen)
i       = seq(0,20,0.01)/100 # Sequenz des Zinses
RBF     = 1/i*(1-1/(1+i)^15) # Rentenbarwertfaktor
C.0     = -8 + RBF - 2/(1+i)^5 - 2/(1+i)^10 + 1.5/(1+i)^15 # Kapitalwert
C.0[1]  = -8 + 4*1-1+4*1-1+4*1+2.5 # Kapitalwert für i = 0

plot(i, C.0, type="l", xlab="Interner Zins",
ylab="Kapitalwert", lwd=4, cex.axis=1.4, cex.lab=1.4, col="blue3")
abline(h=0, lwd=3)

x1 = 0                        # Anfangswert i=0
y1 = C.0[which(i==x1)]
x2 = 0.08                     # Anfangswert i=0.08
y2 = C.0[which(i==x2)]
text(x1,-0.4,"x1",lwd=2,cex=1.3)
abline(v=x1,lty=3,col="darkgrey",lwd=1)
text(x2,-0.4,"x2",lwd=2,cex=1.3)
abline(v=x2,lty=3,col="darkgrey",lwd=1)

a = C.0[1]                    # Schnittpunkt mit y-Achse
b = (C.0[which(i==0.08)]-a)/i[which(i==0.08)] # Steigung
abline(a=a, b=b,col="red3",lwd=3,lty=2)

x3 = x1 - y1*((x2-x1)/(y2-y1)) # 1. Approximation
text(0.15,2,paste("1. Approximation: i* =",round(x3,4)),lwd=2,cex=1.3)
x3
y3 = C.0[which(i==round(x3,4))] # Kapitalwert an der 1. Approximation
text(x3,-0.4,"x3",lwd=2,cex=1.3)
abline(v=x3,lty=3,col="darkgrey",lwd=1)

a = (y3*x2 -y2*x3)/(x2-x3)    # neuer Schnittpunkt mit y-Achse
b = (y2 - a)/x2               # neue Steigung
abline(a=a, b=b,col="red3",lwd=3,lty=2)

x4 = x3 - y3*((x2-x3)/(y2-y3)) # 2. Approximation
text(0.15,1.5,paste("2. Approximation: i* =",round(x4,digits=4)),lwd=2,cex=1.3)
y4 = C.0[which(as.character(i)==as.character(round(x4,digits=4)))]

text(x4,-0.4,"x4",lwd=2,cex=1.3)
abline(v=x4,lty=3,col="darkgrey",lwd=1)

if(abs(y4)>epsilon){ # Teste, ob weitere Approximation notwendig
a = (y4*x3 -y3*x4)/(x3-x4)
b = (y3 - a)/x3#(C.0[which(i==i.star)]-a)/i[which(i==i.star)]
abline(a=a, b=b,col="red3",lwd=3,lty=2)
x5 <- x4 - y4*((x3-x4)/(y3-y4))
text(0.15,1,paste("3. Approximation: i* =",round(multiroot(function(x){a+b*x},start=0.08)$root,digits=4)),lwd=2,cex=1.3)
}
```
