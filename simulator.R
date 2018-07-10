##### Un archivo temporario para hacer una simulacion

#### La entrada de los datos
I <- 322                                 # numero de los individuos
K <- 18                                   # numero de los periodos
nfam <- 72                                 # numero de las familias

## make sure each one is represented at least once
family <- (c(1:nfam, sample(1:nfam, size=I-nfam, replace=TRUE)))
tmp <- rnorm(1:nfam)
## family values are repeated
carez <- tmp[family]
## age of fledling, centered
year <- sort(sample(1:4, size=I, replace=TRUE))
agec <- 1:K-mean(1:K)

## los parametros de usar en la simulacion, por ahora son iguales que el
## adjuste de los datos reales.
sigmayearphi=1.06
sigmaphi=-.77
sigmap=-.25
a=c(1.04, 0.98, 0.41, 0.43, 0.58, 0.97, 0.59, 0.84, 0.38, -0.53,
    1.43, -0.63, -0.57, 0.13, -0.27, -1.44, -1.28)
a1=0.68
b0=c(2.15, 1.93, 3.39, 3.33)
b1=c(0.05, -.3, -.4, -0.06)
## effectos aleatorios
fameffp_raw <-
  c(0.47, -0.86, 1.22, -0.53, -0.36, 1.11, -1.2, -0.81, 0.35, 0.22,
    -0.12, -0.17, -0.88, 1.04, -1.16, 0.22, 1.42, -0.43, 0.16, -0.01,
    0.98, -1.19, -0.02, -0.15, 0.11, 0.28, -1.01, -0.32, 0.17, 0.9,
    1.01, 0.46, -1.35, 0.03, 1.31, -1.36, -0.48, -0.56, -1.4, 0.11,
    0.05, -0.38, 0.1, 1.01, -0.89, 0.92, -0.95, 0.94, 0.74, -0.62,
    0.4, 1.04, 0.42, -0.86, 1.53, -0.83, -0.21, 0.03, -1.01, -0.36,
    -1.02, 0.37, 0.25, -0.15, -0.52, -0.2, 0.21, -0.55, -1.2, 0.11,
    0.6, 0.41)
fameffphi_raw <-
  c(0.43, 0.14, 0.03, -0.08, 0.18, 0.28, -0.03, 0.36, -0.06, -0.49,
    -0.65, 0.35, 1.03, -0.45, 1.06, -0.83, -0.42, -0.54, -1.18, -0.03,
    -0.12, 0.9, 0.99, -0.41, 0.17, -0.61, 0.31, -0.42, 0.06, 0.34,
    -0.77, 0.95, -0.65, 0.5, 0.62, -0.28, -0.2, -0.2, 0.57, 0.24,
    0.52, 0.89, 0.22, -0.19, -0.61, -0.67, -0.38, -0.7, 0.27, 0.14,
    0.48, 0.34, -0.4, 0.34, 0.36, -0.55, 0.2, -0.22, -0.57, -0.19,
    0.07, -0.19, 0.16, -0.51, 0.61, 0.5, -1.31, 0.01, 0.48, -0.36,
    1.31, -0.54)
yeareffphi_raw <- c(0.78, 0.97, 1.1, 0.98)
## componente de observacion
CH <- matrix(NA, I,K)
last <- rep(NA, len=I)

p <- matrix(NA, nrow=I, ncol=K) # capture probability
phi <- matrix(NA, I, K-1);   # survival probability
sigmayearphi2 <- exp(sigmayearphi)
sigmaphi2 <- exp(sigmaphi)
sigmap2 <- exp(sigmap)
inv_logit <- function(x) 1/(1+exp(-x))
alive <- CH <- p
last <- rep(NA, len=I)

## Loop over each individual and calculate the survival, detection
## probabilities and then do the simulated observation
for(i in 1:I){
  ## calculate phi as a function of fixed and random effects
  for(t in 1:(K-1)) {
    x <- a[t]+ a1*carez[i]+
      sigmayearphi2*yeareffphi_raw[year[i]]+
      sigmaphi2*fameffphi_raw[family[i]];
    phi[i,t] <- inv_logit(x);
  }
  ## calculate p as a function of fixed and random effects
  p[i,1] = 1;  ## first occasion is marking occasion
  for(t in 2:K){
    x <- b0[year[i]]+ b1[year[i]]*agec[t]+
      sigmap2*fameffp_raw[family[i]];
    p[i,t] <- inv_logit(x);
  }
  ## Observation component
  CH[i,1] <- alive[i,1] <- 1 ## alive and observed at first marking
  for(t in 2:K){
    ## only alive this period if alive last period and survied last period
    alive[i,t] <- alive[i,t-1]*rbinom(n=1, size=1, prob=phi[i,t-1])
    ## observed?
    CH[i,t] <- alive[i,t]*rbinom(n=1, size=1, prob=p[i,t-1])
  }
  ## last period animal seen
  last[i] <- tail(which(CH[i,]==1), n=1)
}

simdata <- list(I=I, K=K, nfam=nfam, CH=CH, carez=carez,
                year=year, agec=agec, family=family, last=last)

simpars <- list(sigmayearphi=sigmayearphi, sigmaphi=sigmaphi,
                sigmap=sigmap, a=a, a1=a1, b0=b0, b1=b1,
                fameffphi_raw=fameffphi_raw,
                fameffp_raw=fameffp_raw, yeareffphi_raw=yeareffphi_raw)
