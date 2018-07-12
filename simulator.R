##### Un archivo temporario para hacer una simulacion

#### La entrada de los datos
I <- 2000                                 # numero de los individuos
K <- 50                                   # numero de los periodos

## los parametros de usar en la simulacion, por ahora son iguales que el
## adjuste de los datos reales.
p0 <- 2 # constante probabilidad de captura
phi0 <- 2 # constante probabilidad de sobrevivencia
## componente de observacion
alive <- p <- CH <- matrix(NA, I,K)
last <- rep(NA, len=I)
p <- matrix(NA, nrow=I, ncol=K) # capture probability
phi <- matrix(NA, I, K-1);   # survival probability
inv_logit <- function(x) 1/(1+exp(-x))
last <- rep(NA, len=I)

for(i in 1:I){
  ## calcula sobrevivencia por todos los periodos excepto el final porque
  ## no hay nada informacion para saberlo
  for(t in 1:(K-1)) {
    x <- phi0;
    phi[i,t] <- inv_logit(x);
  }
  ## calcula probabilidad de recaptura por todos los periodos except el
  ## primero porque conocido (captura)
  p[i,1] <- 1 # captura
  for(t in 2:K){
    x <- p0
    p[i,t] <- inv_logit(x);
  }
  ## Observation component
  CH[i,1] <- alive[i,1] <- 1 ## alive and observed at first marking
  for(t in 2:K){
    ## vida ahora solo si vida y sobrevivido en el periodo anterior
    alive[i,t] <- alive[i,t-1]*rbinom(n=1, size=1, prob=phi[i,t-1])
    ## recapturado solo si vida
    CH[i,t] <- alive[i,t]*rbinom(n=1, size=1, prob=p[i,t])
  }
  ## last period animal seen
  last[i] <- tail(which(CH[i,]==1), n=1)
}

simdata <- list(I=I, K=K, CH=CH, last=last)
simpars <- list(phi0=phi0, p0=p0)
