##### Un archivo temporario para hacer una simulacion

#### La entrada de los datos
I <- 7000                                 # numero de los individuos
K <- 28                                   # numero de los periodos

## El primero event de las capturas es 1 (2nd octubre hasta 18 octubre
## 2008) con 3000 individuos
## El segundo es en periodo 10 con 2000 individuos
## El ultimo es en 14 con 2000 individuos


## los parametros de usar en la simulacion, por ahora son iguales que el
## adjuste de los datos reales.
M <- .03
r <- .1
r <- .01
k <- .002

## componente de observacion
counts <- alive <- p <- CH <- matrix(NA, I,K)
last <- rep(NA, len=I)
p <- matrix(NA, nrow=I, ncol=K) # capture probability
phi <- matrix(NA, I, K);   # survival probability
last <- rep(NA, len=I)
effort <- rep(1000, len=K)

for(i in 1:I){
  ## inicializacion
  ## estan vivo en periodo uno
  counts[i,1] <- CH[i,1] <- alive[i,1] <- 1
  phi[i,1] <- exp(-M-counts[i,1]*r)
  p[i,1] <- 1 # captura
  for(t in 2:K) {
    ## calcula sobrevivencia
    phi[i,t] <- exp(-M-counts[i,t-1]*r)
    ## calcula tprobabilidad de recaptura por todos los periodos except el
    ## primero porque conocido (captura)
    p[i,t] <- 1*(1-exp(-k*effort[t]))
    ## vida ahora solo si vida y sobrevivido en el periodo anterior
    alive[i,t] <- alive[i,t-1]*rbinom(n=1, size=1, prob=phi[i,t-1])
    ## recapturado solo si vida
    CH[i,t] <- alive[i,t]*rbinom(n=1, size=1, prob=p[i,t])
    counts[i,t] <- counts[i,t-1]+CH[i,t]
  }
  ## el ultimo periodo en que un individuo fue visto
  last[i] <- tail(which(CH[i,]==1), n=1)
}

simdata <- list(I=I, K=K, CH=CH, last=last, counts=counts)
simpars <- list(phi0=phi0, p0=p0)
