message("Cargando los packages...")
library(TMB)
library(ggplot2)
library(reshape2)
library(plyr)


message("Cargando los datos...")
CH <- read.csv("datos/CH_Matrix.csv")
names(CH)[1:2] <- c("numero", "evento")
CH$evento <- as.factor(CH$evento)
individuo <- read.csv("datos/Green_tags.csv")

message("Cargando las funciones...")
##### Una funcion temporaria para hacer una simulacion
simulator <- function(make.plots){
  ## La entrada de los datos
  I <- 7000                                 # numero de los individuos
  K <- 28                                   # numero de los periodos
  ## El primero event de las capturas es 1 (2nd octubre hasta 18 octubre
  ## 2008) con 3000 individuos
  ## El segundo es en periodo 10 con 2000 individuos
  ## El ultimo es en 14 con 2000 individuos
  ##
  ## los parametros de usar en la simulacion, por ahora son iguales que el
  ## adjuste de los datos reales.
  M <- .03
  r <- .01
  k <- .0001
  ## Selectividad
  a <- 20.65
  b <- -.24
  ## componente de observacion
  counts <- matrix(0, I,K)
  alive <- p <- CH <- matrix(NA, I,K)
  first <- last <- rep(NA, len=I)
  p <- matrix(NA, nrow=I, ncol=K) # capture probability
  phi <- matrix(NA, I, K);   # survival probability
  last <- rep(NA, len=I)
  effort <- runif(K, min=500, max=2000)
  effort[14:19] <- 0
  ## las longitudes de los individuos de los datos reales
  lengths <- sample(na.omit(individuo$Size), size=I, replace=TRUE)-20
  for(i in 1:I){
    if(i <= 3000) t0 <- 1
    if(i > 3000 & i <= 5000) t0 <- 10
    if(i > 5000) t0 <- 13
    ## inicializacion
    ## estan vivo en periodo uno
    first[i] <- t0
    counts[i,t0] <- CH[i,t0] <- alive[i,t0] <- 1
    phi[i,t0] <- exp(-M-counts[i,t0]*r)
    p[i,t0] <- 1 # captura
    for(t in (t0+1):K) {
      ## calcula sobrevivencia
      phi[i,t] <- exp(-M-counts[i,t-1]*r)
      ## calcula probabilidad de recaptura para todos los periodos except el
      ## primero porque conocido (captura)
      p[i,t] <- 1*(1-exp(-k*effort[t]))/(1+exp(a+b*lengths[i]))
      ## actualiza esta probabilidad por la longitud
      ## vida ahora solo si vida y sobrevivido en el periodo anterior
      alive[i,t] <- alive[i,t-1]*rbinom(n=1, size=1, prob=phi[i,t-1])
      ## recapturado solo si vida
      CH[i,t] <- alive[i,t]*rbinom(n=1, size=1, prob=p[i,t])
      counts[i,t] <- counts[i,t-1]+ ifelse(!is.na(CH[i,t]), CH[i,t],0)
    }
    ## el ultimo periodo en que un individuo fue visto
    last[i] <- tail(which(CH[i,]==1), n=1)
  }
  CH[, 14:19] <- NA
  simdata <- list(I=I, K=K, CH=CH, last=last, counts=counts, effort=effort,
                  lengths=lengths, first=first)
  simpars <- list(logM=log(M), logr=log(r), logk=log(k))
  if(make.plots){
    par(mfcol=c(2,2))
    x <- seq(50,120, len=1000)
    plot(x, (1-exp(-k*2000))/(1+exp(a+b*x)), ylab='Prob. captura',
         type='l', xlim=c(50,120), ylim=c(0,1), xlab='Longitud')
    legend('topleft', legend=c('Max. Esfuerzo', 'Min. Esfuerzo'), lty=1, col=c(1,2))
    lines(x, (1-exp(-k*500))/(1+exp(a+b*x)), col='red')
    hist(lengths, xlim=c(50,120), xlab='Longitud')
    plot(0,0, xlim=c(1,K), ylim=c(0,1), xlab='Periodo',
         ylab='Pr. captura de 5 individuos',
         type='n')
    trash <- sapply(1:5, function(i) lines(2:K, p[i,-1], col=1))
    plot(0,0, xlim=c(1,K), ylim=c(.8,1), xlab='Periodo',
         ylab='Pr. sobrevivencia de 5 individuos',
         type='n')
    trash <- sapply(1:5, function(i) lines(2:K, phi[i,-1], col=1))
  }
  return(list(simdata=simdata, simpars=simpars, p=p, phi=phi))
}
## test <- simulator(1)
## CH.long <- melt(test$simdata$CH)
## names(CH.long) <- c("numero", "periodo", "recapturas")
## CH.long2 <- droplevels(subset(CH.long, !is.na(recapturas) & recapturas>0))
## g <- ggplot(CH.long2, aes(numero, y=periodo)) +
##   geom_jitter(width=0, alpha=.5, height=.1)
## g
