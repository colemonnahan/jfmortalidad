## Este codigo correr la version original del modelo "swallows" que es un
## CJS modelo. Examina el archivo swallows.stan para ver documentacion del
## modelo. Vamos a usarlo como el punto de partida de construir nuestro
## modelo.
library(TMB)

data <- readRDS('datos/swallows_datos.RDS')
inits <- list(
  sigmayearphi=.7, sigmaphi=.5, sigmap=.9, a=rep(3.5, len=data$K-1), a1=0,
  b0=rep(2, len=4), b1=rep(0, len=4), fameffphi_raw=rep(0, len=data$nfam),
  fameffp_raw=rep(0, len=data$nfam), yeareffphi_raw=rep(0, len=4))

## Compile, link and build model
compile('modelos/swallows.cpp')
dyn.load(dynlib('modelos/swallows'))
obj <- MakeADFun(data=data, parameters=inits,
                 random=c('fameffphi_raw', 'fameffp_raw', 'yeareffphi_raw'))
## Adjustar con efectos aleatorios
opt <- nlminb(obj$par, obj$fn, obj$gr)

## extraer los parametros optimal para usarlos en la simulacion
dput(as.vector(round(opt$par,2)[4:20]))
rep <- obj$report()
dput(round(rep$fameffp_raw,2))
dput(round(rep$fameffphi_raw,2))
dput(round(rep$yeareffphi_raw,2))
