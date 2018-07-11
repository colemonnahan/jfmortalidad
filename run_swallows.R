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


### Ahora simulamos datos similares que los reales y los adjustamos. Tienes
### que correr el codigo abajo primero.
source("simulator.R")
par(mfrow=c(1,3))
## chequea que son similares
plot(data$last, main='Real data')
plot(last, main='Simulated data')
plot(opt$par, as.vector(unlist(simpars[1:7])), main='Parameters',
     xlab='Real', ylab='Simulated')
abline(0,1)
### adjusta el modelo con los datos simulados
compile('modelos/cjs_jf.cpp')
dyn.load(dynlib('modelos/cjs_jf'))
obj <- MakeADFun(data=simdata, parameters=simpars, DLL='cjs_jf',
                 random=c('fameffphi_raw', 'fameffp_raw', 'yeareffphi_raw'))
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)
est <- rep$par.fixed
se <- sqrt(diag(rep$cov.fixed))
true <- unlist(simpars[1:7])
## verifica que sirve
out <- data.frame(par=names(se), true=true, est=est,
                  covered= (true-2*se < est & est < true+2*se))
out
mean(out$covered)
