compile('modelos/cjs_jf.cpp')
dyn.load(dynlib('modelos/cjs_jf'))

### Ahora simulamos datos similares que los reales y los adjustamos. Tienes
### que correr el codigo abajo primero.
nrep <- 10  ## numero de iteraciones de monte carlo
coverage.list <- results.list <- list()
## no se puede usar la variable 'i' porque es usada en el archivo
## simulator.R
for(ii in 1:nrep){
  print(ii); set.seed(ii)
  make.plots <- ii==1
  source("simulator.R")
  ## plot(simdata$last, main='Simulated data')
  ## adjusta el modelo con los datos simulados
  obj <- MakeADFun(data=simdata, parameters=simpars, DLL='cjs_jf')
  obj$env$beSilent()
  opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(trace=0))
  rep <- sdreport(obj)
  est <- rep$par.fixed
  se <- sqrt(diag(rep$cov.fixed))
  true <- unlist(simpars)
  ## verifica que sirve
  results.list[[ii]] <- data.frame(rep=ii,par=names(se), true=true, est=est, se=se,
                    covered= (true-2*se < est & est < true+2*se))
  coverage.list[[ii]] <- c((true-2*se < est & est < true+2*se))
}

## los parametros estimados deberian contener el valor de la simulacion
## aproximadamente 95% de los replicacions
apply(do.call(rbind, coverage.list), 2, mean, na.rm=TRUE)

## lo visualizar
results <- do.call(rbind, results.list)
ggplot(data=results) +
  geom_linerange(aes(x=rep, ymin=est-1.96*se, ymax=est+1.96*se)) +
  geom_point(aes(x=rep, y=est))+
  geom_hline(aes(yintercept=true), col='red')+
  facet_wrap('par', scales='free')

