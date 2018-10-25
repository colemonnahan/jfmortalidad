source("startup.R")

## Construye el modelo
dyn.unload(dynlib('modelos/cjs_jf_sim'))
compile('modelos/cjs_jf_sim.cpp')
dyn.load(dynlib('modelos/cjs_jf_sim'))

## ## Probando que sirve
## set.seed(32)
## out <- simulator(TRUE)
## obj <- MakeADFun(data=out$simdata, parameters=out$simpars, DLL='cjs_jf_sim')
## obj$fn()
## obj$gr()
## rep <- obj$report()
## i <- 1
## cbind(out$p[i,], rep$p[i,])
## cbind(out$phi[i,], rep$phi[i,])


### Ahora simulamos datos similares que los reales y los adjustamos. Tienes
### que correr el codigo abajo primero.
nrep <- 50  ## numero de iteraciones de monte carlo
coverage.list <- results.list <- list()
## no se puede usar la variable 'i' porque es usada en el archivo
## simulator.R
for(ii in 1:nrep){
  print(ii); set.seed(ii)
  out <- simulator(ii==1)
  ## plot(simdata$last, main='Simulated data')
  ## adjusta el modelo con los datos simulados
  obj <- MakeADFun(data=out$simdata, parameters=out$simpars,
                   DLL='cjs_jf_sim')
  obj$env$beSilent()
  opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(trace=0))
  rep <- sdreport(obj)
  est <- rep$par.fixed
  se <- sqrt(diag(rep$cov.fixed))
  true <- unlist(out$simpars)
  ## verifica que sirve
  results.list[[ii]] <- data.frame(replicate=ii,par=names(se), true=true, est=est, se=se,
                    covered= (true-2*se < est & est < true+2*se))
  coverage.list[[ii]] <- c((true-2*se < est & est < true+2*se))
}

## los parametros estimados deberian contener el valor de la simulacion
## aproximadamente 95% de los replicacions
apply(do.call(rbind, coverage.list), 2, mean, na.rm=TRUE)
results <- do.call(rbind, results.list)

## lo visualizar
g <- ggplot(data=results) +
  geom_linerange(aes(x=replicate, ymin=est-1.96*se, ymax=est+1.96*se)) +
  geom_point(aes(x=replicate, y=est))+
  geom_hline(aes(yintercept=true), col='red')+
  facet_wrap('par', scales='free', ncol=2)
ggsave('plots/simulacion_cobertura.png', g, width=7, height=5)
