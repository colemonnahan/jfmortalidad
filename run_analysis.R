source("startup.R")

## Construye el modelo
compile('modelos/cjs_jf.cpp')
dyn.load(dynlib('modelos/cjs_jf'))
pars <- list(logM=log(.03), logr=log(.01), logk=log(.01))
obj <- MakeADFun(data=data, parameters=pars, DLL='cjs_jf')
obj$env$beSilent()
opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(trace=1))
rep <- sdreport(obj)
est <- rep$par.fixed
se <- sqrt(diag(rep$cov.fixed))

cbind(est, se)

