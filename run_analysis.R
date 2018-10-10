source("startup.R")

## Construye el modelo
compile('modelos/cjs_jf.cpp')
dyn.load(dynlib('modelos/cjs_jf'))
pars <- list(logM=log(.03), logr=log(.01), logk=log(.01), a=20, b=-.2)
obj <- MakeADFun(data=data, parameters=pars, DLL='cjs_jf')
obj$env$beSilent()
opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(trace=1))
rep <- sdreport(obj)
est <- data.frame(par=names(rep$value), est=rep$value, se=rep$sd)
selex <- est[est$par=='sel_pred',]
plot(data$lengths_pred, y=selex$est, ylim=c(0,1), pch=16, xlab='Length', ylab='Selex')
segments(x0=data$lengths_pred, y0=selex$est-1.96*selex$se,
         y1=selex$est+1.96*selex$se)
rug(x=data$lengths)
