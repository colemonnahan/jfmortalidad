source("startup.R")

## Construye el modelo
compile('modelos/cjs_jf.cpp')
dyn.load(dynlib('modelos/cjs_jf'))
pars <- list(logM=log(.03), logr=-10, a=3, b=-.2,
             tau=rep(-2,data$K), mu_tau=log(.1), logsigma_tau=1)
obj <- MakeADFun(data=data, parameters=pars, DLL='cjs_jf',
                 random='tau', map=list(logr=factor(NA)))
obj$env$beSilent()
## opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(trace=1))
opt <- TMBhelper::Optimize(obj, control=list(trace=1))
rep <- sdreport(obj)
est <- data.frame(par=names(rep$value), est=rep$value, se=rep$sd)

par(mfrow=c(3,1))
selex <- est[est$par=='sel_pred',]
plot(data$lengths_pred, y=selex$est, ylim=c(0,1), pch=16, xlab='Length', ylab='Selex')
segments(x0=data$lengths_pred, y0=selex$est-1.96*selex$se,
         y1=selex$est+1.96*selex$se)
rug(x=data$lengths)
catchability <- est[est$par=='k',]
x <- 1:data$K
plot(x=x, y=catchability$est, pch=16, xlab='Periodo', ylab='k',
     ylim=c(0, max(catchability$est+1.96*catchability$se)))
segments(x0=x, y0=catchability$est-1.96*catchability$se,
         y1=catchability$est+1.96*catchability$se)
effort <- est[est$par=='catchability_pred',]
x <- data$esfuerzo_pred
plot(x=x, y=effort$est, pch=16, xlab='Esfuerzo (1000s trampas)', ylab='Catchability',
     ylim=c(0, max(effort$est+1.96*effort$se)))
segments(x0=x, y0=effort$est-1.96*effort$se,
         y1=effort$est+1.96*effort$se)

