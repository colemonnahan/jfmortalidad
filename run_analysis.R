source("startup.R")

## Construye el modelo
compile('modelos/cjs_jf.cpp')
dyn.load(dynlib('modelos/cjs_jf'))
pars <- list(logM=log(.03), logr=-10, a=3, b=-.2,
             tauM=rep(-2,data$K), tauH=rep(-2,data$K),
             mu_tauM=-2, mu_tauH=-2, logsigma_tau=1)
obj <- MakeADFun(data=data, parameters=pars, DLL='cjs_jf',
                 random=c('tauM', 'tauH'), map=list(logr=factor(NA)))
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
tmp <- est[est$par=='kM',]
x <- 1:data$K
plot(x=x, y=tmp$est, pch=16, xlab='Periodo', ylab='k',
     ylim=c(0, .4))
segments(x0=x, y0=tmp$est-1.96*tmp$se,
         y1=tmp$est+1.96*tmp$se)
tmp <- est[est$par=='kH',]
x <- 1:data$K+.25
points(x=x, y=tmp$est, pch=16, col=2)
segments(x0=x, y0=tmp$est-1.96*tmp$se, y1=tmp$est+1.96*tmp$se, col=2)
tmp <- est[est$par=='catchabilityM_pred',]
x <- data$esfuerzo_pred
plot(x=x, y=tmp$est, pch=16, xlab='Esfuerzo (1000s trampas)', ylab='Catchability',
     ylim=c(0, .4))
segments(x0=x, y0=tmp$est-1.96*tmp$se,
         y1=tmp$est+1.96*tmp$se)
tmp <- est[est$par=='catchabilityH_pred',]
x <- data$esfuerzo_pred+.03
points(x=x, y=tmp$est, pch=16, col=2)
segments(x0=x, y0=tmp$est-1.96*tmp$se,
         y1=tmp$est+1.96*tmp$se, col=2)

