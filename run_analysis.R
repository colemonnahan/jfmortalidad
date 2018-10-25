source("startup.R")

## Construye el modelo
compile('modelos/cjs_jf.cpp')
dyn.load(dynlib('modelos/cjs_jf'))
pars <- list(logNatM=c(-7,-7), logr=-10, a=.4, b=90,
             tau=array(0, dim=c(data$K, 2,3)),
             mu_tau=matrix(-2, nrow=2,ncol=3), logsigma_tau=1)
map <- list(logr=factor(NA), a=factor(NA), b=factor(NA),
            logNatM=factor(c(1,2)))
obj <- MakeADFun(data=data, parameters=pars, DLL='cjs_jf',
                 random=c('tau'), map=map)
obj$env$beSilent()
## opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(trace=1))
opt <- TMBhelper::Optimize(obj, control=list(trace=1))
rep <- sdreport(obj)
est <- data.frame(par=names(rep$value), est=rep$value, se=rep$sd)

tmp <- cbind(periodo=1:data$K,est[grep('k', est$par),])
tmp$lwr <- pmax(tmp$est-1.96*tmp$se,0)
tmp$upr <- tmp$est+1.96*tmp$se
tmp$genero <- substr(tmp$par,2,2)
tmp$evento <- as.numeric(substr(tmp$par,3,3))
tmp$periodo <- tmp$periodo + ifelse(tmp$genero=='M',0,.25)
g <- ggplot(tmp, aes(periodo, y=est, color=genero)) +
  geom_pointrange(aes(ymin=lwr, ymax=upr), fatten=.5) +
  facet_wrap('evento', ncol=1) + ylab("k")
ggsave('plots/k_by_sex_event.png', g, width=7, height=5)


xx <- droplevels(est[grep('sel_pred', est$par),])
plot(data$lengths_pred, xx$est, ylim=c(0,1), xlim=c(50,130))
lines(data$lengths_pred, xx$est+ xx$se*1.96)
lines(data$lengths_pred, xx$est- xx$se*1.96)
rug(data$lengths)

## par(mfrow=c(3,1))
## x <- 1:data$K
## plot(x=x, pch=16, xlab='Periodo', ylab='k', type='n',
##      ylim=c(0, .4))
## tmp <- est[est$par=='kM1',]
## x <- 1:data$K+.25
## points(x=x, y=tmp$est, pch=16, col=2)
## segments(x0=x, y0=tmp$est-1.96*tmp$se, y1=tmp$est+1.96*tmp$se, col=2)
## tmp <- est[est$par=='catchabilityM_pred',]
## x <- data$esfuerzo_pred
## plot(x=x, y=tmp$est, pch=16, xlab='Esfuerzo (1000s trampas)', ylab='Catchability',
##      ylim=c(0, .4))
## segments(x0=x, y0=tmp$est-1.96*tmp$se,
##          y1=tmp$est+1.96*tmp$se)
## tmp <- est[est$par=='catchabilityH_pred',]
## x <- data$esfuerzo_pred+.03
## points(x=x, y=tmp$est, pch=16, col=2)
## segments(x0=x, y0=tmp$est-1.96*tmp$se,
##          y1=tmp$est+1.96*tmp$se, col=2)

