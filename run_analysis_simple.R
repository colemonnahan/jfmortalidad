source("startup.R")
## data$CH[is.na(data$CH)] <- 0

## Construye el modelo
dyn.unload(dynlib('modelos/cjs_jf_simple'))
compile('modelos/cjs_jf_simple.cpp')
dyn.load(dynlib('modelos/cjs_jf_simple'))
pars <- list(phi2=rep(-3, data$K),
             p2=rep(-5, data$K))
pars$phi2[1] <- 5
## map <- list(logr=factor(NA), a=factor(NA), b=factor(NA),
##             logNatM=factor(c(1,2)))
mapp <- 1:data$K
mapphi <- 1:data$K
mapp[16:24] <- NA
mapphi[c(1,29)] <- NA
map <- list(p2=factor(mapp), phi2=factor(mapphi))
## map <- NULL
obj <- MakeADFun(data=data, parameters=pars, DLL='cjs_jf_simple',
                 map=map)
obj$env$beSilent()
## opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(trace=1))
opt <- TMBhelper::Optimize(obj, control=list(trace=10))
rep <- sdreport(obj)
est <- data.frame(par=names(rep$value), est=rep$value, se=rep$sd)
tmp <- cbind(periodo=1:data$K,est)
tmp$lwr <- pmax(tmp$est-1.96*tmp$se,0)
tmp$upr <- tmp$est+1.96*tmp$se
tmp$par2 <- ifelse(tmp$par=='p3', 'p', 'phi')
g <- ggplot(tmp, aes(periodo, y=est, color=par2)) +
  geom_pointrange(aes(ymin=lwr, ymax=upr), fatten=.5)
g



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

