source('startup.R')
ww <- 6
hh <- 4


CH.long <- melt(df, id.vars=names(df)[1:10],
                variable.name='periodo', value.name='recapturas')
CH.long$periodo <- as.numeric(gsub("X", "", x=CH.long$periodo))
CH.long$numero <- CH.long$Tag_num
CH.long$genero <- CH.long$Sex
CH.long2 <- droplevels(subset(CH.long, !is.na(recapturas) & recapturas>0))


g <- ggplot(CH.long2, aes(numero, y=periodo, color=evento)) +
  geom_jitter(width=0, alpha=.5, height=.1) + facet_wrap('genero')
ggsave('plots/CH_matriz.png', g, width=7, height=5)
g <- ggplot(CH.long2, aes(recapturas)) + geom_bar() + scale_y_log10() +
  xlab("Numero de recapturas en un periodo") + facet_wrap('genero')
ggsave('plots/recapturas_periodo.png', g, width=ww, height=hh)

xx <- ddply(CH.long2, .(numero, genero), summarize,
            recapturas.total=sum(recapturas),
            recapturas.periodo=length(recapturas),
            periodo.rango=max(periodo)-min(periodo),
            periodo.primero=min(periodo))
g <- ggplot(xx, aes(recapturas.total))  + geom_bar() + scale_y_log10() +
  xlab("Numero de recapturas") +facet_wrap('genero')
ggsave('plots/recapturas_total.png', g, width=ww, height=hh)
g <- ggplot(xx, aes(recapturas.periodo))  + geom_bar() + scale_y_log10() +
  xlab("Periodos con un recaptura")+facet_wrap('genero')
ggsave('plots/recapturas_periodo.png', g, width=ww, height=hh)
g <- ggplot(xx, aes(periodo.rango))  + geom_bar() + scale_y_log10() +
  xlab("Rango de los periodos")+facet_wrap('genero')
ggsave('plots/rango_periodo.png', g, width=ww, height=hh)
g <- ggplot(xx, aes(periodo.primero, periodo.rango)) +
  geom_jitter(width=.5, height=0, alpha=.5) +
  xlab('Periodo de la primera captura') +
  ylab("Rango de los periodos de captura") +facet_wrap('genero')
ggsave('plots/rango_primero.png', g, width=ww, height=hh)
