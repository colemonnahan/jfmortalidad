source('startup.R')
ww <- 6
hh <- 4


CH.long <- melt(df, id.vars=names(df)[1:10],
                variable.name='periodo', value.name='recapturas')
CH.long$periodo <- as.numeric(gsub("X", "", x=CH.long$periodo))
CH.long$numero <- CH.long$Tag_num
CH.long$genero <- CH.long$Sex
CH.long2 <- droplevels(subset(CH.long, !is.na(recapturas) & recapturas>0))
CH.long2 <- CH.long2[,c('numero', 'periodo', 'recapturas', 'genero', 'Event')]
CH.long2$data <- "real"
g <- ggplot(CH.long2, aes(numero, y=periodo, color=factor(Event))) +
  geom_jitter(width=0, alpha=.5, height=.1) + facet_grid(genero~data)
ggsave('plots/CH_matriz.png', g, width=7, height=5)


## Simulate data to compare to real
test <- simulator(1)
CH2.long <- melt(test$simdata$CH)
names(CH2.long) <- c("numero", "periodo", "recapturas")
CH2.long$genero <- sample(x=c("M", "F"), size=nrow(CH2.long), repl=TRUE)
CH2.long$Event <- 1
CH2.long2 <- droplevels(subset(CH2.long, !is.na(recapturas) & recapturas>0))
CH2.long2$data <- "simulado"
CH.all <- rbind(CH.long2, CH2.long2)
g <- ggplot(CH.all, aes(numero, y=periodo, color=factor(Event))) +
  geom_jitter(width=0, size=.1, alpha=.5, height=.1) + facet_grid(~data)
ggsave('plots/CH_matriz_all.png', g, width=7, height=5)

## Use both simulated and real or just real?
CH.long2 <- CH.all

g <- ggplot(CH.long2, aes(recapturas)) + geom_bar() + scale_y_log10() +
  xlab("Numero de recapturas en un periodo") + facet_grid(genero~data)
ggsave('plots/recapturas_periodo.png', g, width=ww, height=hh)

xx <- ddply(CH.long2, .(numero, genero, data), summarize,
            recapturas.total=sum(recapturas),
            recapturas.periodo=length(recapturas),
            periodo.rango=max(periodo)-min(periodo),
            periodo.primero=min(periodo))
g <- ggplot(xx, aes(recapturas.total))  + geom_bar() +# scale_y_log10() +
  xlab("Numero de recapturas") +facet_grid(genero~data)
ggsave('plots/recapturas_total.png', g, width=ww, height=hh)
g <- ggplot(xx, aes(recapturas.periodo))  + geom_bar() + scale_y_log10() +
  xlab("Periodos con un recaptura")+facet_grid(genero~data)
ggsave('plots/recapturas_periodo.png', g, width=ww, height=hh)
g <- ggplot(xx, aes(periodo.rango))  + geom_bar() + #scale_y_log10() +
  xlab("Rango de los periodos")+facet_grid(genero~data)
ggsave('plots/rango_periodo.png', g, width=ww, height=hh)
g <- ggplot(xx, aes(periodo.primero, periodo.rango)) +
  geom_jitter(width=.5, height=0, alpha=.5) +
  xlab('Periodo de la primera captura') +
  ylab("Rango de los periodos de captura") +facet_grid(genero~data)
ggsave('plots/rango_primero.png', g, width=ww, height=hh)


