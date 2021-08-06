p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

fun.1= function(x) exp(x^1)

fun.2=function(x) (1+x)^4
fun.3=function(x) (1+x)^1
fun.4=function(x) (1+x)^0.5



plot.range1 <- data.frame(x=c(0, 10), Functions = factor(1))
plot.range2 <- data.frame(x=c(0, 10), Functions = factor(2))
#plot.range3 <- data.frame(x=c(0, 10), Functions = factor(3))
#plot.range4 <- data.frame(x=c(0, 10), Functions = factor(4))

p <- ggplot(NULL, aes(x=x, colour=Functions)) +
  stat_function(data = plot.range1, fun = fun.1, size=2) +
  stat_function(data = plot.range2, fun = fun.2, size=2) +
  #stat_function(data = plot.range3, fun = fun.3) +
  #stat_function(data = plot.range4, fun = fun.4) +
  scale_colour_manual(values = c("royalblue4", "seagreen"))+
  scale_y_log10()+
  scale_x_continuous(breaks=c(0,1, 2,3, 4,5, 6,7, 8,9, 10), limits=c(0,10), expand = c(0,0))+
  theme(panel.grid.minor = element_blank())+
  xlab("Seasonal Score, z")+
  ylab("Fitness, w(z)")+
  theme(legend.position = "none")
ggexport(p, filename = "~/phd_data/fitness_function.pdf", width = 5, height = 4)
