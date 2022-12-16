setwd("/Users/uym2/my_gits/problin/Experiments/ML_brlens")
require(ggplot2)

d_MP = read.table("results_MP_brlen.txt",header=T)
ggplot(d_MP,aes(y=brlen/k,x=nodeDepth,color=as.factor(k))) + 
  stat_summary() + geom_line(stat="summary") +
  theme_classic()


d_MLP = read.table("results_MLP_brlen.txt",header=T)

ggplot(d_MLP,aes(y=brlen,x=nodeDepth,color=as.factor(k))) + 
  stat_summary() + geom_line(stat="summary") + 
  geom_hline(yintercept = 0.2,linetype=2) + 
  scale_y_log10() + theme_classic()

ggplot(d_MLP[d_MLP$rep == "rep019",],aes(y=brlen,x=nodeDepth)) + 
  stat_summary() + geom_line(stat="summary") +
  geom_hline(yintercept = 0.2,linetype=2) +
  facet_wrap(~k,scale="free") + 
  #scale_y_log10() + 
  theme_classic()

ggplot(d_MLP,aes(y=abs(brlen-0.2),x=k)) + 
  stat_summary() + geom_line(stat="summary") + 
  theme_classic()
