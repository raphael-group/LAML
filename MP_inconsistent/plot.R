# setwd("/Users/uym2/my_gits/problin/MP_inconsistent")
setwd("/n/fs/ragr-research/projects/problin/MP_inconsistent")

require(ggplot2)

d = read.table("results.txt",header=T)

ggplot(d,aes(x=k,y=percent_correct*100,color=method)) + geom_line() + geom_point() +
  scale_x_log10() + theme_classic() + xlab("Number of sites") + ylab("Percent correct")

ggsave("MP_inconsistent.pdf")
