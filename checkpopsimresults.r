library(ggplot2)

dat <- at <- read.table(file="/Users/Rob/Dropbox/Current_work/popsim/10_100/results.txt", sep = "\t", header= TRUE)

#let's just look at a single mutation rate
r <- dat[which(dat$u==0.00010),]

#now let's plot it out
p <- ggplot(r, aes(x=S, y=k, group=pN))
p+geom_line(aes(color=pN), size=2, alpha=0.5)+geom_point() + xlab("Overlapping generations (proportion of offspring that survive)") + ylab("substitution rate") + scale_colour_gradient(name="Degree of\npopulation size\nfluctuation")
