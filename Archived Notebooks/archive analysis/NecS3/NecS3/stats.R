# R code. I suggest you install R-Studio. Very nice
library(ggplot2)
library(RJSONIO)

## Control Data

a<-fromJSON('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS1/NecS1_T1_vals_day1.json')
b<-fromJSON('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS1/NecS1_T1_vals_day2.json')
x <- data.frame(T1val=a, cond='day 0')
y <- data.frame(T1val=b, cond='day 1')
combined <- rbind(x,y)
ggplot(combined, aes(x=T1val, fill=cond)) + geom_density(alpha=.7) + xlim(1000, 5000)

ks.test(a,b)
wilcox.test(a,b)
t.test(a,b)

## Treated Data

a<-fromJSON('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS3/NecS3_T1_vals_day1.json')
b<-fromJSON('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS3/NecS3_T1_vals_day2.json')
x <- data.frame(T1val=a, cond='Day 1')
y <- data.frame(T1val=b, cond='Day 2')
combined <- rbind(x,y)
ggplot(combined, aes(x=T1val, fill=cond)) + geom_density(alpha=.7)  + xlim(1000, 5000)

ks.test(a,b)
wilcox.test(a,b)
t.test(a,b)