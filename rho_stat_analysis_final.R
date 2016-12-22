library(gtools)
library(caTools)
library(plyr)
library(car)
library(lme4)
library(zoo)

setwd("/Users/jackdaw/Dropbox/Heterochromatin/prelim_results/final_table")
dat<-read.table("submission_table", header=T)


# replace NA values of rho with ones of adjacent windows

  intermediate_table<-ddply(dat, .(scaffold), 
                       cor2_rho_extended=na.locf(na.locf(cor2_rho, na.rm=F, maxgap = 5), na.rm=F, fromLast=T, maxgap=5),
                       swe_rho_extended=na.locf(na.locf(swe_rho, na.rm=F, maxgap = 5), na.rm=F, fromLast=T, maxgap=5),
                       cor2_theta_extended=na.locf(na.locf(theta_cor2, na.rm=F, maxgap = 5), na.rm=F, fromLast=T, maxgap=5),
                       swe_pol_theta_extended=na.locf(na.locf(theta_SWE_POL, na.rm=F, maxgap = 5), na.rm=F, fromLast=T, maxgap=5),
                       cor2_swe_fst_extended=na.locf(na.locf(cor2_swe_fst, na.rm=F, maxgap = 5), na.rm=F, fromLast=T, maxgap=5),
                       transform)
  dat<-intermediate_table

chr1<-subset(dat, chr == "1")
chr1A<-subset(dat, chr == "1A")
chr2<-subset(dat, chr == "2")
chr3<-subset(dat, chr == "3")
chr4<-subset(dat, chr == "4")
chr4A<-subset(dat, chr == "4A")
chr5<-subset(dat, chr == "5")
chr6<-subset(dat, chr == "6")
chr7<-subset(dat, chr == "7")
chr8<-subset(dat, chr == "8")
chr9<-subset(dat, chr == "9")
chr10<-subset(dat, chr == "10")
chr11<-subset(dat, chr == "11")
chr12<-subset(dat, chr == "12")
chr13<-subset(dat, chr == "13")
chr14<-subset(dat, chr == "14")
chr15<-subset(dat, chr == "15")
chr17<-subset(dat, chr == "17")
chr18<-subset(dat, chr == "18")
chr19<-subset(dat, chr == "19")
chr20<-subset(dat, chr == "20")
chr21<-subset(dat, chr == "21")
chr22<-subset(dat, chr == "22")
chr23<-subset(dat, chr == "23")
chr24<-subset(dat, chr == "24")
chr25<-subset(dat, chr == "25")
chr26<-subset(dat, chr == "26")
chr27<-subset(dat, chr == "27")
chr28<-subset(dat, chr == "28")
chrZ<-subset(dat, chr == "Z")

chromlist<-list(chr1,chr1A,chr2,chr3,chr4,chr4A,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chrZ)
titles<-c("chr1","chr1A","chr2","chr3","chr4","chr4A","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24","chr25","chr26","chr27","chr28","chrZ")

# Rho carrion crow
i<-24
for (i in 1:28){
  mi<-min(chromlist[[i]]$cor2_rho_extended, na.rm=T)
  ma<-max(chromlist[[i]]$cor2_rho_extended, na.rm=T)
  ra<-range(chromlist[[i]]$cor2_rho_extended, na.rm=T)
  plot(chromlist[[i]]$winnumber, chromlist[[i]]$cor2_rho_extended, pch=19, cex=0.1, ylim=c(mi-ra[2]/4, ma), main=titles[i])
  title(main=paste(titles[i]))
  points(chromlist[[i]]$winnumber, chromlist[[i]]$RAM-(1-(ra[2]/4) + (3.5*ra[2]/12)), pch=19, col="red")
  points(chromlist[[i]]$winnumber, chromlist[[i]]$crowSat-(1-(ra[2]/4) + (4*ra[2]/12)), pch=15, col="orange")
  abline(0,1, v=chromlist[[i]]$winnumber[which(chromlist[[i]]$scaf_end==1)], col="#0000ff11")
  i<-i+1
}

# Rho hooded crow
i<-24
for (i in 1:28){
  mi<-min(chromlist[[i]]$swe_rho_extended, na.rm=T)
  ma<-max(chromlist[[i]]$swe_rho_extended, na.rm=T)
  ra<-range(chromlist[[i]]$swe_rho_extended, na.rm=T)
  plot(chromlist[[i]]$winnumber, chromlist[[i]]$swe_rho_extended, pch=19, cex=0.1, ylim=c(mi-ra[2]/4, ma), main=titles[i])
  title(main=paste(titles[i]))
  points(chromlist[[i]]$winnumber, chromlist[[i]]$RAM-(1-(ra[2]/4) + (3.5*ra[2]/12)), pch=19, col="red")
  points(chromlist[[i]]$winnumber, chromlist[[i]]$crowSat-(1-(ra[2]/4) + (4*ra[2]/12)), pch=15, col="orange")
  abline(0,1, v=chromlist[[i]]$winnumber[which(chromlist[[i]]$scaf_end==1)], col="#0000ff11")
  i<-i+1
}

# theta carrion crow
i<-19
for (i in 1:28){
  mi<-min(chromlist[[i]]$theta_cor2, na.rm=T)
  ma<-max(chromlist[[i]]$theta_cor2, na.rm=T)
  ra<-range(chromlist[[i]]$theta_cor2, na.rm=T)
  plot(chromlist[[i]]$winnumber, chromlist[[i]]$theta_cor2, pch=19, cex=0.1, ylim=c(mi-ra[2]/4, ma), main=titles[i])
  title(main=paste(titles[i]))
  points(chromlist[[i]]$winnumber, chromlist[[i]]$RAM-(1-(ra[2]/4) + (3.5*ra[2]/12)), pch=19, col="red")
  points(chromlist[[i]]$winnumber, chromlist[[i]]$crowSat-(1-(ra[2]/4) + (4*ra[2]/12)), pch=15, col="orange")
  abline(0,1, v=chromlist[[i]]$winnumber[which(chromlist[[i]]$scaf_end==1)], col="#0000ff11")
  i<-i+1
}

# theta hooded crow
i<-19
for (i in 1:28){
  mi<-min(chromlist[[i]]$theta_SWE_POL, na.rm=T)
  ma<-max(chromlist[[i]]$theta_SWE_POL, na.rm=T)
  ra<-range(chromlist[[i]]$theta_SWE_POL, na.rm=T)
  plot(chromlist[[i]]$winnumber, chromlist[[i]]$theta_SWE_POL, pch=19, cex=0.1, ylim=c(mi-ra[2]/4, ma), main=titles[i])
  title(main=paste(titles[i]))
  points(chromlist[[i]]$winnumber, chromlist[[i]]$RAM-(1-(ra[2]/4) + (3.5*ra[2]/12)), pch=19, col="red")
  points(chromlist[[i]]$winnumber, chromlist[[i]]$crowSat-(1-(ra[2]/4) + (4*ra[2]/12)), pch=15, col="orange")
  abline(0,1, v=chromlist[[i]]$winnumber[which(chromlist[[i]]$scaf_end==1)], col="#0000ff11")
  i<-i+1
}


# Fst carrion vs hooded crow

i<-19
for (i in 1:30){
  mi<-min(chromlist[[i]]$cor2_swe_fst_extended, na.rm=T)
  ma<-max(chromlist[[i]]$cor2_swe_fst_extended, na.rm=T)
  ra<-range(chromlist[[i]]$cor2_swe_fst_extended, na.rm=T)
  plot(chromlist[[i]]$winnumber, chromlist[[i]]$cor2_swe_fst_extended, pch=19, cex=0.1, ylim=c(mi-ra[2]/4, ma), main=titles[i])
  title(main=paste(titles[i]))
  points(chromlist[[i]]$winnumber, chromlist[[i]]$RAM-(1-(ra[2]/4) + (3.5*ra[2]/12)), pch=19, col="red")
  points(chromlist[[i]]$winnumber, chromlist[[i]]$crowSat-(1-(ra[2]/4) + (5*ra[2]/12)), pch=15, col="orange")
  abline(0,1, v=chromlist[[i]]$winnumber[which(chromlist[[i]]$scaf_end==1)], col="#0000ff11")
  i<-i+1
}

# set occurences of hetchrom, scaffold starts/ends and crowSat as factors

dat$RAM[is.na(dat$RAM)]<-0
dat$scaf_end[is.na(dat$scaf_end)]<-0
dat$crowSat[is.na(dat$crowSat)]<-0

dat$RAM<-factor(dat$RAM)
dat$scaf_end<-factor(dat$scaf_end)
dat$crowSat<-factor(dat$crowSat)

se<-subset(dat, scaf_end=="1")
str(se)


# rho boxplot

  boxplot(
          log(dat$swe_rho_extended[which(dat$scaf_end != 1)]), # whole data without scaffold starts / ends
          log(se$swe_rho_extended[which(se$RAM != 1 | se$crowSat !=1)]), # scaffold starts / ends without RAM pattern
          log(se$swe_rho_extended[which(se$RAM==1 )]),  # 
          log(se$swe_rho_extended[which(se$crowSat==1)]),
          log(dat$cor2_rho_extended[which(dat$scaf_end != 1)]),
          log(se$cor2_rho_extended[which(se$RAM != 1 | se$crowSat !=1)]), 
          log(se$cor2_rho_extended[which(se$RAM==1)]), 
          log(se$cor2_rho_extended[which(se$crowSat==1)]),
          ylab="log(rho)"
          )

#mappability

boxplot(
  dat$map_100[which(dat$scaf_end != 1)], # whole data without scaffold starts / ends
  se$map_100[which(se$RAM != 1 | se$crowSat !=1)], # scaffold starts / ends without RAM pattern
  se$map_100[which(se$RAM==1 )],  # 
  se$map_100[which(se$crowSat==1)],
  ylab="log(GQ)"
)

boxplot(
  dat$map_200[which(dat$scaf_end != 1)], # whole data without scaffold starts / ends
  se$map_200[which(se$RAM != 1 | se$crowSat !=1)], # scaffold starts / ends without RAM pattern
  se$map_200[which(se$RAM==1 )],  # 
  se$map_200[which(se$crowSat==1)],
  ylab="map_200"
)

t.test(dat$map_200[which(dat$scaf_end != 1)], se$map_200[which(se$RAM != 1 | se$crowSat !=1)])
t.test(dat$map_200[which(dat$scaf_end != 1)], se$map_200[which(se$RAM==1 )])  
t.test(dat$cor2_genotype_qual[which(dat$scaf_end != 1)], se$cor2_genotype_qual[which(se$crowSat==1)])


svg(file="mappability_50.svg")
par(mfrow=c(4,2), mar=c(2,2,2,2))
plot(dat$map_50[which(dat$scaf_end != 1)], dat$swe_rho_extended[which(dat$scaf_end != 1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Hooded Crow: Mappability (50 bp) vs. Rho - Genome-wide", cex.main=0.8)
plot(dat$map_50[which(dat$scaf_end != 1)], dat$cor2_rho_extended[which(dat$scaf_end != 1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Carrion crow: Mappability (50 bp) vs. Rho - Genome-wide", cex.main=0.8)
plot(se$map_50[which(se$RAM != 1 | se$crowSat !=1)], se$swe_rho_extended[which(se$RAM != 1 | se$crowSat !=1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Hooded Crow:Mappability (50 bp) vs. Rho - Scaffold ends", cex.main=0.8)
plot(se$map_50[which(se$RAM != 1 | se$crowSat !=1)], se$cor2_rho_extended[which(se$RAM != 1 | se$crowSat !=1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Carrion crow :Mappability (50 bp) vs. Rho - Scaffold ends", cex.main=0.8)
plot(se$map_50[which(se$RAM==1 )], se$swe_rho_extended[which(se$RAM==1 )], pch=19, cex=0.5,ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Hooded Crow:Mappability (50 bp) vs. Rho - RAMs", cex.main=0.8)
plot(se$map_50[which(se$RAM==1 )], se$cor2_rho_extended[which(se$RAM==1 )], pch=19, cex=0.5,ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Carrion crow:Mappability (50 bp) vs. Rho - RAMs", cex.main=0.8)
plot(se$map_50[which(se$crowSat==1)], se$swe_rho_extended[which(se$crowSat==1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Hooded Crow:Mappability (50 bp) vs. Rho - crowSat1", cex.main=0.8)
plot(se$map_50[which(se$crowSat==1)], se$cor2_rho_extended[which(se$crowSat==1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Carrion crow:Mappability (50 bp) vs. Rho - crowSat1", cex.main=0.8)
dev.off()


#svg(file="mappability_100.svg")
par(mfrow=c(4,2), mar=c(2,2,2,2))
plot(dat$map_100[which(dat$scaf_end != 1)], dat$swe_rho_extended[which(dat$scaf_end != 1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Hooded Crow: Mappability (100 bp) vs. Rho - Genome-wide", cex.main=0.8)
plot(dat$map_100[which(dat$scaf_end != 1)], dat$cor2_rho_extended[which(dat$scaf_end != 1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Carrion crow: Mappability (100 bp) vs. Rho - Genome-wide", cex.main=0.8)
plot(se$map_100[which(se$RAM != 1 | se$crowSat !=1)], se$swe_rho_extended[which(se$RAM != 1 | se$crowSat !=1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Hooded Crow:Mappability (100 bp) vs. Rho - Scaffold ends", cex.main=0.8)
plot(se$map_100[which(se$RAM != 1 | se$crowSat !=1)], se$cor2_rho_extended[which(se$RAM != 1 | se$crowSat !=1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Carrion crow :Mappability (100 bp) vs. Rho - Scaffold ends", cex.main=0.8)
plot(se$map_100[which(se$RAM==1 )], se$swe_rho_extended[which(se$RAM==1 )], pch=19, cex=0.5,ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Hooded Crow:Mappability (100 bp) vs. Rho - RAMs", cex.main=0.8)
plot(se$map_100[which(se$RAM==1 )], se$cor2_rho_extended[which(se$RAM==1 )], pch=19, cex=0.5,ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Carrion crow:Mappability (100 bp) vs. Rho - RAMs", cex.main=0.8)
plot(se$map_100[which(se$crowSat==1)], se$swe_rho_extended[which(se$crowSat==1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Hooded Crow:Mappability (100 bp) vs. Rho - crowSat1", cex.main=0.8)
plot(se$map_100[which(se$crowSat==1)], se$cor2_rho_extended[which(se$crowSat==1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Carrion crow:Mappability (100 bp) vs. Rho - crowSat1", cex.main=0.8)
dev.off()

svg(file="mappability_150.svg")
par(mfrow=c(4,2), mar=c(2,2,1,1))
plot(dat$map_150[which(dat$scaf_end != 1)], dat$swe_rho_extended[which(dat$scaf_end != 1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Hooded Crow: Mappability (150 bp) vs. Rho - Genome-wide", cex.main=0.8)
plot(dat$map_150[which(dat$scaf_end != 1)], dat$cor2_rho_extended[which(dat$scaf_end != 1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Carrion crow: Mappability (150 bp) vs. Rho - Genome-wide", cex.main=0.8)
plot(se$map_150[which(se$RAM != 1 | se$crowSat !=1)], se$swe_rho_extended[which(se$RAM != 1 | se$crowSat !=1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Hooded Crow:Mappability (150 bp) vs. Rho - Scaffold ends", cex.main=0.8)
plot(se$map_150[which(se$RAM != 1 | se$crowSat !=1)], se$cor2_rho_extended[which(se$RAM != 1 | se$crowSat !=1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Carrion crow :Mappability (150 bp) vs. Rho - Scaffold ends", cex.main=0.8)
plot(se$map_150[which(se$RAM==1 )], se$swe_rho_extended[which(se$RAM==1 )], pch=19, cex=0.5,ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Hooded Crow:Mappability (150 bp) vs. Rho - RAMs", cex.main=0.8)
plot(se$map_150[which(se$RAM==1 )], se$cor2_rho_extended[which(se$RAM==1 )], pch=19, cex=0.5,ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Carrion crow:Mappability (150 bp) vs. Rho - RAMs", cex.main=0.8)
plot(se$map_150[which(se$crowSat==1)], se$swe_rho_extended[which(se$crowSat==1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Hooded Crow:Mappability (150 bp) vs. Rho - crowSat1", cex.main=0.8)
plot(se$map_150[which(se$crowSat==1)], se$cor2_rho_extended[which(se$crowSat==1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Carrion crow:Mappability (150 bp) vs. Rho - crowSat1", cex.main=0.8)
dev.off()

svg(file="mappability_200.svg")
par(mfrow=c(4,2), mar=c(2,2,1,1))
plot(dat$map_200[which(dat$scaf_end != 1)], dat$swe_rho_extended[which(dat$scaf_end != 1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Hooded Crow: Mappability (200 bp) vs. Rho - Genome-wide", cex.main=0.8)
plot(dat$map_200[which(dat$scaf_end != 1)], dat$cor2_rho_extended[which(dat$scaf_end != 1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Carrion crow: Mappability (200 bp) vs. Rho - Genome-wide", cex.main=0.8)
plot(se$map_200[which(se$RAM != 1 | se$crowSat !=1)], se$swe_rho_extended[which(se$RAM != 1 | se$crowSat !=1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Hooded Crow: Mappability (200 bp) vs. Rho - Scaffold ends", cex.main=0.8)
plot(se$map_200[which(se$RAM != 1 | se$crowSat !=1)], se$cor2_rho_extended[which(se$RAM != 1 | se$crowSat !=1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Carrion crow: Mappability (200 bp) vs. Rho - Scaffold ends", cex.main=0.8)
plot(se$map_200[which(se$RAM==1 )], se$swe_rho_extended[which(se$RAM==1 )], pch=19, cex=0.5,ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Hooded Crow: Mappability (200 bp) vs. Rho - RAMs", cex.main=0.8)
plot(se$map_200[which(se$RAM==1 )], se$cor2_rho_extended[which(se$RAM==1 )], pch=19, cex=0.5,ylim=c(0,1.25),xlim=c(0,0.8), xlab="",ylab="")
title(main="Carrion crow: Mappability (200 bp) vs. Rho - RAMs", cex.main=0.8)
plot(se$map_200[which(se$crowSat==1)], se$swe_rho_extended[which(se$crowSat==1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Hooded Crow: Mappability (200 bp) vs. Rho - crowSat1", cex.main=0.8)
plot(se$map_200[which(se$crowSat==1)], se$cor2_rho_extended[which(se$crowSat==1)], pch=19, cex=0.5, ylim=c(0,1.25),xlim=c(0,0.8), xlab="", ylab="")
title(main="Carrion crow: Mappability (200 bp) vs. Rho - crowSat1", cex.main=0.8)
dev.off()

summary(dat$map_50[which(dat$scaf_end != 1)])
summary(se$map_50[which(se$RAM != 1 | se$crowSat !=1)])
summary(se$map_50[which(se$RAM==1 )])
summary(se$map_50[which(se$crowSat==1)])

summary(dat$map_100[which(dat$scaf_end != 1)])
summary(se$map_100[which(se$RAM != 1 | se$crowSat !=1)])
summary(se$map_100[which(se$RAM==1 )])
summary(se$map_100[which(se$crowSat==1)])

summary(dat$map_150[which(dat$scaf_end != 1)])
summary(se$map_150[which(se$RAM != 1 | se$crowSat !=1)])
summary(se$map_150[which(se$RAM==1 )])
summary(se$map_150[which(se$crowSat==1)])

summary(dat$map_200[which(dat$scaf_end != 1)])
summary(se$map_200[which(se$RAM != 1 | se$crowSat !=1)])
summary(se$map_200[which(se$RAM==1 )])
summary(se$map_200[which(se$crowSat==1)])

# genotype quality 

boxplot(
  log(dat$cor2_genotype_qual[which(dat$scaf_end != 1)]), # whole data without scaffold starts / ends
  log(se$cor2_genotype_qual[which(se$RAM != 1 | se$crowSat !=1)]), # scaffold starts / ends without RAM pattern
  log(se$cor2_genotype_qual[which(se$RAM==1 )]),  # 
  log(se$cor2_genotype_qual[which(se$crowSat==1)]),
  log(dat$cor2_genotype_qual[which(dat$scaf_end != 1)]),
  log(se$cor2_genotype_qual[which(se$RAM != 1 | se$crowSat !=1)]), 
  log(se$cor2_genotype_qual[which(se$RAM==1)]), 
  log(se$cor2_genotype_qual[which(se$crowSat==1)]),
  ylab="log(GQ)"
)

boxplot(
  dat$swe_genotype_qual[which(dat$scaf_end != 1)], # whole data without scaffold starts / ends
  se$swe_genotype_qual[which(se$RAM != 1 | se$crowSat !=1)], # scaffold starts / ends without RAM pattern
  se$swe_genotype_qual[which(se$RAM==1 )],  # 
  se$swe_genotype_qual[which(se$crowSat==1)],
  dat$cor2_genotype_qual[which(dat$scaf_end != 1)], # whole data without scaffold starts / ends
  se$cor2_genotype_qual[which(se$RAM != 1 | se$crowSat !=1)], # scaffold starts / ends without RAM pattern
  se$cor2_genotype_qual[which(se$RAM==1 )],  # 
  se$cor2_genotype_qual[which(se$crowSat==1)])
  

t.test(dat$cor2_genotype_qual[which(dat$scaf_end != 1)], se$cor2_genotype_qual[which(se$RAM != 1 | se$crowSat !=1)])
t.test(dat$cor2_genotype_qual[which(dat$scaf_end != 1)], se$cor2_genotype_qual[which(se$RAM==1 )])  
t.test(dat$cor2_genotype_qual[which(dat$scaf_end != 1)], se$cor2_genotype_qual[which(se$crowSat==1)])


svg(file="GQ.svg")
par(mfrow=c(4,2), mar=c(2,2,2,2))
plot(dat$swe_genotype_qual[which(dat$scaf_end != 1)], dat$swe_rho_extended[which(dat$scaf_end != 1)], pch=19, cex=0.5, ylim=c(0,1.2), xlim=c(15,90),xlab="", ylab="")
title(main=" Hooded Crow: genotype quality vs. Rho - Genome-wide", cex.main=0.8)
plot(dat$cor2_genotype_qual[which(dat$scaf_end != 1)], dat$cor2_rho_extended[which(dat$scaf_end != 1)], pch=19, cex=0.5, ylim=c(0,1.2), xlim=c(15,90),xlab="", ylab="")
title(main=" Carrion Crow: genotype quality vs. Rho - Genome-wide", cex.main=0.8)

plot(se$swe_genotype_qual[which(se$RAM != 1 | se$crowSat !=1)], se$swe_rho_extended[which(se$RAM != 1 | se$crowSat !=1)], pch=19, cex=0.5, ylim=c(0,1.2),xlim=c(15,90),xlab="",ylab="")
title(main=" Hooded Crow: genotype quality vs. Rho - Scaffold ends", cex.main=0.8)
plot(se$cor2_genotype_qual[which(se$RAM != 1 | se$crowSat !=1)], se$cor2_rho_extended[which(se$RAM != 1 | se$crowSat !=1)], pch=19, cex=0.5, ylim=c(0,1.2),xlim=c(15,90),xlab="",ylab="")
title(main=" Carrion Crow: genotype quality vs. Rho - Scaffold ends", cex.main=0.8)

plot(se$swe_genotype_qual[which(se$RAM==1 )], se$swe_rho_extended[which(se$RAM==1 )], pch=19, cex=0.5,ylim=c(0,1.2), xlim=c(15,90),xlab="",ylab="")
title(main=" Hooded Crow: genotype quality vs. Rho - RAMs", cex.main=0.8)
plot(se$cor2_genotype_qual[which(se$RAM==1 )], se$cor2_rho_extended[which(se$RAM==1 )], pch=19, cex=0.5,ylim=c(0,1.2), xlim=c(15,90),xlab="",ylab="")
title(main=" Carrion Crow: genotype quality vs. Rho - RAMs", cex.main=0.8)

plot(se$swe_genotype_qual[which(se$crowSat==1)], se$swe_rho_extended[which(se$crowSat==1)], pch=19, cex=0.5, ylim=c(0,1.2),xlim=c(15,90), xlab="", ylab="")
title(main=" Hooded Crow: genotype quality vs. Rho - crowSat1", cex.main=0.8)
plot(se$cor2_genotype_qual[which(se$crowSat==1)], se$cor2_rho_extended[which(se$crowSat==1)], pch=19, cex=0.5, ylim=c(0,1.2),xlim=c(15,90), xlab="", ylab="")
title(main=" Carrion Crow: genotype quality vs. Rho - crowSat1", cex.main=0.8)
dev.off()

summary(dat$swe_genotype_qual[which(dat$scaf_end != 1)])
summary(se$swe_genotype_qual[which(se$RAM != 1 | se$crowSat !=1)])
summary(se$swe_genotype_qual[which(se$RAM==1 )])
summary(se$swe_genotype_qual[which(se$crowSat==1)])

summary(dat$cor2_genotype_qual[which(dat$scaf_end != 1)])
summary(se$cor2_genotype_qual[which(se$RAM != 1 | se$crowSat !=1)])
summary(se$cor2_genotype_qual[which(se$RAM==1 )])
summary(se$cor2_genotype_qual[which(se$crowSat==1)])

## test for significance:

#rho
a<-t.test(log(dat$swe_rho_extended[which(dat$scaf_end != 1)]), log(se$swe_rho_extended[which(se$RAM != 1 | se$crowSat !=1)]))
b<-t.test(log(dat$swe_rho_extended[which(dat$scaf_end != 1)]), log(se$swe_rho_extended[which(se$RAM==1 )]))  
c<-t.test(log(dat$swe_rho_extended[which(dat$scaf_end != 1)]), log(se$swe_rho_extended[which(se$crowSat==1 )]))
d<-t.test(log(se$swe_rho_extended[which(se$RAM != 1 | se$crowSat !=1)]), log(se$swe_rho_extended[which(se$RAM==1)]))
e<-t.test(log(se$swe_rho_extended[which(se$RAM != 1 | se$crowSat !=1)]), log(se$swe_rho_extended[which(se$crowSat==1)]))

p.adjust(c(a[3], b[3], c[3], d[3], e[3]), method=p.adjust.methods, n=5)       
       
a<-t.test(log(dat$cor2_rho_extended[which(dat$scaf_end != 1)]),  log(se$cor2_rho_extended[which(se$RAM != 1 | se$crowSat !=1)]))
b<-t.test(log(dat$cor2_rho_extended[which(dat$scaf_end != 1)]),  log(se$cor2_rho_extended[which(se$RAM==1)]))
c<-t.test(log(dat$cor2_rho_extended[which(dat$scaf_end != 1)]),   log(se$cor2_rho_extended[which(se$crowSat==1)]))
d<-t.test(log(se$cor2_rho_extended[which(se$RAM != 1 | se$crowSat !=1)]), log(se$cor2_rho_extended[which(se$RAM==1)]))
e<-t.test(log(se$cor2_rho_extended[which(se$RAM != 1 | se$crowSat !=1)]), log(se$cor2_rho_extended[which(se$crowSat==1)]))

c(a[3], b[3], c[3], d[3], e[3])
p.adjust(c(a[3], b[3], c[3], d[3], e[3]), method=p.adjust.methods, n=5)       

######### FINAL stats model #############

## RHO

# first look at rho SWE only in scaf_end windows

mod_rho_swe_scaf_end<-lmer(log(swe_rho_extended) ~ RAM + crowSat + (1|chr), data = se)
Anova(mod_rho_swe_scaf_end, type=3)

mod_rho_swe_genome_wide<-lmer(log(swe_rho_extended) ~ RAM + crowSat + (1|chr), data = dat)
Anova(mod_rho_swe_genome_wide, type=3)

# now the same in cor2

mod_rho_cor2_scaf_end<-lmer(log(cor2_rho_extended) ~ RAM + crowSat + (1|chr), data = se)
Anova(mod_rho_cor2_scaf_end)

mod_rho_cor2_genome_wide<-lmer(log(cor2_rho_extended) ~ RAM + crowSat + (1|chr), data = dat)
Anova(mod_rho_cor2_genome_wide, type=3)

## THETA

# first look at rho SWE only in scaf_end windows

mod_swe_pol_theta_scaf_end<-lmer(log(swe_pol_theta_extended) ~ RAM + crowSat + (1|chr), data = se)
Anova(mod_swe_pol_theta_scaf_end, type=3)

mod_swe_pol_theta_genome_wide<-lmer(log(swe_pol_theta_extended) ~ RAM + crowSat + (1|chr), data = dat)
Anova(mod_swe_pol_theta_genome_wide, type=3)

# now the same in cor2

mod_theta_cor2_scaf_end<-lmer(log(cor2_theta_extended) ~ RAM + crowSat + (1|chr), data = se)
Anova(mod_theta_cor2_scaf_end)

mod_theta_cor2_genome_wide<-lmer(log(cor2_theta_extended) ~ RAM + crowSat + (1|chr), data = dat)
Anova(mod_theta_cor2_genome_wide, type=3)

## FST

# first look at rho SWE only in scaf_end windows

mod_fst_swe_cor2_scaf_end<-lmer(log(cor2_swe_fst_extended) ~ RAM + crowSat + (1|chr), data = se)
Anova(mod_rho_swe_scaf_end, type=3)

mod_fst_swe_cor2_genome_wide<-lmer(log(cor2_swe_fst_extended) ~ RAM + crowSat + (1|chr), data = dat)
Anova(mod_rho_swe_genome_wide, type=3)


