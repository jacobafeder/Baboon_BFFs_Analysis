## LAIKIPIA AGE ANALYSES
setwd("~/Desktop/SPRF/Data/Chacma_Moremi/SourceData")
library(RPostgreSQL)
library(dbplyr)
library(tidyverse)
library(lubridate)
library(reshape2)
library(dplyr)
theme_set(theme_classic(base_size = 16))

## LOAD DATA FOR BOTH POPS
ff.dyads<-read.csv("BFFs_Data_Dyadic.csv")
ff.dyads$kincat2<-as.factor(ff.dyads$kincat2)
ff.dyads$kincat2<-factor(ff.dyads$kincat2, levels=c("Non-kin","Sisters","Mother-daughter"))

ff.dyads$age_diff<-ff.dyads$id2_age-ff.dyads$id1_age

ff.dyads<-ff.dyads %>%
  group_by(year,group)  %>%
  mutate(N.Females=length(unique(id1)))

range(ff.dyads$N.Females[ff.dyads$Population=="Moremi"])
range(ff.dyads$N.Females[ff.dyads$Population=="Laikipia"])
unique(ff.dyads$year[ff.dyads$Population=="Moremi"])

table(ff.dyads$kincat2[!duplicated(ff.dyads$dyad) & ff.dyads$Population=="Moremi"])
table(ff.dyads$kincat2[!duplicated(ff.dyads$dyad) & ff.dyads$Population=="Laikipia"])

range(ff.dyads$id2_proplact[ff.dyads$Population=="Laikipia"])
mean(ff.dyads$id2_proplact[ff.dyads$Population=="Laikipia"])
mean(ff.dyads$id2_proplact[ff.dyads$Population=="Moremi"])
range(ff.dyads$id2_proplact[ff.dyads$Population=="Moremi"])

unique(ff.dyads$id1[ff.dyads$id1_age>18])

## MODEL FROM PERSPECTIVE OF THE *ACTOR*
## MODEL APPROACHES BY AGE
library(brms)
sum(unique(ff.dyads$id1[ff.dyads$year==1992|ff.dyads$year==1993]) %in%
unique(ff.dyads$id1[ff.dyads$year==2001]))

ggplot(data=ff.dyads, aes(y=appr,x=id1_age, color=kincat2)) + geom_point() + geom_smooth(method="gam")
ggplot(data=ff.dyads, aes(y=appr,x=id1_age, color=Population)) + geom_point() + geom_smooth(method="gam")

ggplot(data=ff.dyads, aes(y=aggr,x=id1_age, color=kincat2)) + geom_point() + geom_smooth(method="gam")
sum(ff.dyads$aggr[ff.dyads$id1_age>ff.dyads$id2_age & ff.dyads$kincat2=="Mother-daughter"])/
sum(ff.dyads$aggr[ff.dyads$kincat2=="Mother-daughter"])

ff.dyads$group[ff.dyads$group=="C Group"]<-"C-Late"
ff.dyads$group[ff.dyads$year==1992|ff.dyads$year==1993]<-"C-Early"
ff.dyads$Group.Year<-paste(ff.dyads$group, ff.dyads$year, sep="_")

model.aging.given<-brm(data=ff.dyads, mvbind(appr,groom,aggr)  ~ 
                            scale(id1_age)*kincat2 + Population + scale(N.Females) +
                            scale(id2_proplact) +
                            offset(log(hrsobs)) + (1|p|mm(id1,id2)) + (1|r|dyad), 
                          prior = c(
                            prior(normal(0, 5), "b",resp="appr"),
                            prior(normal(0, 5), "b",resp="groom"),
                            prior(normal(0, 5), "b",resp="aggr")),  
                          chains=4, cores = 2,
                          family="negbinomial")

summary(model.aging.given)
conditional_effects(model.aging.given, "id1_age", resp="appr")
conditional_effects(model.aging.given, "id1_age:kincat2",resp="appr")
conditional_effects(model.aging.given, "id1_age:kincat2",resp="groom")
conditional_effects(model.aging.given, "id1_age:kincat2",resp="aggr")
conditional_effects(model.aging.given, "N.Females",resp="appr")
conditional_effects(model.aging.given, "N.Females",resp="groom")
conditional_effects(model.aging.given, "N.Females",resp="aggr")
conditional_effects(model.aging.given, "id1_age:kincat2",resp="aggr")
conditional_effects(model.aging.given, "kincat2", resp="appr")
conditional_effects(model.aging.given, "kincat2", resp="groom")
conditional_effects(model.aging.given, "kincat2:Population", resp="appr")
conditional_effects(model.aging.given, "id2_proplact", resp="appr")
conditional_effects(model.aging.given, "Population", resp="appr")
conditional_effects(model.aging.given, "Population", resp="groom")
conditional_effects(model.aging.given, "Population", resp="aggr")

tapply(ff.dyads$aggr, ff.dyads$Group.Year, sum)
tapply(ff.dyads$hrsobs, ff.dyads$Group.Year, sum)/(tapply(ff.dyads$id1, ff.dyads$Group.Year, function(x) length(unique(x)))*(tapply(ff.dyads$id1, ff.dyads$Group.Year, function(x) length(unique(x)))-1))

## MODEL AGING RECEIVED
model.aging.received<-brm(data=ff.dyads, mvbind(appr,groom,aggr)  ~ 
                         scale(id2_age)*kincat2 + Population + scale(N.Females) +
                         scale(id2_proplact) +
                         offset(log(hrsobs)) + (1|p|mm(id1,id2)) + (1|r|dyad), 
                         prior = c(
                           prior(normal(0, 5), "b",resp="appr"),
                           prior(normal(0, 5), "b",resp="groom"),
                           prior(normal(0, 5), "b",resp="aggr")),  
                       chains=4, cores = 2,
                       family="negbinomial")
summary(model.aging.received)
conditional_effects(model.aging.received, "id2_age", resp="groom")
conditional_effects(model.aging.received, "id2_age:kincat2","appr")
conditional_effects(model.aging.received, "id2_age:kincat2","groom")
conditional_effects(model.aging.received, "kincat2", resp="appr")
conditional_effects(model.aging.received, "kincat2", resp="groom")
conditional_effects(model.aging.received, "Population", resp="appr")
conditional_effects(model.aging.received, "Population", resp="aggr")
conditional_effects(model.aging.received, "Population", resp="groom")
conditional_effects(model.aging.received, "id2_proplact", resp="appr")

## PLOT EFFECTS
## MOTHER-DAUGHTER
?conditional_effects
plot1a<-conditional_effects(model.aging.given, "id1_age:kincat2", resp="appr")
plot1a$`appr.appr_id1_age:kincat2`$estimate__<-plot1a$`appr.appr_id1_age:kincat2`$estimate__/plot1a$`appr.appr_id1_age:kincat2`$hrsobs
plot1a$`appr.appr_id1_age:kincat2`$upper__<-plot1a$`appr.appr_id1_age:kincat2`$upper__/plot1a$`appr.appr_id1_age:kincat2`$hrsobs
plot1a$`appr.appr_id1_age:kincat2`$lower__<-plot1a$`appr.appr_id1_age:kincat2`$lower__/plot1a$`appr.appr_id1_age:kincat2`$hrsobs
plot1a$`appr.appr_id1_age:kincat2`<- plot1a$`appr.appr_id1_age:kincat2`[plot1a$`appr.appr_id1_age:kincat2`$kincat2 =="Mother-daughter", ]
plot1b<-conditional_effects(model.aging.received, "id2_age:kincat2")
plot1b$`appr.appr_id2_age:kincat2`$estimate__<-plot1b$`appr.appr_id2_age:kincat2`$estimate__/plot1b$`appr.appr_id2_age:kincat2`$hrsobs
plot1b$`appr.appr_id2_age:kincat2`$upper__<-plot1b$`appr.appr_id2_age:kincat2`$upper__/plot1b$`appr.appr_id2_age:kincat2`$hrsobs
plot1b$`appr.appr_id2_age:kincat2`$lower__<-plot1b$`appr.appr_id2_age:kincat2`$lower__/plot1b$`appr.appr_id2_age:kincat2`$hrsobs
plot1b$`appr.appr_id2_age:kincat2`<- plot1b$`appr.appr_id2_age:kincat2`[plot1b$`appr.appr_id2_age:kincat2`$kincat2 =="Mother-daughter", ]

library(RColorBrewer)
brewer.pal(3,"Dark2")
p1<-plot(plot1a, plot = FALSE, prob=0.89, line_args = list(color="black", alpha=0.6))[[1]] +
  scale_x_continuous(limits=c(3,28), breaks=seq(5,25,5)) + xlab("") +
  scale_y_log10(limits=c(0.03,1.25), breaks=c(0.01,0.05,0.25,1.25)) +
  ggtitle("Mother-daughter") + theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ylab("Approaches (events/hour)") + theme(legend.position = c(0.7,0.75)) + scale_fill_manual(values="#1B9E77") +
  theme(legend.title= element_text(size=12), legend.text = element_text(size=12), legend.position = "none") 

plot1b$`appr.appr_id2_age:kincat2`$estimate__
p1a<-p1 + geom_line(data=plot1b$`appr.appr_id2_age:kincat2`, aes(x=id2_age, y = estimate__), color="black", lwd=1.5) +
  geom_ribbon(data=plot1b$`appr.appr_id2_age:kincat2`, aes(x=id2_age, ymin = lower__, ymax=upper__), color=NA, alpha=0.6, fill="#D95F02")

## Sisters
plot1a<-conditional_effects(model.aging.given, "id1_age:kincat2")
plot1a$`appr.appr_id1_age:kincat2`$estimate__<-plot1a$`appr.appr_id1_age:kincat2`$estimate__/plot1a$`appr.appr_id1_age:kincat2`$hrsobs
plot1a$`appr.appr_id1_age:kincat2`$upper__<-plot1a$`appr.appr_id1_age:kincat2`$upper__/plot1a$`appr.appr_id1_age:kincat2`$hrsobs
plot1a$`appr.appr_id1_age:kincat2`$lower__<-plot1a$`appr.appr_id1_age:kincat2`$lower__/plot1a$`appr.appr_id1_age:kincat2`$hrsobs
plot1a$`appr.appr_id1_age:kincat2`<- plot1a$`appr.appr_id1_age:kincat2`[plot1a$`appr.appr_id1_age:kincat2`$kincat2 =="Sisters", ]

plot1b<-conditional_effects(model.aging.received, "id2_age:kincat2")
plot1b$`appr.appr_id2_age:kincat2`$estimate__<-plot1b$`appr.appr_id2_age:kincat2`$estimate__/plot1b$`appr.appr_id2_age:kincat2`$hrsobs
plot1b$`appr.appr_id2_age:kincat2`$upper__<-plot1b$`appr.appr_id2_age:kincat2`$upper__/plot1b$`appr.appr_id2_age:kincat2`$hrsobs
plot1b$`appr.appr_id2_age:kincat2`$lower__<-plot1b$`appr.appr_id2_age:kincat2`$lower__/plot1b$`appr.appr_id2_age:kincat2`$hrsobs
plot1b$`appr.appr_id2_age:kincat2`<- plot1b$`appr.appr_id2_age:kincat2`[plot1b$`appr.appr_id2_age:kincat2`$kincat2 =="Sisters", ]

brewer.pal(3,"Dark2")
p1<-plot(plot1a, plot = FALSE, prob=0.89, line_args = list(color="black", alpha=0.6))[[1]] + ylim(0,1.0) +
  scale_x_continuous(limits=c(3,28), breaks=seq(5,25,5)) + xlab("") +
  scale_y_log10(limits=c(0.03,1.25), breaks=c(0.01,0.05,0.25,1.25)) +
  ggtitle("Sisters") + theme(plot.title = element_text(hjust=0.5,face="bold")) +
  ylab("") + theme(legend.position = c(0.7,0.75)) + scale_fill_manual(values="#1B9E77") +
  theme(legend.title= element_text(size=12), legend.text = element_text(size=12), legend.position = "none") 

p1b<-p1 + geom_line(data=plot1b$`appr.appr_id2_age:kincat2`, aes(x=id2_age, y = estimate__), color="black", lwd=1.5) +
  geom_ribbon(data=plot1b$`appr.appr_id2_age:kincat2`, aes(x=id2_age, ymin = lower__, ymax=upper__), color=NA, alpha=0.6, fill="#D95F02")

## NON-KIN
plot1a<-conditional_effects(model.aging.given, "id1_age:kincat2")
plot1a$`appr.appr_id1_age:kincat2`$estimate__<-plot1a$`appr.appr_id1_age:kincat2`$estimate__/plot1a$`appr.appr_id1_age:kincat2`$hrsobs
plot1a$`appr.appr_id1_age:kincat2`$upper__<-plot1a$`appr.appr_id1_age:kincat2`$upper__/plot1a$`appr.appr_id1_age:kincat2`$hrsobs
plot1a$`appr.appr_id1_age:kincat2`$lower__<-plot1a$`appr.appr_id1_age:kincat2`$lower__/plot1a$`appr.appr_id1_age:kincat2`$hrsobs
plot1a$`appr.appr_id1_age:kincat2`<- plot1a$`appr.appr_id1_age:kincat2`[plot1a$`appr.appr_id1_age:kincat2`$kincat2 =="Non-kin", ]
plot1b<-conditional_effects(model.aging.received, "id2_age:kincat2")
plot1b$`appr.appr_id2_age:kincat2`$estimate__<-plot1b$`appr.appr_id2_age:kincat2`$estimate__/plot1b$`appr.appr_id2_age:kincat2`$hrsobs
plot1b$`appr.appr_id2_age:kincat2`$upper__<-plot1b$`appr.appr_id2_age:kincat2`$upper__/plot1b$`appr.appr_id2_age:kincat2`$hrsobs
plot1b$`appr.appr_id2_age:kincat2`$lower__<-plot1b$`appr.appr_id2_age:kincat2`$lower__/plot1b$`appr.appr_id2_age:kincat2`$hrsobs

plot1b$`appr.appr_id2_age:kincat2`<- plot1b$`appr.appr_id2_age:kincat2`[plot1b$`appr.appr_id2_age:kincat2`$kincat2 =="Non-kin", ]

brewer.pal(3,"Dark2")
p1<-plot(plot1a, plot = FALSE, prob=0.89, line_args = list(color="black", alpha=0.6))[[1]] + ylim(0,0.5) +
  ggtitle("Non-kin") + theme(plot.title=element_text(hjust=0.5,face="bold")) +
  scale_y_log10(limits=c(0.03,1.25), breaks=c(0.01,0.05,0.25,1.25)) +
  scale_x_continuous(limits=c(3,28), breaks=seq(5,25,5)) + xlab("") + 
  ylab("") + theme(legend.position = c(0.7,0.75)) + scale_fill_manual(values="#1B9E77") + 
  theme(legend.title= element_text(size=12), legend.text = element_text(size=12), legend.position = "none") 

p1c<-p1 + geom_line(data=plot1b$`appr.appr_id2_age:kincat2`, aes(x=id2_age, y = estimate__), color="black", lwd=1.5) +
  geom_ribbon(data=plot1b$`appr.appr_id2_age:kincat2`, aes(x=id2_age, ymin = lower__, ymax=upper__), color=NA, alpha=0.6, fill="#D95F02")

library(patchwork)
p1a+p1b+p1c + plot_annotation(tag_levels = "a")

## PLOT EFFECTS
## MOTHER-DAUGHTER
?conditional_effects
plot2a<-conditional_effects(model.aging.given, "id1_age:kincat2", resp="groom")
plot2a$`groom.groom_id1_age:kincat2`$estimate__<-plot2a$`groom.groom_id1_age:kincat2`$estimate__/plot2a$`groom.groom_id1_age:kincat2`$hrsobs
plot2a$`groom.groom_id1_age:kincat2`$upper__<-plot2a$`groom.groom_id1_age:kincat2`$upper__/plot2a$`groom.groom_id1_age:kincat2`$hrsobs
plot2a$`groom.groom_id1_age:kincat2`$lower__<-plot2a$`groom.groom_id1_age:kincat2`$lower__/plot2a$`groom.groom_id1_age:kincat2`$hrsobs
plot2a$`groom.groom_id1_age:kincat2`<- plot2a$`groom.groom_id1_age:kincat2`[plot2a$`groom.groom_id1_age:kincat2`$kincat2 =="Mother-daughter", ]
plot2b<-conditional_effects(model.aging.received, "id2_age:kincat2", resp="groom")
plot2b$`groom.groom_id2_age:kincat2`$estimate__<-plot2b$`groom.groom_id2_age:kincat2`$estimate__/plot2b$`groom.groom_id2_age:kincat2`$hrsobs
plot2b$`groom.groom_id2_age:kincat2`$upper__<-plot2b$`groom.groom_id2_age:kincat2`$upper__/plot2b$`groom.groom_id2_age:kincat2`$hrsobs
plot2b$`groom.groom_id2_age:kincat2`$lower__<-plot2b$`groom.groom_id2_age:kincat2`$lower__/plot2b$`groom.groom_id2_age:kincat2`$hrsobs
plot2b$`groom.groom_id2_age:kincat2`<- plot2b$`groom.groom_id2_age:kincat2`[plot2b$`groom.groom_id2_age:kincat2`$kincat2 =="Mother-daughter", ]

brewer.pal(3,"Dark2")
p2<-plot(plot2a, plot = FALSE, prob=0.89, line_args = list(color="black", alpha=0.6))[[1]] + ylim(0,0.5) +
  scale_x_continuous(limits=c(3,28), breaks=seq(5,25,5)) + xlab("") +
  scale_y_log10(limits=c(0.001,1), breaks=c(0.001,0.01,0.10,1.00)) +
  ylab("Grooming (events/hour)") + theme(legend.position = c(0.7,0.75)) + scale_fill_manual(values="#1B9E77") +
  theme(legend.title= element_text(size=12), legend.text = element_text(size=12), legend.position = "none") 

plot2b$`groom.groom_id2_age:kincat2`$estimate__
p2a<-p2 + geom_line(data=plot2b$`groom.groom_id2_age:kincat2`, aes(x=id2_age, y = estimate__), color="black", lwd=1.5) +
  geom_ribbon(data=plot2b$`groom.groom_id2_age:kincat2`, aes(x=id2_age, ymin = lower__, ymax=upper__), color=NA, alpha=0.6, fill="#D95F02")

## Sisters
plot2a<-conditional_effects(model.aging.given, "id1_age:kincat2", resp="groom")
plot2a$`groom.groom_id1_age:kincat2`$estimate__<-plot2a$`groom.groom_id1_age:kincat2`$estimate__/plot2a$`groom.groom_id1_age:kincat2`$hrsobs
plot2a$`groom.groom_id1_age:kincat2`$upper__<-plot2a$`groom.groom_id1_age:kincat2`$upper__/plot2a$`groom.groom_id1_age:kincat2`$hrsobs
plot2a$`groom.groom_id1_age:kincat2`$lower__<-plot2a$`groom.groom_id1_age:kincat2`$lower__/plot2a$`groom.groom_id1_age:kincat2`$hrsobs
plot2a$`groom.groom_id1_age:kincat2`<- plot2a$`groom.groom_id1_age:kincat2`[plot2a$`groom.groom_id1_age:kincat2`$kincat2 =="Sisters", ]

plot2b<-conditional_effects(model.aging.received, "id2_age:kincat2", resp="groom")
plot2b$`groom.groom_id2_age:kincat2`$estimate__<-plot2b$`groom.groom_id2_age:kincat2`$estimate__/plot2b$`groom.groom_id2_age:kincat2`$hrsobs
plot2b$`groom.groom_id2_age:kincat2`$upper__<-plot2b$`groom.groom_id2_age:kincat2`$upper__/plot2b$`groom.groom_id2_age:kincat2`$hrsobs
plot2b$`groom.groom_id2_age:kincat2`$lower__<-plot2b$`groom.groom_id2_age:kincat2`$lower__/plot2b$`groom.groom_id2_age:kincat2`$hrsobs
plot2b$`groom.groom_id2_age:kincat2`<- plot2b$`groom.groom_id2_age:kincat2`[plot2b$`groom.groom_id2_age:kincat2`$kincat2 =="Sisters", ]

brewer.pal(3,"Dark2")
p2<-plot(plot2a, plot = FALSE, prob=0.89, line_args = list(color="black", alpha=0.6))[[1]] + ylim(0,0.2) +
  scale_x_continuous(limits=c(3,28), breaks=seq(5,25,5)) + xlab("") +
  scale_y_log10(limits=c(0.001,1), breaks=c(0.001,0.01,0.10,1.00)) +
  ylab("") + theme(legend.position = c(0.7,0.75)) + scale_fill_manual(values="#1B9E77") +
  theme(legend.title= element_text(size=12), legend.text = element_text(size=12), legend.position = "none") 

p2b<-p2 + geom_line(data=plot2b$`groom.groom_id2_age:kincat2`, aes(x=id2_age, y = estimate__), color="black", lwd=1.5) +
  geom_ribbon(data=plot2b$`groom.groom_id2_age:kincat2`, aes(x=id2_age, ymin = lower__, ymax=upper__), color=NA, alpha=0.6, fill="#D95F02")

## NON-KIN
plot2a<-conditional_effects(model.aging.given, "id1_age:kincat2", resp="groom")
plot2a$`groom.groom_id1_age:kincat2`$estimate__<-plot2a$`groom.groom_id1_age:kincat2`$estimate__/plot2a$`groom.groom_id1_age:kincat2`$hrsobs
plot2a$`groom.groom_id1_age:kincat2`$upper__<-plot2a$`groom.groom_id1_age:kincat2`$upper__/plot2a$`groom.groom_id1_age:kincat2`$hrsobs
plot2a$`groom.groom_id1_age:kincat2`$lower__<-plot2a$`groom.groom_id1_age:kincat2`$lower__/plot2a$`groom.groom_id1_age:kincat2`$hrsobs
plot2a$`groom.groom_id1_age:kincat2`<- plot2a$`groom.groom_id1_age:kincat2`[plot2a$`groom.groom_id1_age:kincat2`$kincat2 =="Non-kin", ]
plot2b<-conditional_effects(model.aging.received, "id2_age:kincat2", resp="groom")
plot2b$`groom.groom_id2_age:kincat2`$estimate__<-plot2b$`groom.groom_id2_age:kincat2`$estimate__/plot2b$`groom.groom_id2_age:kincat2`$hrsobs
plot2b$`groom.groom_id2_age:kincat2`$upper__<-plot2b$`groom.groom_id2_age:kincat2`$upper__/plot2b$`groom.groom_id2_age:kincat2`$hrsobs
plot2b$`groom.groom_id2_age:kincat2`$lower__<-plot2b$`groom.groom_id2_age:kincat2`$lower__/plot2b$`groom.groom_id2_age:kincat2`$hrsobs

plot2b$`groom.groom_id2_age:kincat2`<- plot2b$`groom.groom_id2_age:kincat2`[plot2b$`groom.groom_id2_age:kincat2`$kincat2 =="Non-kin", ]

brewer.pal(3,"Dark2")
p2<-plot(plot2a, plot = FALSE, prob=0.89, line_args = list(color="black", alpha=0.6))[[1]] + ylim(0,0.05) +
  scale_x_continuous(limits=c(3,28), breaks=seq(5,25,5)) + xlab("") + 
  scale_y_log10(limits=c(0.001,1), breaks=c(0.001,0.01,0.10,1.00)) +
  ylab("") + theme(legend.position = c(0.7,0.75)) + scale_fill_manual(values="#1B9E77") + 
  theme(legend.title= element_text(size=12), legend.text = element_text(size=12), legend.position = "none") 

p2c<-p2 + geom_line(data=plot2b$`groom.groom_id2_age:kincat2`, aes(x=id2_age, y = estimate__), color="black", lwd=1.5) +
  geom_ribbon(data=plot2b$`groom.groom_id2_age:kincat2`, aes(x=id2_age, ymin = lower__, ymax=upper__), color=NA, alpha=0.6, fill="#D95F02")

library(patchwork)
p2a+p2b+p2c + plot_annotation(tag_levels = "a")

## PLOT EFFECTS
## MOTHER-DAUGHTER
?conditional_effects
plot3a<-conditional_effects(model.aging.given, "id1_age:kincat2", resp="aggr")
plot3a$`aggr.aggr_id1_age:kincat2`$estimate__<-plot3a$`aggr.aggr_id1_age:kincat2`$estimate__/plot3a$`aggr.aggr_id1_age:kincat2`$hrsobs
plot3a$`aggr.aggr_id1_age:kincat2`$upper__<-plot3a$`aggr.aggr_id1_age:kincat2`$upper__/plot3a$`aggr.aggr_id1_age:kincat2`$hrsobs
plot3a$`aggr.aggr_id1_age:kincat2`$lower__<-plot3a$`aggr.aggr_id1_age:kincat2`$lower__/plot3a$`aggr.aggr_id1_age:kincat2`$hrsobs
plot3a$`aggr.aggr_id1_age:kincat2`<- plot3a$`aggr.aggr_id1_age:kincat2`[plot3a$`aggr.aggr_id1_age:kincat2`$kincat2 =="Mother-daughter", ]
plot3b<-conditional_effects(model.aging.received, "id2_age:kincat2", resp="aggr")
plot3b$`aggr.aggr_id2_age:kincat2`$estimate__<-plot3b$`aggr.aggr_id2_age:kincat2`$estimate__/plot3b$`aggr.aggr_id2_age:kincat2`$hrsobs
plot3b$`aggr.aggr_id2_age:kincat2`$upper__<-plot3b$`aggr.aggr_id2_age:kincat2`$upper__/plot3b$`aggr.aggr_id2_age:kincat2`$hrsobs
plot3b$`aggr.aggr_id2_age:kincat2`$lower__<-plot3b$`aggr.aggr_id2_age:kincat2`$lower__/plot3b$`aggr.aggr_id2_age:kincat2`$hrsobs
plot3b$`aggr.aggr_id2_age:kincat2`<- plot3b$`aggr.aggr_id2_age:kincat2`[plot3b$`aggr.aggr_id2_age:kincat2`$kincat2 =="Mother-daughter", ]

brewer.pal(3,"Dark2")
p3<-plot(plot3a, plot = FALSE, prob=0.89, line_args = list(color="black", alpha=0.6))[[1]] + ylim(0,0.25) +
  scale_x_continuous(limits=c(3,28), breaks=seq(5,25,5)) + xlab("Age (years)") +
  scale_y_log10(limits=c(0.0003,1), breaks=c(0.001,0.01,0.10,1.00)) +
  ylab("Aggression (events/hour)") + theme(legend.position = c(0.7,0.75)) + scale_fill_manual(values="#1B9E77") +
  theme(legend.title= element_text(size=12), legend.text = element_text(size=12), legend.position = "none") 

plot3b$`aggr.aggr_id2_age:kincat2`$estimate__
p3a<-p3 + geom_line(data=plot3b$`aggr.aggr_id2_age:kincat2`, aes(x=id2_age, y = estimate__), color="black", lwd=1.5) +
  geom_ribbon(data=plot3b$`aggr.aggr_id2_age:kincat2`, aes(x=id2_age, ymin = lower__, ymax=upper__), color=NA, alpha=0.6, fill="#D95F02")

## Sisters
plot3a<-conditional_effects(model.aging.given, "id1_age:kincat2", resp="aggr")
plot3a$`aggr.aggr_id1_age:kincat2`$estimate__<-plot3a$`aggr.aggr_id1_age:kincat2`$estimate__/plot3a$`aggr.aggr_id1_age:kincat2`$hrsobs
plot3a$`aggr.aggr_id1_age:kincat2`$upper__<-plot3a$`aggr.aggr_id1_age:kincat2`$upper__/plot3a$`aggr.aggr_id1_age:kincat2`$hrsobs
plot3a$`aggr.aggr_id1_age:kincat2`$lower__<-plot3a$`aggr.aggr_id1_age:kincat2`$lower__/plot3a$`aggr.aggr_id1_age:kincat2`$hrsobs
plot3a$`aggr.aggr_id1_age:kincat2`<- plot3a$`aggr.aggr_id1_age:kincat2`[plot3a$`aggr.aggr_id1_age:kincat2`$kincat2 =="Sisters", ]

plot3b<-conditional_effects(model.aging.received, "id2_age:kincat2", resp="aggr")
plot3b$`aggr.aggr_id2_age:kincat2`$estimate__<-plot3b$`aggr.aggr_id2_age:kincat2`$estimate__/plot3b$`aggr.aggr_id2_age:kincat2`$hrsobs
plot3b$`aggr.aggr_id2_age:kincat2`$upper__<-plot3b$`aggr.aggr_id2_age:kincat2`$upper__/plot3b$`aggr.aggr_id2_age:kincat2`$hrsobs
plot3b$`aggr.aggr_id2_age:kincat2`$lower__<-plot3b$`aggr.aggr_id2_age:kincat2`$lower__/plot3b$`aggr.aggr_id2_age:kincat2`$hrsobs
plot3b$`aggr.aggr_id2_age:kincat2`<- plot3b$`aggr.aggr_id2_age:kincat2`[plot3b$`aggr.aggr_id2_age:kincat2`$kincat2 =="Sisters", ]

brewer.pal(3,"Dark2")
p3<-plot(plot3a, plot = FALSE, prob=0.89, line_args = list(color="black", alpha=0.6))[[1]] + ylim(0,0.15) +
  scale_x_continuous(limits=c(3,28), breaks=seq(5,25,5)) + xlab("Age (years)") +
  scale_y_log10(limits=c(0.0003,1), breaks=c(0.001,0.01,0.10,1.00)) +
  ylab("") + theme(legend.position = c(0.7,0.75)) + scale_fill_manual(values="#1B9E77") +
  theme(legend.title= element_text(size=12), legend.text = element_text(size=12), legend.position = "none") 

p3b<-p3 + geom_line(data=plot3b$`aggr.aggr_id2_age:kincat2`, aes(x=id2_age, y = estimate__), color="black", lwd=1.5) +
  geom_ribbon(data=plot3b$`aggr.aggr_id2_age:kincat2`, aes(x=id2_age, ymin = lower__, ymax=upper__), color=NA, alpha=0.6, fill="#D95F02")

## NON-KIN
plot3a<-conditional_effects(model.aging.given, "id1_age:kincat2", resp="aggr")
plot3a$`aggr.aggr_id1_age:kincat2`$estimate__<-plot3a$`aggr.aggr_id1_age:kincat2`$estimate__/plot3a$`aggr.aggr_id1_age:kincat2`$hrsobs
plot3a$`aggr.aggr_id1_age:kincat2`$upper__<-plot3a$`aggr.aggr_id1_age:kincat2`$upper__/plot3a$`aggr.aggr_id1_age:kincat2`$hrsobs
plot3a$`aggr.aggr_id1_age:kincat2`$lower__<-plot3a$`aggr.aggr_id1_age:kincat2`$lower__/plot3a$`aggr.aggr_id1_age:kincat2`$hrsobs
plot3a$`aggr.aggr_id1_age:kincat2`<- plot3a$`aggr.aggr_id1_age:kincat2`[plot3a$`aggr.aggr_id1_age:kincat2`$kincat2 =="Non-kin", ]
plot3b<-conditional_effects(model.aging.received, "id2_age:kincat2", resp="aggr")
plot3b$`aggr.aggr_id2_age:kincat2`$estimate__<-plot3b$`aggr.aggr_id2_age:kincat2`$estimate__/plot3b$`aggr.aggr_id2_age:kincat2`$hrsobs
plot3b$`aggr.aggr_id2_age:kincat2`$upper__<-plot3b$`aggr.aggr_id2_age:kincat2`$upper__/plot3b$`aggr.aggr_id2_age:kincat2`$hrsobs
plot3b$`aggr.aggr_id2_age:kincat2`$lower__<-plot3b$`aggr.aggr_id2_age:kincat2`$lower__/plot3b$`aggr.aggr_id2_age:kincat2`$hrsobs

plot3b$`aggr.aggr_id2_age:kincat2`<- plot3b$`aggr.aggr_id2_age:kincat2`[plot3b$`aggr.aggr_id2_age:kincat2`$kincat2 =="Non-kin", ]

brewer.pal(3,"Dark2")
p3<-plot(plot3a, plot = FALSE, prob=0.89, line_args = list(color="black", alpha=0.6))[[1]] + ylim(0,0.08) +
  scale_x_continuous(limits=c(3,28), breaks=seq(5,25,5)) + xlab("Age (years)") + 
  scale_y_log10(limits=c(0.0003,1), breaks=c(0.001,0.01,0.10,1.00)) +
  ylab("") + scale_fill_manual(values="#1B9E77", labels=c("Given","Received"))  +
  labs(fill="Direction") + guides(fill = guide_legend(override.aes = list(linetype = 0))) +
  theme(legend.title= element_text(size=14, face="bold"), legend.text = element_text(size=12), legend.position=c(0.3,0.7)) 

p3c<-p3 + geom_line(data=plot3b$`aggr.aggr_id2_age:kincat2`, aes(x=id2_age, y = estimate__, group=kincat2), color="black", lwd=1.5) +
  ggnewscale::new_scale_fill() +  ggnewscale::new_scale_color() +
  theme(legend.title= element_blank(), legend.text = element_text(size=14), legend.position="right",
        legend.box = "vertical", legend.spacing = unit(0, "mm"), legend.key.size = unit(7,"mm")) +
  geom_ribbon(data=plot3b$`aggr.aggr_id2_age:kincat2`, aes(x=id2_age, ymin = lower__, ymax=upper__, fill=kincat2, color=kincat2), alpha=0.6, color=NA) +
  scale_color_manual(values="#D95F02", labels=c("Received")) + labs(fill="", color="") +
  scale_fill_manual(values="#D95F02", labels=c("Received")) 

library(patchwork)
p3a+p3b+p3c + plot_annotation(tag_levels = "a")

(p1a+p1b+p1c)/
  (p2a+p2b+p2c)/
  (p3a+p3b+p3c) +plot_annotation(tag_levels = "a")
ggsave(file="BFFs_Figure 2.jpg", units="cm", width=30, height=27, dpi=300)

ggplot(data=ff.dyads, aes(x=id1_age, y=aggr, color=kincat2)) + geom_smooth(method = "glm", family="poisson") + geom_point()
ggplot(data=ff.dyads, aes(x=id2_age, y=aggr, color=kincat2)) + geom_smooth(method = "glm", family="poisson") + geom_point()

## MERGE INTO PROPORTION DATASET
ff.props<- ff.dyads

ff.props$focal<-NA
ff.props$partner<-NA
ff.props$partner_lact<-NA
ff.props$focal_age<-NA
ff.props$partner_age<-NA

i<-1
for (i in 1:nrow(ff.props)) {
  names<-c(ff.props$id1[i], ff.props$id2[i])
  first<-sample(names,1)
  second<-names[names!=first]
  ff.props$dyad[ff.props$dyad==ff.props$dyad[i] & ff.props$year==ff.props$year[i] & ff.props$group==ff.props$group[i]]<-paste(first,second,sep="_")
  ff.props$focal[ff.props$dyad==ff.props$dyad[i] & ff.props$year==ff.props$year[i] & ff.props$group==ff.props$group[i]]<-first
  ff.props$partner[ff.props$dyad==ff.props$dyad[i] & ff.props$year==ff.props$year[i] & ff.props$group==ff.props$group[i]]<-second
  ff.props$focal_rank[ff.props$dyad==ff.props$dyad[i] & ff.props$year==ff.props$year[i] & ff.props$group==ff.props$group[i]]<-ff.props$rank1[first==ff.props$id1 & second==ff.props$id2 & ff.props$group==ff.props$group[i] & ff.props$year[i]==ff.props$year]
  ff.props$partner_rank[ff.props$dyad==ff.props$dyad[i] & ff.props$year==ff.props$year[i] & ff.props$group==ff.props$group[i]]<-ff.props$rank2[first==ff.props$id1 & second==ff.props$id2 & ff.props$group==ff.props$group[i] & ff.props$year[i]==ff.props$year]
  ff.props$partner_lact[ff.props$dyad==ff.props$dyad[i] & ff.props$year==ff.props$year[i] & ff.props$group==ff.props$group[i]]<-ff.props$id2_proplact[first==ff.props$id1 & second==ff.props$id2 & ff.props$group==ff.props$group[i] & ff.props$year[i]==ff.props$year]
  ff.props$focal_age[ff.props$dyad==ff.props$dyad[i] & ff.props$year==ff.props$year[i] & ff.props$group==ff.props$group[i]]<-ff.props$id1_age[first==ff.props$id1 & second==ff.props$id2 & ff.props$group==ff.props$group[i] & ff.props$year[i]==ff.props$year]
  ff.props$partner_age[ff.props$dyad==ff.props$dyad[i] & ff.props$year==ff.props$year[i] & ff.props$group==ff.props$group[i]]<-ff.props$id2_age[ff.props$focal[i]==ff.props$id1 & ff.props$partner[i]==ff.props$id2 & ff.props$group==ff.props$group[i] & ff.props$year[i]==ff.props$year]
  }

warnings()
colnames(ff.props)
ff.props2<-ff.props %>%
  group_by(dyad, focal, partner, kincat2, Population,group,year, focal_age, partner_age, focal_rank, partner_rank, partner_lact) %>%
  summarise(a1=sum(appr[id1==focal]), a.total=sum(appr),
            g1=sum(groom[id1==focal]), g.total=sum(groom),
            agg1=sum(aggr[id1==focal]), agg.total=sum(aggr))

sum(ff.props2$a1)/sum(ff.props2$a.total)
sum(ff.props2$a1[ff.props2$focal_rank>ff.props2$partner_rank])/sum(ff.props2$a.total[ff.props2$focal_rank>ff.props2$partner_rank])
sum(ff.props2$a1[ff.props2$partner_lact>0.5])/sum(ff.props2$a.total[ff.props2$partner_lact>0.5])

ff.props2$rankdiff<-ff.props2$focal_rank-ff.props2$partner_rank
ff.props2$agediff<-ff.props2$focal_age-ff.props2$partner_age
ff.props2$agecat<-"Younger"
ff.props2$agecat[ff.props2$focal_age>ff.props2$partner_age]<-"Older"

sum(ff.props2$g1[ff.props2$kincat2=="Sisters" & ff.props2$focal_age>ff.props2$partner_age])/sum(ff.props2$g.total[ff.props2$kincat2=="Sisters" & ff.props2$focal_age>ff.props2$partner_age])
sum(ff.props2$g1[ff.props2$kincat2=="Sisters" & ff.props2$focal_age<ff.props2$partner_age])/sum(ff.props2$g.total[ff.props2$kincat2=="Sisters" & ff.props2$focal_age<ff.props2$partner_age])

sum(ff.props2$g1[ff.props2$kincat2=="Mother-daughter" & ff.props2$focal_age>ff.props2$partner_age])/sum(ff.props2$g.total[ff.props2$kincat2=="Mother-daughter" & ff.props2$focal_age>ff.props2$partner_age])
sum(ff.props2$g1[ff.props2$kincat2=="Mother-daughter" & ff.props2$focal_age<ff.props2$partner_age])/sum(ff.props2$g.total[ff.props2$kincat2=="Mother-daughter" & ff.props2$focal_age<ff.props2$partner_age])

sum(ff.props2$g1[ff.props2$kincat2=="Non-kin" & ff.props2$focal_age>ff.props2$partner_age])/sum(ff.props2$g.total[ff.props2$kincat2=="Non-kin" & ff.props2$focal_age>ff.props2$partner_age])
sum(ff.props2$g1[ff.props2$kincat2=="Non-kin" & ff.props2$focal_age<ff.props2$partner_age])/sum(ff.props2$g.total[ff.props2$kincat2=="Non-kin" & ff.props2$focal_age<ff.props2$partner_age])

ff.props2$agecat<-as.factor(ff.props2$agecat)
ff.props2$agecat<-factor(ff.props2$agecat, levels = c("Younger","Older"))

## INFER PROPORTIONS FROM GIVEN AND RECEIVED
## INSTEAD OF MODELING DIRECTLY

palette<-c("gray30", "#21908CFF","#440154FF")

## APPROACHES
newdata<-data.frame(id1_age=rep(seq(min(ff.dyads$id1_age), max(ff.dyads$id1_age)),3),
                    kincat2=c(rep("Non-kin",23), rep("Sisters",23), rep("Mother-daughter",23)),
                    N.Females=mean(ff.dyads$N.Females),
                    id2_proplact=mean(ff.dyads$id2_proplact),
                    Population="Moremi",
                    hrsobs=1)
predicted.given<-fitted(model.aging.given, newdata=newdata, re_formula = NA, resp="appr", summary = F)
predicted.given

newdata2<-data.frame(id2_age=rep(seq(min(ff.dyads$id2_age), max(ff.dyads$id2_age)),3),
                    kincat2=c(rep("Non-kin",23), rep("Sisters",23), rep("Mother-daughter",23)),
                    N.Females=mean(ff.dyads$N.Females),
                    id2_proplact=mean(ff.dyads$id2_proplact),
                    Population="Moremi",
                    hrsobs=1)

predicted.received<-fitted(model.aging.received, newdata=newdata2, re_formula = NA, resp="appr", summary = F)
predicted.received

predicted.props<-predicted.given/(predicted.given+predicted.received)

newdata$Prop.Given<-colMeans(predicted.props)
newdata$Prop.Given.LL<-apply(predicted.props, 2, quantile, probs = c(0.025))
newdata$Prop.Given.UL<-apply(predicted.props, 2, quantile, probs = c(0.975))

newdata$kincat2<-as.factor(newdata$kincat2)
newdata$kincat2<-factor(newdata$kincat2,levels=c("Non-kin","Sisters","Mother-daughter"))

prop.plot.1a<-ggplot(data=newdata[!newdata$kincat2=="Sisters",], aes(x=id1_age, y=Prop.Given, fill=kincat2)) + ylim(0,1) +
  geom_ribbon(aes(x=id1_age, ymin=Prop.Given.LL, ymax=Prop.Given.UL), alpha=0.5) + geom_line() +
  ylab("Proportion by focal female") + xlab("Age (years)") + 
  ggtitle("Approaches") + 
  scale_fill_manual(values = palette[c(1,3)])  +
  geom_hline(yintercept=0.5, lty=2) +
  labs(fill="Kin category") + theme(legend.position="none", plot.title =element_text(face="bold", hjust=0.5)) +
  scale_x_continuous(limits=c(3.5,27), breaks=seq(5,25,5))
prop.plot.1a

prop.plot.1b<-ggplot(data=newdata[!newdata$kincat2=="Mother-daughter",], aes(x=id1_age, y=Prop.Given, fill=kincat2)) + ylim(0,1) +
  geom_ribbon(aes(x=id1_age, ymin=Prop.Given.LL, ymax=Prop.Given.UL), alpha=0.5) + geom_line() +
  ylab("Proportion by focal female") + xlab("Age (years)") + 
  scale_fill_manual(values = palette[c(1,2)])  +
  geom_hline(yintercept=0.5, lty=2) +
  labs(fill="Kin category") + theme(legend.position="none") +
  scale_x_continuous(limits=c(3.5,27), breaks=seq(5,25,5))
prop.plot.1b

## GROOMING
newdata<-data.frame(id1_age=rep(seq(min(ff.dyads$id1_age), max(ff.dyads$id1_age)),3),
                    kincat2=c(rep("Non-kin",23), rep("Sisters",23), rep("Mother-daughter",23)),
                    N.Females=mean(ff.dyads$N.Females),
                    id2_proplact=mean(ff.dyads$id2_proplact),
                    Population="Moremi",
                    hrsobs=1)
predicted.given<-fitted(model.aging.given, newdata=newdata, re_formula = NA, resp="groom", summary = F)
predicted.given

newdata2<-data.frame(id2_age=rep(seq(min(ff.dyads$id2_age), max(ff.dyads$id2_age)),3),
                     kincat2=c(rep("Non-kin",23), rep("Sisters",23), rep("Mother-daughter",23)),
                     N.Females=mean(ff.dyads$N.Females),
                     id2_proplact=mean(ff.dyads$id2_proplact),
                     Population="Moremi",
                     hrsobs=1)

predicted.received<-fitted(model.aging.received, newdata=newdata2, re_formula = NA, resp="groom", summary = F)
predicted.received

predicted.props<-predicted.given/(predicted.given+predicted.received)

newdata$Prop.Given<-colMeans(predicted.props)
newdata$Prop.Given.LL<-apply(predicted.props, 2, quantile, probs = c(0.025))
newdata$Prop.Given.UL<-apply(predicted.props, 2, quantile, probs = c(0.975))

newdata$kincat2<-as.factor(newdata$kincat2)
newdata$kincat2<-factor(newdata$kincat2,levels=c("Non-kin","Sisters","Mother-daughter"))

prop.plot.2a<-ggplot(data=newdata[!newdata$kincat2=="Sisters",], aes(x=id1_age, y=Prop.Given, fill=kincat2)) + ylim(0,1) +
  geom_ribbon(aes(x=id1_age, ymin=Prop.Given.LL, ymax=Prop.Given.UL), alpha=0.5) + geom_line() +
  ylab("") + xlab("Age (years)") + 
  ggtitle("Grooming") + 
  scale_fill_manual(values = palette[c(1,3)])  +
  geom_hline(yintercept=0.5, lty=2) +
  labs(fill="Kin category") + theme(legend.position="none", plot.title =element_text(face="bold", hjust=0.5)) +
  scale_x_continuous(limits=c(3.5,27), breaks=seq(5,25,5))
prop.plot.2a

prop.plot.2b<-ggplot(data=newdata[!newdata$kincat2=="Mother-daughter",], aes(x=id1_age, y=Prop.Given, fill=kincat2)) + ylim(0,1) +
  geom_ribbon(aes(x=id1_age, ymin=Prop.Given.LL, ymax=Prop.Given.UL), alpha=0.5) + geom_line() +
  ylab("") + xlab("Age (years)") + 
  scale_fill_manual(values = palette[c(1,2)])  +
  geom_hline(yintercept=0.5, lty=2) +
  labs(fill="Kin category") + theme(legend.position="none") +
  scale_x_continuous(limits=c(3.5,27), breaks=seq(5,25,5))
prop.plot.2b

## AGGRESSION
newdata<-data.frame(id1_age=rep(seq(min(ff.dyads$id1_age), max(ff.dyads$id1_age)),3),
                    kincat2=c(rep("Non-kin",23), rep("Sisters",23), rep("Mother-daughter",23)),
                    N.Females=mean(ff.dyads$N.Females),
                    id2_proplact=mean(ff.dyads$id2_proplact),
                    Population="Moremi",
                    hrsobs=1)
predicted.given<-fitted(model.aging.given, newdata=newdata, re_formula = NA, resp="aggr", summary = F)
predicted.given

newdata2<-data.frame(id2_age=rep(seq(min(ff.dyads$id2_age), max(ff.dyads$id2_age)),3),
                     kincat2=c(rep("Non-kin",23), rep("Sisters",23), rep("Mother-daughter",23)),
                     N.Females=mean(ff.dyads$N.Females),
                     id2_proplact=mean(ff.dyads$id2_proplact),
                     Population="Moremi",
                     hrsobs=1)

predicted.received<-fitted(model.aging.received, newdata=newdata2, re_formula = NA, resp="aggr", summary = F)
predicted.received

predicted.props<-predicted.given/(predicted.given+predicted.received)

newdata$Prop.Given<-colMeans(predicted.props)
newdata$Prop.Given.LL<-apply(predicted.props, 2, quantile, probs = c(0.025))
newdata$Prop.Given.UL<-apply(predicted.props, 2, quantile, probs = c(0.975))

newdata$kincat2<-as.factor(newdata$kincat2)
newdata$kincat2<-factor(newdata$kincat2,levels=c("Non-kin","Sisters","Mother-daughter"))

prop.plot3a<-ggplot(data=newdata[!newdata$kincat2=="Sisters",], aes(x=id1_age, y=Prop.Given, fill=kincat2)) + ylim(0,1) +
  geom_ribbon(aes(x=id1_age, ymin=Prop.Given.LL, ymax=Prop.Given.UL), alpha=0.5, show.legend = T) + geom_line() +
  ylab("") + xlab("Age (years)") + 
  ggtitle("Aggression") + 
  scale_fill_manual(values = palette,drop=F)  +
  geom_hline(yintercept=0.5, lty=2) +
  labs(fill="Kin category") + theme(plot.title =element_text(face="bold", hjust=0.5)) +
  scale_x_continuous(limits=c(3.5,27), breaks=seq(5,25,5))
prop.plot3a

prop.plot3b<-ggplot(data=newdata[!newdata$kincat2=="Mother-daughter",], aes(x=id1_age, y=Prop.Given, fill=kincat2)) + ylim(0,1) +
  geom_ribbon(aes(x=id1_age, ymin=Prop.Given.LL, ymax=Prop.Given.UL), alpha=0.5) + geom_line() +
  ylab("") + xlab("Age (years)") + 
  scale_fill_manual(values = palette[c(1,2)])  +
  geom_hline(yintercept=0.5, lty=2) +
  labs(fill="Kin category") + theme(legend.position = "none") +
  scale_x_continuous(limits=c(3.5,27), breaks=seq(5,25,5))
prop.plot3b

prop.plot.1a + prop.plot.2a + prop.plot3a + 
prop.plot.1b + prop.plot.2b + prop.plot3b + plot_annotation(tag_levels = "a")
ggsave(file="BFFs_Supp_Plot.jpg", units="cm", width=34, height=18, dpi=300)

## JUST LAIKIPIA DATA
model.aging.given.anubis<-brm(data=ff.dyads[ff.dyads$Population=="Laikipia",], mvbind(appr,groom,aggr)  ~ 
                                scale(id1_age)*kincat2 + group +
                                scale(id2_proplact) +(1|p|mm(id1,id2))+ (1|r|dyad), 
                              chains=4, cores = 2, save_all_pars = T,
                              family="negbinomial")
summary(model.aging.given.anubis)

## JUST CHACMA DATA
model.aging.given.chacma<-brm(data=ff.dyads[ff.dyads$Population=="Moremi",], mvbind(appr,groom,aggr) ~ 
                                scale(id1_age)*kincat2 + group +
                                scale(id2_proplact) +
                                offset(log(hrsobs)) + (1|p|mm(id1,id2)) + (1|r|dyad), 
                              chains=4, cores = 2, save_all_pars = T,
                              family="negbinomial")
summary(model.aging.given.chacma)

## JUST LAIKIPIA DATA
model.aging.received.anubis<-brm(data=ff.dyads[ff.dyads$Population=="Laikipia",], mvbind(appr,groom,aggr)  ~ 
                                scale(id2_age)*kincat2 + group +
                                scale(id2_proplact) +(1|p|mm(id1,id2))+ (1|r|dyad), 
                              chains=4, cores = 2, save_all_pars = T,
                              family="negbinomial")
summary(model.aging.received.anubis)

## JUST CHACMA DATA
model.aging.received.chacma<-brm(data=ff.dyads[ff.dyads$Population=="Moremi",], mvbind(appr,groom,aggr) ~ 
                                scale(id2_age)*kincat2 + group +
                                scale(id2_proplact) +
                                offset(log(hrsobs)) + (1|p|mm(id1,id2)) + (1|r|dyad), 
                              chains=4, cores = 2, save_all_pars = T,
                              family="negbinomial")
summary(model.aging.received.chacma)

## COMPARE EFFECTS
conditional_effects(model.aging.given.anubis, "id1_age:kincat2",resp="appr")
conditional_effects(model.aging.given.chacma, "id1_age:kincat2",resp="appr")
conditional_effects(model.aging.given.anubis, "kincat2",resp="appr")
conditional_effects(model.aging.given.chacma, "kincat2",resp="appr")
conditional_effects(model.aging.given.anubis, "id2_proplact",resp="appr")
conditional_effects(model.aging.given.chacma, "id2_proplact",resp="appr")

##
conditional_effects(model.aging.given.anubis, "id1_age:kincat2",resp="groom")
conditional_effects(model.aging.given.chacma, "id1_age:kincat2",resp="groom")
conditional_effects(model.aging.given.anubis, "kincat2",resp="groom")
conditional_effects(model.aging.given.chacma, "kincat2",resp="groom")

##
conditional_effects(model.aging.given.anubis, "id1_age:kincat2",resp="aggr")
conditional_effects(model.aging.given.chacma, "id1_age:kincat2",resp="aggr")
conditional_effects(model.aging.given.anubis, "kincat2",resp="aggr")
conditional_effects(model.aging.given.chacma, "kincat2",resp="aggr")

## APPROACH EFFECTS
fixef(model.aging.given.anubis)[4,]
fixef(model.aging.given.chacma)[4,]
fixef(model.aging.given)[4,]

fixef(model.aging.given.anubis)[5,]
fixef(model.aging.given.chacma)[5,]
fixef(model.aging.given)[5,]

fixef(model.aging.given.anubis)[6,]
fixef(model.aging.given.chacma)[6,]
fixef(model.aging.given)[6,]

fixef(model.aging.given.anubis)[9,]
fixef(model.aging.given.chacma)[8,]
fixef(model.aging.given)[9,]

## MAKE FOREST PLOT
beta.ests1<-as.data.frame(fixef(model.aging.given, probs=c(0.025,0.25,0.75,0.975)))
beta.ests1$behav<-sub("^(([^_]*_){0}[^_]*).*", "\\1", rownames(beta.ests1))
beta.ests1$effect<-c(rep("Intercept",3),
                     rep(c("Actor age","Sister","Mother-daughter","Population (Moremi)","# of females",
                         "Recipient lactating", "Sister x Actor age", "Mother-daughter x Actor age"),3))


beta.ests1<-beta.ests1[!beta.ests1$effect=="Intercept",]
beta.ests1$effect<-as.factor(beta.ests1$effect)
beta.ests1$effect<-factor(beta.ests1$effect, levels=c("Mother-daughter x Actor age","Sister x Actor age","Recipient lactating","# of females","Population (Moremi)", "Mother-daughter","Sister", "Actor age"
                                                       ))
beta.ests1$behav<-as.factor(beta.ests1$behav)
beta.ests1$behav<-factor(beta.ests1$behav, levels=c("appr","groom","aggr"))
levels(beta.ests1$behav)<-c("Approaches","Grooming","Aggression")

forest.plota<-ggplot(data=beta.ests1,aes(y=effect,x=Estimate)) + 
  geom_errorbarh(aes(xmin = Q2.5, xmax=Q97.5),width=0, cex=2.5,color="#58e3b9") + 
  theme_linedraw() +
  geom_errorbarh(aes(xmin = Q25, xmax=Q75),width=0, cex=2.5,color="#1B9E77") + 
  facet_wrap(vars(behav), nrow=3) + geom_vline(xintercept=0, lty=2) + ggtitle("Given effects") +
  theme(plot.title = element_text(face="bold", hjust=0.5, size=16), axis.title = element_text(face="bold", hjust=0.5, size=14)) + ylab("Covariate")

## MAKE FOREST PLOT
## RECEIVED
beta.ests2<-as.data.frame(fixef(model.aging.received, probs=c(0.025,0.25,0.75,0.975)))
beta.ests2$behav<-sub("^(([^_]*_){0}[^_]*).*", "\\1", rownames(beta.ests2))
beta.ests2$effect<-c(rep("Intercept",3),
                     rep(c("Recipient age","Sister","Mother-daughter","Population (Moremi)","# of females",
                           "Recipient lactating", "Sister x Recipient age", "Mother-daughter x Recipient age"),3))


beta.ests2<-beta.ests2[!beta.ests2$effect=="Intercept",]
beta.ests2$effect<-as.factor(beta.ests2$effect)
beta.ests2$effect<-factor(beta.ests2$effect, levels=c("Mother-daughter x Recipient age","Sister x Recipient age","Recipient lactating","# of females","Population (Moremi)", "Mother-daughter","Sister", "Recipient age"
))

beta.ests2$behav<-as.factor(beta.ests2$behav)
beta.ests2$behav<-factor(beta.ests2$behav, levels=c("appr","groom","aggr"))
levels(beta.ests2$behav)<-c("Approaches","Grooming","Aggression")

forest.plotb<-ggplot(data=beta.ests2,aes(y=effect,x=Estimate)) + 
  geom_errorbarh(aes(xmin = Q2.5, xmax=Q97.5),width=0, cex=2.5,color="#fd9c52") + 
  theme_linedraw() +
  geom_errorbarh(aes(xmin = Q25, xmax=Q75),width=0, cex=2.5,color="#D95F02") + 
  facet_wrap(vars(behav), nrow=3) + geom_vline(xintercept=0, lty=2) + ggtitle("Received effects") +
  theme(plot.title = element_text(face="bold", hjust=0.5, size=16), axis.title = element_text(face="bold", hjust=0.5, size=14)) + ylab("")

forest.plota + forest.plotb
ggsave(file="BFFs_Figure 1.jpg", units="cm", width=21, height=16, dpi=300)

## PRINT MODEL OUTPUTS
sink("AGING_MODELS.txt")
print(summary(model.aging.given, prob=0.89), digits=3)
print(summary(model.aging.received, prob=0.89), digits=3)
print(summary(model.aging.given.anubis, prob=0.89), digits=3)
print(summary(model.aging.given.chacma, prob=0.89), digits=3)
print(summary(model.aging.received.anubis, prob=0.89), digits=3)
print(summary(model.aging.received.chacma, prob=0.89), digits=3)
sink()
