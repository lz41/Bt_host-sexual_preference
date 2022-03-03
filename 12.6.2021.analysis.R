#################################################################################
############ host and sexual preference, June 2021 by Linyi Zhang #################
###################################################################################
########################### package to load ##########################
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(bbmle)
library(plyr)
library(patchwork)
library(grid)
library(cowplot)
library(stringr)
library(car)
library(ggpattern)
####################################################################################
################# read host preference data ############
host.preference<-read.csv("Bt host preference new.csv")
######### curate host preference data ##########
head(host.preference)
host.preference$host.plant<-"Qv"
host.preference$host.plant[host.preference$HP=="Qg"]<-"Qg"
host.preference$Block<-paste(host.preference$Year,host.preference$Block,sep="_")
host.preference$species<-"B. treatae"
host.preference$species[host.preference$HP=="Qv-T"]<-"B. kinseyi"
host.preference$species[host.preference$HP=="Qg"]<-"B. fossoria"
head(host.preference)
#######################################################################################################
#### define a function to calculate the amount of time spend on each host plant ####
H.pre1<-function(data,n){
  for(i in 1:length(data[,1])){
    data$Qv.Times[i]<-length(which(data[i,(14+n):34]=="Qv"))
    data$Qg.Times[i]<-length(which(data[i,(14+n):34]=="Qg"))
  }
  data
}  

host.pre<-H.pre1(host.preference,1)
#### calculate the preference = #times on native host/(times on native host+times on non-native host) ####
host.pre$preference[host.pre$host.plant=="Qv"]<-
  host.pre$Qv.Times[host.pre$host.plant=="Qv"]/(host.pre$Qv.Times[host.pre$host.plant=="Qv"]+host.pre$Qg.Times[host.pre$host.plant=="Qv"])

host.pre$preference[host.pre$host.plant=="Qg"]<-
  host.pre$Qg.Times[host.pre$host.plant=="Qg"]/(host.pre$Qv.Times[host.pre$host.plant=="Qg"]+host.pre$Qg.Times[host.pre$host.plant=="Qg"])

host.pre$Total.times<-host.pre$Qv.Times+host.pre$Qg.Times
#### calculate the latency = time takes to land on one host plant ####
for(i in 1:length(host.pre[,1])){
  host.pre$latency1[i]<-which(host.pre[i,15:34]!="C")[1]
}

for(i in 1:length(host.pre[,1])){
  host.pre$latency2[i]<-which(host.pre[i,15:34]==host.pre[i,"host.plant"])[1]
}
###############################################################################
########## overview of the distribution of host preference ####################
ggplot(data=host.pre,aes(x=preference))+
  facet_grid(vars(Sex),vars(HP))+
  geom_histogram()
ggplot(data=host.pre,aes(x=preference))+
  facet_grid(vars(Method),vars(HP))+
  geom_histogram()
################################################################################
################# host preference, only females #################
host.pre$HP<-factor(host.pre$HP,levels=c("Qv-T","Qv-F","Qg"))
host.pre$species[host.pre$HP=="Qv-T"]<-"B. kinseyi"
host.pre$species[host.pre$HP=="Qv-F"]<-"B. treatae"
host.pre$species[host.pre$HP=="Qg"]<-"B. fossoria"
host.pre$Year<-as.factor(host.pre$Year)
host.pre$count<-1

host.pre.new<-host.pre[host.pre$Total.times!=0,]
host.pre.new$preference[host.pre.new$preference==0]<-host.pre.new$preference[host.pre.new$preference==0]+0.0000001
host.pre.new$preference[host.pre.new$preference==1]<-host.pre.new$preference[host.pre.new$preference==1]-0.0000001
host.pre.new<-host.pre.new[!is.na(host.pre.new$Source),]
host.pre.new$pre.times[host.pre.new$host.plant=="Qv"]<-host.pre.new$Qv.Times[host.pre.new$host.plant=="Qv"]
host.pre.new$pre.times[host.pre.new$host.plant=="Qg"]<-host.pre.new$Qg.Times[host.pre.new$host.plant=="Qg"]

head(host.pre.new)
ddply(host.pre.new,"Source",summarize,N=sum(count))
host.pre.new<-host.pre.new[host.pre.new$Source!="Lake Griffin",]
hp.sum<-ddply(host.pre.new,c("HP","Source"),summarize,N=sum(count))
###### write.csv(hp.sum,"hp.sum.csv")
#########################################################
unique(host.pre.new$Source[host.pre.new$species=="B. kinseyi"])
t.test(host.pre.new$preference[host.pre.new$species=="B. kinseyi" &
                                 host.pre.new$Source=="Rice"],mu=0.5,alternative = 'greater')
### mean = 0.6161 ####
t.test(host.pre.new$preference[host.pre.new$species=="B. kinseyi" &
                                 host.pre.new$Source=="Golden Meadow"],mu=0.5,alternative = 'greater')
### mean = 0.6597 ####
t.test(host.pre.new$preference[host.pre.new$species=="B. kinseyi" &
                                 host.pre.new$Source=="Picayune "],mu=0.5,alternative = 'greater')
### mean = 0.6273 ####
t.test(host.pre.new$preference[host.pre.new$species=="B. kinseyi"],mu=0.5,alternative = 'greater')
### mean = 0.6289 ####
hp.pop<-c(0.6161,0.6597,0.6273)
hp.pop.name<-c("Rice","GM","PY")
hp.pop<-data.frame(hp.pop.name,hp.pop)
hp.mean<-0.6289
hp.species<-"B. kinseyi"
hp.mean<-data.frame(hp.species,hp.mean)
hp.pop$species<-"B. kinseyi"
##########################################################################################
hp.sum<-ddply(host.pre.new,c("species","Source"),summarize,N=sum(count),pref=mean(preference),
              pre.SE=sd(preference)/sqrt(N))

ddply(host.pre.new[host.pre.new$species=="B. kinseyi",],c("species"),summarize,N=sum(count),pref=mean(preference),
      pre.SE=sd(preference)/sqrt(N))

ddply(host.pre.new[host.pre.new$species=="B. treatae",],c("species"),summarize,N=sum(count),pref=mean(preference),
      pre.SE=sd(preference)/sqrt(N))

ddply(host.pre.new[host.pre.new$species=="B. fossoria",],c("species"),summarize,N=sum(count),pref=mean(preference),
      pre.SE=sd(preference)/sqrt(N))


hp.sum$species<-factor(hp.sum$species,levels=c("B. treatae","B. fossoria","B. kinseyi"))
hp.sum<-hp.sum[order(hp.sum$species),]
hp.sum$site<-c("CC","KRE","ABS","DCK","LL","GM","PY","RU")
ddply(host.pre.new,c("species"),summarize,N=sum(count),pref=mean(preference),
      pre.SE=sd(preference)/sqrt(N))

ddply(hp.sum,c("species"),summarize,N=sum(N))
##########################################################################################
### host fidelity of allopatric species ####
HP.Bk<-ggplot(data=hp.sum[hp.sum$species=="B. kinseyi",])+
  geom_bar(aes(y=pref,x=site),fill="mediumspringgreen",color="black",
             stat="identity", position=position_dodge2(width=0.1),width=0.6)+
  geom_errorbar(aes(x=site,ymin=pref-pre.SE,ymax=pref+pre.SE),stat="identity", width=.2,
                position=position_dodge2(width=0.1))+
  scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=seq(0,1,0.2))+
  ylab("Host fidelity\n (% time spent on native host)")+xlab("")+
  geom_hline(yintercept = 0.5,linetype="dashed",size=1.5)+
  annotate("text",x=c(1,2,3),y=hp.sum[hp.sum$species=="B. kinseyi","pref"]+hp.sum[hp.sum$species=="B. kinseyi","pre.SE"]+0.04,
           label=paste("n= ",hp.sum[hp.sum$species=="B. kinseyi","N"],sep=""),size=8)+
  annotate("text",x=c(1,2,3),y=c(0.8,0.8,0.78),label=c("***","*","***"),size=18)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),  
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(colour="black"),
        axis.text.x=element_text(size=28),
        axis.title.x = element_text(size=28),
        axis.text.y = element_text(size=24),
        axis.ticks = element_blank(),  
        axis.title.y = element_text(size=26,margin=margin(0,20,0,0)),
        plot.margin=unit(c(1,0,3,0),"lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.title= element_blank(),
        legend.position = c(0.25,0.96),
        legend.text = element_text(size=20),
        plot.title = element_text(size=34,hjust = 0.5))+xlab(expression(paste(italic("B. kinseyi")," (allopatry)")))
HP.Bk

boxes <- data.frame(x=0.28,y=0.1)
HP.Bk1<-ggdraw(HP.Bk) + 
  geom_rect(data = boxes, aes(xmin = x, xmax = x +0.64, ymin = y, ymax = y +0.0018),colour = "black", fill = "black")
HP.Bk1
################################################################
##### compare host preference of sympatry with allopatry #####
hist(host.pre.new[host.pre.new$HP!="Qg","preference"])
hist(host.pre.new[host.pre.new$HP!="Qg","preference1"])

hp1<-glmmTMB(pre.times/Total.times~Method*HP+(1|Source),family="binomial",weights =Total.times,
             data=host.pre.new[host.pre.new$HP!="Qg",])
simulationOutput1<-simulateResiduals(fittedModel=hp1)
testDispersion(simulationOutput1)

hp1<-glmmTMB(pre.times/Total.times~Method*HP+(1|Source),family=betabinomial(link = "logit"),weights =Total.times,
             data=host.pre.new[host.pre.new$HP!="Qg",])
summary(hp1)
hp1<-glmmTMB(pre.times/Total.times~Method+HP+(1|Source),family=betabinomial(link = "logit"),weights =Total.times,
             data=host.pre.new[host.pre.new$HP!="Qg",])
summary(hp1)

hp1.sum<-data.frame(lsmeans(hp1,~HP,type="response"))
hp1.sum$species<-c("B. kinyesi","B. treatae")
hp1.sum$species<-factor(hp1.sum$species,levels=c("B. treatae","B. kinyesi"))
hp1.sum<-hp1.sum[order(hp1.sum$species),]
hp1.anova<-data.frame(Anova(hp1,type="III"))
#### figure of host preference sympatry vs. allopatry ####
hp.sum1<-hp.sum[hp.sum$species!="B. fossoria",]
hp.sum1[6,]<-c("B. treatae","BB",0,-1,0,"BB")
hp.sum1$N<-as.numeric(hp.sum1$N)
hp.sum1$pref<-as.numeric(hp.sum1$pref)
hp.sum1$pre.SE<-as.numeric(hp.sum1$pre.SE)
hp.sum1<-hp.sum1[order(hp.sum1$species),]
HP.geo<-ggplot(data=hp.sum1,aes(x=species,group=site))+
  geom_bar(aes(y=pref,fill=species),width=0.6,color="black",stat="identity", position=position_dodge(width=0.95))+
  geom_errorbar(aes(ymin=pref-pre.SE,ymax=pref+pre.SE),stat="identity", width=.2,
                position=position_dodge(width=0.95))+
  scale_fill_manual(values=c("blue3","mediumspringgreen"))+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.03),breaks=seq(0,1,0.2))+
  ylab("Host fidelity")+xlab("")+
  annotate("text",x=1.5,y=0.96,label="**",size=20)+
  annotate("text",x=c(1,1.34,0.7,1.7,2.02,2.34),y=hp.sum1$pref+hp.sum1$pre.SE+0.03,
           label=paste("n= ",hp.sum1$N,sep=""),size=7)+
  geom_hline(yintercept = 0.5, linetype="dashed",size=0.5)+
  geom_segment(aes(x=0.95,y=0.93,xend=1.35,yend=0.93),size=0.8)+
  geom_segment(aes(x=1.6,y=0.82,xend=2.45,yend=0.82),size=0.8)+
  geom_segment(aes(x=1.15,y=0.93,xend=1.15,yend=0.96),size=0.8)+
  geom_segment(aes(x=2.02,y=0.82,xend=2.02,yend=0.96),size=0.8)+
  geom_segment(aes(x=1.15,y=0.958,xend=2.02,yend=0.958),size=0.8)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),  
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(colour="black"),
        axis.text.x=element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size=24),
        axis.title.y = element_text(size=26,margin=margin(0,20,0,0)),
        plot.margin=unit(c(1,0,3,0),"lines"),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
HP.geo
HP.geo1<-ggdraw(add_sub(HP.geo,hp.sum[hp.sum$species!="B. fossoria","site"], 
                        x=c(0.275,0.42,0.59,0.73,0.88),y=0.5,size=18,color="black"))
HP.geo1
boxes <- data.frame( x =c(0.35,0.62),y =c(0.07,0.07))
HP.geo2<-ggdraw(HP.geo1) + 
  geom_rect(data = boxes, aes(xmin = x, xmax = x +c(0.18,0.29), ymin = y, ymax = y +0.0015),colour = "black", fill = "black")
HP.geo2
HP.geo3<-ggdraw(add_sub(HP.geo2,c(expression(paste(italic("B. treatae"))),
                                                   expression(paste(italic("B. kinseyi")))), 
                        x=c(0.43,0.76),y=1.26,size=22,color="black"))
HP.geo3
HP.geo4<-ggdraw(add_sub(HP.geo3,c("Sympatry","Allopatry"), x=c(0.43,0.76),y=1.4,size=22,color="black"))
HP.geo4
#### compare host preference of Qv with Qg ####
hp2<-glmmTMB(pre.times/Total.times~Method*HP+(1|Source),family=betabinomial(link = "logit"),weights=Total.times,
             data=host.pre.new[host.pre.new$HP!="Qv-T",])
summary(hp2)

hp2<-glmmTMB(pre.times/Total.times~Method+HP+(1|Source),family=betabinomial(link = "logit"),weights=Total.times,
             data=host.pre.new[host.pre.new$HP!="Qv-T",])
summary(hp2)
hp2.sum<-data.frame(lsmeans(hp2,~HP,type="response"))
hp2.sum$species<-c("B. treatae","B. fossoria")
hp2.sum$species<-factor(hp2.sum$species,levels=c("B. treatae","B. fossoria"))
hp2.sum<-hp2.sum[order(hp2.sum$species),]
hp2.anova<-data.frame(Anova(hp2,type="III"))

hp.anova<-list(hp1.anova,hp2.anova)

write.csv(hp.anova,"host.preference.anova.csv")
#### figure of host preference Qv vs. Qg in sympatry ####
hp.sum2<-hp.sum[hp.sum$species!="B. kinseyi",]
hp.sum2[6,]<-c("B. treatae","BB",0,-1,0,"BB")
hp.sum2$N<-as.numeric(hp.sum2$N)
hp.sum2$pref<-as.numeric(hp.sum2$pref)
hp.sum2$pre.SE<-as.numeric(hp.sum2$pre.SE)
hp.sum2<-hp.sum2[order(hp.sum2$species),]
HP.sym<-ggplot(data=hp.sum2,aes(x=species,group=site))+
  geom_bar(aes(y=pref,fill=species),width=0.6,color="black",stat="identity", position=position_dodge(width=0.95))+
  geom_errorbar(aes(ymin=pref-pre.SE,ymax=pref+pre.SE),stat="identity", width=.2,
                position=position_dodge(width=0.95))+
  scale_fill_manual(values=c("blue3","orange2"))+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.03),breaks=seq(0,1,0.2))+
  ylab("Host fidelity")+xlab("")+
  annotate("text",x=1.63,y=0.96,label="*",size=20)+
  annotate("text",x=c(1,1.3,0.7,1.7,2,2.3),y=hp.sum2$pref+hp.sum2$pre.SE+0.03,
           label=paste("n= ",hp.sum2$N,sep=""),size=7)+
  geom_hline(yintercept = 0.5, linetype="dashed",size=0.5)+
  geom_segment(aes(x=0.95,y=0.93,xend=1.35,yend=0.93),size=0.8)+
  geom_segment(aes(x=1.6,y=0.82,xend=2.45,yend=0.82),size=0.8)+
  geom_segment(aes(x=1.15,y=0.93,xend=1.15,yend=0.96),size=0.8)+
  geom_segment(aes(x=2.02,y=0.82,xend=2.02,yend=0.96),size=0.8)+
  geom_segment(aes(x=1.15,y=0.958,xend=2.02,yend=0.958),size=0.8)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),  
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(colour="black"),
        axis.text.x=element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size=24),
        axis.title.y = element_text(size=26,margin=margin(0,20,0,0)),
        plot.margin=unit(c(1,0,3,0),"lines"),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")

HP.sym1<-ggdraw(add_sub(HP.sym,hp.sum[hp.sum$species!="B. kinseyi","site"], 
                        x=c(0.275,0.42,0.59,0.73,0.88),y=0.5,size=18,color="black"))
HP.sym1
boxes <- data.frame( x =c(0.34,0.62),y =c(0.07,0.07))
HP.sym2<-ggdraw(HP.sym1) + 
  geom_rect(data = boxes, aes(xmin = x, xmax = x +c(0.18,0.29), ymin = y, ymax = y +0.0015),colour = "black", fill = "black")
HP.sym2
HP.sym3<-ggdraw(add_sub(HP.sym2,c(expression(paste(italic("B. treatae"))),
                                  expression(paste(italic("B. fossoria")))), 
                        x=c(0.43,0.76),y=1.4,size=22,color="black"))
HP.sym3
HP.sym4<-ggdraw(add_sub(HP.sym3,c("Higher dispersal","Lower dispersal"), x=c(0.44,0.76),y=1.4,size=22,color="black"))
HP.sym4
####################################################################
###### testing host preference from 0.5 ######
t.test(host.pre$preference[host.pre$HP=="Qv-T"& host.pre$Total.times!=0],mu=0.5,alternative = 'greater')
t.test(host.pre$preference[host.pre$HP=="Qv-F"& host.pre$Total.times!=0],mu=0.5,alternative = 'greater')
t.test(host.pre$preference[host.pre$HP=="Qg" & host.pre$Total.times!=0],mu=0.5,alternative = 'greater')
####################################################################
#####################################################################################################
library(ggplot2)
library(multcomp)
library(gridExtra)
library(cowplot)
library(brms)
library(plyr)
library(emmeans)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(patchwork)
#######################################################################################
############################ sexual isolation  ########################################
#######################################################################################
library(stringr)
Mating.T<-read.csv("mate.preference.7.27.csv")
Mating.T<-Mating.T[,-1]
unique(Mating.T$Source..M.)
unique(Mating.T$Source..F.)
levels(Mating.T$Source..F.)
sf.data<-list()
Mating.T$pair5<-Mating.T$pair4
Mating.T$pair5<-as.character(Mating.T$pair5)
Mating.T$Source..F.<-as.character(Mating.T$Source..F.)
Mating.T$Source..M.<-as.character(Mating.T$Source..M.)
Mating.T$pair5[Mating.T$Source..F.!=Mating.T$Source..M. & Mating.T$pair4 =="same"]<-"different.pop"
Mating.T$pair6<-"same"
Mating.T$pair6[Mating.T$Source..M.!=Mating.T$Source..F.]<-"different"
Mating.T$pairs<-paste(Mating.T$Source..M.,Mating.T$Source..F.,sep="_")
#####################################################################################################
### sexual isolation between Qv-T, Qv-F  ###################
##########################################################################
#### 1. including Qv-T, Qv-F ########
Mating.T3<-Mating.T[Mating.T$HP.M.%in% c("Qv-T","Qv-F") & Mating.T$HP.F.%in% c("Qv-T","Qv-F"),]
head(Mating.T3)
summary(factor(Mating.T3$pair2))
Mating.T3$pair4[Mating.T3$pair3=="Qv-T_Qv-F"]<-"different"
ddply(Mating.T3,c("pair4"),summarize,Mt=mean(Mt),count=sum(count))
Mating.T3.sum<-ddply(Mating.T3,c("HP.M.","HP.F."),summarize,Mt=mean(Mt),count=sum(count))

Mating.T3.sum$Mtcount<-Mating.T3.sum$Mt*Mating.T3.sum$count
### analysis model:  pair4 (same/different) ###
Mt2.m1<-glmer(Mt~pair4+(1|Source..F.)+(1|Source..M.),data=Mating.T3,
              family="binomial")
summary(Mt2.m1)
Mt2.ls<-lsmeans(Mt2.m1,~pair4)
Mt2.lsm<-data.frame(summary(Mt2.ls,type="response"))

Mt2.lsm$type<-"different.geo"
Mt2.lsm$N<-c(84,238)
MtQvT.lsm<-rbind(Mt1.lsm,Mt2.lsm)
MtQvT.lsm$pair4<-factor(MtQvT.lsm$pair4,levels=c("same","different"))
MtQvT.lsm<-MtQvT.lsm[order(MtQvT.lsm$type,MtQvT.lsm$pair4),]

1-2*(heterQv.fre/(heterQv.fre+homoQv.fre))

### Qv-T and Qv-F SI with bootstrap ###
Qvt<-c(rep(0,104),rep(1,51))
Qvf<-c(rep(0,38),rep(1,9))
Qvtf<-c(rep(0,61),rep(1,21))
Qvft<-c(rep(0,30),rep(1,7))
home.fre<-(sum(Qvt)/length(Qvt)+sum(Qvf)/length(Qvf))/2
heto.fre<-(sum(Qvtf)/length(Qvtf)+sum(Qvft)/length(Qvft))/2
TFSI.mean<-1-2*heto.fre/(home.fre+heto.fre)
TF.SI<-c()
for (i in 1:10000){
  Qvt.new<-sample(Qvt,replace=TRUE)
  Qvf.new<-sample(Qvf,replace=TRUE)
  Qvtf.new<-sample(Qvtf,replace=TRUE)
  Qvft.new<-sample(Qvft,replace=TRUE)
  home.fre<-(sum(Qvt.new)/length(Qvt.new)+sum(Qvf.new)/length(Qvf.new))/2
  heto.fre<-(sum(Qvtf.new)/length(Qvtf.new)+sum(Qvft.new)/length(Qvft.new))/2
  TF.SI[i]<-1-2*heto.fre/(home.fre+heto.fre)
}

TF.SI
quantile(TF.SI, c(0.025,0.975)) 

### Qv-T and Qg SI with bootstrap ###
#### 1. including Qv-T, Qg ########
Mating.T2<-Mating.T[Mating.T$HP.M.%in% c("Qv-T","Qg") & Mating.T$HP.F.%in% c("Qv-T","Qg"),]
head(Mating.T2)
summary(factor(Mating.T2$pair2))
ddply(Mating.T2,c("pair4"),summarize,Mt=mean(Mt),count=sum(count))
ddply(Mating.T2,c("HP.M.","HP.F."),summarize,Mt=mean(Mt),count=sum(count),Mt.count=Mt*count)
### analysis model:  pair4 (same/different) ###
Mt1.m1<-glmer(Mt~pair4+(1|Source..F.)+(1|Source..M.),data=Mating.T2,
              family="binomial")
summary(Mt1.m1)
Mt1.ls<-lsmeans(Mt1.m1,~pair4)
Mt1.lsm<-data.frame(summary(Mt1.ls,type="response"))
Mt1.lsm$type<-"different.host"
Mt1.lsm$N<-c(427,370)

Qg<-c(rep(0,178),rep(1,38))
Qvt<-c(rep(0,104),rep(1,51))
Qvtg<-c(rep(0,191),rep(1,87))
Qvgt<-c(rep(0,133),rep(1,10))
home.fre<-(sum(Qvt)/length(Qvt)+sum(Qg)/length(Qg))/2
heto.fre<-(sum(Qvtg)/length(Qvtg)+sum(Qvgt)/length(Qvgt))/2
TGSI.mean<-1-2*heto.fre/(home.fre+heto.fre)
TG.SI<-c()
for (i in 1:10000){
  Qg.new<-sample(Qg,replace=TRUE)
  Qvt.new<-sample(Qvt,replace=TRUE)
  Qvtg.new<-sample(Qvtg,replace=TRUE)
  Qvgt.new<-sample(Qvgt,replace=TRUE)
  home.fre<-(sum(Qvt.new)/length(Qvt.new)+sum(Qg.new)/length(Qg.new))/2
  heto.fre<-(sum(Qvtg.new)/length(Qvtg.new)+sum(Qvgt.new)/length(Qvgt.new))/2
  TG.SI[i]<-1-2*heto.fre/(home.fre+heto.fre)
}

TG.SI
c(quantile(TG.SI, c(0.025,0.975)))
t.test(TF.SI,TG.SI)

Qv_T.SI<-data.frame(matrix(nrow=2,ncol=4))
colnames(Qv_T.SI)<-c("pair","mean.SI","lower.SI","upper.SI")
Qv_T.SI[,1]<-c("Qv-T_Qg","Qv-T_Qv-F")
Qv_T.SI[,2]<-c(TGSI.mean,TFSI.mean)
Qv_T.SI[1,3:4]<-c(quantile(TG.SI, c(0.025,0.975)))
Qv_T.SI[2,3:4]<-c(quantile(TF.SI, c(0.025,0.975)))

Qv_T.SInew1<-data.frame(pair=rep("B. kinseyi x B. fossoria",10000),SI=TG.SI)
Qv_T.SInew2<-data.frame(pair=rep("B. kinseyi x B. treatae",10000),SI=TF.SI)
Qv_T.SInew<-rbind(Qv_T.SInew1,Qv_T.SInew2)
##################################################################
QvT_RI_plot<-ggplot(data=Qv_T.SInew,aes(x=pair,y=SI,fill=pair),pattern_fill=
                      "stripe")+
  geom_boxplot()+
  scale_fill_manual(values=c("mediumspringgreen","white"))+
  scale_y_continuous(expand=c(0,0),limits=c(-0.4,0.8),breaks=seq(-0.4,0.8,0.2))+
  geom_hline(yintercept =0,linetype="dashed",size=1.1)+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),  
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=24,colour="black"),
        axis.title=element_text(size=26),
        axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=26,face="italic"),
        plot.margin=unit(c(1,0,3,0),"lines"),
        plot.title = element_text(size = 28, face = "bold",hjust = 0.5))+
  labs(x=" ",y="Strength of sexual isolation")+
  annotate("text",x=1,y=0.43,label="*",size=18)
QvT_RI_plot
QvT_RI_plot1<-ggdraw(add_sub(QvT_RI_plot,c("different host","same host"),x=c(0.27,0.72),y=0.4,size=28))
########## save as 9.5*7.5 ###################
plot_grid(HP.Bk1,QvT_RI_plot1,nrow=1,labels=c("A","B"),label_size = 38) ## save as 20*11 ##
###########################################################################
######## male and female from Qv and Qg, test the asymmetric SI#######
###########################################################################
Mating.T$Block.id<-paste(Mating.T$Year,Mating.T$Block.ID,sep="_")
Mating.T$Year<-as.factor(Mating.T$Year)
Mating.T.sym<-Mating.T[Mating.T$pair3 %in% c("Qg_Qg","Qv-F_Qg","Qv-F_Qv-F"),]
#########################################################################
### male wing buzz ######
### statistical analysis ###
Wp.m1<-glmmTMB(W~Year+HP.M.*pair4+(1|Source..M.)+(1|Source..F.),
               data=Mating.T.sym,
               family="binomial")
summary(Wp.m1)
coef(summary(Wp.m1))$cond[6,4]
fixef(Wp.m1)
simulationOutput <- simulateResiduals(fittedModel = Wp.m1)
plot(simulationOutput)
plotResiduals(simulationOutput,Mating.T.sym$Year )
plotResiduals(simulationOutput,Mating.T.sym$HP.M. )
plotResiduals(simulationOutput,Mating.T.sym$pair4)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)
lsmeans(Wp.m1,pairwise~HP.M.*pair4,adjust="tukey")

Wp.m1.anova<-data.frame(Anova(Wp.m1,type="III"))
### power analysis ####
Mating.T.sym.Qv<-Mating.T.sym[Mating.T.sym$HP.M.!="Qg",]

Qg.same<-Mating.T.sym[Mating.T.sym$HP.M.=="Qg" & Mating.T.sym$pair4=="same",]
Qg.diff<-Mating.T.sym[Mating.T.sym$HP.M.=="Qg" & Mating.T.sym$pair4=="different",]
same.no<-length(Qg.same[,1])
diff.no<-length(Qg.diff[,1])
z.list<-c()
p.list<-c()
p.list2<-c()
for (i in 1:1000){
  Qg.same.temp<-Qg.same[sample(1:same.no,47),]
  Qg.diff.temp<-Qg.diff[sample(1:diff.no,71),]
  Mating.T.sym.temp<-rbind(Mating.T.sym.Qv,Qg.same.temp,Qg.diff.temp)
  Wp.mnew<-glmmTMB(W~Year+HP.M.*pair4+(1|Source..M.)+(1|Source..F.),
                   data=Mating.T.sym.temp,family="binomial")
  p.list2[i]<-data.frame(lsmeans(Wp.mnew,pairwise~HP.M.*pair4,adjust="tukey")$contrast)[5,6]
  p.list[i]<-coef(summary(Wp.mnew))$cond[6,4]
  z.list[i]<-coef(summary(Wp.mnew))$cond[6,3]
}
hist(z.list)
hist(p.list2)
length(p.list[p.list<0.05])
length(p.list[p.list<0.1])
length(p.list2[p.list2<0.05])
length(z.list[z.list>0])
##### Total Wing buzz period m2 ####
Wp.m2<-glmmTMB(T.W~Year+HP.M.*pair4+(1|Source..M.)+(1|Source..F.),
               data=Mating.T.sym,
               family="nbinom2")
summary(Wp.m2)
simulationOutput<- simulateResiduals(fittedModel = Wp.m2)
plot(simulationOutput)
plotResiduals(simulationOutput,Mating.T.sym$Year )
plotResiduals(simulationOutput,Mating.T.sym$HP.M. )
plotResiduals(simulationOutput,Mating.T.sym$pair4)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)
lsmeans(Wp.m2,pairwise~HP.M.*pair4,adjust="tukey")

Wp.m2.anova<-data.frame(Anova(Wp.m2,type="III"))
### Qv wing buzz plot1 ######
### least square means of wing buzz frequency ####
W.lsm<-lsmeans(Wp.m1,~HP.M.*pair4)
W.lsm<-data.frame(summary(W.lsm,type="response"))
W.lsm$HP.M.<-c("B. fossoria","B. treatae","B. fossoria","B. treatae")
W.lsm$HP.M.<-factor(W.lsm$HP.M.,levels=c("B. treatae","B. fossoria"))
W.lsm$pair4<-factor(W.lsm$pair4,levels=c("same","different"))
W.lsm<-W.lsm[order(W.lsm$HP.M.,W.lsm$pair4,decreasing=FALSE),]
### adding number of individuals per group ###
W.count<-ddply(Mating.T.sym,c("HP.M.","pair4"),summarize,count=sum(count))
W.count$HP.M.<-factor(W.count$HP.M.,levels=c("Qv-F","Qg"))
W.count$pair4<-factor(W.count$pair4,levels=c("same","different"))
W.count<-W.count[order(W.count$HP.M.,W.count$pair4,decreasing = FALSE),]
W.lsm$count<-W.count$count
#########################################
W.lsm$pair<-"conspecific"
W.lsm$pair[W.lsm$pair4=="different"]<-"heterospecific"
W.lsm$comp<-paste(W.lsm$HP.M.,W.lsm$pair4,sep="_")
W.lsm$HP.F.<-c("B. treatae","B. fossoria", "B. fossoria", "B. treatae")
W.lsm$HP.F.<-factor(W.lsm$HP.F.,levels=c("B. treatae","B. fossoria"))
W.lsm<-W.lsm[order(W.lsm$HP.M.,W.lsm$HP.F.),]
#################################################
Wing.plot<-ggplot(data=W.lsm,aes(x=HP.M.,group=HP.F.))+
stat_summary(aes(pattern=pair,y=prob,fill=HP.M.,pattern_fill=HP.M.),
             fun = "mean", position=position_dodge(width=0.8),width=0.6,pattern_density=0.5,
             geom = "bar_pattern", colour="black")+
  scale_pattern_manual(values = c(conspecific = "none", heterospecific = "stripe"))+
  scale_fill_manual(values=c("blue3","orange2"))+
  scale_pattern_fill_manual(values=c("orange2","blue3"))+
 geom_errorbar(aes(ymin=prob-SE,ymax=prob+SE),position=position_dodge(width=0.8),width=0.2)+
  scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=seq(0,1,0.2))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),  
    panel.background=element_blank(),
    axis.line=element_line(colour="black"),
    axis.text=element_text(size=22,colour="black"),
    axis.title=element_text(size=26),
    axis.title.y=element_text(margin=margin(0,20,0,0)),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(1,0,3,0),"lines"))+
  labs(x=" ",y="Male mate preference")+
  annotate("text",x=c(0.78,1.23,1.79,2.23),y=W.lsm$prob+W.lsm$SE+0.04,label=paste("n= ",W.lsm$count,sep=""),size=9)+
  geom_segment(aes(x=1.8,y=0.53,xend=2.2,yend=0.53),size=1)+annotate("text",x=2,y=0.56,label="**",size=11)

Wing.plot
### annotate("text",x=1.5,y=0.7,label="Male species x Species pair: ",size=9)+
### annotate("text",x=1.5,y=0.64,label=expression(paste(italic("z"),"= 2.664, ",italic("p")," =0.008")),size=9)+
### geom_segment(aes(x=1.8,y=0.56,xend=2.2,yend=0.56),size=1)+annotate("text",x=2,y=0.58,label="**",size=11)

Wing.plot1<-ggdraw(add_sub(Wing.plot,c("Female:",expression(italic("B. treatae")),expression(italic("B. fossoria")),
            expression(italic("B. treatae")),expression(italic("B. fossoria"))),x=c(0.03,0.19,0.37,0.63,0.81),y=0.4,size=20))
Wing.plot1
Wing.plot2<-ggdraw(add_sub(Wing.plot1,c("Male:",expression(italic("B. treatae")),expression(italic("B. fossoria"))),
x=c(0.17,0.38,0.77),y=1.6,size=20))
Wing.plot2
boxes <- data.frame( x = c(0.24,0.62),y = c(0.1,0.10))
Wing.plot3<-ggdraw(Wing.plot2) + 
 geom_rect(data = boxes, aes(xmin = x, xmax = x +0.28, ymin = y, ymax = y +0.0018),colour = "black", fill = "black")
Wing.plot3 
Wing.plot4<-ggdraw(add_sub(Wing.plot3,c("Higher dispersal","Lower dispersal"),
                           x=c(0.38,0.77),y=1.6,size=20))
Wing.plot4
### save as 10*8 ###
###########################
###########################################################################
#### female copulation frequency ####
F.Mating.sym<-Mating.T.sym[Mating.T.sym$W==1,]
F.Mating.sym$Fa<-F.Mating.sym$Mo+F.Mating.sym$Mt
F.Mating.sym$Fa[F.Mating.sym$Fa!=0]<-1
######## male courtship behavior: including wing buzz and mounting ####
####### female mating behavior: copulation frequency when  the male courtship behavior presents ###
#### stastistical analysis ####
str(F.Mating.sym$Year)
Fsym.m1<-glmmTMB(Fa~Year+HP.F.*pair4+(1|Source..F.)+(1|Source..F.),data=F.Mating.sym,
                 family=binomial)
summary(Fsym.m1)

simulationOutput <- simulateResiduals(fittedModel = Fsym.m1)
plot(simulationOutput)
plotResiduals(simulationOutput,F.Mating.sym$Year )
plotResiduals(simulationOutput,F.Mating.sym$HP.F. )
plotResiduals(simulationOutput,F.Mating.sym$pair4)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

lsmeans(Fsym.m1,pairwise~HP.F.*pair4,adjust="tukey")

Fsym.m1.anova<-data.frame(Anova(Fsym.m1,type="III"))
##### least square means of copulation frequency ####
F.lsm<-lsmeans(Fsym.m1,~HP.F.*pair4)
F.lsm<-data.frame(summary(F.lsm,type="response"))
F.lsm$HP.F.<-factor(F.lsm$HP.F.,levels=c("Qv-F","Qg"))
F.lsm$pair4<-factor(F.lsm$pair4,levels=c("same","different"))
F.lsm<-F.lsm[order(F.lsm$HP.F.,F.lsm$pair4,decreasing = FALSE),]
F.lsm$HP.F.<-c("B. treatae","B. treatae","B. fossoria","B. fossoria")
F.lsm$HP.F.<-factor(F.lsm$HP.F.,levels=c("B. treatae","B. fossoria"))
F.lsm$pair<-c("conspecific","heterospecific","conspecific","heterospecific")
##### add the number of individuals #####
F.count<-ddply(F.Mating.sym,c("HP.F.","pair4"),summarize,count=sum(count))

ddply(Mating.T,c("HP.F.","pair4"),summarize,count=sum(count))

F.count$HP.F.<-factor(F.count$HP.F.,levels=c("Qv-F","Qg"))
F.count$pair4<-factor(F.count$pair4,levels=c("same","different"))
F.count<-F.count[order(F.count$HP.F.,F.count$pair4,decreasing = FALSE),]
F.lsm$count<-F.count$count
F.lsm$HP.M.<-c("B. treatae","B. fossoria","B. fossoria","B. treatae")
F.lsm$HP.M.<-factor(F.lsm$HP.M.,levels=c("B. treatae","B. fossoria"))
F.lsm<-F.lsm[order(F.lsm$HP.F.,F.lsm$HP.M.),]

F.plot<-ggplot(data=F.lsm,aes(y=prob,x=HP.F.,group=HP.M.))+
  stat_summary(aes(pattern=pair,y=prob,fill=HP.F.,pattern_fill=HP.F.),
               fun = "mean", position=position_dodge(width=0.8),width=0.6,pattern_density=0.5,
               geom = "bar_pattern", colour="black")+
  scale_pattern_manual(values = c(conspecific = "none", heterospecific = "stripe"))+
  scale_fill_manual(values=c("blue3","orange2"))+
  scale_pattern_fill_manual(values=c("orange2","blue3"))+
  geom_errorbar(aes(ymin=prob-SE,ymax=prob+SE),position=position_dodge(width=0.8),width=0.2)+
  scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=seq(0,1,0.2))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),  
    panel.background=element_blank(),
    axis.line=element_line(colour="black"),
    axis.text=element_text(size=22,colour="black"),
    axis.title=element_text(size=26),
    axis.title.y=element_text(margin=margin(0,20,0,0)),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(1,0,3,0),"lines"))+
  labs(x=" ",y="Female mate preference")+
  annotate("text",x=c(0.79,1.24,1.78,2.24),y=F.lsm$prob+F.lsm$SE+0.04,label=paste("n= ",F.lsm$count,sep=""),size=9)
F.plot

F.plot1<-ggdraw(add_sub(F.plot,c("Male:",expression(italic("B. treatae")),expression(italic("B. fossoria")),
                                       expression(italic("B. treatae")),expression(italic("B. fossoria"))),x=c(0.03,0.19,0.37,0.63,0.81),y=0.4,size=20))
F.plot1
F.plot2<-ggdraw(add_sub(F.plot1,c("Female:",expression(italic("B. treatae")),expression(italic("B. fossoria"))),
                           x=c(0.16,0.4,0.76),y=1.6,size=20))
F.plot2
boxes <- data.frame( x = c(0.24,0.62),y = c(0.1,0.1))
F.plot3<-ggdraw(F.plot2) + 
  geom_rect(data = boxes, aes(xmin = x, xmax = x +0.28, ymin = y, ymax = y +0.0018),colour = "black", fill = "black")
F.plot3 
F.plot4<-ggdraw(add_sub(F.plot3,c("Higher dispersal","Lower dispersal"),
                           x=c(0.42,0.77),y=1.6,size=20))
F.plot4

###############################
### combine male and female plot ####
plot_grid(Wing.plot4,F.plot4,ncol=2,labels=c("A","B"),label_size = 38)
#######################################################################
########## test of reinforcement #######################
Mating.T$Block.id<-paste(Mating.T$Year,Mating.T$Block.ID,sep="_")
Mating.T$Year<-as.factor(Mating.T$Year)
#########################################################################
### reinforcement for male mate preference
Mating.Tnew<-Mating.T[Mating.T$pair2 %in% c("Qv-T_Qv-T","Qg_Qg",
                                            "Qv-T_Qg","Qg_Qv-T"),]

Wp.Qv.1<-glmmTMB(W~Year+HP.M.*pair4+(1|Source..M.),
                 data=Mating.Tnew,
                 family="binomial")
summary(Wp.Qv.1)
simulationOutput <- simulateResiduals(fittedModel = Wp.Qv.1)
plot(simulationOutput)
lsmeans(Wp.Qv.1,pairwise~HP.M.*pair4,adjust="tukey")
Wp.Qv.1.anova<-data.frame(Anova(Wp.Qv.1,type="III"))
WQv<-lsmeans(Wp.Qv.1,~HP.M.*pair4)
WQv<-data.frame(summary(WQv,type="response"))
WQv$pair4<-factor(WQv$pair4,levels=c("same","different"))
WQv<-WQv[order(WQv$HP.M.,WQv$pair4,decreasing=FALSE),]
WQv$species<-c("B. fossoria","B. fossoria","B. kinyesi","B. kinyesi")
WQv$species<-factor(WQv$species,levels=c("B. fossoria","B. kinyesi"))
WQv$HP.F.<-c("B. fossoria","B. kinyesi","B. kinyesi","B. fossoria")
WQv$HP.F.<-factor(WQv$HP.F.,levels=c("B. fossoria","B. kinyesi"))


### adding number of individuals per group ###
WQv.count<-ddply(Mating.Tnew,c("HP.M.","pair4"),summarize,count=sum(count))
WQv.count$HP.M.<-factor(WQv.count$HP.M.,levels=c("Qv-F","Qg","Qv-T"))
WQv.count$pair4<-factor(WQv.count$pair4,levels=c("same","different"))
WQv.count<-WQv.count[order(WQv.count$HP.M.,WQv.count$pair4,decreasing = FALSE),]
WQv$count<-WQv.count$count

WQv<-WQv[order(WQv$species,WQv$HP.F.,decreasing=FALSE),]
WQv
#########################################
Wing.Qv.plot<-ggplot(data=WQv,aes(y=prob,x=HP.M.,group=HP.F.))+
  stat_summary(aes(pattern=pair4,y=prob,fill=HP.M.,pattern_fill=HP.M.),
               fun = "mean", position=position_dodge(width=0.8),width=0.6,pattern_density=0.5,
               geom = "bar_pattern", colour="black")+
  scale_pattern_manual(values = c(same = "none", different = "stripe"))+
  scale_pattern_fill_manual(values=c("mediumspringgreen","orange2"))+
  geom_errorbar(aes(ymin=prob-SE,ymax=prob+SE),position=position_dodge(width=0.8),width=0.2)+
  scale_fill_manual(values=c("orange2","mediumspringgreen"))+
  scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=seq(0,1,0.2))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),  
    panel.background=element_blank(),
    axis.line=element_line(colour="black"),
    axis.text=element_text(size=22,colour="black"),
    axis.title=element_text(size=26),
    axis.title.y=element_text(margin=margin(0,20,0,0)),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(1,0,3,0),"lines"))+
  labs(x=" ",y="Male mate preference")+
  annotate("text",x=c(0.76,1.24,1.76,2.24),y=WQv$prob+WQv$SE+0.06,label=paste("n= ",WQv$count,sep=""),size=7)+
  geom_segment(aes(x=0.8,y=0.56,xend=1.2,yend=0.56),size=1)+annotate("text",x=1,y=0.58,label="***",size=12)
Wing.Qv.plot

Wing.Qv.plot1<-ggdraw(add_sub(Wing.Qv.plot,c("Female:",expression(italic("B. fossoria")),expression(italic("B. kinseyi")),
                              expression(italic("B. fossoria")),expression(italic("B. kinseyi"))),x=c(0.03,0.19,0.37,0.63,0.81),y=0.4,size=20))
Wing.Qv.plot1
Wing.Qv.plot2<-ggdraw(add_sub(Wing.Qv.plot1,c("Male:",expression(paste(italic("B. fossoria")," (sympatry)")),expression(paste(italic("B. kinseyi")," (allopatry)"))),
                        x=c(0.14,0.37,0.76),y=1.6,size=20))
Wing.Qv.plot2
boxes <- data.frame( x = c(0.22,0.61),y = c(0.12,0.12))
Wing.Qv.plot3<-ggdraw(Wing.Qv.plot2) + 
  geom_rect(data = boxes, aes(xmin = x, xmax = x +0.29, ymin = y, ymax = y +0.0018),colour = "black", fill = "black")
Wing.Qv.plot3


### save image 12.89*8.7 #######
###########################################################################
