setwd("H:\\1New_paper\\12GBD_thyroid_cancer_smoking")
library(readxl)
library(ggplot2)
library(tidyverse)
library(scatterplot3d)
library(ggpubr)
library(RColorBrewer)
library(popPyramid)


  #####绘制年龄段的4个指标的金字塔图######
  ####DALY####

data <- read_excel("0427F1DALYs.xlsx")
# data <- read_excel("DALY_pyramid.xlsx")
  
data$age <- factor(data$age,ordered=TRUE,levels=c('5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79','80-84','85-89','90-94','95+'))
  
p1 <- plotPyramid(df=data, age="age", sex="sex", pop="val",twocolors = c('#82E7E7', '#E66060'), labx="Age-standarzied DALYs, per 100,000",laby= "Age, years", value.labels=FALSE) + 
      theme(legend.position=c(0.95,0.05), legend.justification=c(1,0))+
      theme(panel.grid=element_blank(),
      panel.background=element_blank(),
      axis.line = element_line(color="black"))
p1 <- p1+ theme(axis.text.x = element_text(face = "bold",
                                           size = 9, angle = 360),
                axis.text.y = element_text(face = "bold",
                                           size = 8))
####Incidence####
data <- read_excel("0427F1incidence.xlsx")

data$age <- factor(data$age,ordered=TRUE,levels=c('5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79','80-84','85-89','90-94','95+'))

p2 <- plotPyramid(df=data, age="age", sex="sex", pop="val",twocolors = c('#82E7E7', '#E66060'), labx="ASIR, per 100,000",laby= "Age, years", value.labels=FALSE) + 
  theme(legend.position=c(0.95,0.05), legend.justification=c(1,0))+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line = element_line(color="black"))
p2 <- p2+ theme(axis.text.x = element_text(face = "bold",
                                           size = 9, angle = 360),
                axis.text.y = element_text(face = "bold",
                                           size = 8))

####Deths####
data <- read_excel("0427F1deaths.xlsx")

data$age <- factor(data$age,ordered=TRUE,levels=c('5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79','80-84','85-89','90-94','95+'))

p3 <- plotPyramid(df=data, age="age", sex="sex", pop="val",twocolors = c('#82E7E7', '#E66060'), labx="ASDR, per 100,000",laby= "Age, years", value.labels=FALSE) + 
  theme(legend.position=c(0.95,0.05), legend.justification=c(1,0))+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line = element_line(color="black"))
p3 <- p3+ theme(axis.text.x = element_text(face = "bold",
                                           size = 9, angle = 360),
                axis.text.y = element_text(face = "bold",
                                           size = 8))

####Prevalence####
data <- read_excel("0427F1prevalence.xlsx")

data$age <- factor(data$age,ordered=TRUE,levels=c('5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79','80-84','85-89','90-94','95+'))

p4 <- plotPyramid(df=data, age="age", sex="sex", pop="val",twocolors = c('#82E7E7', '#E66060'), labx="ASPR, per 100,000",laby= "Age, years", value.labels=FALSE) + 
  theme(legend.position=c(0.95,0.05), legend.justification=c(1,0))+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line = element_line(color="black"))
p4<- p4+ theme(axis.text.x = element_text(face = "bold",
                                           size = 9, angle = 360),
                axis.text.y = element_text(face = "bold",
                                           size = 8))
p4

#######Figure1绘制结束#######
#拼图
library(ggpubr)
ggarrange(p3,p1,p4,p2,ncol=2,nrow=2,common.legend = T,
          legend= "bottom",
          labels = c("A","B","C","D"))


#######Figure2双轴图形绘制#######


####F2-DALY
data <- read_excel("f2_DALY.xlsx")

p <- ggplot(data, aes(x = year, y=Number, fill=sex))+
    scale_fill_manual(values=c('#E66060','#82E7E7'))

p <- p + geom_errorbar(aes(ymin = nlower, ymax = nupper),
                width = .3, position = position_dodge(0.7))+
                geom_bar(stat = "identity",width = 0.8,
                position = position_dodge()) 

p <- p + 
  geom_path(aes(x=year,y=ASR*10000,group=sex,color=sex),size=1)+
  geom_ribbon(aes(ymin=alower*10000,ymax=aupper*10000,group=sex),alpha=0.5)
p <- p + scale_y_continuous(sec.axis = sec_axis(~./10000 , name = "Age-standardized DALYs, per 100,000"))  #geom_line()+
p<-p+ theme(panel.grid=element_blank(),
            panel.background=element_blank(),
            axis.line = element_line(color="black"))+theme(legend.position=c('top')) +
            scale_x_continuous(breaks=c(1990, 1995, 2000, 2005, 2010,  2015, 2019))+
            labs(x='Year',y="Number of DALYs")
p1 <- p

####F2-deaths
data <- read_excel("f2_deaths.xlsx")

p <- ggplot(data, aes(x = year, y=Number, fill=sex))+
  scale_fill_manual(values=c('#E66060','#82E7E7'))

p <- p + geom_errorbar(aes(ymin = nlower, ymax = nupper),
                       width = .3, position = position_dodge(0.7))+
  geom_bar(stat = "identity",width = 0.8,
           position = position_dodge()) 

p <- p + 
  geom_path(aes(x=year,y=ASR*10000,group=sex,color=sex),size=1)+
  geom_ribbon(aes(ymin=alower*10000,ymax=aupper*10000,group=sex),alpha=0.5)
p <- p + scale_y_continuous(sec.axis = sec_axis(~./10000 , name = "ASDR, per 100,000"))  #geom_line()+
p<-p+ theme(panel.grid=element_blank(),
            panel.background=element_blank(),
            axis.line = element_line(color="black"))+theme(legend.position=c('top')) +
  scale_x_continuous(breaks=c(1990, 1995, 2000, 2005, 2010,  2015, 2019))+
  labs(x='Year',y="Number of deaths")
p2 <- p

####f2-prevalence
data <- read_excel("f2_prevalence.xlsx")

p <- ggplot(data, aes(x = year, y=Number, fill=sex))+
  scale_fill_manual(values=c('#E66060','#82E7E7'))

p <- p + geom_errorbar(aes(ymin = nlower, ymax = nupper),
                       width = .3, position = position_dodge(0.7))+
  geom_bar(stat = "identity",width = 0.8,
           position = position_dodge()) 

p <- p + 
  geom_path(aes(x=year,y=ASR*10000,group=sex,color=sex),size=1)+
  geom_ribbon(aes(ymin=alower*10000,ymax=aupper*10000,group=sex),alpha=0.5)
p <- p + scale_y_continuous(sec.axis = sec_axis(~./10000 , name = "ASPR, per 100,000"))  #geom_line()+
p<-p+ theme(panel.grid=element_blank(),
            panel.background=element_blank(),
            axis.line = element_line(color="black"))+theme(legend.position=c('top')) +
  scale_x_continuous(breaks=c(1990, 1995, 2000, 2005, 2010,  2015, 2019))+
  labs(x='Year',y="Number of prevalence")
p3 <- p

####f2-incidence

data <- read_excel("f2_incidence.xlsx")

p <- ggplot(data, aes(x = year, y=Number, fill=sex))+
  scale_fill_manual(values=c('#E66060','#82E7E7'))

p <- p + geom_errorbar(aes(ymin = nlower, ymax = nupper),
                       width = .3, position = position_dodge(0.7))+
  geom_bar(stat = "identity",width = 0.8,
           position = position_dodge()) 

p <- p + 
  geom_path(aes(x=year,y=ASR*10000,group=sex,color=sex),size=1)+
  geom_ribbon(aes(ymin=alower*10000,ymax=aupper*10000,group=sex),alpha=0.5)
p <- p + scale_y_continuous(sec.axis = sec_axis(~./10000 , name = "ASIR, per 100,000"))  #geom_line()+
p<-p+ theme(panel.grid=element_blank(),
            panel.background=element_blank(),
            axis.line = element_line(color="black"))+theme(legend.position=c('top')) +
  scale_x_continuous(breaks=c(1990, 1995, 2000, 2005, 2010,  2015, 2019))+
  labs(x='Year',y="Number of incidence")
p4 <- p

###
ggarrange(p2,p1,p3,p4,ncol=2,nrow=2,common.legend = T,
          legend= "bottom",
          labels = c("A","B","C","D"))
############尝试气泡图展示年龄段和年份对rate的关系#######


data <- read_excel("BMIandTC1.xlsx")
data$measure <- ifelse(data$measure=="DALYs (Disability-Adjusted Life)","DALYs","ASDR")

###ASDR
data1 <- subset(data,data$measure=="ASDR")
data1 <- subset(data1,data1$sex!="Both")
data1$ASDR <- data1$val

par(mar=c(5,3,3,2))
ggballoonplot(data1,x="year",y="age",
              size="ASDR",fill="ASDR",
              facet.by = "sex")+
              scale_fill_gradient(low= '#82E7E7', high='#E66060')


###DALYs
data1 <- subset(data,data$measure=="DALYs")
data1 <- subset(data1,data1$sex!="Both")
data1$DALYs <- data1$val

ggballoonplot(data1,x="year",y="age",
              size="DALYs",fill="DALYs",
              facet.by = "sex")+
  scale_fill_gradient(low= '#82E7E7', high='#E66060')

###Both
data1 <- subset(data,data$sex=="Both")
ggballoonplot(data1,x="year",y="age",
              size="val",fill="val",
              facet.by = "measure")+
  scale_fill_gradient(low= '#82E7E7', high='#E66060')


###############%%%%%%%%%%%%%%%%%%%%%%%%%%##################
data <- read_excel("BMIandTC2.xlsx")
data$measure <- ifelse(data$measure=="YLDs (Years Lived with Disability)","YLDs","YLLs")

###YLDs
data1 <- subset(data,data$measure=="YLDs")
data1 <- subset(data1,data1$sex!="Both")
data1$YLDs <- data1$val

par(mar=c(5,3,3,2))
ggballoonplot(data1,x="year",y="age",
              size="YLDs",fill="YLDs",
              facet.by = "sex")+
  scale_fill_gradient(low= '#82E7E7', high='#E66060')

###YLLs
data1 <- subset(data,data$measure=="YLLs")
data1 <- subset(data1,data1$sex!="Both")
data1$YLLs <- data1$val

ggballoonplot(data1,x="year",y="age",
              size="YLLs",fill="YLLs",
              facet.by = "sex")+
  scale_fill_gradient(low= '#82E7E7', high='#E66060')

###Both
data1 <- subset(data,data$sex=="Both")
ggballoonplot(data1,x="year",y="age",
              size="val",fill="val",
              facet.by = "measure")+
  scale_fill_gradient(low= '#82E7E7', high='#E66060')



data <- GBDread('20230413EAPC.zip')

data <- subset(data,data$age=="Age-standardized")
data <- subset(data,data$metric=="Rate")

result <- GBDeapc(data, rei=F,EAPC_95CI = T,digits = 2,sep=", ")
write.csv(result,"EAPC_easyGBD.csv")

data <- data[,c("measure","sex","year","val")]


####以下代码功能由GBDeapc函数代替####

# ##############
# data1 <- subset(data,data$measure=="DALYs (Disability-Adjusted Life Years)"&data$sex=="Both")
# 
# data1$lograte<-log(data1$val)
# fit <- lm(lograte ~ year, data1)
# b <- coefficients(fit)
# CI <- confint(fit)
# EAPC0.5 <- format(round(100 * (exp(b[2]) - 1), digits = 2), nsmall = 2)
# EAPC0.975 <- format(round(100 * (exp(CI[4]) - 1), digits = 2), nsmall = 2)
# EAPC0.025 <- format(round(100 * (exp(CI[2]) - 1), digits = 2), nsmall = 2)
# paste(EAPC0.5,EAPC0.025,EAPC0.975)
# 
# data1 <- subset(data,data$measure=="DALYs (Disability-Adjusted Life Years)"&data$sex=="Male")
# 
# data1$lograte<-log(data1$val)
# fit <- lm(lograte ~ year, data1)
# b <- coefficients(fit)
# CI <- confint(fit)
# EAPC0.5 <- format(round(100 * (exp(b[2]) - 1), digits = 2), nsmall = 2)
# EAPC0.975 <- format(round(100 * (exp(CI[4]) - 1), digits = 2), nsmall = 2)
# EAPC0.025 <- format(round(100 * (exp(CI[2]) - 1), digits = 2), nsmall = 2)
# paste(EAPC0.5,EAPC0.025,EAPC0.975)
# 
# data1 <- subset(data,data$measure=="DALYs (Disability-Adjusted Life Years)"&data$sex=="Female")
# 
# data1$lograte<-log(data1$val)
# fit <- lm(lograte ~ year, data1)
# b <- coefficients(fit)
# CI <- confint(fit)
# EAPC0.5 <- format(round(100 * (exp(b[2]) - 1), digits = 2), nsmall = 2)
# EAPC0.975 <- format(round(100 * (exp(CI[4]) - 1), digits = 2), nsmall = 2)
# EAPC0.025 <- format(round(100 * (exp(CI[2]) - 1), digits = 2), nsmall = 2)
# paste(EAPC0.5,EAPC0.025,EAPC0.975)
# 
# ####################
# data1 <- subset(data,data$measure=="Deaths"&data$sex=="Both")
# 
# data1$lograte<-log(data1$val)
# fit <- lm(lograte ~ year, data1)
# b <- coefficients(fit)
# CI <- confint(fit)
# EAPC0.5 <- format(round(100 * (exp(b[2]) - 1), digits = 2), nsmall = 2)
# EAPC0.975 <- format(round(100 * (exp(CI[4]) - 1), digits = 2), nsmall = 2)
# EAPC0.025 <- format(round(100 * (exp(CI[2]) - 1), digits = 2), nsmall = 2)
# paste(EAPC0.5,EAPC0.025,EAPC0.975)
# 
# data1 <- subset(data,data$measure=="Deaths"&data$sex=="Male")
# 
# data1$lograte<-log(data1$val)
# fit <- lm(lograte ~ year, data1)
# b <- coefficients(fit)
# CI <- confint(fit)
# EAPC0.5 <- format(round(100 * (exp(b[2]) - 1), digits = 2), nsmall = 2)
# EAPC0.975 <- format(round(100 * (exp(CI[4]) - 1), digits = 2), nsmall = 2)
# EAPC0.025 <- format(round(100 * (exp(CI[2]) - 1), digits = 2), nsmall = 2)
# paste(EAPC0.5,EAPC0.025,EAPC0.975)
# 
# data1 <- subset(data,data$measure=="Deaths"&data$sex=="Female")
# 
# data1$lograte<-log(data1$val)
# fit <- lm(lograte ~ year, data1)
# b <- coefficients(fit)
# CI <- confint(fit)
# EAPC0.5 <- format(round(100 * (exp(b[2]) - 1), digits = 2), nsmall = 2)
# EAPC0.975 <- format(round(100 * (exp(CI[4]) - 1), digits = 2), nsmall = 2)
# EAPC0.025 <- format(round(100 * (exp(CI[2]) - 1), digits = 2), nsmall = 2)
# paste(EAPC0.5,EAPC0.025,EAPC0.975)
# 
# #######################
# data1 <- subset(data,data$measure=="Prevalence"&data$sex=="Both")
# 
# data1$lograte<-log(data1$val)
# fit <- lm(lograte ~ year, data1)
# b <- coefficients(fit)
# CI <- confint(fit)
# EAPC0.5 <- format(round(100 * (exp(b[2]) - 1), digits = 2), nsmall = 2)
# EAPC0.975 <- format(round(100 * (exp(CI[4]) - 1), digits = 2), nsmall = 2)
# EAPC0.025 <- format(round(100 * (exp(CI[2]) - 1), digits = 2), nsmall = 2)
# paste(EAPC0.5,EAPC0.025,EAPC0.975)
# 
# data1 <- subset(data,data$measure=="Prevalence"&data$sex=="Male")
# 
# data1$lograte<-log(data1$val)
# fit <- lm(lograte ~ year, data1)
# b <- coefficients(fit)
# CI <- confint(fit)
# EAPC0.5 <- format(round(100 * (exp(b[2]) - 1), digits = 2), nsmall = 2)
# EAPC0.975 <- format(round(100 * (exp(CI[4]) - 1), digits = 2), nsmall = 2)
# EAPC0.025 <- format(round(100 * (exp(CI[2]) - 1), digits = 2), nsmall = 2)
# paste(EAPC0.5,EAPC0.025,EAPC0.975)
# 
# data1 <- subset(data,data$measure=="Prevalence"&data$sex=="Female")
# 
# data1$lograte<-log(data1$val)
# fit <- lm(lograte ~ year, data1)
# b <- coefficients(fit)
# CI <- confint(fit)
# EAPC0.5 <- format(round(100 * (exp(b[2]) - 1), digits = 2), nsmall = 2)
# EAPC0.975 <- format(round(100 * (exp(CI[4]) - 1), digits = 2), nsmall = 2)
# EAPC0.025 <- format(round(100 * (exp(CI[2]) - 1), digits = 2), nsmall = 2)
# paste(EAPC0.5,EAPC0.025,EAPC0.975)
# 
# 
# #######################
# data1 <- subset(data,data$measure=="Incidence"&data$sex=="Both")
# 
# data1$lograte<-log(data1$val)
# fit <- lm(lograte ~ year, data1)
# b <- coefficients(fit)
# CI <- confint(fit)
# EAPC0.5 <- format(round(100 * (exp(b[2]) - 1), digits = 2), nsmall = 2)
# EAPC0.975 <- format(round(100 * (exp(CI[4]) - 1), digits = 2), nsmall = 2)
# EAPC0.025 <- format(round(100 * (exp(CI[2]) - 1), digits = 2), nsmall = 2)
# paste(EAPC0.5,EAPC0.025,EAPC0.975)
# 
# data1 <- subset(data,data$measure=="Incidence"&data$sex=="Male")
# 
# data1$lograte<-log(data1$val)
# fit <- lm(lograte ~ year, data1)
# b <- coefficients(fit)
# CI <- confint(fit)
# EAPC0.5 <- format(round(100 * (exp(b[2]) - 1), digits = 2), nsmall = 2)
# EAPC0.975 <- format(round(100 * (exp(CI[4]) - 1), digits = 2), nsmall = 2)
# EAPC0.025 <- format(round(100 * (exp(CI[2]) - 1), digits = 2), nsmall = 2)
# paste(EAPC0.5,EAPC0.025,EAPC0.975)
# 
# data1 <- subset(data,data$measure=="Incidence"&data$sex=="Female")
# 
# data1$lograte<-log(data1$val)
# fit <- lm(lograte ~ year, data1)
# b <- coefficients(fit)
# CI <- confint(fit)
# EAPC0.5 <- format(round(100 * (exp(b[2]) - 1), digits = 2), nsmall = 2)
# EAPC0.975 <- format(round(100 * (exp(CI[4]) - 1), digits = 2), nsmall = 2)
# EAPC0.025 <- format(round(100 * (exp(CI[2]) - 1), digits = 2), nsmall = 2)
# paste(EAPC0.5,EAPC0.025,EAPC0.975)



































