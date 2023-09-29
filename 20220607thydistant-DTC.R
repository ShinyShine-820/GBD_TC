
rm(list = ls())

library(car)
library(ggplot2)
library(glmnet)
library(MASS)
library(survival)
library(rms)
library(survminer)
library(ggridges)
library(pROC)
library(plotROC)
library(riskRegression)
library(magrittr)
library(DynNom)
library(packrat)
library(rsconnect)
library(readxl)
library(plyr)
library(survivalROC)

setwd("H:\\Frontiers in Genetics\\0 20220606\\修改稿")
thy <- read.csv("20220607_raw_data_10.csv")
thy <- thy[,1:11]
str(thy)

thy$Age1 <- thy$Age
thy$Tumor_size1 <- thy$Tumor_size

thy$Age <- ifelse(thy$Age>5&thy$Age<=50,"Age_Low",ifelse(thy$Age>50&thy$Age<77,"Age_Middle","Age_High"))
thy$Tumor_size <- ifelse(thy$Tumor_size<28,"Size_Small",ifelse(thy$Tumor_size>27&thy$Tumor_size<66,"Size_Middle","Size_Big"))

fit <- survfit(Surv(Survivalmonths,status) ~ Chemotherapy, data = thy)
ggsurvplot(fit, data = thy,
           pval = TRUE,
           pval.coord = c(0, 0.03), 
           surv.median.line = "hv", 
           legend.title="Type",
           xlim=c(1,180),
           palette=c("#005824", "#E41A1C","#377EB8","#984EA3"),
           # title="Cancer-Specific Survival",
           ylab="Cumulative survival (percentage)",xlab = " Survival time (Months)", #???ĺ???????
           censor.shape = 124,censor.size = 2,conf.int = FALSE,
           break.x.by = 30,
           risk.table = TRUE, # 添加风险表
           risk.table.col = "strata", # 根据分层更改风险表颜色
           # xlab = "PDTC CSS(m)", # 指定x轴标签
           ggtheme = theme_bw())


aa <- thy
#1.批量单因素cox的y
y<- Surv(time = aa$Survivalmonths,event = aa$status==1)
#2.循环函数
Uni_cox_model<- function(x){
        FML <- as.formula(paste0 ("y~",x))
        cox<- coxph(FML,data=aa)
        cox1<-summary(cox)
        HR <- round(cox1$coefficients[,2],2)
        PValue <- round(cox1$coefficients[,5],3)
        CI5 <-round(cox1$conf.int[,3],2)
        CI95 <-round(cox1$conf.int[,4],2)
        Uni_cox_model<- data.frame('Characteristics' = x,
                                   'HR' = HR,
                                   'CI5' = CI5,
                                   'CI95' = CI95,
                                   'Uni_P' = PValue)
        return(Uni_cox_model)}  
#3.需要挑选需要进行单因素分析的变量。
variable.names<- colnames(aa)[c(1:9)];variable.names
#4.进行循环并调整最终结果
Uni_cox <- lapply(variable.names, Uni_cox_model)
Uni_cox<- ldply(Uni_cox,data.frame)
Uni_cox$HR.CI95<-paste0(Uni_cox$HR," (",Uni_cox$CI5,'-',Uni_cox$CI95,")")
Uni_cox<-Uni_cox[,-2:-4] 
View(Uni_cox)
write.csv(Uni_cox,"Uni_cox.csv")

Outcome <- "Surv(Survivalmonths,status)"
CandidateVariables <- c("Age","Sex","Race","Tumor_size" ,"Histologic",
                        "RN_positive","Surgery","Radiation",
                       "Chemotherapy") 
#,"Grade","Year","Race","Sex",
Formula <- formula(paste(paste(Outcome,"~", collapse=" "), 
                         paste(CandidateVariables, collapse=" + ")))
model.full <- coxph(Formula,thy,x=TRUE)
summary(model.full)
model.step <- stepAIC(model.full, direction="both")
summary(model.step)

#除了Race外其他因素都是预后因素
#################################################################################
# View(thy)
thy00 <- subset(thy,thy$Survivalmonths<120&thy$status==0)
thy01 <- subset(thy,thy$Survivalmonths>=120&thy$status==0)
thy02 <- subset(thy,thy$status==1)
thy1 <- rbind(thy01,thy02)

str(thy1)

set.seed(123456)
prop_train <- 0.7   # proportion of training data
train <- sample(nrow(thy1), nrow(thy1) * prop_train)
train_data <- thy1[train, ]
test_data <- thy1[-train, ]
nrow(train_data)
nrow(test_data)
# write.csv(train_data,"train_data.csv")
# write.csv(test_data,"test_data.csv")
################################################################################
#生成哑变量

dataset2 = train_data
library(plyr)
into_factor = function(x){
  
  if(class(x) == "character"){
    n = length(x)
    data.fac = data.frame(x = x,y = 1:n)
    output = model.matrix(y~x,data.fac)[,-1]
    ## Convert factor into dummy variable matrix
  }else{
    output = x
    ## if x is numeric, output is x
  }
  output
}
#into_factor(dataset$V4)[1:5,]
dataset2 = colwise(into_factor)(dataset2)
dataset2 = do.call(cbind,dataset2)
dataset2 = as.data.frame(dataset2);head(dataset2)

dataset3 = test_data
library(plyr)
into_factor = function(x){
  
  if(class(x) == "character"){
    n = length(x)
    data.fac = data.frame(x = x,y = 1:n)
    output = model.matrix(y~x,data.fac)[,-1]
    ## Convert factor into dummy variable matrix
  }else{
    output = x
    ## if x is numeric, output is x
  }
  output
}
#into_factor(dataset$V4)[1:5,]
dataset3 = colwise(into_factor)(dataset3)
dataset3 = do.call(cbind,dataset3)
dataset3 = as.data.frame(dataset3);head(dataset3)

train_data <- dataset2
test_data <- dataset3
################################################################################

str(thy11)
###############################################################  ROC
#XGBoost 模型

library(xgboost)
library(Matrix)
#########################################################################################################
#20220225调整代码
traindata1 <- data.matrix(train_data[,1:14]) # 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
traindata2 <- Matrix(traindata1,sparse=T) 
traindata3 <- data.matrix(train_data[,16])# 将自变量和因变量拼接为list
traindata4 <- list(data=traindata2,label=traindata3) # 构造模型需要的xgb.DMatrix对象，处理对象为稀疏矩阵
dtrain <- xgb.DMatrix(data = traindata4$data, label = traindata4$label) 

testset1 <- data.matrix(test_data[,1:14]) # 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
testset2 <- Matrix(testset1,sparse=T) # 将因变量转化为numeric
testset3 <- data.matrix(test_data[,16])# 将自变量和因变量拼接为list
testset4 <- list(data=testset2,label=testset3) # 构造模型需要的xgb.DMatrix对象，处理对象为稀疏矩阵
dtest <- xgb.DMatrix(data = testset4$data, label = testset4$label) 

xgb <- xgboost(data = dtrain,max_depth=100, eta=0.6,  
               objective='binary:logistic', nround=100)

XGBoost_pred. = predict(xgb,newdata = dtrain)
XGBoost_ROC <- roc(train_data$status,as.numeric(XGBoost_pred.))
plot(XGBoost_ROC, type="l",col="#E41A1C", ##线条设置
     xlim=c(1,0), ylim=c(0,1),
     # xlab="1-Specificity" , 
     ylab="Sensitivity")

XGBoost_pred. = predict(xgb,newdata = dtest)
XGBoost_ROC <- roc(test_data$status,as.numeric(XGBoost_pred.))
plot(XGBoost_ROC, type="l",col="#E41A1C", ##线条设置
     xlim=c(1,0), ylim=c(0,1),
     # xlab="1-Specificity" , 
     ylab="Sensitivity")
XGBoost_ROC

importance_matrix <- xgb.importance(colnames(dtrain), model = xgb)
xgb.plot.importance(importance_matrix, rel_to_first = TRUE, 
                    xlab = "Relative importance",
                    col = "#3182BD",
                    main = "XGBoost Feature Importance")

library(kernlab)
SVM <- ksvm(train_data$status ~ xAge_Low+xAge_Middle+Sex+xOther+xWhite+
              xSize_Middle+xSize_Small+Histologic+xTotal_thyroidectomy+xOther
                +Radiation+xNot_examined+xPositive+Chemotherapy,
                data = train_data, kernel = "rbfdot")

SVM_pred. <- predict(SVM, train_data)
SVM_ROC <- roc(train_data$status,as.numeric(SVM_pred.))
lines(SVM_ROC, type="l",col="#005824")

SVM_pred. <- predict(SVM, test_data)
SVM_ROC <- roc(test_data$status,as.numeric(SVM_pred.))
lines(SVM_ROC, type="l",col="#005824")


logit <- glm(status ~ xAge_Low+xAge_Middle+Sex+xOther+xWhite+
               xSize_Middle+xSize_Small+Histologic+xTotal_thyroidectomy+xOther
             +Radiation+xNot_examined+xPositive+Chemotherapy,
             data = train_data,family=binomial)

LR_pred. <- predict(logit, train_data)
LR_ROC <- roc(train_data$status,LR_pred.)
lines(LR_ROC, type="l",col="#377EB8")

LR_pred. <- predict(logit, test_data)
LR_ROC <- roc(test_data$status,LR_pred.)
lines(LR_ROC, type="l",col="#377EB8")

library(randomForest)
mod <- randomForest(as.factor(status) ~ xAge_Low+xAge_Middle+Sex+xOther+xWhite+
                      xSize_Middle+xSize_Small+Histologic+xTotal_thyroidectomy+xOther
                    +Radiation+xNot_examined+xPositive+Chemotherapy,
                    data = train_data,
                    importance = T,
                    na.action = na.omit)

RF_pred. <- predict(mod, newdata=train_data,type="prob")
str(RF_pred.)
RF_pred. <- as.data.frame(RF_pred.)
a <- RF_pred.[,1]
RF_ROC <- roc(train_data$status,a)
lines(RF_ROC, type="l",col="#984EA3")

RF_pred. <- predict(mod, newdata=test_data,type="prob")
str(RF_pred.)
RF_pred. <- as.data.frame(RF_pred.)
a <- RF_pred.[,1]
RF_ROC <- roc(test_data$status,a )
lines(RF_ROC, type="l",col="#984EA3")


legend("bottomright",c(paste("AUC of XGBoost_ROC = 0.948"),
                  paste("AUC of SVM_ROC = 0.888"),
                  paste("AUC of LR_ROC = 0.873"),
                  paste("AUC of RF_ROC = 0.881")),
       x.intersp=1, y.intersp=0.9,
       lty= 1 ,lwd= 2,col=c("#E41A1C","#005824","#377EB8","#984EA3"),
       bty = "n",# bty框的类型
       seg.len=1,cex=0.8)# 

legend("bottomright",c(paste("AUC of XGBoost_ROC = 0.864"),
                  paste("AUC of SVM_ROC = 0.871"),
                  paste("AUC of LR_ROC = 0.889"),
                  paste("AUC of RF_ROC = 0.858")),
       x.intersp=1, y.intersp=0.9,
       lty= 1 ,lwd= 2,col=c("#E41A1C","#005824","#377EB8","#984EA3"),
       bty = "n",# bty框的类型
       seg.len=1,cex=0.8)# 


###############################################################################################################

#####################################################################    DCA
dca <- function(data, outcome, predictors, xstart=0.01, xstop=0.99, xby=0.01, 
                ymin=-0.05, probability=NULL, harm=NULL,graph=TRUE, intervention=FALSE, 
                interventionper=100, smooth=FALSE,loess.span=0.10) {
  
  # LOADING REQUIRED LIBRARIES
  require(stats)
  
  # data MUST BE A DATA FRAME
  if (class(data)!="data.frame") {
    stop("Input data must be class data.frame")
  }
  
  #ONLY KEEPING COMPLETE CASES
  data=data[complete.cases(data[append(outcome,predictors)]),append(outcome,predictors)]
  
  # outcome MUST BE CODED AS 0 AND 1
  if (max(data[[outcome]])>1 | min(data[[outcome]])<0) {
    stop("outcome cannot be less than 0 or greater than 1")
  }
  # xstart IS BETWEEN 0 AND 1
  if (xstart<0 | xstart>1) {
    stop("xstart must lie between 0 and 1")
  }
  
  # xstop IS BETWEEN 0 AND 1
  if (xstop<0 | xstop>1) {
    stop("xstop must lie between 0 and 1")
  }
  
  # xby IS BETWEEN 0 AND 1
  if (xby<=0 | xby>=1) {
    stop("xby must lie between 0 and 1")
  }
  
  # xstart IS BEFORE xstop
  if (xstart>=xstop) {
    stop("xstop must be larger than xstart")
  }
  
  #STORING THE NUMBER OF PREDICTORS SPECIFIED
  pred.n=length(predictors)
  
  #IF probability SPECIFIED ENSURING THAT EACH PREDICTOR IS INDICATED AS A YES OR NO
  if (length(probability)>0 & pred.n!=length(probability)) {
    stop("Number of probabilities specified must be the same as the number of predictors being checked.")
  }
  
  #IF harm SPECIFIED ENSURING THAT EACH PREDICTOR HAS A SPECIFIED HARM
  if (length(harm)>0 & pred.n!=length(harm)) {
    stop("Number of harms specified must be the same as the number of predictors being checked.")
  }
  
  #INITIALIZING DEFAULT VALUES FOR PROBABILITES AND HARMS IF NOT SPECIFIED
  if (length(harm)==0) {
    harm=rep(0,pred.n)
  }
  if (length(probability)==0) {
    probability=rep(TRUE,pred.n)
  }
  
  
  #CHECKING THAT EACH probability ELEMENT IS EQUAL TO YES OR NO, 
  #AND CHECKING THAT PROBABILITIES ARE BETWEEN 0 and 1
  #IF NOT A PROB THEN CONVERTING WITH A LOGISTIC REGRESSION
  for(m in 1:pred.n) { 
    if (probability[m]!=TRUE & probability[m]!=FALSE) {
      stop("Each element of probability vector must be TRUE or FALSE")
    }
    if (probability[m]==TRUE & (max(data[predictors[m]])>1 | min(data[predictors[m]])<0)) {
      stop(paste(predictors[m],"must be between 0 and 1 OR sepcified as a non-probability in the probability option",sep=" "))  
    }
    if(probability[m]==FALSE) {
      model=NULL
      pred=NULL
      model=glm(data.matrix(data[outcome]) ~ data.matrix(data[predictors[m]]), family=binomial("logit"))
      pred=data.frame(model$fitted.values)
      pred=data.frame(pred)
      names(pred)=predictors[m]
      data=cbind(data[names(data)!=predictors[m]],pred)
      print(paste(predictors[m],"converted to a probability with logistic regression. Due to linearity assumption, miscalibration may occur.",sep=" "))
    }
  }
  
  # THE PREDICTOR NAMES CANNOT BE EQUAL TO all OR none.
  if (length(predictors[predictors=="all" | predictors=="none"])) {
    stop("Prediction names cannot be equal to all or none.")
  }  
  
  #########  CALCULATING NET BENEFIT   #########
  N=dim(data)[1]
  event.rate=colMeans(data[outcome])
  
  # CREATING DATAFRAME THAT IS ONE LINE PER THRESHOLD PER all AND none STRATEGY
  nb=data.frame(seq(from=xstart, to=xstop, by=xby))
  names(nb)="threshold"
  interv=nb
  
  nb["all"]=event.rate - (1-event.rate)*nb$threshold/(1-nb$threshold)
  nb["none"]=0
  
  # CYCLING THROUGH EACH PREDICTOR AND CALCULATING NET BENEFIT
  for(m in 1:pred.n){
    for(t in 1:length(nb$threshold)){
      # COUNTING TRUE POSITIVES AT EACH THRESHOLD
      tp=mean(data[data[[predictors[m]]]>=nb$threshold[t],outcome])*sum(data[[predictors[m]]]>=nb$threshold[t])
      # COUNTING FALSE POSITIVES AT EACH THRESHOLD
      fp=(1-mean(data[data[[predictors[m]]]>=nb$threshold[t],outcome]))*sum(data[[predictors[m]]]>=nb$threshold[t])
      #setting TP and FP to 0 if no observations meet threshold prob.
      if (sum(data[[predictors[m]]]>=nb$threshold[t])==0) {
        tp=0
        fp=0
      }
      
      # CALCULATING NET BENEFIT
      nb[t,predictors[m]]=tp/N - fp/N*(nb$threshold[t]/(1-nb$threshold[t])) - harm[m]
    }
    interv[predictors[m]]=(nb[predictors[m]] - nb["all"])*interventionper/(interv$threshold/(1-interv$threshold))
  }
  
  # CYCLING THROUGH EACH PREDICTOR AND SMOOTH NET BENEFIT AND INTERVENTIONS AVOIDED 
  for(m in 1:pred.n) {
    if (smooth==TRUE){
      lws=loess(data.matrix(nb[!is.na(nb[[predictors[m]]]),predictors[m]]) ~ data.matrix(nb[!is.na(nb[[predictors[m]]]),"threshold"]),span=loess.span)
      nb[!is.na(nb[[predictors[m]]]),paste(predictors[m],"_sm",sep="")]=lws$fitted
      
      lws=loess(data.matrix(interv[!is.na(nb[[predictors[m]]]),predictors[m]]) ~ data.matrix(interv[!is.na(nb[[predictors[m]]]),"threshold"]),span=loess.span)
      interv[!is.na(nb[[predictors[m]]]),paste(predictors[m],"_sm",sep="")]=lws$fitted
    }
  }
  
  # PLOTTING GRAPH IF REQUESTED
  if (graph==TRUE) {
    require(graphics)
    
    # PLOTTING INTERVENTIONS AVOIDED IF REQUESTED
    if(intervention==TRUE) {
      # initialize the legend label, color, and width using the standard specs of the none and all lines
      legendlabel <- NULL
      legendcolor <- NULL
      legendwidth <- NULL
      legendpattern <- NULL
      
      #getting maximum number of avoided interventions
      ymax=max(interv[predictors],na.rm = TRUE)
      
      #INITIALIZING EMPTY PLOT WITH LABELS
      plot(x=nb$threshold, y=nb$all, type="n" ,xlim=c(xstart, xstop), ylim=c(ymin, ymax), xlab="Threshold probability", ylab=paste("Net reduction in interventions per",interventionper,"patients"))
      
      #PLOTTING INTERVENTIONS AVOIDED FOR EACH PREDICTOR
      for(m in 1:pred.n) {
        if (smooth==TRUE){
          lines(interv$threshold,data.matrix(interv[paste(predictors[m],"_sm",sep="")]),col=m,lty=2)
        } else {
          lines(interv$threshold,data.matrix(interv[predictors[m]]),col=m,lty=2)
        }
        
        # adding each model to the legend
        legendlabel <- c(legendlabel, predictors[m])
        legendcolor <- c(legendcolor, m)
        legendwidth <- c(legendwidth, 1)
        legendpattern <- c(legendpattern, 2)
      }
    } else {
      # PLOTTING NET BENEFIT IF REQUESTED
      
      # initialize the legend label, color, and width using the standard specs of the none and all lines
      legendlabel <- c("None", "All")
      legendcolor <- c(17, 8)
      legendwidth <- c(2, 2)
      legendpattern <- c(1, 1)
      
      #getting maximum net benefit
      ymax=max(nb[names(nb)!="threshold"],na.rm = TRUE)
      
      # inializing new benfit plot with treat all option
      plot(x=nb$threshold, y=nb$all, type="l", col=8, lwd=2 ,xlim=c(xstart, xstop), ylim=c(ymin, ymax), xlab="Threshold probability", ylab="Net benefit")
      # adding treat none option
      lines(x=nb$threshold, y=nb$none,lwd=2)
      #PLOTTING net benefit FOR EACH PREDICTOR
      for(m in 1:pred.n) {
        if (smooth==TRUE){
          lines(nb$threshold,data.matrix(nb[paste(predictors[m],"_sm",sep="")]),col=m,lty=2) 
        } else {
          lines(nb$threshold,data.matrix(nb[predictors[m]]),col=m,lty=2)
        }
        # adding each model to the legend
        legendlabel <- c(legendlabel, predictors[m])
        legendcolor <- c(legendcolor, m)
        legendwidth <- c(legendwidth, 1)
        legendpattern <- c(legendpattern, 2)
      }
    }
    # then add the legend
    legend("topright", legendlabel, cex=0.8, col=legendcolor, lwd=legendwidth, lty=legendpattern)
    
  }
  
  #RETURNING RESULTS
  results=list() 
  results$N=N
  results$predictors=data.frame(cbind(predictors,harm,probability))
  names(results$predictors)=c("predictor","harm.applied","probability")
  results$interventions.avoided.per=interventionper
  results$net.benefit=nb
  results$interventions.avoided=interv
  
  return(results)
  
}  


library(nricens)
library(rms)
library(foreign)

#################  train
xgb <- xgboost(data = dtrain,max_depth=6, eta=0.5,  
               objective='binary:logistic', nround=25)

SVM <- ksvm(as.numeric(train_data$status) ~ xAge_Low+xAge_Middle+Sex+xOther+xWhite+
              xSize_Middle+xSize_Small+Histologic+xTotal_thyroidectomy+xOther
            +Radiation+xNot_examined+xPositive+Chemotherapy,
            data = train_data, kernel = "rbfdot")

logit <- glm(status ~  xAge_Low+xAge_Middle+Sex+xOther+xWhite+
               xSize_Middle+xSize_Small+Histologic+xTotal_thyroidectomy+xOther
             +Radiation+xNot_examined+xPositive+Chemotherapy,
             data = train_data,family=binomial(link="logit"))

RF <- randomForest(as.factor(status) ~ xAge_Low+xAge_Middle+Sex+xOther+xWhite+
                     xSize_Middle+xSize_Small+Histologic+xTotal_thyroidectomy+xOther
                   +Radiation+xNot_examined+xPositive+Chemotherapy,
                    data = train_data,
                    importance = T,
                    na.action = na.omit)


train_data$XGBoost_pred. = predict(xgb,newdata = dtrain,"response")
train_data$LR_pred. <- predict(logit,train_data,"response")
train_data$SVM_pred. <- predict(SVM, train_data,"response")
data.scale <- (train_data$SVM_pred.-min(train_data$SVM_pred.))/(max(train_data$SVM_pred.)-min(train_data$SVM_pred.))
train_data$SVM_pred. <- data.scale

RF_pred. <- predict(mod, newdata=train_data,type="prob")
str(RF_pred.)
RF_pred. <- as.data.frame(RF_pred.)
a <- RF_pred.[,2]
train_data <- cbind(train_data,a)
train_data$RF_pred. <- train_data$a


library(tidyverse)
library(lubridate)
train_data <- as.data.frame(train_data)
dca(data=train_data, outcome="status", predictors=c("XGBoost_pred.","LR_pred.","SVM_pred.","RF_pred."))

##############           test
xgb <- xgboost(data = dtrain,max_depth=6, eta=0.5,  
               objective='binary:logistic', nround=25)

SVM <- ksvm(as.numeric(train_data$status) ~ xMage+xSage+Sex+xL+xM+xFTC+xMTC+xPTC+
              xNot_examined+xPositive+xOther+xThyroidectomy
            +Radiation+Chemotherapy,
            data = train_data, kernel = "rbfdot")

logit <- glm(status ~  xMage+xSage+Sex++xL+xM+xFTC+xMTC+xPTC+
               xNot_examined+xPositive+xOther+xThyroidectomy
             +Radiation+Chemotherapy,
             data = train_data,family=binomial(link="logit"))
RF <- randomForest(as.factor(status) ~ xMage+xSage+Sex+xL+xM+xFTC+xMTC+xPTC+
                     xNot_examined+xPositive+xOther+xThyroidectomy
                   +Radiation+Chemotherapy,
                   data = train_data,
                   importance = T,
                   na.action = na.omit)

test_data$XGBoost_pred. = predict(xgb,newdata = dtest,"response")
test_data$LR_pred. <- predict(logit,test_data,"response")
test_data$SVM_pred. <- predict(SVM, test_data,"response")
data.scale <- (test_data$SVM_pred.-min(test_data$SVM_pred.))/(max(test_data$SVM_pred.)-min(test_data$SVM_pred.))
test_data$SVM_pred. <- data.scale

RF_pred. <- predict(RF, newdata=test_data,type="prob")
str(RF_pred.)
RF_pred. <- as.data.frame(RF_pred.)
a <- RF_pred.[,2]
test_data <- cbind(test_data,a)
test_data$RF_pred. <- test_data$a


test_data <- as.data.frame(test_data)
dca(data=test_data, outcome="status", predictors=c("XGBoost_pred.","LR_pred.","SVM_pred.","RF_pred."))#, smooth=F, xstop=0.50)"SVM__pred.","LR_pred."

##############################################################################校准曲线
plotCalibration <-
  function (data, cOutcome, predRisk, groups, rangeaxis, plottitle, 
            xlabel, ylabel, filename, fileplot, plottype) 
  {
    if (missing(groups)) {
      groups <- 5
    }
    if (missing(plottitle)) {
      plottitle <- " "
    }
    if (missing(xlabel)) {
      xlabel <- "Predicted risk"
    }
    if (missing(ylabel)) {
      ylabel <- "Observed risk"
    }
    if (missing(rangeaxis)) {
      rangeaxis <- c(0, 1)
    }
    p = predRisk
    y = data[, cOutcome]
    if (length(unique(y)) != 2) {
      stop(" The specified outcome is not a binary variable.\n")
    }
    else {
      matres <- matrix(NA, nrow = groups, ncol = 5)
      sor <- order(p)
      p <- p[sor]
      y <- y[sor]
      groep <- cut2(p, g = groups)
      total <- tapply(y, groep, length)
      predicted <- round(tapply(p, groep, sum), 2)
      observed <- tapply(y, groep, sum)
      meanpred <- round(tapply(p, groep, mean), 3)
      meanobs <- round(tapply(y, groep, mean), 3)
      matres <- cbind(total, meanpred, meanobs, predicted, 
                      observed)
      plot(matres[, 2], matres[, 3], main = plottitle, xlab = xlabel, 
           ylab = ylabel, pch = 16,type = "o", ps = 2, xlim = rangeaxis,
           ylim = rangeaxis, cex.lab = 1.2, cex.axis = 1.1, 
           las = 1)
      contr <- ((observed - predicted)^2)/(total * meanpred * 
                                             (1 - meanpred))
      chisqr <- sum(contr)
      df <- (groups - 2)
      pval <- 1 - pchisq(chisqr, df)
      lines(x = c(0, 1), y = c(0, 1))
      if (missing(plottype)) {
        plottype <- "jpg"
      }
      if (!missing(fileplot)) 
        savePlot(filename = fileplot, type = plottype, device = dev.cur())
      if (!missing(filename)) 
        write.table(matres, file = filename, row.names = TRUE, 
                    sep = "\t", dec = ",")
      out <- list(Table_HLtest = matres, Chi_square = round(chisqr, 
                                                            3), df = df, p_value = round(pval, 4))
      return(out)
    }
  }

cOutcome <- 16
groups <- 5
rangeaxis <- c(0,1)
#######################         train                 
plotCalibration(data=train_data, cOutcome=cOutcome, 
                predRisk=train_data$LR_pred.,groups = groups)
legend(0.7,0.1,c(paste("Calibration plot of LR_pred." )),bty = "n",seg.len=1,cex=0.8)

plotCalibration(data=train_data, cOutcome=cOutcome, 
                predRisk=train_data$XGBoost_pred.,groups = groups)
legend(0.7,0.1,c(paste("Calibration plot of XGBoost_pre." )),bty = "n",seg.len=1,cex=0.8)

plotCalibration(data=train_data, cOutcome=cOutcome, 
                predRisk=train_data$SVM_pred.,groups = groups)
legend(0.7,0.1,c(paste("Calibration plot of SVM_pred." )),bty = "n",seg.len=1,cex=0.8)

plotCalibration(data=train_data, cOutcome=cOutcome, 
                predRisk=train_data$RF_pred.,groups = groups)
legend(0.7,0.1,c(paste("Calibration plot of RF_pred." )),bty = "n",seg.len=1,cex=0.8)


#####################   test

plotCalibration(data=test_data, cOutcome=cOutcome, 
                predRisk=test_data$LR_pred.,groups = groups)
legend(0.7,0.1,c(paste("Calibration plot of LR_pred." )),bty = "n",seg.len=1,cex=0.8)

plotCalibration(data=test_data, cOutcome=cOutcome, 
                predRisk=test_data$XGBoost_pred.,groups = groups)
legend(0.7,0.1,c(paste("Calibration plot of XGBoost_pre." )),bty = "n",seg.len=1,cex=0.8)

plotCalibration(data=test_data, cOutcome=cOutcome, 
                predRisk=test_data$SVM_pred.,groups = groups)
legend(0.7,0.1,c(paste("Calibration plot of SVM_pred." )),bty = "n",seg.len=1,cex=0.8)

plotCalibration(data=test_data, cOutcome=cOutcome, 
                predRisk=test_data$RF_pred.,groups = groups)
legend(0.7,0.1,c(paste("Calibration plot of RF_pred." )),bty = "n",seg.len=1,cex=0.8)




































































