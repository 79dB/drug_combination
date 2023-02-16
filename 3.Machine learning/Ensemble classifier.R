#######github项目
##3.Machine learning

library(data.table)
library(randomForest)
library(varSelRF)
library(pROC)
library(ROCR)
library(caret)
library(modEvA)
library(e1071)
library(nnet)

##整合投票相关信息
setwd("G:\\毕业设计R程序\\结果整理\\github\\3.Machine learning\\Raw data")
##设置种子
seed <- 431
set.seed(seed)
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
seed <- 431
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##数据
set <- rbind(pos,neg)

##SVM
seed <- 3
set.seed(seed)
svm_train <- set[unique(sample(nrow(set),nrow(set),replace = T)),]
svm_test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(svm_train)))),]
set.seed(seed)
svm_model <- svm(class~.,data = svm_train[,-c(1,2)],kernel="radial",
                 cost=1,gamma=1/ncol(svm_train[,-c(1,2)]),probability = TRUE)
svm_pred <- predict(svm_model,newdata = svm_test[,-c(1,2)])

##Logistic分类器
seed <- 3
set.seed(seed)
log_train <- set[unique(sample(nrow(set),nrow(set),replace = T)),]
log_test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(log_train)))),]
set.seed(seed)
log_model <- glm(class~.,data = log_train[,-c(1,2)],family = "binomial")
log_prob <- predict(object =log_model,newdata=log_test[,-c(1,2)],type = "response")
log_pred<-ifelse(log_prob>=0.5,"synergistic","antagonistic")
log_pred <- factor(log_pred)

##ANN分类器
seed <- 6
set.seed(seed)
ann_train <- set[unique(sample(nrow(set),nrow(set),replace = T)),]
ann_test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(ann_train)))),]
set.seed(seed)
ann_model <- nnet(class~.,data = ann_train[,-c(1,2)],size=10)
ann_pred <- predict(ann_model,newdata = ann_test[,-c(1,2)],type="class")
ann_pred <- factor(ann_pred)
names(ann_pred) <- rownames(ann_test)

##NB分类器
seed <- 6
set.seed(seed)
bys_train <- set[unique(sample(nrow(set),nrow(set),replace = T)),]
bys_test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(bys_train)))),]
set.seed(seed)
bys_model <- naiveBayes(class~. ,data = bys_train[,-c(1,2)])
bys_pred <- predict(bys_model,newdata = bys_test[,-c(1,2)])
names(bys_pred) <- rownames(bys_test)

##交集信息
test_ins <- c(rownames(svm_test),rownames(bys_test),
              rownames(log_test),rownames(ann_test))
test_ins <- as.data.frame(table(test_ins))
test_ins$test_ins <- as.character(test_ins$test_ins)
test_ins <- test_ins[test_ins$Freq>=3,]

##添加类别
for (i in 1:nrow(test_ins)) {
  if(length(svm_pred[names(svm_pred)%in%test_ins$test_ins[i]])!=0)
  {test_ins$SVM[i] <- svm_pred[names(svm_pred)%in%test_ins$test_ins[i]]
  }else
    test_ins$SVM[i] <- NA
  
}

for (i in 1:nrow(test_ins)) {
  if(length(bys_pred[names(bys_pred)%in%test_ins$test_ins[i]])!=0)
  {test_ins$Bys[i] <- bys_pred[names(bys_pred)%in%test_ins$test_ins[i]]
  }else
    test_ins$Bys[i] <- NA
  
}

for (i in 1:nrow(test_ins)) {
  if(length(log_pred[names(log_pred)%in%test_ins$test_ins[i]])!=0)
  {test_ins$Log[i] <- log_pred[names(log_pred)%in%test_ins$test_ins[i]]
  }else
    test_ins$Log[i] <- NA
  
}

for (i in 1:nrow(test_ins)) {
  if(length(ann_pred[names(ann_pred)%in%test_ins$test_ins[i]])!=0)
  {test_ins$ANN[i] <- ann_pred[names(ann_pred)%in%test_ins$test_ins[i]]
  }else
    test_ins$ANN[i] <- NA
  
}
##去除频率信息
test_ins <- test_ins[,-2]

##结果数据框
result <- data.frame(classification=c("EL1","EL2","EL3","EL4","EL5"))

##不同类别进行投票判断
##Ensemble Model 1 :SVM + NB + LR
el1 <- test_ins[,-5]
for (i in 1:nrow(el1)) {
  if(length(which(el1[i,-1]%in%2))>=2){
    el1$pred_class[i] <- "synergistic"
  }else
    el1$pred_class[i] <- "antagonistic"
}
##添加正式类别
for (i in 1:nrow(el1)) {
  el1$actual_class[i] <- as.character(set$class[which(rownames(set)%in%el1$test_ins[i])])
}
el1$pred_class <- factor(el1$pred_class)
el1$actual_class <- factor(el1$actual_class)
##四格表
el1_perf <- table(el1$actual_class,el1$pred_class,
                  dnn=c("Actual","Predicted"))
el1_perf
##整合ROC曲线
el1_roc <- roc(as.numeric(el1$actual_class),
               as.numeric(el1$pred_class))
plot(el1_roc,print.auc=T)
##计算效能得分
##灵敏度,又叫真阳性率
result$Sensitivity[1] <- el1_perf[2,2]/(el1_perf[2,2]+el1_perf[2,1])
##特异度,又叫真阴性率
result$Specificity[1] <- el1_perf[1,1]/(el1_perf[1,1]+el1_perf[1,2])
##伪阳性率
result$FPR[1] <- 1-result$Specificity[1]
##假发现率
result$FDR[1] <- 1-result$Sensitivity[1]
##准确度
result$Accuracy[1] <- (el1_perf[1,1]+el1_perf[2,2])/sum(el1_perf)
##阳性预测率
result$PPV[1] <- el1_perf[2,2]/(el1_perf[2,2]+el1_perf[1,2])
##阴性预测率
result$NPV[1] <- el1_perf[1,1]/(el1_perf[1,1]+el1_perf[2,1])
##F1评分
result$F1[1] <- (2*el1_perf[2,2])/(2*el1_perf[2,2]+el1_perf[1,2]+el1_perf[2,1])
##Mathews系数(MCC)，即Phi相关系数
result$MCC[1] <- (el1_perf[2,2]*el1_perf[1,1]-el1_perf[1,2]*el1_perf[2,1])/sqrt((el1_perf[2,2]+el1_perf[1,2])*(el1_perf[2,2]+el1_perf[2,1])*(el1_perf[1,1]+el1_perf[1,2])*(el1_perf[1,1]+el1_perf[2,1]))
##基尼指数
result$Gini[1] <- 2*result$Accuracy[1] - 1
##Kappa,通常0.75以上表示一致性较满意，0.4以下一致性不好
result$Kappa[1] <- ((el1_perf[1,1]+el1_perf[1,2])*(el1_perf[1,1]+el1_perf[2,1])+(el1_perf[2,1]+el1_perf[2,2])*(el1_perf[1,2]+el1_perf[2,2]))/(sum(el1_perf)*sum(el1_perf))
result$AUC[1] <- el1_roc$auc

##Ensemble Model 2 :SVM + NB + ANN
el2 <- test_ins[,-4]
for (i in 1:nrow(el2)) {
  if(length(which(el2[i,-1]%in%2))>=2){
    el2$pred_class[i] <- "synergistic"
  }else
    el2$pred_class[i] <- "antagonistic"
}
##添加正式类别
for (i in 1:nrow(el2)) {
  el2$actual_class[i] <- as.character(set$class[which(rownames(set)%in%el2$test_ins[i])])
}
el2$pred_class <- factor(el2$pred_class)
el2$actual_class <- factor(el2$actual_class)
##四格表
el2_perf <- table(el2$actual_class,el2$pred_class,
                  dnn=c("Actual","Predicted"))
el2_perf
##整合ROC曲线
el2_roc <- roc(as.numeric(el2$actual_class),
               as.numeric(el2$pred_class))
plot(el2_roc,print.auc=T)
##计算效能得分
##灵敏度,又叫真阳性率
result$Sensitivity[2] <- el2_perf[2,2]/(el2_perf[2,2]+el2_perf[2,1])
##特异度,又叫真阴性率
result$Specificity[2] <- el2_perf[1,1]/(el2_perf[1,1]+el2_perf[1,2])
##伪阳性率
result$FPR[2] <- 1-result$Specificity[2]
##假发现率
result$FDR[2] <- 1-result$Sensitivity[2]
##准确度
result$Accuracy[2] <- (el2_perf[1,1]+el2_perf[2,2])/sum(el2_perf)
##阳性预测率
result$PPV[2] <- el2_perf[2,2]/(el2_perf[2,2]+el2_perf[1,2])
##阴性预测率
result$NPV[2] <- el2_perf[1,1]/(el2_perf[1,1]+el2_perf[2,1])
##F1评分
result$F1[2] <- (2*el2_perf[2,2])/(2*el2_perf[2,2]+el2_perf[1,2]+el2_perf[2,1])
##Mathews系数(MCC)，即Phi相关系数
result$MCC[2] <- (el2_perf[2,2]*el2_perf[1,1]-el2_perf[1,2]*el2_perf[2,1])/sqrt((el2_perf[2,2]+el2_perf[1,2])*(el2_perf[2,2]+el2_perf[2,1])*(el2_perf[1,1]+el2_perf[1,2])*(el2_perf[1,1]+el2_perf[2,1]))
##基尼指数
result$Gini[2] <- 2*result$Accuracy[2] - 1
##Kappa,通常0.75以上表示一致性较满意，0.4以下一致性不好
result$Kappa[2] <- ((el2_perf[1,1]+el2_perf[1,2])*(el2_perf[1,1]+el2_perf[2,1])+(el2_perf[2,1]+el2_perf[2,2])*(el2_perf[1,2]+el2_perf[2,2]))/(sum(el2_perf)*sum(el2_perf))
result$AUC[2] <- el2_roc$auc

##Ensemble Model 3 :SVM + LR + ANN
el3 <- test_ins[,-3]
for (i in 1:nrow(el3)) {
  if(length(which(el3[i,-1]%in%2))>=2){
    el3$pred_class[i] <- "synergistic"
  }else
    el3$pred_class[i] <- "antagonistic"
}
##添加正式类别
for (i in 1:nrow(el3)) {
  el3$actual_class[i] <- as.character(set$class[which(rownames(set)%in%el3$test_ins[i])])
}
el3$pred_class <- factor(el3$pred_class)
el3$actual_class <- factor(el3$actual_class)
##四格表
el3_perf <- table(el3$actual_class,el3$pred_class,
                  dnn=c("Actual","Predicted"))
el3_perf
##整合ROC曲线
el3_roc <- roc(as.numeric(el3$actual_class),
               as.numeric(el3$pred_class))
plot(el3_roc,print.auc=T)
##计算效能得分
##灵敏度,又叫真阳性率
result$Sensitivity[3] <- el3_perf[2,2]/(el3_perf[2,2]+el3_perf[2,1])
##特异度,又叫真阴性率
result$Specificity[3] <- el3_perf[1,1]/(el3_perf[1,1]+el3_perf[1,2])
##伪阳性率
result$FPR[3] <- 1-result$Specificity[3]
##假发现率
result$FDR[3] <- 1-result$Sensitivity[3]
##准确度
result$Accuracy[3] <- (el3_perf[1,1]+el3_perf[2,2])/sum(el3_perf)
##阳性预测率
result$PPV[3] <- el3_perf[2,2]/(el3_perf[2,2]+el3_perf[1,2])
##阴性预测率
result$NPV[3] <- el3_perf[1,1]/(el3_perf[1,1]+el3_perf[2,1])
##F1评分
result$F1[3] <- (2*el3_perf[2,2])/(2*el3_perf[2,2]+el3_perf[1,2]+el3_perf[2,1])
##Mathews系数(MCC)，即Phi相关系数
result$MCC[3] <- (el3_perf[2,2]*el3_perf[1,1]-el3_perf[1,2]*el3_perf[2,1])/sqrt((el3_perf[2,2]+el3_perf[1,2])*(el3_perf[2,2]+el3_perf[2,1])*(el3_perf[1,1]+el3_perf[1,2])*(el3_perf[1,1]+el3_perf[2,1]))
##基尼指数
result$Gini[3] <- 2*result$Accuracy[3] - 1
##Kappa,通常0.75以上表示一致性较满意，0.4以下一致性不好
result$Kappa[3] <- ((el3_perf[1,1]+el3_perf[1,2])*(el3_perf[1,1]+el3_perf[2,1])+(el3_perf[2,1]+el3_perf[2,2])*(el3_perf[1,2]+el3_perf[2,2]))/(sum(el3_perf)*sum(el3_perf))
result$AUC[3] <- el3_roc$auc

##Ensemble Model 4 :NB + LR + ANN
el4 <- test_ins[,-2]
for (i in 1:nrow(el4)) {
  if(length(which(el4[i,-1]%in%2))>=2){
    el4$pred_class[i] <- "synergistic"
  }else
    el4$pred_class[i] <- "antagonistic"
}
##添加正式类别
for (i in 1:nrow(el4)) {
  el4$actual_class[i] <- as.character(set$class[which(rownames(set)%in%el4$test_ins[i])])
}
el4$pred_class <- factor(el4$pred_class)
el4$actual_class <- factor(el4$actual_class)
##四格表
el4_perf <- table(el4$actual_class,el4$pred_class,
                  dnn=c("Actual","Predicted"))
el4_perf
##整合ROC曲线
el4_roc <- roc(as.numeric(el4$actual_class),
               as.numeric(el4$pred_class))
plot(el4_roc,print.auc=T)
##计算效能得分
##灵敏度,又叫真阳性率
result$Sensitivity[4] <- el4_perf[2,2]/(el4_perf[2,2]+el4_perf[2,1])
##特异度,又叫真阴性率
result$Specificity[4] <- el4_perf[1,1]/(el4_perf[1,1]+el4_perf[1,2])
##伪阳性率
result$FPR[4] <- 1-result$Specificity[4]
##假发现率
result$FDR[4] <- 1-result$Sensitivity[4]
##准确度
result$Accuracy[4] <- (el4_perf[1,1]+el4_perf[2,2])/sum(el4_perf)
##阳性预测率
result$PPV[4] <- el4_perf[2,2]/(el4_perf[2,2]+el4_perf[1,2])
##阴性预测率
result$NPV[4] <- el4_perf[1,1]/(el4_perf[1,1]+el4_perf[2,1])
##F1评分
result$F1[4] <- (2*el4_perf[2,2])/(2*el4_perf[2,2]+el4_perf[1,2]+el4_perf[2,1])
##Mathews系数(MCC)，即Phi相关系数
result$MCC[4] <- (el4_perf[2,2]*el4_perf[1,1]-el4_perf[1,2]*el4_perf[2,1])/sqrt((el4_perf[2,2]+el4_perf[1,2])*(el4_perf[2,2]+el4_perf[2,1])*(el4_perf[1,1]+el4_perf[1,2])*(el4_perf[1,1]+el4_perf[2,1]))
##基尼指数
result$Gini[4] <- 2*result$Accuracy[4] - 1
##Kappa,通常0.75以上表示一致性较满意，0.4以下一致性不好
result$Kappa[4] <- ((el4_perf[1,1]+el4_perf[1,2])*(el4_perf[1,1]+el4_perf[2,1])+(el4_perf[2,1]+el4_perf[2,2])*(el4_perf[1,2]+el4_perf[2,2]))/(sum(el4_perf)*sum(el4_perf))
result$AUC[4] <- el4_roc$auc

##Ensemble Model 5 :SVM + NB + LR + ANN
el5 <- test_ins
for (i in 1:nrow(el5)) {
  if(length(which(el5[i,-1]%in%2))>=3){
    el5$pred_class[i] <- "synergistic"
  }else
    el5$pred_class[i] <- "antagonistic"
}
##添加正式类别
for (i in 1:nrow(el5)) {
  el5$actual_class[i] <- as.character(set$class[which(rownames(set)%in%el5$test_ins[i])])
}
el5$pred_class <- factor(el5$pred_class)
el5$actual_class <- factor(el5$actual_class)
##四格表
el5_perf <- table(el5$actual_class,el5$pred_class,
                  dnn=c("Actual","Predicted"))
el5_perf
test_perf <- el5_perf
##整合ROC曲线
el5_roc <- roc(as.numeric(el5$actual_class),
               as.numeric(el5$pred_class))
plot(el5_roc,print.auc=T)
##计算效能得分
##灵敏度,又叫真阳性率
result$Sensitivity[5] <- test_perf[2,2]/(test_perf[2,2]+test_perf[2,1])
##特异度,又叫真阴性率
result$Specificity[5] <- test_perf[1,1]/(test_perf[1,1]+test_perf[1,2])
##伪阳性率
result$FPR[5] <- 1-result$Specificity[5]
##假发现率
result$FDR[5] <- 1-result$Sensitivity[5]
##准确度
result$Accuracy[5] <- (test_perf[1,1]+test_perf[2,2])/sum(test_perf)
##阳性预测率
result$PPV[5] <- test_perf[2,2]/(test_perf[2,2]+test_perf[1,2])
##阴性预测率
result$NPV[5] <- test_perf[1,1]/(test_perf[1,1]+test_perf[2,1])
##F1评分
result$F1[5] <- (2*test_perf[2,2])/(2*test_perf[2,2]+test_perf[1,2]+test_perf[2,1])
##Mathews系数(MCC)，即Phi相关系数
result$MCC[5] <- (test_perf[2,2]*test_perf[1,1]-test_perf[1,2]*test_perf[2,1])/sqrt((test_perf[2,2]+test_perf[1,2])*(test_perf[2,2]+test_perf[2,1])*(test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1]))
##基尼指数
result$Gini[5] <- 2*result$Accuracy[5] - 1
##Kappa,通常0.75以上表示一致性较满意，0.4以下一致性不好
result$Kappa[5] <- ((test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1])+(test_perf[2,1]+test_perf[2,2])*(test_perf[1,2]+test_perf[2,2]))/(sum(test_perf)*sum(test_perf))
result$AUC[5] <- el5_roc$auc
