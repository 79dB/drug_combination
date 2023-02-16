#######github项目
##3.Machine learning

library(data.table)

result <- data.frame(seed=seq(1,1000))
setwd("G:\\毕业设计R程序\\结果整理\\github\\3.Machine learning\\Raw data")
##SVM种子选择为322
for (i in 1:1000) {
  std <- as.data.frame(fread("Gold standard for drug combination.csv"))
  std$class <- factor(std$class)
  set.seed(i)
  pos <- std[std$class=="synergistic",]
  neg <- std[std$class=="antagonistic",]
  ##抽取1：1
  neg <- neg[sample(nrow(neg),nrow(pos)),]
  ##设置种子
  set.seed(i)
  ##抽取train和test
  train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
  set <- rbind(pos,neg)
  test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(train)))),]
  svm_model <- svm(class~.,data = train[,-c(1,2)],kernel="radial",
                   cost=1,gamma=1/ncol(train[,-c(1,2)]),probability = TRUE)
  ##SVM预测
  svm_pred <- predict(svm_model,newdata = test[,-c(1,2)])
  svm_prob <- attr(predict(svm_model,newdata = test[,-c(1,2)],probability = T), "probabilities")
  ##生成四格表
  svm_perf<-table(test$class,svm_pred,dnn=c("Actual","Predicted"))
  test_perf <- svm_perf
  svm_roc <- roc(test$class,svm_prob[,1])
  ##灵敏度,又叫真阳性率
  result$Sensitivity[i] <- test_perf[2,2]/(test_perf[2,2]+test_perf[2,1])
  ##特异度,又叫真阴性率
  result$Specificity[i] <- test_perf[1,1]/(test_perf[1,1]+test_perf[1,2])
  ##伪阳性率
  result$FPR[i] <- 1-result$Specificity[i]
  ##假发现率
  result$FDR[i] <- 1-result$Sensitivity[i]
  ##准确度
  result$Accuracy[i] <- (test_perf[1,1]+test_perf[2,2])/sum(test_perf)
  ##阳性预测率
  result$PPV[i] <- test_perf[2,2]/(test_perf[2,2]+test_perf[1,2])
  ##阴性预测率
  result$NPV[i] <- test_perf[1,1]/(test_perf[1,1]+test_perf[2,1])
  ##F1评分
  result$F1[i] <- (2*test_perf[2,2])/(2*test_perf[2,2]+test_perf[1,2]+test_perf[2,1])
  ##Mathews系数(MCC)，即Phi相关系数
  result$MCC[i] <- (test_perf[2,2]*test_perf[1,1]-test_perf[1,2]*test_perf[2,1])/sqrt((test_perf[2,2]+test_perf[1,2])*(test_perf[2,2]+test_perf[2,1])*(test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1]))
  ##基尼指数
  result$Gini[i] <- 2*result$Accuracy[i] - 1
  ##Kappa,通常0.75以上表示一致性较满意，0.4以下一致性不好
  result$Kappa[i] <- ((test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1])+(test_perf[2,1]+test_perf[2,2])*(test_perf[1,2]+test_perf[2,2]))/(sum(test_perf)*sum(test_perf))
  result$AUC[i] <- round(svm_roc$auc,3)
  
  print(i)
}


##Bys种子选择为623
for (i in 1:1000) {
  std <- as.data.frame(fread("Gold standard for drug combination.csv"))
  std$class <- factor(std$class)
  set.seed(i)
  pos <- std[std$class=="synergistic",]
  neg <- std[std$class=="antagonistic",]
  ##抽取1：1
  neg <- neg[sample(nrow(neg),nrow(pos)),]
  ##设置种子
  set.seed(i)
  ##抽取train和test
  train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
  set <- rbind(pos,neg)
  test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(train)))),]
  bys_model <- naiveBayes(class~. ,data = train[,-c(1,2)])
  ##贝叶斯预测
  bys_pred <- predict(bys_model,newdata = test[,-c(1,2)])
  ##提取分类概率
  bys_prob <- predict(bys_model,newdata = test[,-c(1,2)],type = "raw")
  ##生成四格表
  bys_perf<-table(test$class,bys_pred,dnn=c("Actual","Predicted"))
  test_perf <- bys_perf
  bys_roc <- roc(test$class,bys_prob[,2])
  ##灵敏度,又叫真阳性率
  result$Sensitivity[i] <- test_perf[2,2]/(test_perf[2,2]+test_perf[2,1])
  ##特异度,又叫真阴性率
  result$Specificity[i] <- test_perf[1,1]/(test_perf[1,1]+test_perf[1,2])
  ##伪阳性率
  result$FPR[i] <- 1-result$Specificity[i]
  ##假发现率
  result$FDR[i] <- 1-result$Sensitivity[i]
  ##准确度
  result$Accuracy[i] <- (test_perf[1,1]+test_perf[2,2])/sum(test_perf)
  ##阳性预测率
  result$PPV[i] <- test_perf[2,2]/(test_perf[2,2]+test_perf[1,2])
  ##阴性预测率
  result$NPV[i] <- test_perf[1,1]/(test_perf[1,1]+test_perf[2,1])
  ##F1评分
  result$F1[i] <- (2*test_perf[2,2])/(2*test_perf[2,2]+test_perf[1,2]+test_perf[2,1])
  ##Mathews系数(MCC)，即Phi相关系数
  result$MCC[i] <- (test_perf[2,2]*test_perf[1,1]-test_perf[1,2]*test_perf[2,1])/sqrt((test_perf[2,2]+test_perf[1,2])*(test_perf[2,2]+test_perf[2,1])*(test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1]))
  ##基尼指数
  result$Gini[i] <- 2*result$Accuracy[i] - 1
  ##Kappa,通常0.75以上表示一致性较满意，0.4以下一致性不好
  result$Kappa[i] <- ((test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1])+(test_perf[2,1]+test_perf[2,2])*(test_perf[1,2]+test_perf[2,2]))/(sum(test_perf)*sum(test_perf))
  result$AUC[i] <- round(bys_roc$auc,3)
  
  print(i)
}

##Log种子选择为37
for (i in 1:1000) {
  std <- as.data.frame(fread("Gold standard for drug combination.csv"))
  std$class <- factor(std$class)
  set.seed(i)
  pos <- std[std$class=="synergistic",]
  neg <- std[std$class=="antagonistic",]
  ##抽取1：1
  neg <- neg[sample(nrow(neg),nrow(pos)),]
  ##设置种子
  set.seed(i)
  ##抽取train和test
  train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
  set <- rbind(pos,neg)
  test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(train)))),]
  log_model <- glm(class~.,data = train[,-c(1,2)],family = "binomial")
  log_prob <- predict(object =log_model,newdata=test[,-c(1,2)],type = "response")
  log_pred<-ifelse(log_prob>=0.5,"synergistic","antagonistic")
  log_pred <- factor(log_pred)
  log_perf<-table(test$class,log_pred,dnn=c("Actual","Predicted"))
  test_perf <- log_perf
  log_roc <- roc(test$class,log_prob)
  ##灵敏度,又叫真阳性率
  result$Sensitivity[i] <- test_perf[2,2]/(test_perf[2,2]+test_perf[2,1])
  ##特异度,又叫真阴性率
  result$Specificity[i] <- test_perf[1,1]/(test_perf[1,1]+test_perf[1,2])
  ##伪阳性率
  result$FPR[i] <- 1-result$Specificity[i]
  ##假发现率
  result$FDR[i] <- 1-result$Sensitivity[i]
  ##准确度
  result$Accuracy[i] <- (test_perf[1,1]+test_perf[2,2])/sum(test_perf)
  ##阳性预测率
  result$PPV[i] <- test_perf[2,2]/(test_perf[2,2]+test_perf[1,2])
  ##阴性预测率
  result$NPV[i] <- test_perf[1,1]/(test_perf[1,1]+test_perf[2,1])
  ##F1评分
  result$F1[i] <- (2*test_perf[2,2])/(2*test_perf[2,2]+test_perf[1,2]+test_perf[2,1])
  ##Mathews系数(MCC)，即Phi相关系数
  result$MCC[i] <- (test_perf[2,2]*test_perf[1,1]-test_perf[1,2]*test_perf[2,1])/sqrt((test_perf[2,2]+test_perf[1,2])*(test_perf[2,2]+test_perf[2,1])*(test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1]))
  ##基尼指数
  result$Gini[i] <- 2*result$Accuracy[i] - 1
  ##Kappa,通常0.75以上表示一致性较满意，0.4以下一致性不好
  result$Kappa[i] <- ((test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1])+(test_perf[2,1]+test_perf[2,2])*(test_perf[1,2]+test_perf[2,2]))/(sum(test_perf)*sum(test_perf))
  result$AUC[i] <- round(log_roc$auc,3)
  
  print(i)
}

##ANN种子选择为150
for (i in 1:1000) {
  std <- as.data.frame(fread("Gold standard for drug combination.csv"))
  std$class <- factor(std$class)
  set.seed(i)
  pos <- std[std$class=="synergistic",]
  neg <- std[std$class=="antagonistic",]
  ##抽取1：1
  neg <- neg[sample(nrow(neg),nrow(pos)),]
  ##设置种子
  set.seed(i)
  ##抽取train和test
  train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
  set <- rbind(pos,neg)
  test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(train)))),]
  ann_model <- nnet(class~.,data = train[,-c(1,2)],size=10)
  #测试集
  ann_pred <- predict(ann_model,newdata = test[,-c(1,2)],type="class")
  ann_prob <- predict(ann_model,newdata = test[,-c(1,2)])
  ann_perf<-table(test$class,ann_pred,dnn=c("Actual","Predicted"))
  test_perf <- ann_perf
  log_roc <- roc(test$class,ann_prob[,1])
  ##灵敏度,又叫真阳性率
  result$Sensitivity[i] <- test_perf[2,2]/(test_perf[2,2]+test_perf[2,1])
  ##特异度,又叫真阴性率
  result$Specificity[i] <- test_perf[1,1]/(test_perf[1,1]+test_perf[1,2])
  ##伪阳性率
  result$FPR[i] <- 1-result$Specificity[i]
  ##假发现率
  result$FDR[i] <- 1-result$Sensitivity[i]
  ##准确度
  result$Accuracy[i] <- (test_perf[1,1]+test_perf[2,2])/sum(test_perf)
  ##阳性预测率
  result$PPV[i] <- test_perf[2,2]/(test_perf[2,2]+test_perf[1,2])
  ##阴性预测率
  result$NPV[i] <- test_perf[1,1]/(test_perf[1,1]+test_perf[2,1])
  ##F1评分
  result$F1[i] <- (2*test_perf[2,2])/(2*test_perf[2,2]+test_perf[1,2]+test_perf[2,1])
  ##Mathews系数(MCC)，即Phi相关系数
  result$MCC[i] <- (test_perf[2,2]*test_perf[1,1]-test_perf[1,2]*test_perf[2,1])/sqrt((test_perf[2,2]+test_perf[1,2])*(test_perf[2,2]+test_perf[2,1])*(test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1]))
  ##基尼指数
  result$Gini[i] <- 2*result$Accuracy[i] - 1
  ##Kappa,通常0.75以上表示一致性较满意，0.4以下一致性不好
  result$Kappa[i] <- ((test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1])+(test_perf[2,1]+test_perf[2,2])*(test_perf[1,2]+test_perf[2,2]))/(sum(test_perf)*sum(test_perf))
  result$AUC[i] <- round(log_roc$auc,3)
  
  print(i)
}
