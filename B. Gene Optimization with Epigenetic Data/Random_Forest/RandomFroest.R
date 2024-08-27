#https://www.bilibili.com/read/cv16304221

library(randomForest)
library(caret)
library(pROC)
set.seed(1234)
colnames(TSS)<-c("6mA","H3K4me3","H2A.Z","nucleosome","type")

TSS<-read.csv("Z:/UTR_annotation/机器学习模型/随机森林/model2/TSS-input.txt",sep = "\t",header = 1, row.names = 1)
TES<-read.csv("Z:/UTR_annotation/机器学习模型/随机森林/model2/TES-input.txt",sep = "\t",header = 1, row.names = 1)
TSS_TES<-rbind(TSS,TES)
trains <- createDataPartition(y = TSS_TES$type,p = 0.75,list = F)
traindata <- TSS_TES[trains,]
testdata <- TSS_TES[-trains,]
rf.train <- randomForest(as.factor(type)~.,data = traindata,na.action = na.roughfix,importance = TRUE, ntree=500, mtry=1)

plot(rf.train,main="ERROR&TREES")

importance(rf.train)

trainpredprob<-predict(rf.train, newdata = traindata, type="prob")

trainroc<-roc(response=traindata$type, predictor=trainpredprob[,2])

plot(trainroc, print.auc=TRUE, auc.polygon=TRUE, grid=T, max.auc.polygon=T, auc.polygon.col="skyblue", print.thres=T, legacy.axes=T, bty="l")

bestp <- trainroc$threshold[which.max(trainroc$sensitivities + trainroc$specificities -1)]

trainpredlab <- as.factor(ifelse(trainpredprob[,2] > bestp ,1,0))

confusionMatrix(data = as.factor(trainpredlab),  reference = as.factor(traindata$type), positive = '1', mode = 'everything')

testpredprob <- predict(rf.train, newdata = testdata,type = 'prob')

testpredlab <- as.factor(ifelse(testpredprob[,2] > bestp ,1,0))

confusionMatrix(data = as.factor(testpredlab),  reference = as.factor(testdata$type),  positive = '1', mode = 'everything') 

testroc <- roc(response = testdata$type,predictor = testpredprob[,2]) 

plot(trainroc,print.auc = T,grid = c(0.1,0.2),auc.polygon = F,max.auc.polygon = T,main = "random Forest ROC", grid.col = c("green","red"))
plot(testroc,print.auc = T, print.auc.y = 0.4,add = T, col = 'red')
legend("bottomright",legend = c("traindata","testdata"), col = c(par('fg'),'red'), lwd = 2, cex = 0.9)





rf.test <- predict(rf.train,traindata,type="response")
trn_pred<-ifelse(predict(rf.train, type = "response")> 0.3, 1, 0)
trn_tab <- table(predicted = trn_pred, actual = traindata$type)
tst_pred <- ifelse(predict(rf.train, newdata = testdata, type = "response") > 0.3, 1, 0)
tst_tab <- table(predicted = tst_pred, actual = testdata$type)
calc_class_err(actual = teat$type, predicted = tst_pred)
confusionMatrix(trn_tab, positive = "1")
test_roc <- roc(testdata$type , test_prob, plot = TRUE, print.auc = TRUE)


