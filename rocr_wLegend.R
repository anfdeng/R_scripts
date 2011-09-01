# ROC curve(s) with ROCR
# Chupakhin Vladimir (chupvl@gmail.com)
# loading ROCR library
library("ROCR")
# loading active compounds, or compounds with label1
active <- read.table("sample.active", sep=",", header=FALSE)
# loading inactive compounds, or compounds with label2
inactive <- read.table("sample.inactive", sep=",", header=FALSE)
# binding them and converting to matrix because ROCR works with matrix data
target_pred <- as.matrix(rbind(active,inactive))
# because number of the colums should be the same - making additional param
ncol <- ncol(inactive)
# generating classes (1 for active, 0 for inactive, but it can be -1 and 1 - there is no difference)
class.active <- matrix(sample(1, (ncol(active)*nrow(active)), replace=T), ncol=ncol)
class.inactive <- matrix(sample(0, (ncol(inactive)*nrow(inactive)), replace=T), ncol=ncol)
# binding the classes
target_class <- rbind(class.active,class.inactive)
#target_class1 <- target_class[,1]
# calculating the values for ROC curve
pred <- prediction(target_pred, target_class)
perf <- performance(pred,"tpr","fpr")
# changing params for the ROC plot - width, etc
par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=1.3,cex.lab=1.4)
# plotting the ROC curve
plot(perf,col="black",lty=3, lwd=3)
# calculating AUC
auc <- performance(pred,"auc")
# now converting S4 class to vector
auc <- unlist(slot(auc, "y.values"))
# adding min and max ROC AUC to the center of the plot
minauc<-min(round(auc, digits = 2))
maxauc<-max(round(auc, digits = 2))
minauct <- paste(c("min(AUC)  = "),minauc,sep="")
maxauct <- paste(c("max(AUC) = "),maxauc,sep="")
legend(0.3,0.6,c(minauct,maxauct,"\n"),border="white",cex=1.7,box.col = "white")
#