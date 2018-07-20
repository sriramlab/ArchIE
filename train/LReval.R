library(data.table)
setwd("~/Documents/archaic_net/neanderthal/two_percent_explore/")
train <- data.frame(fread("./sim_data/train_10K_clean_16288.txt"))
SS.java <- read.table("../predictions/ss_prob_class_8338.txt.pr")

lr <- train[,1:209]
model <- glm(V209 ~ .,family=binomial(link='logit'),data=lr)

call_nn <- function(df, threshold){
  arch <- subset(df, df$V2 == 1)
  na <- subset(df, df$V2 == 0)
  tp <- nrow(subset(arch, arch$V1 >= threshold))#/nrow(df)
  fp <- nrow(subset(na, na$V1 >= threshold))#/nrow(df)
  tn <- nrow(subset(na, na$V1 <= threshold))#/nrow(df)
  fn <- nrow(subset(arch, arch$V1 <= threshold))#/nrow(df)
  return(c(tp, fp, tn, fn))
}
f <- function(dat, label){
  predicted <- plogis(predict.glm(model, dat))
  pred <- data.frame(predicted, test$V209)
  names(pred) <- c("V1", "V2")
  thresh <- seq(0,1,0.001)
  out <- sapply(thresh, call_nn, df=pred)
  prec <- out[1,]/(out[1,]+out[2,])
  recall <-  out[1,]/(out[1,]+out[4,])
  pdf(paste0(label,".pdf"))
  par(cex=1.5)
  plot(recall, prec, pch=4, ylim=c(0,1), xlim=c(0,1), xlab="Recall", ylab="Precision")
  points(SS.java$V1, SS.java$V2, pch=5, col="red")
  dev.off()
}


test <- data.frame(fread("../sim_data/test_clean_8338.txt"))
f(test, "LR_8338")
test <- data.frame(fread("../sim_data/test_clean_15161.txt"))
f(test, "LR_15161")
test <- data.frame(fread("../sim_data/test_clean_17091.txt"))
f(test, "LR_17091")
test <- data.frame(fread("../sim_data/test_clean_26769.txt"))
f(test, "LR_26769")

### get the feature importance
# pdf("LR_weights.pdf")
# barplot(model$coefficients[2:209], ylab="Weights")
# coefs.norm <-  model$coefficients[c(203,204,205,206,207,208,209)]
# barplot(coefs.norm, ylab="weights", las=2,names=c('mean_dist','var',"skew", "kurtosis","min_d","S*","n_priv"))
# dev.off()

pdf("weights-IFS.pdf")
barplot(model$coefficients[2:101], ylab="Weights", names.arg = 1:100, main="Individual frequency spectrum", xlab = "Frequency")
dev.off()
pdf("weights-d_vec.pdf")
barplot(model$coefficients[103:202], ylab="Weights", names.arg = 1:100, main="Distance vector", xlab = "Index")
dev.off()
pdf("weights-stats.pdf")
coefs.norm <-  model$coefficients[c(203,204,205,206,207,208,209)]
barplot(coefs.norm, ylab="Weights", las=2, xlab="Features",xaxt="n")
labs<-c('mean_dist','var',"skew", "kurtosis","Min dist","S*","Private SNPs")
text(cex=1, x=bp-0, y=-1.5, labs, xpd=TRUE, srt=20)
dev.off()

predicted <- plogis(predict.glm(model, test))
pred <- data.frame(predicted, test$V209)
names(pred) <- c("V1", "V2")
admix <- subset(pred, pred$V2 == 1)
none <- subset(pred, pred$V2 == 0)
none2 <- none[sample(nrow(none), nrow(admix)), ]
plot(density(none2$V1), col="red", lwd=2, xlim=c(-0.1,1.1), ylim=c(0,2))
lines(density(admix$V1), col="blue", lwd=2)

pdf("LR_hist.pdf")
par(cex=1.5)
hist(none2$V1, xlab="P(admix)", col=rgb(0.75,0.75,0, 0.4), main="", border=NA)
hist(admix$V1, add=T, col=rgb(0,0.4,1, 0.4), border=NA)
legend("topright", c("Archaic","Not archaic"), fill = c(rgb(0,0.4,1, 0.4), rgb(0.75,0.75,0, 0.4)), border=NA, bty="n")
dev.off()


####
predicted <- plogis(predict.glm(model, data.frame(fread("../sim_data/test_clean_26769.txt"))))
pred <- data.frame(predicted, test$V209)
names(pred) <- c("V1", "V2")
thresh <- seq(0.8,1,0.001)
out <- sapply(thresh, call_nn, df=pred)
prec <- out[1,]/(out[1,]+out[2,])
recall <-  out[1,]/(out[1,]+out[4,])
cbind(prec, recall, thresh)


### save predictions for AUC.jar
predicted <- plogis(predict.glm(model, data.frame(fread("../sim_data/test_clean_26769.txt"))))
pred <- data.frame(predicted, test$V209)
names(pred) <- c("V1", "V2")
write.table(pred, "LR-predictions-26769.txt", row.names=F, col.names = F, quote=F, sep=" ")

test.8338 <- data.frame(fread("../sim_data/test_clean_8338.txt"))
ss <- data.frame(test.8338$V207, test.8338$V209)
names(ss) <- c("ss", "true")
ss_prob <- ecdf(ss$ss)(ss$ss)
ss_pred <- data.frame(ss_prob, ss$true)
write.table(ss_pred, "SS-predictions-8338.txt", row.names = F, col.names = F)

#### once ROC stuff is done
# java -jar ../1000G_crflabel/auc.jar SS-predictions-8338.txt list
# java -jar ../1000G_crflabel/auc.jar LR-predictions-26769.txt list

ss.8338.roc <- read.table("SS-predictions-8338.txt.roc")
lr.26769.roc <- read.table("LR-predictions-26769.txt.roc")

pdf("ROC.pdf")
par(cex=1.5)
plot(ss.8338.roc$V1, ss.8338.roc$V2, col="red", ylab="True Positive Rate", xlab="False Positive Rate", pch=5)
points(lr.26769.roc$V1, lr.26769.roc$V2, pch=4)
dev.off()
