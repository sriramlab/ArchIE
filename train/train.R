library(data.table)
train <- data.frame(fread("../simulations/training_data.txt")) ## read in the training data

lr <- train[,1:209] # cut down the extra columns -- 209 is the label (0=not archaic, 1=archaic)
model <- glm(V209 ~ .,family=binomial(link='logit'),data=lr) # train the model
save.image("trained_model.Rdata") # save the trained model so we don't have to train it again. can load it with load("trained_model.Rdata")
test <- read.table("../simulations/test_data.txt") # read in test data
predicted <- plogis(predict.glm(model, test[,1:209])) # predict!
