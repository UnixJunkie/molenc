# Train a Multiple Linear Regression model given a dataset (input CSV file)
#
# background on MLR: https://www.investopedia.com/terms/m/mlr.asp
# R example: https://www.tutorialspoint.com/r/r_multiple_regression.htm

# FBR: TODO
# center and scale all dependant variables, store the scaling parameters
# with the model

# read data in
train <- read.csv("data/gold_3clpro_std_forR.csv", header = T, sep = ",")
# only keep interesting columns
train <- train[,c("gold.3clpro.score","MolW","cLogP","RotB")]

# train model
model <- lm(gold.3clpro.score ~ MolW + cLogP + RotB, data = train)

# Show the model.
print(model)

a <- coef(model)[1]
b <- coef(model)[2]
c <- coef(model)[3]
d <- coef(model)[4]

print(c(a, b, c, d))
