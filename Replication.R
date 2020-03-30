library(rbounds)
#Loads library containing all functions and datasets necessary

d <- dw_data_1_
#assigns the Lalonde Observational data
X <- cbind(d$age, d$education, d$black, d$hispanic, d$married, d$nodegree, d$re74, d$re75)
Y <- d$re78
treat <- d$treat

genout <- GenMatchRose(Tr=treat, X=X, Y=Y)
genout2 <- GenMatch(Tr = treat, X = X)
mout <- Match(Y=Y, Tr=treat, X=X, estimand="ATT", Weight.matrix=genout)
mout2 <- Match(Y=Y, Tr=treat, X=X, estimand = "ATT", Weight.matrix = genout2)
summary(mout)
summary(mout2)
sens <- psens(mout, Gamma = 5, GammaInc = 0.25)
sens2 <- psens(mout2, Gamma = 5, GammaInc = 0.25)
sens$bounds
sens2$bounds

data(Lalonde)
#assigns the Lalonde RCT data
X <- cbind(lalonde$age, lalonde$education, lalonde$black, lalonde$hispanic, lalonde$married, lalonde$nodegree, lalonde$re74, lalonde$re75)
Y <- lalonde$re78
treat <- lalonde$treat

genout <- GenMatchRose(Tr=treat, X=X, Y=Y)
genout2 <- GenMatch(Tr = treat, X = X)
mout <- Match(Y=Y, Tr=treat, X=X, estimand="ATT", Weight.matrix=genout)
mout2 <- Match(Y=Y, Tr=treat, X=X, estimand = "ATT", Weight.matrix = genout2)
summary(mout)
summary(mout2)
sens <- psens(mout, Gamma = 5, GammaInc = 0.25)
sens2 <- psens(mout2, Gamma = 5, GammaInc = 0.25)
sens$bounds
sens2$bounds

gencon <- genConfoundSearch(X = X, Y = Y, Tr = treat, Gamma = 3, search.pop.size = 25, match.pop.size = 25)
gencon$correlations
gencon$fitness
gencon$kappa