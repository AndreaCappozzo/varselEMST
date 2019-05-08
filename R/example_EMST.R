# library(MASS)
# library(varselEMST)
# library(mclust)
# library(clustvarsel)
#
# set.seed(73)
# n <- 200
# pro <- 0.5
# mu1 <- c(0,0)
# mu2 <- c(3,3)
# sigma1 <- matrix(c(1,0.5,0.5,1),2,2,byrow=TRUE)
# sigma2 <- matrix(c(1.5,-0.7,-0.7,1.5),2,2,byrow=TRUE)
# XX <- matrix(0, n, 5)
# colnames(XX) <- paste("XX", 1:ncol(XX), sep ="")
# # generate the grouping variables
# u <- runif(n)
# Class <- ifelse(u < pro, 1, 2)
# XX[u < pro, 1:2]  <- mvrnorm(sum(u < pro), mu = mu1, Sigma = sigma1)
# XX[u >= pro, 1:2] <- mvrnorm(sum(u >= pro), mu = mu2, Sigma = sigma2)
# # generate the non-grouping variables
# XX[,3] <- XX[,1] + rnorm(n)
# XX[,4] <- rnorm(n, mean = 1.5, sd = 2)
# XX[,5] <- rnorm(n, mean = 2, sd = 1)
#
# # I add 5% outliers
#
# # XX <- rbind(XX, mvrnorm(10, mu = rep(8,5), Sigma = diag(1,5)))
#
# # plot the data
# clPairs(XX, Class, gap = 0)
# # clPairs(XX, c(Class,rep(3,10)), gap = 0)
#
# X_train = XX
# class_train <- Class
# n_relevant_variables <- 2
# alpha_train <- .05
#
# clustvarsel(X_train, G = 2, emModels2 = "EII", direction = "forward")
#
# # model_name <- "VVV"
# # ctrl_init <- raedda::control_init()
#
# test_0 <- varselEMST::var_sel_EMST_raedda_l(
#   X_train = X_train,
#   class_train = Class,
#   n_relevant_variables = 2,
#   alpha_train = .1,
#   model_names = "EVE",
#   ctrl_init = control_init_EMST(n_samp = 1),
#   swap_step = "exhaustive",
#   ctrl_GA = control_GA(ngen = 10),
#   verbose = F
# )
#
# test_0$Best$relevant_variables
#
# test_all <-
#   lapply(mclust.options("emModelNames"), raedda_l_EMST_model, X_train = X_train,
#          class_train = Class,
#          n_relevant_variables = 2,
#          alpha_train = .05,
#          ctrl_init = control_init_EMST(),
#          swap_step = "ga",
#          ctrl_GA = control_GA(ngen = 10))
# names(test_all) <- mclust.options("emModelNames")
#
