#' Function for prediction on testing data and get misclassification rate
#'
#' @param xtest Predictors of testing dataset
#' @param ytest Binary response vector of testing dataset
#' @param rf.model A random forest predictive model
#' @returns err A scalar of misclassification error of the rf.model on given testing dataset
#' @examples
#' sim.data <- gen_data(num.sample = 500, num.var = 1000, scenario = 2, AveEffect = 2)
#' res.one.rep <- networkRF(x = sim.data$rna[[1]]$x_train,
#'                          y = sim.data$rna[[1]]$y_train,
#'                          network = sim.data$network,
#'                          method = "topology",
#'                          min.num.gene = 20,
#'                          num.trees = 2000)
#' pred.accuracy <- pred_accuracy(x = sim.data$rna[[1]]$x_test,
#'                                y = sim.data$rna[[1]]$y_test,
#'                                rf.model = res.one.rep$networkRF_model)
#'
#' @export
pred_accuracy <- function(xtest, ytest, rf.model){
  xtest <- scale(xtest, center = TRUE, scale = TRUE)
  xtest <- data.frame(xtest)

  ypred <- predict(rf.model, xtest)$predictions
  ypred <- as.numeric(levels(ypred))[ypred]
  err <- mean(abs(ypred - ytest))
  err
}
