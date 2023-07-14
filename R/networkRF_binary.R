#' Construction of network-guided random forests for binary classification problem
#'
#' @param x Matrix of predictors
#' @param y Binary response vector
#' @param Wt Sampling probability of predictors used in the RF construction at each node
#' @param min.num.gene The number of genes retained at the end of the variable selection procedure
#' @param num.trees Number of trees in the RF construction
#' @param ... Any other arguments allowed in the ranger function
#' @returns A list of objects giving the results of variable selection based on network-guided RF. The list consists of the following subjects:
#'
#' * networkRF_model The final random forest predictive model constructed with top min.num.gene number of genes
#' * geneSelected A matrix recording the top selected genes at each step of the feature elimination procedure
#' * OOB.err A vector recording the OOB prediction error estimates for each RF predictive model constructed during the feature elimination procedure
#' * Wt Sampling probability of all predictors used in the RF construction at each node
#'
#' @examples
#'\dontrun{
#' sim.data <- gen_data(num.sample = 500, num.var = 1000, scenario = 2, AveEffect = 2)
#' res.one.rep <- networkRF_binary(x = sim.data$rna[[1]]$x_train,
#'                                 y = sim.data$rna[[1]]$y_train,
#'                                 Wt = rep(1, num.var),
#'                                 min.num.gene = 20,
#'                                 num.trees = 2000)
#'}
#' @export
networkRF_binary <- function(x, y, Wt, min.num.gene, num.trees, ...){
  # library(ranger)

  x <- scale(x, center = TRUE, scale = TRUE)
  y <- as.factor(y)
  data.all <- data.frame(x, y, check.names = TRUE)
  colnames(data.all)[ncol(data.all)] <- "group"
  x.names <- colnames(data.all)[1:(ncol(data.all)-1)]

  Wt.temp <- Wt
  names(Wt.temp) <- x.names

  # create a grid of number of genes selected from the networkRF
  ladder <- gen_ladder(Ntotal = ncol(x), ratio = 0.9, Nmin = min.num.gene)
  geneSelected <- matrix(0, nrow = ncol(x), ncol = length(ladder) - 1)
  rownames(geneSelected) <- x.names
  colnames(geneSelected) <- ladder[-1]
  genes <- x.names
  oob.err <- rep(NA, length(ladder))
  # an iterative procedure for gene selection
  # the procedure is similar to feature elimination
  # at each step, a network-guided RF is constructed with top genes from previous step
  # then top genes from the RF is retained and used to construct RF in the next step
  # OOB prediction error estimates are also recorded
  for (i in 1:length(ladder)) {
    x.temp <- data.all[, c(genes, "group")]
    x.grow <- ranger(group ~ ., x.temp,
                     num.trees = num.trees,
                     split.select.weights = Wt.temp[genes],
                     importance = "permutation", ...)
    vimp.i <- sort(x.grow$variable.importance, decreasing = TRUE)
    if (i < length(ladder)) {
      genes <- names(vimp.i)[1:ladder[i + 1]]
      geneSelected[genes, i] <- geneSelected[genes, i] + 1
    }
    oob.err[i] <- x.grow$prediction.err
  }
  rownames(geneSelected) <- colnames(x)

  x.temp <- data.all[, c(genes, "group")]
  x.grow <- ranger(group ~ ., x.temp, num.trees = num.trees, ...)

  list(networkRF_model = x.grow,
       OOB.err = oob.err,
       Wt = Wt,
       geneSelected = geneSelected)
}
