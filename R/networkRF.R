#' Construction of network-guided random forests
#'
#' @param x Matrix of predictors
#' @param y Binary response vector
#' @param network Network to reflect the interactions among all predictors
#' @param method Method for constructing sampling probability of predictors, currently supported methods are
#' * \"uniform\": all predictors have equal probability to be sampled which leads to the standard RF, default option
#' * \"marginal-p\": sampling probability derived from p-values of marginal two-sample t-test per gene
#' * \"topology\": sampling probability derived from network structure to only reflect the network topology
#' * \"network-p\": sampling probability derived from both network structure and marginal p-values
#' * \"network-q\": sampling probability derived from both network structure and FDR-adjusted marginal p-values
#' * NULL: if no method is supplied, then a sampling probabibility must be supplied in the argument Wt
#' @param Wt Sampling probability of predictors used in the RF construction at each node
#' @param min.num.gene The number of genes retained at the end of the variable selection procedure
#' @param num.trees Number of trees in the RF construction
#' @param gamma The restart parameter in the directed random walk algorithm, the default value is 0.3
#' @param ... Any other arguments allowed in the ranger function for RF construction
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
#' res.one.rep <- networkRF(x = sim.data$rna[[1]]$x_train,
#'                          y = sim.data$rna[[1]]$y_train,
#'                          network = sim.data$network,
#'                          method = "topology",
#'                          min.num.gene = 20,
#'                          num.trees = 2000)
#'}
#' @export
networkRF <- function(x, y, network,
                      method = "uniform", Wt = NULL,
                      min.num.gene = 10, num.trees = 1000,
                      gamma = 0.3, ...){

  if(is.null(method)){
    if(is.null(Wt)) stop("method and Wt cannot be empty at the same time!")
  }else{
    if(any(is.na(Wt)) | any(Wt<0)) stop("Wt supplied should not have missing value and should all be non-negative!")
  }
  # construct the sampling probability
  if(!is.null(method)){
    Wt <- get_weight(method = method, x, y, network, gamma = gamma)
  }else{
    Wt <- Wt
  }
  # construct RF and conduct variable selection
  res <- networkRF_binary(x, y, Wt, min.num.gene, num.trees, ...)
  res

}
