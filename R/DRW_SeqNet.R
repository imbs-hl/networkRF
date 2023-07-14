#' Function for directed random walk algorithm
#'
#' @param network Network to reflect the interactions among all predictors
#' @param p0 The initial probability distribution over all nodes within the network
#' @param gamma The restart parameter in the directed random walk algorithm, the default value is 0
#' @param tol The tolerance parameter to reflect the precision of the estimated equilibrium probability distribution
#' @returns A vectror of estimated equilibrium probability distribution over all nodes in the network
#'
#' @import igraph
#' @examples
#'\dontrun{
#' sim.data <- gen_data(num.sample = 500, num.var = 1000, scenario = 2, AveEffect = 2)
#' prob.topology <- DRW_network(network = sim.data$network)
#'}
#' @export
DRW_network <- function (network, p0 = NULL, gamma = 0, tol = 1e-10){
  if(class(network) == "igraph"){
    W <- igraph::as_adjacency_matrix(network, type = "both", sparse = FALSE)
  }else{
    W <- get_adjacency_matrix(network)
  }
  W.col.normed <- apply(W, 2, function(x) x / max(1, sum(x)))

  if(is.null(p0)){
    p0 <- rep(1, ncol(W))
  }
  p0 <- t(as.matrix(p0/sum(p0))) # self-normalization to a prob vector

  PT <- p0
  k <- 0
  delta <- 1
  while (delta > tol) {
    PT1 <- (1 - gamma) * W.col.normed
    PT2 <- PT1 %*% t(PT)
    PT3 <- (gamma * p0)
    PT4 <- t(PT2) + PT3
    delta <- sum(abs(PT4 - PT)) # measured in l_1 norm
    PT <- PT4
    k <- k + 1
  }
  PT <- t(PT)
  rownames(PT) <- NULL
  res <- drop(PT[1:dim(PT)[1]])

  res
}
