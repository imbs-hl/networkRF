#' Function for constructing sampling probability of predictors used in the RF construction
#'
#' @param method Method for constructing sampling probability of predictors, currently supported methods are
#' * \"uniform\": all predictors have equal probability to be sampled which leads to the standard RF, default option
#' * \"marginal-p\": sampling probability derived from p-values of marginal two-sample t-test per gene
#' * \"topology\": sampling probability derived from network structure to only reflect the network topology
#' * \"network-p\": sampling probability derived from both network structure and marginal p-values
#' * \"network-q\": sampling probability derived from both network structure and FDR-adjusted marginal p-values
#' @param x Matrix of predictors
#' @param y Binary response vector
#' @param network Network to reflect the interactions among all predictors
#' @param gamma The restart parameter in the directed random walk algorithm, the default value is 0.3
#' @returns A vectror of sampling probability of all predictors used in the RF construction at each node
#'
#' @import qvalue
#' @examples
#'\dontrun{
#' sim.data <- gen_data(num.sample = 500, num.var = 1000, scenario = 2, AveEffect = 2)
#' Wt.marginal.p <- get_weight(method = "marginal-p",
#'                  x = sim.data$rna[[1]]$x_train,
#'                  y = sim.data$rna[[1]]$y_train)
#' Wt.topology <- get_weight(method = "topology", network = sim.data$network)
#'}
#' @export
get_weight <- function(method = "uniform", x, y, network, gamma = 0.3){
  if(is.null(method)) stop("Need to supply the method to construct sampling probability")

  `%notin%` <- Negate(`%in%`)
  if(method %notin% c("uniform", "marginal-p", "topology", "network-p", "network-q"))
    stop("The specified method is currently not supported!")

  if(method == "uniform"){
    Wt <- rep(1, ncol(x)) / ncol(x)
  }else if(method == "topology"){
    if(is.null(network)) stop("Need to supply network for topology method!")
    A <- get_adjacency_matrix(network)
    w.init <- rep(1, ncol(A))
    Wt <- DRW_network(network = network, p0 = w.init, gamma = 0)
  }else{
    x <- scale(x, center = TRUE, scale = TRUE)
    y <- as.factor(y)

    data.all <- data.frame(x, y, check.names = TRUE)
    colnames(data.all)[ncol(data.all)] <- "group"
    gene.ZP <- matrix(NA, nrow = ncol(x), ncol = 2)
    rownames(gene.ZP) <- colnames(x)
    colnames(gene.ZP) <- c("ZScore", "p-value")
    for (i in 1:ncol(x)) {
      res.marginal <- glm(as.formula(paste("group~",
                                           colnames(data.all)[i])), data.all, family = "binomial")
      gene.ZP[i, ] <- summary(res.marginal)$coefficients[2, c(3, 4)]
    }

    if(method == "network-q"){
      marginal.q <- qvalue::qvalue(gene.ZP[, 2])
      gene.weight <- 1 / (marginal.q$qvalues) - 1 # reference: enriched RF (bioinfor. 2008)
      w.init <- gene.weight / sum(gene.weight)
      Wt <- DRW_network(network = network, p0 = w.init, gamma = gamma)
    }else{
      gene.weight <- -log(gene.ZP[, 2] + 2.2e-16)
      gene.weight[which(is.na(gene.weight))] <- 0
      gene.weight <- (gene.weight - min(gene.weight))/(max(gene.weight) - min(gene.weight))
      gene.weight <- gene.weight / sum(gene.weight)
      if(method == "marginal-p"){
        Wt <- gene.weight
      }else{
        Wt <- DRW_network(network = network, p0 = gene.weight, gamma = gamma)
      }
    }
  }
  names(Wt) <- colnames(x)
  Wt
}
