#' Generate synthetic gene network and RNA-Seq data for a given scenario
#'
#' @param num.sample Number of samples of training set. The testing set has the same sample size.
#' @param num.var Number of genes simulated.
#' @param scenario A numeric ID to specify the simulation scenario. 1 for Null case,
#' 2 for randomly distributed disease genes, 3 for one disease module without main disease gene,
#' 4 for one disease module with one main disease gene, 5 for two disease modules without main disease genes,
#' 6 for two disease modules with one main disease gene in each module
#' @param AveEffect The average effect size used in the logistic regression model to generate disease status
#' @param gamma The probability of re-start in the directed random walk algorithm
#' @param num.rep Number of repetitions
#' @param seed An integer number serving as the seed of the (pseudo-)random number generation
#' @returns A list of objects giving the synthetic datasets for a given scenario. The list consists of the following subjects:
#'
#' * scenario The index of scenario
#' * num_sample Number of samples in training set and testing set
#' * num_var Number of genes simulated
#' * network The synthetic underlying gene network
#' * rna A list of 100 replications of RNA-Seq data with each replication consisting of both training and testing datasets
#' * disease_gene A vector of indices of genes corresponding to the disease genes of the simulation scenario
#' * main_disease_gene A vector of indices of genes corresponding to the main disease gene(s) of the scenario. When all disease genes having the same effect size, then all genes will be main disease genes. Otherwise, the index(es) of main disease gene(s) from disease module(s) is returned.
#'
#' @examples
#'\dontrun{
#' sim.data <- gen_data(num.sample = 500, num.var = 1000, scenario = 2, AveEffect = 2)
#'}
#' @export
gen_data <- function(num.sample = 1000, num.var = 1000,
                     scenario = 2, AveEffect = 2,
                     gamma = 0.3, num.rep = 100, seed = 11223344)
{
  # load dependent library and function
  # library(SeqNet)

  # ----------------------------------------------------------------------------
  # this chunk generates fixed underlying network
  # and based on the scenario, select signal genes
  set.seed(seed)

  p <- num.var # number of nodes in the network
  network <- random_network(p)
  network <- gen_partial_correlations(network)
  module.length <- sapply(network$modules, function(module) length(module$nodes))
  mod.len.q1 <- floor(quantile(module.length, probs = 0.25))
  # reference for signal module size

  scenario <- scenario # scenario is a parameter -------------------------------

  isNull <- (scenario == 1)
  isModule <- (scenario > 2)
  isTwoMod <- (scenario >  4)
  isEqual <- (scenario %in% c(3, 5))

  effect.size <- AveEffect # average effect size is a parameter ----------------

  if(isNull){
    seed.node <- NA
    gene.sig <- NA
    beta <- 0
  }else if(!isModule){
    # randomly distributed disease genes
    num.sig.gene <- sample(module.length[module.length <= mod.len.q1], 1)
    gene.sig <- sample(seq(p), num.sig.gene, replace = FALSE)
    seed.node <- sort(gene.sig) # idx of disease genes

    beta <- rep(effect.size, length(gene.sig))
  }else if(!isTwoMod){
    mod.candidate <- which(module.length <= mod.len.q1)
    mod.signal <- sample(mod.candidate, 1)
    seed.node <- network$modules[[mod.signal]]$nodes
    gene.sig <- seed.node

    if(isEqual){
      beta <- rep(effect.size, length(gene.sig))
    }else{
      gene.sig <- network$modules[[mod.signal]]$nodes
      seed.node <- sample(gene.sig, 1)

      gene.similarity.mod <- rep(0, length(gene.sig))
      gene.similarity.mod[which(gene.sig == seed.node)] <- 1
      gene.similarity.mod <- DRW_network(network = network$modules[[mod.signal]],
                                         p0 = gene.similarity.mod, gamma = gamma)
      beta <- gene.similarity.mod * effect.size * length(gene.sig)
    }
  }else{
    mod.candidate <- which(module.length <= mod.len.q1)

    mod.signal.1st <- sample(mod.candidate, 1)
    seed.node.1st <- network$modules[[mod.signal.1st]]$nodes

    mod.intersect <- sapply(mod.candidate,
                            function(mod){length(intersect(seed.node.1st, network$modules[[mod]]$nodes))>0})
    mod.mutual <- mod.candidate[!mod.intersect]
    mod.signal.2nd <- sample(mod.mutual, 1)
    seed.node.2nd <- network$modules[[mod.signal.2nd]]$nodes

    gene.sig <- c(seed.node.1st, seed.node.2nd)

    if(isEqual){
      beta <- rep(effect.size, length(gene.sig))
      seed.node <- gene.sig
    }else{
      mod.sig.beta <- function(mod, size){
        gene.sig <- network$modules[[mod]]$nodes
        seed.node <- sample(gene.sig, 1)

        gene.similarity.mod <- rep(0, length(gene.sig))
        gene.similarity.mod[which(gene.sig == seed.node)] <- 1
        gene.similarity.mod <- DRW_network(network = network$modules[[mod]],
                                           p0 = gene.similarity.mod, gamma = gamma)
        beta <- gene.similarity.mod * size * length(gene.sig)
        list(seed.node = seed.node, similarity = gene.similarity.mod, beta = beta)
      }

      mod.1st <- mod.sig.beta(mod.signal.1st, effect.size)
      mod.2nd <- mod.sig.beta(mod.signal.2nd, effect.size)

      beta <- c(mod.1st$beta, mod.2nd$beta)
      seed.node <- c(mod.1st$seed.node, mod.2nd$seed.node)
    }
  }

  # ----------------------------------------------------------------------------
  # this chunk generates class label
  # and completes the data sets generation

  # total number of both training and testing samples
  # with each set having num.sample number of samples
  num.sample.total <- num.sample * 2
  # initialize the RNA-Seq data output
  rna.data <- vector(mode = "list", length = num.rep)
  for(iter in seq(num.rep)){
    set.seed(9876 * iter + 1234) # iter is a parameter ---------------------------

    # generate RNASeq data and take log-transformation
    x.total <- gen_rnaseq(num.sample.total, network)
    x.total <- log2(x.total$x + 1)
    x.total <- scale(x.total, center = TRUE, scale = TRUE)

    # generate label
    if(isNull){
      log.odds <- rep(0, nrow(x.total))
    }else{
      log.odds <- x.total[, gene.sig] %*% beta
    }
    prob <- as.vector(exp(log.odds) / (1 + exp(log.odds)))
    Y <- sapply(prob, function(p) rbinom(1, 1, p))

    # training and testing samples for the repetition
    trainSmpl.Idx <- sample(seq(num.sample.total), num.sample.total / 2, replace = FALSE)
    testSmpl.Idx <- setdiff(seq(num.sample.total), trainSmpl.Idx)
    x.trainSmpl <- x.total[trainSmpl.Idx, ]
    x.testSmpl <- x.total[testSmpl.Idx, ]
    Y.trainSmpl <- Y[trainSmpl.Idx]
    Y.testSmpl <- Y[testSmpl.Idx]

    rna.data[[iter]] <- list(x_train = x.trainSmpl,
                             y_train = Y.trainSmpl,
                             x_test = x.testSmpl,
                             y_test = Y.testSmpl)
  }
  list(scenario = scenario, num_sample = num.sample, num_var = num.var,
       network = network, rna = rna.data,
       disease_gene = gene.sig, main_disease_gene = seed.node)
}
