# networkRF: network-guided random forest
Jianchang Hu

## Introduction
networkRF implements network-guided random forests (RFs) and a couple of other variants of random forests based on ranger. The sampling probability of predictors in the decision tree construction is not necessarily uniform, and can be constrcuted based on external network information and other source of information about the predictors. The package also provides functions to conduct simulation studies with synthetic RNA-Seq gene expression data and underlying gene network to demonstrate the effects of incorporating network information into RF construction for disease gene discovery. Two TCGA breast cancer datasets are also included in the package, the one from microarray and the other from RNA-Seq technology. The phenotype is PR status and all expression data has been log-transformed and standardized to have mean 0 and standard deviation 1.

## Installation
To install the networkRF R package, run
```
install.packages("devtools")
library(devtools)
devtools::install_github("imbs-hl/networkRF")
```

## Examples
The package includes functions to generate synthetic RNA-Seq data along with underlying gene network. The following codes provide a simple example to generate synthetic datasets and conduct analysis for disease gene identification and predictive model construction. More details on the functions and their arguments, please refer to the help pages of the package.
```
library(networkRF)

# scenario setting for synthetic datasets
seed = 1111
num.sample = 200
num.var = 1000
scenario = 1
AveEffect = 0.5
method = "uniform"

# simulate synthetic datasets
sim.data <- gen_data(num.sample = num.sample,
                     num.var = num.var,
                     scenario = scenario, 
                     AveEffect = AveEffect)

# analyze one data set
res.one.rep <- networkRF(x = sim.data$rna[[1]]$x_train,
                         y = sim.data$rna[[1]]$y_train,
                         network = sim.data$network,
                         method = method,
                         min.num.gene = 20,
                         num.trees = 2000, 
                         # importance = "permutation",
                         seed = seed)
```

If you find any bugs, or if you experience any crashes, please report to us.
