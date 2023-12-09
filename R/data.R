#' Microarray and RNA-Seq breast cancer datasets from The Cancer Genome Atlas (TCGA) for classification of progesterone receptor status
#'
#' @description Two preprocessed and curated TCGA breast cancer gene expression datasets
#' generated with two different technologies (microarray and RNA-sequencing)
#' for prediction of progesterone receptor (PR) status. Only data from primary tumors are used.
#' We retain only white female patients and create two non-overlapping and balanced datasets.
#' The corresponding protein-protein interaction (PPI) network information from the STRING database is also included.
#'
#' @format A list with three elements corresponding to PPI network, microarray dataset (with phenotype and genotype),
#' and RNA-Seq dataset (with phenotype and genotype). The PPI network is saved as an igraph object with 14167 nodes.
#' The gene expression datasets are saved as lists with two elements, one character vector giving the PR status
#' and the other numeric matrix with rows corresponding to samples and columns corresponding to genes.
#'
#' @source Bioconductor package curatedTCGAData and the STRING database
#'
#' @references Marcel Ramos, Ludwig Geistlinger, Sehyun Oh, Lucas Schiffer, Rimsha Azhar, Hanish Kodali, Ino de Bruijn, Jianjiong Gao, Vincent J Carey, Martin Morgan, et al. Multiomic integration of public oncology databases in bioconductor. JCO Clinical Cancer Informatics, 1:958–971, 2020.
#'
#' Damian Szklarczyk, Rebecca Kirsch, Mikaela Koutrouli, Katerina Nastou, Farrokh Mehryary, Radja Hachilif, Annika L Gable, Tao Fang, Nadezhda T Doncheva, Sampo Pyysalo, et al. The STRING database in 2023: protein–protein association networks and functional enrichment analyses for any sequenced genome of interest. Nucleic Acids Research, 51(D1):D638–D646, 2023.
#'
#' @examples
#' dim(tcga_breast_pr$microarray$geno)
#' length(tcga_breast_pr$microarray$pheno)
"tcga_breast_pr"
