#' @title COVID-19 multi-omics data
#'
#' @description A dataset including metabolome and proteome data of 127 patients in Albany Medical Center, NY from 6 April 2020 to 1 May 2020
#'
#'
#' @format A list containing element metabolite (metabolite data matrix with 127 samples, 111 variables), protein (protein data matrix with 127 samples, 516 variables) and covariates (meta information).
#'
#' @references Overmyer K A, Shishkova E, Miller I J, et al. Large-scale multi-omic analysis of COVID-19 severity[J]. Cell systems, 2021, 12(1): 23-40. e7.
#'
#' @keywords covid
#'

"covid"


#' @title Barley microarray data
#'
#' @description This dataset contains gene expression measurements for barley (Hordeum vulgare) in response to the powdery mildew pathogen, as described in Caldo et al. (2004). The study investigates the interaction-dependent gene expression in the Mla-specified response to barley powdery mildew.
#'
#' @format A list containing raw_data, raw_meta, data (microarray data matrix with 104 samples, 22840 genes) and covariates (information of 104 samples).
#'
#' @references Caldo R A, Nettleton D, Wise R P. Interaction-dependent gene expression in Mla-specified response to barley powdery mildew[J]. The Plant Cell, 2004, 16(9): 2514-2528.
#'
#' @keywords barley
#'

"barley"

#' @title Read me data examples
#'
#' @description Read me data examples including simulated data for cPCA, 3CA and real barley and covid data for scPCA and s3CA.
#'
#'
#' @format A list containing four sublists, each including data and covariates.
#'
#' @references Refer references for data(barley) and data(covid)
#'
#' @keywords intro
#'

"intro"

#' @title Contrastive spectral clustering simulated data
#'
#' @description Contrastive spectral clustering simulated data using n11 = 50, n12 = 25, n21 = 55, n22 = 30, m1tilde = 20, m2tilde = 30,
#'                      mu1 = 5, mu2 = 2, theta1 = -5, theta2 = 6, p1 = 5, p2 = 10, Sigma_t = Identity, Sigma_a = Identity, and set.seed(111)
#'
#' @format A list containing target sample Xt, ancillary sample Xa, target covariates covariates_t and ancillary covariates covariates_a.
#'
#'
#' @keywords sim_data_CSC

"sim_data_CSC"


#' @title Acute myelogenous leukemia single cell sequencing data
#'
#' @description scRNAseq for two healthy individuals (healthy 1 and 2), and two AML patients (AML027, 035) with top 1000 most variable genes. Pre-processing and gene selection is performed following Boileau 2020.
#'
#' @format A list containing four single cell experiments object (sce)
#'
#' @references Boileau P, Hejazi N S, Dudoit S. Exploring high-dimensional biological data with sparse contrastive principal component analysis[J]. Bioinformatics, 2020, 36(11): 3422-3430. Zheng G X Y, Terry J M, Belgrader P, et al. Massively parallel digital transcriptional profiling of single cells[J]. Nature communications, 2017, 8(1): 14049.
#' @keywords leukemia

"leukemia"

