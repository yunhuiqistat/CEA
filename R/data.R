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

