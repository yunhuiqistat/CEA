#' Simulation function for contrastive spectral clustering (CSC)
#'
#' This function performs simulation to prove the interesting of hidden clusters by using cPCA and eta tuning algorithm.
#' Same as in the example of Section 2, we consider two target clusters and two nuisance clusters
#'
#' @import PMA
#' @import MASS
#' @param n11 number of target samples in interesting cluster 1 and nuisance cluster 1.
#' @param n12 number of target samples in interesting cluster 1 and nuisance cluster 2.
#' @param n21 number of target samples in interesting cluster 2 and nuisance cluster 1.
#' @param n22 number of target samples in interesting cluster 2 and nuisance cluster 2.
#' @param m1tilde number of ancillary samples in nuisance cluster 1.
#' @param m2tilde number of ancillary samples in nuisance cluster 2.
#' @param mu1 non-zero values in mu_t1.
#' @param mu2 non-zero values in mu_t2.
#' @param theta1 non-zero values in theta_1.
#' @param theta2 non-zero values in theta_2
#' @param p1 number of non-zero elements in mu_t1 and mu_t2.
#' @param p2 number of non-zero elements in theta_1 and theta_2.
#' @param sigma_t standard deviation of the uncorrelated errors in target group.
#' @param sigma_a standard deviation of the uncorrelated errors in ancillary group.
#' @param Sigma_t the error covariance matrix in target group, if is NULL, use standard deviation specified and assume uncorrelation.
#' @param Sigma_a the error covariance matrix in ancillary group, if is NULL, use standard deviation specified and assume uncorrelation.
#' @param eta if eta is a vector, eta_tuning algorithm will be called to select optimal eta; if a number, eta will be used.
#' @param repetition number of replications of the simulation.
#' @param save_plot logical, whether save the score plot and loading plot for the first replication.
#' @export
#' @return A data frame containing the selected optimal eta, the four metrics values.
#' @examples
#' # library(CEA)
#' # library(ggplot2)
#' # library(tidyverse)
#' # library(ggsci)
#' # library(ggh4x)
#' # set.seed(111)
#' # mu2_candidates <- c(1, 2)
#' # p3_candidates <- c(300, 500)
#' # eta_range <-  seq(0.1, 50, 0.2)
#' # Sigmat_candidates <- c("Identity", "Toep(0.4)", "Toep(0.5)")
#' # results <- setNames(data.frame(matrix(ncol = 8, nrow = 0)),
#' #                     c("ARI_ctst_interesting", "ARI_ctst_nuisance",
#' #                      "ARI_trt_interesting", "ARI_trt_nuisance","eta",
#' #                       "mu2", "p1", "Sigmat"))
#' # for (mu2 in mu2_candidates) {
#' #   for (p3 in p3_candidates) {
#' #     for(Sigmat in Sigmat_candidates){
#' #       p1 <- 5
#' #       p2 <- 10
#' #       p <- p1+p2+p3
#' #       Sigma_t <- case_when(Sigmat == "Identity" ~ diag(1, p),
#' #                            Sigmat == "Toep(0.4)" ~ toeplitz(0.4^(0:(p-1))),
#' #                            Sigmat == "Toep(0.5)" ~ toeplitz(0.5^(0:(p-1))))
#' #       res <- sim_sCSC(n11 = 50, n12 = 25, n21 = 55, n22 = 30, m1tilde = 20, m2tilde = 30,
#' #                      mu1 = 5, mu2 = mu2, theta1 = -5, theta2 = 6, p1 = p1, p2 = p2, p3 = p3,
#' #                      sigma_t = 1, sigma_a = 1,
#' #                      Sigma_t = Sigma_t, Sigma_a = NULL,
#' #                      eta = eta_range, repetition = 10)
#' #       res$mu2 <- mu2
#' #       res$mu <- abs(mu2-5)
#' #       res$p3 <- p3
#' #       res$Sigmat <- Sigmat
#' #       results <- rbind(results, res)
#' #     }
#' #   }
#' # }
#' # plot_results <- pivot_longer(results, 1:4, names_to = "Type", values_to = "ARI")
#' # plot_results$Analysis <- sapply(plot_results$Type, function(x){strsplit(x,"_")[[1]][2]})
#' # plot_results$Cluster <- sapply(plot_results$Type, function(x){strsplit(x,"_")[[1]][3]})
#' # plot_results$Analysis <- ifelse(plot_results$Analysis=="ctst","Contrastive","Target-only")
#' # plot_results$Cluster <- ifelse(plot_results$Cluster=="interesting","Interesting","Nuisance")
#' # plot_results$eta <- paste("eta==", plot_results$eta, sep = "")
#' # plot_results$mu <- paste("mu==", plot_results$mu, sep = "")
#' # plot_results$p3 <- paste("p[3]==", plot_results$p1, sep = "")
#' # plot_results <- plot_results %>% group_by(mu, p3, Sigmat, Cluster, Analysis) %>%  summarise(Average_ARI = mean(ARI), SE_ARI = sd(ARI))
#' # print(knitr::kable(plot_results, format = "markdown"))
sim_sCSC <- function(n11, n12, n21, n22, m1tilde, m2tilde,
                    mu1, mu2, theta1, theta2, p1, p2, p3,
                    sigma_t = 1, sigma_a = 1,
                    Sigma_t = NULL, Sigma_a = NULL,
                    eta, repetition = 10){

  if(length(eta) > 1){
    eta_range <- eta
  }
  else{eta_use <- eta}

  ARI_trt_nuisance <- c()
  ARI_trt_interesting <- c()
  ARI_ctst_nuisance <- c()
  ARI_ctst_interesting <- c()
  eta_opt_set <- c()

  p <- p1+p2+p3

  mu_t1 <- matrix(rep(c(mu1,0,0), c(p1,p2,p3)), ncol=1)
  mu_t2 <- matrix(rep(c(mu2,0,0), c(p1,p2,p3)), ncol=1)
  theta_1 <- matrix(rep(c(0,theta1,0), c(p1,p2,p3)), ncol=1)
  theta_2 <- matrix(rep(c(0,theta2,0), c(p1,p2,p3)), ncol=1)

  if(is.null(Sigma_t)){
    Sigma_t <- diag(sigma_t, p)
  }
  if(is.null(Sigma_a)){
    Sigma_a <- diag(sigma_a, p)
  }

  covariates_t <- data.frame(dataset = rep("target", n11+n12+n21+n22),
                             interesting = rep(c("int=1","int=2"), c(n11+n12, n21+n22)),
                             nuisance = rep(c("nui=1","nui=2","nui=1","nui=2"), c(n11, n12, n21, n22)))
  covariates_a <- data.frame(dataset = rep("ancillary", m1tilde+m2tilde),
                             nuisance = rep(c("nui=1","nui=2"), c(m1tilde, m2tilde)))
  for(repe in 1:repetition){
    # generate data
    Xt <- rbind(matrix(mu_t1, nrow = n11, ncol = p, byrow = TRUE) + matrix(theta_1, nrow = n11, ncol = p, byrow = TRUE) + mvrnorm(n = n11, mu = rep(0, p), Sigma = Sigma_t),
                matrix(mu_t1, nrow = n12, ncol = p, byrow = TRUE) + matrix(theta_2, nrow = n12, ncol = p, byrow = TRUE) + mvrnorm(n = n12, mu = rep(0, p), Sigma = Sigma_t),
                matrix(mu_t2, nrow = n21, ncol = p, byrow = TRUE) + matrix(theta_1, nrow = n21, ncol = p, byrow = TRUE) + mvrnorm(n = n21, mu = rep(0, p), Sigma = Sigma_t),
                matrix(mu_t2, nrow = n22, ncol = p, byrow = TRUE) + matrix(theta_2, nrow = n22, ncol = p, byrow = TRUE) + mvrnorm(n = n22, mu = rep(0, p), Sigma = Sigma_t))

    Xa <- rbind(matrix(theta_1, nrow = m1tilde, ncol = p, byrow = TRUE) + mvrnorm(n = m1tilde, mu = rep(0, p), Sigma = Sigma_a),
                matrix(theta_2, nrow = m2tilde, ncol = p, byrow = TRUE) + mvrnorm(n = m2tilde, mu = rep(0, p), Sigma = Sigma_a))


    ## pre-processing
    Xt <- scale(Xt, center = TRUE, scale = FALSE)
    Xa <- scale(Xa, center = TRUE, scale = FALSE)

    if(length(eta) > 1){
      eta_use <- eta_tuning_general(Xt = Xt, Xa = Xa, eta_range = eta_range,
                                    sparse_trt = sqrt((p2)/p),
                                    sparse_ctst = sqrt((p1)/p),
                                    num_comps = 20,
                                    km_cluster = 2,
                                    ckm_cluster = 2,
                                    plot = TRUE)$eta_opt
      eta_opt_set <- c(eta_opt_set, eta_use)
      cat(paste("Use optimal eta:", eta_use, "\n"))
    }



    # scPCA note using sumabsv and sumabsu is equivalent to using sumabs = sumabsv/sqrt(dimension) which is in cPCA()
    cV <- PMD(cov(Xt)-eta_use*cov(Xa), K=20, type = "standard", sumabsv = sqrt(p1), sumabsu = sqrt(p1), center = FALSE, trace = FALSE)$v
    cU <- PMD(cov(Xt)-eta_use*cov(Xa), K=20, type = "standard", sumabsv = sqrt(p1), sumabsu = sqrt(p1), center = FALSE, trace = FALSE)$u
    cPC <- Xt%*%cU[,which.max(diag(t(cU)%*%cV))]

    Ut <- PMD(cov(Xt), K=2, type = "standard", sumabsv = sqrt(p2), sumabsu = sqrt(p2), center = FALSE, trace = FALSE)$v
    PCt <- Xt%*%Ut[,1:2]

    # kmeans clustering using PCs
    ARI_trt_nuisance <- c(ARI_trt_nuisance, spec_clust(PCt[,1], group = covariates_t$nuisance))
    ARI_trt_interesting <- c(ARI_trt_interesting, spec_clust(PCt[,1], group = covariates_t$interesting))
    ARI_ctst_nuisance <- c(ARI_ctst_nuisance, spec_clust(cPC[,1], group = covariates_t$nuisance))
    ARI_ctst_interesting <- c(ARI_ctst_interesting, spec_clust(cPC[,1], group = covariates_t$interesting))

  }
  return(data.frame(ARI_ctst_interesting = ARI_ctst_interesting,
                    ARI_ctst_nuisance = ARI_ctst_nuisance,
                    ARI_trt_interesting = ARI_trt_interesting,
                    ARI_trt_nuisance = ARI_trt_nuisance,
                    eta = eta_opt_set))
}
