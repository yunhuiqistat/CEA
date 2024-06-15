#' Simulation function for scPCA
#'
#' This function performs simulation to prove the removal of nuisance factor and discovery of hidden factor by using scPCA and eta tuning algorithm.
#'
#' @import PMA
#' @import MASS
#' @param gnt sample size for treatment group.
#' @param gnc sample size for control group.
#' @param p number of features that do not differs with respect any subgroups.
#' @param q number of features that differs with respect to interesting subgroups.
#' @param t number of features that differs with respect to nuisance subgroups.
#' @param sigma_1t standard deviation of the first p features in treatment group.
#' @param sigma_2t standard deviation of the second q features in treatment group.
#' @param sigma_3t standard deviation of the third t features in treatment group.
#' @param sigma_1c standard deviation of the first p features in control group.
#' @param sigma_2c standard deviation of the second q features in control group.
#' @param sigma_3c standard deviation of the third t features in control group.
#' @param Sigma_t the error covariance matrix in treatment group, if is NULL, use standard deviation specified and assume uncorrelation.
#' @param Sigma_c the error covariance matrix in control group, if is NULL, use standard deviation specified and assume uncorrelation.
#' @param beta_11 mean of the first interesting subgroup in treatment group.
#' @param beta_12 mean of the second interesting subgroup in treatment group.
#' @param beta_21 mean of the first interesting subgroup in control group.
#' @param beta_22 mean of the second interesting subgroup in control group.
#' @param gamma_11 mean of the first nuisance subgroup in treatment group.
#' @param gamma_12 mean of the second nuisance subgroup in treatment group.
#' @param gamma_21 mean of the first nuisance subgroup in control group.
#' @param gamma_22 mean of the second nuisance subgroup in control group.
#' @param eta if eta is a vector, eta_tuning algorithm will be called to select optimal eta; if a number, eta will be used.
#' @param repetition number of replications of the simulation.
#' @param save_plot logical, whether save the score plot and loading plot for the first replication.
#' @export
#' @return A data frame containing the selected optimal eta, the four metrics values.
#' @examples
#' library(CEA)
#' library(ggplot2)
#' library(ggsci)
#' library(ggh4x)
#' library(tidyverse)
#' set.seed(123)
#' beta_12_candidates <- c(11,12)
#' n_candidates <- c(100)
#' p_candidates <- c(630, 830, 1230)
#' eta_range <-  seq(0.1, 5, 0.05)
#' results <- setNames(data.frame(matrix(ncol = 8, nrow = 0)),
#'                     c("ARI_ctst_discovery", "ARI_ctst_nuisance",
#'                       "ARI_trt_discovery", "ARI_trt_nuisance",
#'                        "eta", "beta12", "SampleGroupSize", "Dimension"))
#' for (beta_12 in beta_12_candidates) {
#'   for (n in n_candidates) {
#'     for(p in p_candidates){
#'       res <- sim_scPCA (gnt = n, gnc = n, p = p, q = 10, t = 10,
#'                        sigma_1t = 1, sigma_1c = 1, sigma_2t = 1,
#'                         sigma_2c = 1, sigma_3t = 1, sigma_3c = 1,
#'                         Sigma_t = NULL, Sigma_c = NULL,
#'                         beta_11 = 7, beta_12 = beta_12, beta_21 = 3, beta_22 = 6,
#'                         gamma_11 = -3, gamma_12 = 3, gamma_21 = -3, gamma_22 = 3,
#'                         eta = eta_range, repetition = 10,
#'                         save_plot = FALSE)
#'       res$beta_12 <- beta_12
#'       res$SampleGroupSize <- n
#'       res$Dimension <- p + 20
#'       results <- rbind(results, res)
#'    }
#'   }
#' }
#' plot_results <- pivot_longer(results, 1:4, names_to = "Type", values_to = "ARI")
#' plot_results$Analysis <- sapply(plot_results$Type, function(x){strsplit(x,"_")[[1]][2]})
#' plot_results$Factor <- sapply(plot_results$Type, function(x){strsplit(x,"_")[[1]][3]})
#' plot_results$Analysis <- ifelse(plot_results$Analysis=="ctst","Contrastive","Treatment")
#' plot_results$Factor <- ifelse(plot_results$Factor=="discovery","Interesting","Nuisance")
#' plot_results$eta <- paste("eta==", plot_results$eta, sep = "")
#' plot_results$beta_12 <- paste("beta[12]==", plot_results$beta_12, sep = "")
#' plot_results$SampleGroupSize <- paste("n=", plot_results$SampleGroupSize, sep = "")
#' plot_results$Dimension <- paste("Dimension==", plot_results$Dimension, sep = "")
#' print(knitr::kable(plot_results %>% group_by(beta_12, Dimension,Factor, Analysis) %>%  summarise(Average_ARI = mean(ARI), SE_ARI = sd(ARI))), format = "markdown")
#' plot_results <- plot_results %>% group_by(beta_12, Factor, Dimension, Analysis) %>%  summarise(Average_ARI = mean(ARI), SE_ARI = sd(ARI))
#' ## barplot
#' plot_results$Interaction <- ifelse(plot_results$beta_12=="beta[12]==11", "Interaction==1","Interaction==2")
#' ggplot(plot_results)+
#'   geom_bar(aes(x = Factor, fill = Analysis, y = Average_ARI+0.5), stat="identity", position ="dodge")+
#'   facet_nested( ~ Interaction + factor(Dimension, levels = c("Dimension==650","Dimension==850","Dimension==1250")), labeller = label_parsed)+
#'   labs(x = "", y="Average ARI")+
#'   scale_y_continuous(breaks = c(0.5,1,1.5), limits = c(0,1.5), labels=c("0","0.5","1"))+
#'   theme_bw()+
#'   theme(
#'     panel.grid = element_blank(),
#'     strip.background = element_rect(fill = "lightblue", color = "black"),
#'     strip.text = element_text(color = "black", face = "bold", size = 14),
#'     text = element_text(size = 14),  # Adjust text size
#'     axis.title = element_text(size = 16),  # Adjust axis title size
#'     axis.text = element_text(size = 12),   # Adjust axis text size
#'     legend.position = "bottom")+
#'   scale_color_aaas()

sim_scPCA <- function(gnt, gnc, p, q, t,
                      sigma_1t, sigma_2t, sigma_3t,
                      sigma_1c, sigma_2c, sigma_3c,
                      Sigma_t = NULL, Sigma_c = NULL,
                      beta_11, beta_12, beta_21, beta_22,
                      gamma_11, gamma_12, gamma_21, gamma_22,
                      eta, repetition = 10, save_plot = FALSE){
  if(length(eta) > 1){
    eta_range <- eta
    eta_true <- c()
    if(is.null(Sigma_t) & is.null(Sigma_c)){
      # feasible eta values can be computed using Inequalities
      for (eta_candi in eta_range) {
        eta_true <- c(eta_true, eta_truth(gnt = gnt, gnc = gnc, p = p, q = q, t = t,
                                          sigma_1t = sigma_1t, sigma_2t = sigma_2t, sigma_3t = sigma_3t,
                                          sigma_1c = sigma_1c, sigma_2c = sigma_2c, sigma_3c = sigma_3c,
                                          beta_11 = beta_11, beta_12 = beta_12, beta_21 = beta_21, beta_22 = beta_22,
                                          gamma_11 = gamma_11, gamma_12 = gamma_12, gamma_21 = gamma_21, gamma_22 = gamma_22,
                                          eta = eta_candi))
      }
    }
    else{
      # feasible eta values can only be computed by eigen decomposition of the expected ctst matrix
      for (eta_candi in eta_range) {
        eta_true <- c(eta_true, eta_truth_general(gnt = gnt, gnc = gnc, p = p, q = q, t = t,
                                                  sigma_1t = sigma_1t, sigma_2t = sigma_2t, sigma_3t = sigma_3t,
                                                  sigma_1c = sigma_1c, sigma_2c = sigma_2c, sigma_3c = sigma_3c,
                                                  Sigma_t = Sigma_t, Sigma_c = Sigma_c,
                                                  beta_11 = beta_11, beta_12 = beta_12, beta_21 = beta_21, beta_22 = beta_22,
                                                  gamma_11 = gamma_11, gamma_12 = gamma_12, gamma_21 = gamma_21, gamma_22 = gamma_22,
                                                  eta = eta_candi))
      }
    }
  }
  else{eta_range <- NULL}

  ARI_trt_nuisance <- c()
  ARI_trt_discovery <- c()
  ARI_ctst_nuisance <- c()
  ARI_ctst_discovery <- c()
  for(repe in 1:repetition){
    # generate data
    if(is.null(Sigma_t)){
      Sigma_t <- diag(c(rep(sigma_1t^2, p), rep(sigma_2t^2, q), rep(sigma_3t^2, t)))
    }
    if(is.null(Sigma_c)){
      Sigma_c <- diag(c(rep(sigma_1c^2, p), rep(sigma_2c^2, q), rep(sigma_3c^2, t)))
    }
    Xt <- cbind(matrix(0, nrow = 4*gnt, ncol=p), # first p features: assumed to have mean 0
                rbind(matrix(beta_11, nrow = 2*gnt, ncol = q), # second q features: have mean beta_11/12 in two consequtive groups
                      matrix(beta_12, nrow = 2*gnt, ncol = q)),
                rbind(matrix(gamma_11, nrow = gnt, ncol = t), # thrid t features: have mean gamma_11/12 in each nt row
                      matrix(gamma_12, nrow = gnt, ncol = t),
                      matrix(gamma_11, nrow = gnt, ncol = t),
                      matrix(gamma_12, nrow = gnt, ncol = t))) +
      mvrnorm(n = 4*gnt, mu = rep(0, p+q+t), Sigma = Sigma_t)

    Xc <- cbind(matrix(0, nrow = 4*gnc, ncol=p), # first p features: assumed to have mean 0
                rbind(matrix(beta_21, nrow = 2*gnc, ncol = q), # second q features: have mean beta_11/12 in two consequtive groups
                      matrix(beta_22, nrow = 2*gnc, ncol = q)),
                rbind(matrix(gamma_21, nrow = gnc, ncol = t), # thrid t features: have mean gamma_11/12 in each nt row
                      matrix(gamma_22, nrow = gnc, ncol = t),
                      matrix(gamma_21, nrow = gnc, ncol = t),
                      matrix(gamma_22, nrow = gnc, ncol = t))) +
      mvrnorm(n = 4*gnc, mu = rep(0, p+q+t), Sigma = Sigma_c)

    covariates_t <- data.frame(A = rep("A=low", 4*gnt),
                               B = rep(c("B=low","B=high"), each = 2*gnt),
                               C = rep(rep(c("C=low","C=high"), each = gnt), 2))
    covariates_c <- data.frame(A = rep("A=high", 4*gnt),
                               B = rep(c("B=low","B=high"), each = 2*gnt),
                               C = rep(rep(c("C=low","C=high"), each = gnt), 2))



    if(!is.null(eta_range)){
      eta <- eta_tuning(Xt = Xt, Xc = Xc, eta_range = eta_range,
                        sparse_trt = sqrt(t/(p+q+t)),
                        sparse_ctst = sqrt(q/(p+q+t)),
                        km_cluster = 2,
                        ckm_cluster = 2,
                        plot = FALSE)$eta_opt
      cat(paste("Use optimal eta:", eta, "\n"))
      # check whether this is a good eta
      if(eta_true[eta_range == eta]=='GoodEta'){
        cat("Seleted optimal eta is GoodEta \n")
      }
      else{
        cat("Seleted optimal eta is BadEta \n")
      }
    }

    ## pre-processing
    Xt <- scale(Xt, center = TRUE, scale = FALSE)
    Xc <- scale(Xc, center = TRUE, scale = FALSE)

    # scPCA note using sumabsv and sumabsu is equalivalent to using sumabs = sumabsv/sqrt(dimension) which is in cPCA()
    Ut <- PMD(cov(Xt), K=2, type = "standard", sumabsv = sqrt(t), sumabsu = sqrt(t), center = FALSE, trace = FALSE)$v
    cU <- PMD(cov(Xt)-eta*cov(Xc), K=2, type = "standard", sumabsv = sqrt(q), sumabsu = sqrt(q), center = FALSE, trace = FALSE)$v

    PCt <- Xt%*%Ut[,1:2]
    cPC <- Xt%*%cU[,1:2]

    # kmeans clustering using PCs
    ARI_trt_nuisance <- c(ARI_trt_nuisance, spec_clust(PCt[,1], group = covariates_t$C))
    ARI_trt_discovery <- c(ARI_trt_discovery, spec_clust(PCt[,1], group = covariates_t$B))
    ARI_ctst_nuisance <- c(ARI_ctst_nuisance, spec_clust(cPC[,1], group = covariates_t$C))
    ARI_ctst_discovery <- c(ARI_ctst_discovery, spec_clust(cPC[,1], group = covariates_t$B))

    if(repe == 1 & save_plot){

      plot_data <- data.frame(cPC1 = cPC[,1], cPC2 = cPC[,2],
                              PC1 = PCt[,1], PC2 = PCt[,2],
                              covariates_t)
      cP <- ggplot(plot_data)+
        geom_point(aes(x = cPC1, y = cPC2, color = B, shape = C))+
        labs(title = "Contrastive")+
        theme_bw()+
        theme(panel.grid = element_blank(),
              text = element_text(size = 14),  # Adjust text size
              axis.title = element_text(size = 16),  # Adjust axis title size
              axis.text = element_text(size = 12))
      Pt <- ggplot(plot_data)+
        geom_point(aes(x=PC1, y=PC2, color= B, shape = C), show.legend = FALSE)+
        labs(title = "Treatment")+
        theme_bw()+
        theme(panel.grid = element_blank(),
              text = element_text(size = 14),  # Adjust text size
              axis.title = element_text(size = 16),  # Adjust axis title size
              axis.text = element_text(size = 12))
      PUt1 <- ggplot(data.frame(index = 1:nrow(Ut), Ut1 = Ut[,1]))+
        geom_point(aes(x=index, y=Ut1))+
        ylim(-0.8,0.8)+
        labs(title="Treatment PCA eigenvector 1", x = "Variable index")+
        theme_bw()+
        theme(text = element_text(size = 14),  # Adjust text size
              axis.title = element_text(size = 16),  # Adjust axis title size
              axis.text = element_text(size = 12))
      PUt2 <- ggplot(data.frame(index = 1:nrow(Ut), Ut2 = Ut[,2]))+
        geom_point(aes(x=index, y=Ut2))+
        ylim(-0.8,0.8)+
        labs(title = "Treatment PCA eigenvector 2", x = "Variable index")+
        theme_bw()+
        theme(text = element_text(size = 14),  # Adjust text size
              axis.title = element_text(size = 16),  # Adjust axis title size
              axis.text = element_text(size = 12))
      PcU1 <- ggplot(data.frame(index = 1:nrow(cU), cU1 = cU[,1]))+
        geom_point(aes(x=index, y=cU1))+
        ylim(-0.8,0.8)+
        labs(title = "Contrastive PCA eigenvector 1", x = "Variable index")+
        theme_bw()+
        theme(text = element_text(size = 14),  # Adjust text size
              axis.title = element_text(size = 16),  # Adjust axis title size
              axis.text = element_text(size = 12))
      PcU2 <- ggplot(data.frame(index = 1:nrow(cU), cU2 = cU[,2]))+
        geom_point(aes(x=index, y=cU2))+
        ylim(-0.8,0.8)+
        labs(title = "Contrastive PCA eigenvector 2", x = "Variable index")+
        theme_bw()+
        theme(text = element_text(size = 14),  # Adjust text size
              axis.title = element_text(size = 16),  # Adjust axis title size
              axis.text = element_text(size = 12))
      p_all <- ggarrange(Pt,cP, PUt1, PcU1, PUt2, PcU2, nrow=3, ncol=2,
                         common.legend = TRUE, legend = "bottom", align = "v")
      print(p_all)
      ggsave(paste("../inst/SavedFigures/", "scPCA","_dim",p+20,"_beta",beta_12,"_etaopt",".pdf", sep = ""), plot = p_all, width = 8, height = 10, units = "in")
    }
  }
  return(data.frame(ARI_ctst_discovery = ARI_ctst_discovery,
                    ARI_ctst_nuisance = ARI_ctst_nuisance,
                    ARI_trt_discovery = ARI_trt_discovery,
                    ARI_trt_nuisance = ARI_trt_nuisance,
                    eta = eta))
}
