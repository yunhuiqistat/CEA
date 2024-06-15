#' Simulation function for \eqn{\eta} tuning algorithm validation
#'
#' This function perform the simulation for validation of eta tuning algorithm.
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
#' @param eta_range candidate range of \eqn{\eta}.
#' @param sparse_trt parameter controlling the sparsity of treatment sparse PCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @param sparse_ctst parameter controlling the sparsity of scPCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @param km_cluster number showing number of nuisance clusters defined from kmeans clustering on treatment PCA.
#' @param ckm_cluster number of interesting clusters defined from kmeans clustering on cPCA/scPCA.
#' @param repetition number of replications of the simulation.
#' @export
#' @return A list containing the selected optimal eta, the two metrics values.
#' @examples
#' library(CEA)
#' library(ggplot2)
#' library(tidyverse)
#' library(ggsci)
#' beta_12_candidates <- c(10, 11, 12)
#' n_candidates <- c(300, 200, 100)
#' p_candidates <- c(10)
#' eta_range <- seq(0.1,5,0.05)
#' results <- setNames(data.frame(matrix(ncol = 7, nrow = 0)),
#'                    c("eta", "eta_true", "nuisance_ARI", "silhouette",
#'                      "beta12", "SampleGroupSize", "FeatureGroupSize"))

#' eta_opt_res <- setNames(data.frame(matrix(ncol = 4, nrow = 0)),
#'                        c("eta_opt", "beta12",  "SampleGroupSize", "FeatureGroupSize"))
#' nuisance_ARI_repe_list <- list()
#' silhouette_repe_list <- list()
#' eta_opt_repe_list <- list()
#' for (beta_12 in beta_12_candidates) {
#'  for (n in n_candidates) {
#'    for(p in p_candidates){
#'      res <- sim_eta_tuning(gnt = n, gnc = n, p = p, q = p, t = p,
#'                            sigma_1t = 1, sigma_2t = 1, sigma_3t = 1,
#'                            sigma_1c = 1, sigma_2c = 1, sigma_3c = 1,
#'                            Sigma_t = NULL, Sigma_c = NULL,
#'                            beta_11 = 7, beta_12 = beta_12, beta_21 = 3, beta_22 = 6,
#'                            gamma_11 = -3, gamma_12 = 3, gamma_21 = -3, gamma_22 = 3,
#'                            eta_range = eta_range, repetition = 10)
#'      res_df <- res$res_df
#'      res_df$beta_12 <- beta_12
#'      res_df$SampleGroupSize <- n
#'      res_df$FeatureGroupSize <- p
#'      results <- rbind(results, res_df)
#'      eta_opt_res_df <- data.frame(eta_opt=res$eta_opt,
#'                                   beta_12=beta_12,
#'                                   SampleGroupSize=n,
#'                                   FeatureGroupSize=p)
#'      eta_opt_res <- rbind(eta_opt_res, eta_opt_res_df)
#'      nuisance_ARI_repe_list <- append(nuisance_ARI_repe_list, list(res$nuisance_ARI_repe))
#'      silhouette_repe_list <- append(silhouette_repe_list, list(res$silhouette_repe))
#'      eta_opt_repe_list <- append(eta_opt_repe_list, list(res$eta_opt_repe))
#'    }
#'  }
#' }
#' eta_opt_res$beta_12 <- paste("beta[12]==", eta_opt_res$beta_12, sep = "")
#' eta_opt_res$SampleGroupSize <- paste("n==", eta_opt_res$SampleGroupSize, sep = "")
#' # evaluate variability of eta_opt, ARI and silhouette
#' unlist(lapply(eta_opt_repe_list, mean))
#' unlist(lapply(eta_opt_repe_list, min))
#' unlist(lapply(eta_opt_repe_list, max))
#' # average SE across candidate eta of the metric
#' round(unlist(lapply(nuisance_ARI_repe_list, function(x){mean(apply(x,2,sd))})),3)
#' round(unlist(lapply(silhouette_repe_list, function(x){mean(apply(x,2,sd))})),3)
#' plot_df <- pivot_longer(results, 3:4, names_to = "Metric", values_to = "Values")
#' plot_df$Interaction <- plot_df$beta_12
#' plot_df$Interaction[plot_df$beta_12 == 10] <- "Interaction==0"
#' plot_df$Interaction[plot_df$beta_12 == 11] <- "Interaction==1"
#' plot_df$Interaction[plot_df$beta_12 == 12] <- "Interaction==2"
#' plot_df$beta_12 <- paste("beta[12]==", plot_df$beta_12, sep = "")
#' plot_df$SampleGroupSize <- paste("n==", plot_df$SampleGroupSize, sep = "")
#' plot_df$eta_true <- factor(plot_df$eta_true, levels = c("BadEta","GoodEta"), labels = c( "Outside range","In range"))
#' plot_df$Metric <- factor(plot_df$Metric, levels = c("nuisance_ARI","silhouette"), labels = c("ARI(nuisance)", "Silhouette score"))
#' eta_opt_res$Interaction <- rep(c("Interaction==0","Interaction==1","Interaction==2"), each=3)
#' p_eta <- ggplot(data = plot_df)+
#'  geom_point(aes(x = eta, y = Values, shape = Metric, color = eta_true), size = 2)+
#'  geom_vline(data = eta_opt_res, aes(xintercept = eta_opt), linetype = "dashed")+
#'  facet_grid(cols = vars(Interaction), rows =vars(SampleGroupSize), labeller = label_parsed)+
#'  theme_bw()+
#'  scale_shape_manual(values=c(5,16))+
#'  theme(
#'   panel.grid = element_blank(),
#'    strip.background = element_rect(fill = "lightblue", color = "black"),
#'    strip.text = element_text(color = "black", face = "bold", size = 14),
#'    text = element_text(size = 14),  # Adjust text size
#'    axis.title = element_text(size = 16),  # Adjust axis title size
#'    axis.text = element_text(size = 12),   # Adjust axis text size
#'    legend.position = "bottom")+
#'  labs(x = expression(eta), y = "", color = "", shape = "")+
#'  scale_color_aaas()
#' p_eta



sim_eta_tuning <- function(gnt, gnc, p, q, t,
                           sigma_1t, sigma_2t, sigma_3t,
                           sigma_1c, sigma_2c, sigma_3c,
                           Sigma_t = NULL, Sigma_c = NULL,
                           beta_11, beta_12, beta_21, beta_22,
                           gamma_11, gamma_12, gamma_21, gamma_22,
                           eta_range, sparse_trt = NULL, sparse_ctst = NULL,
                           km_cluster = 2, ckm_cluster = 2,
                           repetition = 10){

  nuisance_ARI <- matrix(NA, ncol = length(eta_range), nrow = repetition)
  silhouette <- matrix(NA, ncol = length(eta_range), nrow = repetition)
  eta_opt <- c()
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
    tuning_res <- eta_tuning(Xt = Xt, Xc = Xc, eta_range = eta_range,
                             sparse_trt = sparse_trt, sparse_ctst = sparse_ctst,
                             km_cluster = km_cluster, ckm_cluster = ckm_cluster,
                             plot = FALSE)
    nuisance_ARI[repe,] <- tuning_res$nuisance_ARI
    silhouette[repe,] <- tuning_res$silhouette
    eta_opt <- c(eta_opt, tuning_res$eta_opt)
  }
  res_df <- data.frame(eta = eta_range, eta_true = eta_true,
                       nuisance_ARI = apply(nuisance_ARI, 2, mean),
                       silhouette = apply(silhouette, 2, mean))
  return(list(res_df = res_df, eta_opt = mean(eta_opt), nuisance_ARI_repe = nuisance_ARI, silhouette_repe = silhouette, eta_opt_repe = eta_opt))
}
