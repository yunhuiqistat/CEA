% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_eta_tuning.R
\name{sim_eta_tuning}
\alias{sim_eta_tuning}
\title{Simulation function for \eqn{\eta} tuning algorithm validation}
\usage{
sim_eta_tuning(
  gnt,
  gnc,
  p,
  q,
  t,
  sigma_1t,
  sigma_2t,
  sigma_3t,
  sigma_1c,
  sigma_2c,
  sigma_3c,
  Sigma_t = NULL,
  Sigma_c = NULL,
  beta_11,
  beta_12,
  beta_21,
  beta_22,
  gamma_11,
  gamma_12,
  gamma_21,
  gamma_22,
  eta_range,
  sparse_trt = NULL,
  sparse_ctst = NULL,
  km_cluster = 2,
  ckm_cluster = 2,
  repetition = 10
)
}
\arguments{
\item{gnt}{sample size for treatment group.}

\item{gnc}{sample size for control group.}

\item{p}{number of features that do not differs with respect any subgroups.}

\item{q}{number of features that differs with respect to interesting subgroups.}

\item{t}{number of features that differs with respect to nuisance subgroups.}

\item{sigma_1t}{standard deviation of the first p features in treatment group.}

\item{sigma_2t}{standard deviation of the second q features in treatment group.}

\item{sigma_3t}{standard deviation of the third t features in treatment group.}

\item{sigma_1c}{standard deviation of the first p features in control group.}

\item{sigma_2c}{standard deviation of the second q features in control group.}

\item{sigma_3c}{standard deviation of the third t features in control group.}

\item{Sigma_t}{the error covariance matrix in treatment group, if is NULL, use standard deviation specified and assume uncorrelation.}

\item{Sigma_c}{the error covariance matrix in control group, if is NULL, use standard deviation specified and assume uncorrelation.}

\item{beta_11}{mean of the first interesting subgroup in treatment group.}

\item{beta_12}{mean of the second interesting subgroup in treatment group.}

\item{beta_21}{mean of the first interesting subgroup in control group.}

\item{beta_22}{mean of the second interesting subgroup in control group.}

\item{gamma_11}{mean of the first nuisance subgroup in treatment group.}

\item{gamma_12}{mean of the second nuisance subgroup in treatment group.}

\item{gamma_21}{mean of the first nuisance subgroup in control group.}

\item{gamma_22}{mean of the second nuisance subgroup in control group.}

\item{eta_range}{candidate range of \eqn{\eta}.}

\item{sparse_trt}{parameter controlling the sparsity of treatment sparse PCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.}

\item{sparse_ctst}{parameter controlling the sparsity of scPCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.}

\item{km_cluster}{number showing number of nuisance clusters defined from kmeans clustering on treatment PCA.}

\item{ckm_cluster}{number of interesting clusters defined from kmeans clustering on cPCA/scPCA.}

\item{repetition}{number of replications of the simulation.}
}
\value{
A list containing the selected optimal eta, the two metrics values.
}
\description{
This function perform the simulation for validation of eta tuning algorithm.
}
\examples{
library(CEA)
library(ggplot2)
library(tidyverse)
library(ggsci)
beta_12_candidates <- c(10, 11, 12)
n_candidates <- c(300, 200, 100)
p_candidates <- c(10)
eta_range <- seq(0.1,5,0.05)
results <- setNames(data.frame(matrix(ncol = 7, nrow = 0)),
                   c("eta", "eta_true", "nuisance_ARI", "silhouette",
                     "beta12", "SampleGroupSize", "FeatureGroupSize"))
eta_opt_res <- setNames(data.frame(matrix(ncol = 4, nrow = 0)),
                       c("eta_opt", "beta12",  "SampleGroupSize", "FeatureGroupSize"))
nuisance_ARI_repe_list <- list()
silhouette_repe_list <- list()
eta_opt_repe_list <- list()
for (beta_12 in beta_12_candidates) {
 for (n in n_candidates) {
   for(p in p_candidates){
     res <- sim_eta_tuning(gnt = n, gnc = n, p = p, q = p, t = p,
                           sigma_1t = 1, sigma_2t = 1, sigma_3t = 1,
                           sigma_1c = 1, sigma_2c = 1, sigma_3c = 1,
                           Sigma_t = NULL, Sigma_c = NULL,
                           beta_11 = 7, beta_12 = beta_12, beta_21 = 3, beta_22 = 6,
                           gamma_11 = -3, gamma_12 = 3, gamma_21 = -3, gamma_22 = 3,
                           eta_range = eta_range, repetition = 10)
     res_df <- res$res_df
     res_df$beta_12 <- beta_12
     res_df$SampleGroupSize <- n
     res_df$FeatureGroupSize <- p
     results <- rbind(results, res_df)
     eta_opt_res_df <- data.frame(eta_opt=res$eta_opt,
                                  beta_12=beta_12,
                                  SampleGroupSize=n,
                                  FeatureGroupSize=p)
     eta_opt_res <- rbind(eta_opt_res, eta_opt_res_df)
     nuisance_ARI_repe_list <- append(nuisance_ARI_repe_list, list(res$nuisance_ARI_repe))
     silhouette_repe_list <- append(silhouette_repe_list, list(res$silhouette_repe))
     eta_opt_repe_list <- append(eta_opt_repe_list, list(res$eta_opt_repe))
   }
 }
}
eta_opt_res$beta_12 <- paste("beta[12]==", eta_opt_res$beta_12, sep = "")
eta_opt_res$SampleGroupSize <- paste("n==", eta_opt_res$SampleGroupSize, sep = "")
# evaluate variability of eta_opt, ARI and silhouette
unlist(lapply(eta_opt_repe_list, mean))
unlist(lapply(eta_opt_repe_list, min))
unlist(lapply(eta_opt_repe_list, max))
# average SE across candidate eta of the metric
round(unlist(lapply(nuisance_ARI_repe_list, function(x){mean(apply(x,2,sd))})),3)
round(unlist(lapply(silhouette_repe_list, function(x){mean(apply(x,2,sd))})),3)
plot_df <- pivot_longer(results[results$eta \%in\% seq(0.1,5,0.25),], 3:4, names_to = "Metric", values_to = "Values")
plot_df$Interaction <- plot_df$beta_12
plot_df$Interaction[plot_df$beta_12 == 10] <- "Interaction==0"
plot_df$Interaction[plot_df$beta_12 == 11] <- "Interaction==1"
plot_df$Interaction[plot_df$beta_12 == 12] <- "Interaction==2"
plot_df$beta_12 <- paste("beta[12]==", plot_df$beta_12, sep = "")
plot_df$SampleGroupSize <- paste("n==", plot_df$SampleGroupSize, sep = "")
plot_df$eta_true <- factor(plot_df$eta_true, levels = c("BadEta","GoodEta"), labels = c( "Outside range","In range"))
plot_df$Metric <- factor(plot_df$Metric, levels = c("nuisance_ARI","silhouette"), labels = c("ARI(nuisance)", "Silhouette score"))
eta_opt_res$Interaction <- rep(c("Interaction==0","Interaction==1","Interaction==2"), each=3)
p_eta <- ggplot(data = plot_df)+
 geom_point(aes(x = eta, y = Values, shape = Metric, color = eta_true), size = 2)+
 geom_vline(data = eta_opt_res, aes(xintercept = eta_opt), linetype = "dashed")+
 facet_grid(cols = vars(Interaction), rows =vars(SampleGroupSize), labeller = label_parsed)+
 theme_bw()+
 scale_shape_manual(values=c(5,16))+
 theme(
  panel.grid = element_blank(),
   strip.background = element_rect(fill = "lightblue", color = "black"),
   strip.text = element_text(color = "black", face = "bold", size = 14),
   text = element_text(size = 14),  # Adjust text size
   axis.title = element_text(size = 16),  # Adjust axis title size
   axis.text = element_text(size = 12),   # Adjust axis text size
   legend.position = "bottom")+
 labs(x = expression(eta), y = "", color = "", shape = "")+
 scale_color_aaas()
p_eta
}
