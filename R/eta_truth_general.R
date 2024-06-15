#' Evaluation of candidate \eqn{\eta} for general error covariances
#'
#' This function evaluates whether an given \eqn{\eta} value can identify the interesting subgroups by checking whether the first eigenvector of expected contrast matrix is equal to \eqn{(0_p, \pm 1_q/\sqrt{q}, 0_t)} with tolerance of difference smaller than 0.05.
#' Diagonal error covariances can be specified by sigmas and complex structure can be specified by Sigmas.
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
#' @param eta the value to be evaluated.
#' @export
#' @return A character GoodEta or BadEta indicating whether this \eqn{\eta} value gives correct first population eigenvector in expected contrast matrix.
#' @examples
#' # A case when eta_truth_general == eta_truth
#' library(CEA)
#' eta_range <- seq(0.1,5,0.05)
#' eta_res <- data.frame(eta = eta_range, truth = NA, truth_general = NA)
#' for(eta in eta_range){
#' eta_res[eta_res$eta == eta, "truth_general"] <- eta_truth_general(gnt = 100, gnc = 100,
#'  p = 10, q = 10, t = 10,
#'  sigma_1t = 1, sigma_1c = 1, sigma_2t = 1,
#'  sigma_2c = 1, sigma_3t = 1, sigma_3c = 1,
#'  Sigma_t = NULL, Sigma_c = NULL,
#'  beta_11 = 7, beta_12 = 12, beta_21 = 3, beta_22 = 6,
#'  gamma_11 = -3, gamma_12 = 3, gamma_21 = -3, gamma_22 = 3,
#'  eta = eta)
#'  eta_res[eta_res$eta == eta, "truth"] <- eta_truth(gnt = 100, gnc = 100,
#'  p = 10, q = 10, t = 10,
#'  sigma_1t = 1, sigma_1c = 1, sigma_2t = 1,
#'  sigma_2c = 1, sigma_3t = 1, sigma_3c = 1,
#'  beta_11 = 7, beta_12 = 12, beta_21 = 3, beta_22 = 6,
#'  gamma_11 = -3, gamma_12 = 3, gamma_21 = -3, gamma_22 = 3,
#'  eta = eta)
#'  }
#'  mean(eta_res$truth == eta_res$truth_general)

eta_truth_general <- function(gnt, gnc, p, q, t,
                             sigma_1t, sigma_2t, sigma_3t,
                             sigma_1c, sigma_2c, sigma_3c,
                             Sigma_t = NULL, Sigma_c = NULL,
                             beta_11, beta_12, beta_21, beta_22,
                             gamma_11, gamma_12, gamma_21, gamma_22,
                             eta){
  if(is.null(Sigma_t)){
    Sigma_t <- diag(c(rep(sigma_1t^2, p), rep(sigma_2t^2, q), rep(sigma_3t^2, t)))
  }
  if(is.null(Sigma_c)){
    Sigma_c <- diag(c(rep(sigma_1c^2, p), rep(sigma_2c^2, q), rep(sigma_3c^2, t)))
  }
  B_quant <- (gnt/(4*gnt-1)*(beta_11-beta_12)^2-eta*gnc/(4*gnc-1)*(beta_21-beta_22)^2)
  C_quant <- (gnt/(4*gnt-1)*(gamma_11-gamma_12)^2-eta*gnc/(4*gnc-1)*(gamma_21-gamma_22)^2)
  exp_ctst <- rbind(matrix(0, nrow = p, ncol = p+q+t),
                    cbind(matrix(0, nrow = q, ncol = p), B_quant*matrix(1,nrow = q, ncol = q), matrix(0, nrow = q, ncol = t)),
                    cbind(matrix(0, nrow = t, ncol = p), matrix(0, nrow = t, ncol = q), C_quant*matrix(1,nrow = t, ncol = t)))+
    (Sigma_t - eta*Sigma_c)
  cU1 <- eigen(exp_ctst)$vectors[,1]
  if(all(abs(cU1-c(rep(0,p),rep(1/sqrt(q),q), rep(0, t)))<0.05) | all(abs(cU1-c(rep(0,p),rep(-1/sqrt(q),q), rep(0, t)))<0.05)){return("GoodEta")}
  else{return("BadEta")}
}
