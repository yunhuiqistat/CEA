#' Select \eqn{\eta} for subgroup identification by cPCA or scPCA
#' @import PMA
#' @import ggplot2
#' @importFrom aricode ARI
#' @importFrom cluster silhouette
#' @param Xt data frame for treatment group, samples in rows, variables in columns. If PCA is applied to correlation matrix, Xt should be scaled and centered.
#' @param Xc data frame for control group, samples in rows, variables in columns. If cPCA is applied to correlation matrix, Xc should also be scaled and centered.
#' @param eta_range a vector containing candidate eta range to be tuning from.
#' @param plot a logical quantity indicating whether plot the ARI(nuisance) and silhouette score.
#' @param sparse_trt parameter controlling the sparsity of treatment sparse PCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @param sparse_ctst parameter controlling the sparsity of scPCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @param km_cluster number showing number of nuisance clusters defined from kmeans clustering on treatment PCA
#' @param ckm_cluster number of interesting clusters defined from kmeans clustering on cPCA/scPCA
#' @export
#' @return A list containing selected value of \eqn{\eta}  \code{eta_opt}, the vector of ARI(nuisance) \code{nuisance_ARI} and the vector of Silhouette \code{silhouette}.
#' @examples
#' library(CEA)
#' data(intro)
#' eta_range <- seq(0.1, 5, 0.05)
#' eta_opt <- eta_tuning(Xt = intro$cPCA$Xt, Xc = intro$cPCA$Xc, eta_range = eta_range, plot = TRUE, sparse_trt = NULL, sparse_ctst = NULL)$eta_opt

eta_tuning <- function(Xt, Xc, eta_range = seq(0.5, 10, 0.5), plot = TRUE,
                       sparse_trt = NULL, sparse_ctst = NULL,
                       km_cluster = 2, ckm_cluster = 2){

  ## pre-processing
  Xt <- scale(Xt, center = TRUE, scale = FALSE)
  Xc <- scale(Xc, center = TRUE, scale = FALSE)

  ## trtPCA/sPCA to define nuisance labels
  if(!is.null(sparse_trt)){
    # use sPCA
    if(length(sparse_trt) == 1){sumabs_trt <- sparse_trt} # with provided sparsity
    else{# cv to choose sparsity
      sumabs_trt <- PMD.cv(cov(Xt), type="standard", sumabss = sparse_trt, center = FALSE, trace = FALSE)$bestsumabs
    }
    Ut <- PMD(cov(Xt), K=2, type = "standard", sumabs = sumabs_trt, center = FALSE, trace = FALSE)$v}
  else{Ut <- eigen(cov(Xt))$vectors[,1:2]} # use PCA

  PCt <- Xt%*%Ut
  km <- kmeans(PCt[,1], centers = km_cluster)
  nuisance_labels <- km$cluster

  ## cPCA/scPCA for a range of eta
  nuisance_ARI <- c()
  silhouette <- c()
  for (eta in eta_range){
    if(!is.null(sparse_ctst)){ # use scPCA
      if(length(sparse_ctst) == 1){sumabs_ctst <- sparse_ctst}
      else{# cv to choose best sparsity
        sumabs_ctst <- PMD.cv(cov(Xt)-eta*cov(Xc), type="standard",
                         sumabss = sparse_ctst, center = FALSE, trace = FALSE)$bestsumabs
      }
      cU <- PMD(cov(Xt)-eta*cov(Xc), type="standard", sumabs=sumabs_ctst, K=2,
                center = FALSE, trace = FALSE)$v
    }
    else{cU <- eigen(cov(Xt)-eta*cov(Xc))$vectors[,1:2]} # use cPCA
    cPC <- Xt%*%cU
    ckm <- kmeans(cPC[,1], centers = ckm_cluster)
    interesting_labels <- ckm$cluster
    silhouette_score <- silhouette(interesting_labels, dist(cPC[,1]))
    silhouette <- c(silhouette, mean(silhouette_score[, "sil_width"]))
    nuisance_ARI <- c(nuisance_ARI, ARI(nuisance_labels, interesting_labels))
  }
  eta_opt <- eta_range[which.max(silhouette - nuisance_ARI)]
  if(plot){
    plot_df <- data.frame(eta = rep(eta_range,2), metrics = c(nuisance_ARI, silhouette), type = rep(c("ARI(nuisance)","Silhouette"), each = length(eta_range)))
    p <- ggplot(plot_df)+
      geom_point(aes(x = eta, y = metrics, color = type, shape = type))+
      theme_bw()+
      geom_vline(xintercept = eta_opt, linetype = "dashed")+
      theme(legend.position = "bottom")+
      labs(color = "", shape = "", title = "eta tuning results")
    print(p)
  }
  return(list(eta_opt = eta_opt, eta_range = eta_range, nuisance_ARI = nuisance_ARI, silhouette = silhouette))
}
