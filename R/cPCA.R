#' Contrastive PCA
#'
#' This function performs contrastive PCA (cPCA) or sparse contrastive PCA (scPCA), in comparison, also perform treatment/control group only PCA/sPCA
#' @import PMA
#' @param Xt data frame for treatment group, samples in rows, variables in columns.
#' @param Xc data frame for control group, samples in rows, variables in columns.
#' @param eta tuning parameter controlling how much control variation should be contrasted off from the treatment variation. It can be theo \code{eta_opt} from function eta_tuning()
#' @param sparse_ctst parameter controlling the sparsity of scPCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @param sparse_trt parameter controlling the sparsity of treatment sparse PCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @param sparse_ctrl parameter controlling the sparsity of control sparse PCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @export
#' @return A list containing contrastive scores \code{cPC} and loadings \code{cU}, treatment scores \code{PCt} and loadings \code{Ut}, and control scores \code{PCc} and loadings \code{Uc}.
#' @examples
#' library(CEA)
#' data(intro)
#' cpca_res <- cPCA(Xt = intro$cPCA$Xt, Xc = intro$cPCA$Xc, eta = 1.2, sparse_ctst = NULL, sparse_trt = NULL, sparse_ctrl = NULL)


cPCA <- function(Xt, Xc, eta, sparse_ctst = NULL, sparse_trt = NULL, sparse_ctrl = NULL){

  ## pre-processing
  Xt <- scale(Xt, center = TRUE, scale = FALSE)
  Xc <- scale(Xc, center = TRUE, scale = FALSE)

  ## contrastive
  if(is.null(sparse_ctst)){
    # cPCA
    cU <- eigen(cov(Xt)-eta*cov(Xc))$vectors
  }
  else{
    # scPCA
    if(length(sparse_ctst) > 1){
      # use CV to find penalty parameters
      sumabs_ctst <- PMD.cv(cov(Xt)-eta*cov(Xc), type="standard",
                       sumabss = sparse_ctst, center = FALSE, trace = FALSE)$bestsumabs
    }
    else{sumabs_ctst <- sparse_ctst}
    cU <- PMD(cov(Xt)-eta*cov(Xc), K=5, type = "standard", sumabs = sumabs_ctst, center = FALSE, trace = FALSE)$v
  }
  cPC <- Xt%*%cU

  ## treatment only
  if(!is.null(sparse_trt)){
    # use sPCA for treatment group
    if(length(sparse_trt) == 1){sumabs_trt <- sparse_trt} # with provided sparsity
    else{# cv to choose sparsity
      sumabs_trt <- PMD.cv(cov(Xt), type="standard", sumabss = sparse_trt, center = FALSE, trace = FALSE)$bestsumabs
    }
    Ut <- PMD(cov(Xt), K=5, type = "standard", sumabs = sumabs_trt, center = FALSE, trace = FALSE)$v}
  else{Ut <- eigen(cov(Xt))$vectors}  # use PCA for treatment

  PCt <- Xt%*%Ut

  ## control only
  if(!is.null(sparse_ctrl)){
    # use sPCA for control group
    if(length(sparse_ctrl) == 1){sumabs_ctrl <- sparse_ctrl} # with provided sparsity
    else{# cv to choose sparsity
      sumabs_ctrl <- PMD.cv(cov(Xc), type="standard", sumabss = sparse_ctrl, center = FALSE, trace = FALSE)$bestsumabs
    }
    Uc <- PMD(cov(Xc), K=5, type = "standard", sumabs = sumabs_ctrl, center = FALSE, trace = FALSE)$v}
  else{Uc <- eigen(cov(Xc))$vectors}  # use PCA for control

  PCc <- Xc%*%Uc

  return(list(PCt = PCt, Ut = Ut, PCc = PCc, Uc = Uc, cPC = cPC, cU = cU))
}
