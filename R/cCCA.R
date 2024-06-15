#' Contrastive Cross Covariance/Correlation Analysis
#'
#' This function performs contrastive cross covariance/correlation analysis (3CA) or sparse contrastive cross covariance/correlation (s3CA), in comparison, also perform treatment/control group only (sparse) cross covariance analysis.
#' @import PMA
#' @param Xt data frame for one data type in treatment group, samples in rows, variables in columns. To conduct analysis based on cross corrlation matrix, each feature should be centered and standardized.
#' @param Xc data frame for one data type in control group, samples in rows, variables in columns. To conduct analysis based on cross corrlation matrix, each feature should be centered and standardized.
#' @param Yt data frame for another data type in treatment group, samples in rows, variables in columns. To conduct analysis based on cross corrlation matrix, each feature should be centered and standardized.
#' @param Yc data frame for another data type in control group, samples in rows, variables in columns. To conduct analysis based on cross corrlation matrix, each feature should be centered and standardized.
#' @param eta tuning parameter controlling how much control variation should be contrasted off from the treatment variation. If NULL, the estimated eta introduced in paper will be used.
#' @param sparse_ctst parameter controlling the sparsity of s3CA, if NULL, use 3CA, else use s3CA. For use of s3CA, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @param sparse_trt parameter controlling the sparsity of treatment sparse cross covariance/correlation analysis, if NULL, use cross covariance/correlation analysis, else use sparse cross covariance/correlation analysis. For use of sparse version, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @param sparse_ctrl parameter controlling the sparsity of control sparse cross covariance/correlation analysis, if NULL, use cross covariance/correlation analysis, else use sparse cross covariance/correlation analysis. For use of sparse version, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @export
#' @return A list containing contrastive SVD object \code{CCA_ctst} including singular values \code{d}, loadings for X \code{u}, loadings for Y \code{v}; treatment only SVD object \code{CCA_trt} and control only SVD object \code{CCA_ctrl} .
#' @examples
#' library(CEA)
#' data(intro)
#' res_3CA <- cCCA(Xc = intro$cCCA$Xc, Yc = intro$cCCA$Yc,
#' Xt = intro$cCCA$Xt, Yt = intro$cCCA$Yt, eta  = NULL,
#' sparse_ctst = NULL, sparse_trt = NULL, sparse_ctrl = NULL)

cCCA <- function(Xc, Yc, Xt, Yt, eta  = NULL,
                    sparse_ctst = NULL, sparse_trt = NULL, sparse_ctrl = NULL){
  set.seed(111)
  # empirical cross-covariance / cross-correlation
  p <- ncol(Xc)
  q <- ncol(Yc)
  nc <- nrow(Xc)
  nt <- nrow(Xt)

  Xt <- scale(Xt, center = TRUE, scale = FALSE)
  Yt <- scale(Yt, center = TRUE, scale = FALSE)
  Xc <- scale(Xc, center = TRUE, scale = FALSE)
  Yc <- scale(Yc, center = TRUE, scale = FALSE)

  ## Treatment group cross covariance analysis
  if(is.null(sparse_trt)){ # non-sparse version
    CCA_trt <- svd(cov(Xt, Yt))
  }
  else{ # sparse version
    if(length(sparse_trt) > 1){ # use cross validation to choose sparse parameters
      cv_pmd_trt <- PMD.cv(cov(Xt, Yt), type="standard",sumabss = sparse_trt, center = FALSE, trace = FALSE)
      sumabs_trt <- cv_pmd_trt$bestsumabs
      v_init_trt <- cv_pmd_trt$v.init
    }
    else{
      sumabs_trt <- sparse_trt
      v_init_trt <- NULL
    }
    CCA_trt <- PMD(cov(Xt, Yt), type="standard", sumabs = sumabs_trt,
                   K = min(c(p, q, nt-1)), v= v_init_trt, center = FALSE, trace = FALSE)
  }


  ## Control group cross covariance analysis
  if(is.null(sparse_ctrl)){ # non-sparse version
    CCA_ctrl <- svd(cov(Xc, Yc))
  }
  else{ # sparse version
    if(length(sparse_ctrl) > 1){ # use cross validation to choose sparse parameters
      cv_pmd_ctrl <- PMD.cv(cov(Xc, Yc), type="standard",sumabss = sparse_ctrl, center = FALSE, trace = FALSE)
      sumabs_ctrl <- cv_pmd_ctrl$bestsumabs
      v_init_ctrl <- cv_pmd_ctrl$v.init
    }
    else{
      sumabs_ctrl <- sparse_ctrl
      v_init_ctrl <- NULL
    }
    CCA_ctrl <- PMD(cov(Xc, Yc), type="standard", sumabs = sumabs_ctrl,
                   K = min(c(p, q, nc-1)), v = v_init_ctrl, center = FALSE, trace = FALSE)
  }

  # eta choice for contrastive analysis
  if(is.null(eta)){
    eta <- CCA_trt$d[1]/CCA_ctrl$d[1]
    cat(paste("\n optimal eta is ", eta, sep = ""))
  }
  else{
    cat("Use provided eta value. \n")
  }

  # (s)3CA
  contrast <- cov(Xt, Yt) - eta * cov(Xc, Yc)
  if(is.null(sparse_ctst)){ # non-sparse version
    CCA_ctst <- svd(contrast)
  }
  else{ # sparse version
    if(length(sparse_ctst) > 1){ # use cross validation to choose sparse parameters
      cv_pmd_ctst <- PMD.cv(contrast, type="standard",sumabss = sparse_ctst, center = FALSE, trace = FALSE)
      sumabs_ctst <- cv_pmd_ctst$bestsumabs
      v_init_ctst <- cv_pmd_ctst$v.init
    }
    else{
      sumabs_ctst <- sparse_ctst
      v_init_ctst <- NULL
    }
    CCA_ctst <- PMD(contrast, type="standard", sumabs = sumabs_ctst,
                    K = min(c(p, q, nt-1)), v = v_init_ctst, center = FALSE, trace = FALSE)
  }
  return(list(CCA_trt = CCA_trt, CCA_ctrl = CCA_ctrl, CCA_ctst = CCA_ctst))
}
