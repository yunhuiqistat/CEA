#' Compute ARI of k means clustering and give groups
#'
#' This function performs K means clustering on the input score and computes the ARI between the resulting clusters and the given groups.
#' @importFrom aricode ARI
#' @param score_mat score matrix from multivariate analysis like cPCA, samples in rows.
#' @param group the group member ship of the samples.
#' @export
#' @return ARI between the resulting clusters and the given groups.

spec_clust <- function(score_mat, group){
  km <- kmeans(score_mat, centers = length(unique(group)))
  cut <- km$cluster
  return(round(ARI(cut, group),4))
}
