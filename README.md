
# CEA

R package for contrastive eigen analysis.

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

### Authors

Yunhui Qi, Peng Liu, Yumou Qiu

### Contact

<qyh@iastate.edu> (Yunhui Qi)

## Installation

You can install the development version of CEA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("yunhuiqistat/CEA")
```

## Website

Check the vignette at <https://yunhuiqistat.github.io/CEA/>.

## Introduction

The goal of CEA is to conduct contrastive eigen analysis including
contrastive principal component (cPCA) and its sparse version (scPCA),
contrastive cross-covariance analysis (3CA) and its sparse version
(3CA). (s)cPCA can be implemented by function `cPCA`. (s)3CA can be
implemented by function `cCCA`.

``` r
library(CEA)
library(ggpubr)
library(ggplot2)
library(ggsci)
```

### cPCA

`eta_tuning_general()`

- Input

  - `Xt` data frame for target group, samples in rows, variables in
    columns. If PCA is applied to correlation matrix, Xt should be
    scaled and centered.  
  - `Xa` data frame for ancillary group, samples in rows, variables in
    columns. If cPCA is applied to correlation matrix, Xc should also be
    scaled and centered.  
  - `eta_range` a vector containing candidate eta range to be tuning
    from.  
  - `plot` a logical quantity indicating whether plot the ARI(nuisance)
    and silhouette score.  
  - `sparse_trt` parameter controlling the sparsity of target sparse
    PCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector,
    refer sumabss in PMD.cv from package PMA, if single number, refer
    sumabs in PMD in package PMA.  
  - `sparse_ctst` parameter controlling the sparsity of scPCA, if NULL,
    use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss
    in PMD.cv from package PMA, if single number, refer sumabs in PMD in
    package PMA.  
  - `num_comps` if scPCA, specify the number of contrastive components
    that will be estimated by PMD, default is 20.
  - `km_cluster` number showing number of nuisance clusters defined from
    kmeans clustering on target PCA.
  - `ckm_cluster` number of interesting clusters defined from kmeans
    clustering on cPCA/scPCA.

- Output

  - A list containing selected value of $\eta$ `eta_opt`, the vector of
    ARI(nuisance) `nuisance_ARI` and the vector of Silhouette
    `silhouette`, and the ggplot object `p`.

`cPCA()`

- Input

  - `Xt` data frame for target group, samples in rows, variables in
    columns.
  - `Xa` data frame for ancillary group, samples in rows, variables in
    columns.
  - `eta` tuning parameter controlling how much ancillary variation
    should be contrasted off from the target variation. It can be theo
    `eta_opt` from function eta_tuning().
  - `sparse_ctst` parameter controlling the sparsity of scPCA, if NULL,
    use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss
    in PMD.cv from package PMA, if single number, refer sumabs in PMD in
    package PMA.  
  - `sparse_trt` parameter controlling the sparsity of target sparse
    PCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector,
    refer sumabss in PMD.cv from package PMA, if single number, refer
    sumabs in PMD in package PMA.  
  - `sparse_ctrl` parameter controlling the sparsity of ancillary sparse
    PCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector,
    refer sumabss in PMD.cv from package PMA, if single number, refer
    sumabs in PMD in package PMA.

- Output

  - A list containing contrastive scores `cPC` and loadings `cU`, target
    scores `PCt` and loadings `Ut`, and ancillary scores `PCa` and
    loadings `Ua`.

1.  We load the simulated data from the `data(intro)`, and use
    `eta_tuning()` to get optimal value of tuning parameter $\eta$.

``` r
eta_range <- seq(0.1, 5, 0.05)
eta_opt <- eta_tuning_general(Xt = intro$cPCA$Xt, Xa = intro$cPCA$Xa, eta_range = eta_range, plot = TRUE, sparse_trt = NULL, sparse_ctst = NULL, num_comps = 20, km_cluster = 2, ckm_cluster = 2)$eta_opt
```

<img src="man/figures/README-eta tuning for cPCA-1.png" style="display: block; margin: auto;" />

2.  We then input the data and the resulting $\eta$ value to cPCA.

``` r
# speccify sparsity parameters to be NULL for cPCA
cpca_res <- cPCA(Xt = intro$cPCA$Xt, Xa = intro$cPCA$Xa, eta = eta_opt, sparse_ctst = NULL, sparse_trt = NULL, sparse_ctrl = NULL)
```

3.  We get the score plot of cPCA. In comparison, we also get the score
    plot of target PCA. Here, as in the paper, we have nuisance subgroup
    factor $C$ and interesting subgroup factor $B$. In the score plot,
    we see the first target score is associated with the nuisance
    subgroups, while the first contrastive score is associated with the
    interesting subgroups. This validates the removal of nuisance effect
    and discovery of hidden subgroups by cPCA.

``` r
plot_data <- data.frame(cPC1 = cpca_res$cPC[,1], cPC2 = cpca_res$cPC[,2],
                        PC1 = cpca_res$PCt[,1], PC2 = cpca_res$PCt[,2],
                        intro$cPCA$covariates_t)
cP <- ggplot(plot_data)+
        geom_point(aes(x = cPC1, y = cPC2, color = B, shape = C))+
        labs(title = "Contrastive", color = "Interesting \n subgroups", shape = "Nuisance \nsubgroups")+
        theme_bw()+
        theme(panel.grid = element_blank(),
              text = element_text(size = 14),  # Adjust text size
              axis.title = element_text(size = 16),  # Adjust axis title size
              axis.text = element_text(size = 12))
Pt <- ggplot(plot_data)+
        geom_point(aes(x=PC1, y=PC2, color= B, shape = C), show.legend = FALSE)+
        labs(title = "Target", color = "Interesting \n subgroups", shape = "Nuisance \nsubgroups")+
        theme_bw()+
        theme(panel.grid = element_blank(),
              text = element_text(size = 14),  # Adjust text size
              axis.title = element_text(size = 16),  # Adjust axis title size
              axis.text = element_text(size = 12))
p_score <- ggarrange(Pt,cP, nrow = 1, common.legend = TRUE, legend = "bottom", align = "v")
print(p_score)
```

<img src="man/figures/README-cPCA usage score plot-1.png" style="display: block; margin: auto;" />

### scPCA

To illustrate the use of function `cPCA()` for sparse cPCA, we use the
barley expression data example. For details of data processing and
references, please refer to `data(barley)`.

1.  We load the processed data from the `data(intro)`, and use assigned
    $\eta = 1$, we also specify the sparsity parameters to be 0.3.

``` r
set.seed(111)
# specify the sparsity parameters for scPCA
scPCA_res <- cPCA(Xt = intro$scPCA$Xt, Xa = intro$scPCA$Xa, eta = 1, sparse_ctst = 0.3, sparse_trt = 0.3, sparse_ctrl = 0.3)
```

2.  We get the score plot of scPCA. In comparison, we also get the score
    plot of target PCA. As seen from the figure, Genotype effect is
    removed by scPCA, and the fungus subgroups is revealed by the first
    scPCA score.

``` r
PCt <- scPCA_res$PCt
PCa <- scPCA_res$PCa
cPC <- scPCA_res$cPC
Ut <- scPCA_res$Ut
Ua <- scPCA_res$Ua
cU <- scPCA_res$cU

plot_data_t <- data.frame(cPC1 = cPC[,1], cPC2 = cPC[,2],
                          PC1 = PCt[,1], PC2 = PCt[,2],
                          Genotype = intro$scPCA$covariates_t$Genotype,
                          Isolate = intro$scPCA$covariates_t$Isolate,
                          Time_bi = factor(intro$scPCA$covariates_t$Time_bi, levels = c("Early","Late")))
plot_data_a <- data.frame(PC1 = PCa[,1], PC2 = PCa[,2], 
                          Genotype = intro$scPCA$covariates_a$Genotype,
                          Isolate = intro$scPCA$covariates_a$Isolate,
                          Time = intro$scPCA$covariates_a$Time,
                          Time_bi = factor(intro$scPCA$covariates_a$Time_bi, levels = c("Early","Late")))

cP <- ggplot(plot_data_t)+
  geom_point(aes(x = cPC1, y = cPC2, color = Isolate, shape = Genotype), size = 2)+
  theme_bw()+
  labs(title = "Contrastive", color = "Fungus")+
  theme(panel.grid = element_blank(),
        text = element_text(size = 14),  # Adjust text size
        axis.title = element_text(size = 16),  # Adjust axis title size
        axis.text = element_text(size = 12),
        legend.position = "right")+
  scale_color_d3()
Pt <- ggplot(plot_data_t)+
  geom_point(aes(x=PC1, y=PC2, color= Isolate, shape = Genotype ), size = 2)+
  labs(title = "Late", color = "Fungus")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size = 14),  # Adjust text size
        axis.title = element_text(size = 16),  # Adjust axis title size
        axis.text = element_text(size = 12),
        legend.position = "right")+
  scale_color_d3()
Pa <- ggplot(plot_data_a)+
  geom_point(aes(x=PC1, y=PC2, color= Isolate, shape = Genotype), size = 2)+
  labs(title = "Early", color = "Fungus")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size = 14),  # Adjust text size
        axis.title = element_text(size = 16),  # Adjust axis title size
        axis.text = element_text(size = 12),
        legend.position = "right")+
  scale_color_d3()

p_score <- ggarrange(Pa, Pt, cP,  nrow=1, ncol=3,  common.legend = TRUE, legend = "bottom")
print(p_score)
```

<img src="man/figures/README-scPCA usage score plot-1.png" style="display: block; margin: auto;" />

### 3CA

`cCCA()`

- Input

  - `Xt` data frame for one data type in target group, samples in rows,
    variables in columns. To conduct analysis based on cross corrlation
    matrix, each feature should be centered and standardized.  
  - `Xa` data frame for one data type in ancillary group, samples in
    rows, variables in columns. To conduct analysis based on cross
    corrlation matrix, each feature should be centered and
    standardized.  
  - `Yt` data frame for another data type in target group, samples in
    rows, variables in columns. To conduct analysis based on cross
    corrlation matrix, each feature should be centered and
    standardized.  
  - `Ya` data frame for another data type in ancillary group, samples in
    rows, variables in columns. To conduct analysis based on cross
    corrlation matrix, each feature should be centered and
    standardized.  
  - `eta` tuning parameter controlling how much ancillary variation
    should be contrasted off from the target variation. If NULL, the
    estimated eta introduced in paper will be used.
  - `sparse_ctst` parameter controlling the sparsity of s3CA, if NULL,
    use 3CA, else use s3CA. For use of s3CA, if vector, refer sumabss in
    PMD.cv from package PMA, if single number, refer sumabs in PMD in
    package PMA.  
  - `sparse_trt` parameter controlling the sparsity of target sparse
    cross covariance/correlation analysis, if NULL, use cross
    covariance/correlation analysis, else use sparse cross
    covariance/correlation analysis. For use of sparse version, if
    vector, refer sumabss in PMD.cv from package PMA, if single number,
    refer sumabs in PMD in package PMA.  
  - `sparse_ctrl` parameter controlling the sparsity of ancillary sparse
    cross covariance/correlation analysis, if NULL, use cross
    covariance/correlation analysis, else use sparse cross
    covariance/correlation analysis. For use of sparse version, if
    vector, refer sumabss in PMD.cv from package PMA, if single number,
    refer sumabs in PMD in package PMA.

- Output

  - A list containing contrastive SVD object `CCA_ctst` including
    singular values `d`, loadings for X `u`, loadings for Y `v`; target
    only SVD object `CCA_trt` and ancillary only SVD object `CCA_ctrl`.

1.  We load the simulated data from the `data(intro)`, and get the top
    estimated latent factor (score) from 3CA. As in the paper, we also
    simulated the external continuous response from linear model using
    the target population latent factor with coefficient 1.5. Here we
    conduct 3CA and get the estimated top latent factor and regress the
    response on the score to recover the coefficient.

``` r
res_3CA <- cCCA(Xa = intro$cCCA$Xa, Ya = intro$cCCA$Ya, 
                Xt = intro$cCCA$Xt, Yt = intro$cCCA$Yt, eta  = NULL,
                sparse_ctst = NULL, sparse_trt = NULL, sparse_ctrl = NULL)
#> 
#>  optimal eta is 3.14944030744854
```

2.  The coefficient 1.5 is recovered by the estimated score from 3CA,
    validating that the target latent factor is estimated by 3CA.

``` r
coefficients(lm(intro$cCCA$response ~ intro$cCCA$Xt %*% res_3CA$CCA_ctst$u[,1]))[2]
#> intro$cCCA$Xt %*% res_3CA$CCA_ctst$u[, 1] 
#>                                  1.506891
```

### s3CA

To illustrate the use of function `s3CA()` for sparse 3CA, we use the
covid metabolome and proteome data example. For details of data
processing and references, please refer to the paper and `data(covid)`.

1.  We load the processed data from the `data(intro)`, and use assigned
    $\eta = 1$, we also specify the sparsity parameters range to use
    cross validation to select the sparsity parameters. s3CA is
    implemented to get the score plot.

``` r
s3CA_res <- cCCA(Xa = intro$s3CA$Xa, Ya = intro$s3CA$Ya, 
                 Xt = intro$s3CA$Xt, Yt = intro$s3CA$Yt, eta  = 1,
                 sparse_ctst = seq(0.2, 1, len=20), 
                 sparse_trt = seq(0.2, 1, len=40), 
                 sparse_ctrl = seq(0.2, 1, len=40))
#> Use provided eta value.
```

2.  We get the score plot from s3CA. In comparison, we also get the
    score plot of covid group and non-covid group cross-covariance
    analysis respectively. As described in the paper, ICU effect is
    removed by contrasting, and the estimated top latent factor by s3CA
    is associated with the external continuous variable - platelet
    counts.

``` r
CCA_ctst <- s3CA_res$CCA_ctst
CCA_ctrl <- s3CA_res$CCA_ctrl
CCA_trt <- s3CA_res$CCA_trt


df_score <- data.frame(cCCX = intro$s3CA$Xt%*%as.matrix(CCA_ctst$u[,1]), 
                       cCCY = intro$s3CA$Yt%*%as.matrix(CCA_ctst$v[,1]),
                       intro$s3CA$covariates_t)
df_score_t <- data.frame(CCtX = intro$s3CA$Xt%*%as.matrix(CCA_trt$u[,1]), 
                         CCtY = intro$s3CA$Yt%*%as.matrix(CCA_trt$v[,1]), 
                         intro$s3CA$covariates_t)

df_score_a <- data.frame(CCcX = intro$s3CA$Xa%*%as.matrix(CCA_ctrl$u[,1]), 
                         CCcY = intro$s3CA$Ya%*%as.matrix(CCA_ctrl$v[,1]), 
                         intro$s3CA$covariates_a)

interested_y <- "Platelet_K.uL"
p_reg <- ggplot(data = data.frame(df_score[,c("cCCY","icu",interested_y)]),
                aes_string(x = "cCCY",y = interested_y))+
  geom_point(aes_string(color = "icu" ))+ 
  geom_smooth(method = "lm", se = TRUE)+
  theme_bw()+
  labs(title = "Contrastive", x = "Yv", color = "")+
  theme(panel.grid = element_blank(),
        text = element_text(size = 14),  # Adjust text size
        axis.title = element_text(size = 16),  # Adjust axis title size
        axis.text = element_text(size = 12),
        legend.position = "right")+
  scale_color_d3()
p_reg_t <- ggplot(data = data.frame(df_score_t[, c("CCtY", "icu",interested_y)]),
                aes_string(x = "CCtY",y = interested_y))+
  geom_point(aes_string(color = "icu" ))+  
  geom_smooth(method = "lm", se = TRUE)+
  theme_bw()+
  labs(title = "COVID", x = "Yv", color = "")+
  theme(panel.grid = element_blank(),
        text = element_text(size = 14),  # Adjust text size
        axis.title = element_text(size = 16),  # Adjust axis title size
        axis.text = element_text(size = 12),
        legend.position = "right")+
  scale_color_d3()
p_reg_a <- ggplot(data = data.frame(df_score_a[, c("CCcY", "icu",interested_y)]),
                aes_string(x = "CCcY",y = interested_y))+
  geom_point(aes_string(color = "icu" ))+ 
  geom_smooth(method = "lm", se = TRUE)+
  theme_bw()+
  labs(title = "Non COVID", x = "Yv", color = "")+
  theme(panel.grid = element_blank(),
        text = element_text(size = 14),  # Adjust text size
        axis.title = element_text(size = 16),  # Adjust axis title size
        axis.text = element_text(size = 12),
        legend.position = "right")+
  scale_color_d3()

p_score <- ggarrange(p_reg_t, p_reg_a, p_reg, common.legend = T, 
                     legend = "bottom", nrow=1, align = "hv")
print(p_score)
```

<img src="man/figures/README-s3CA usage score plot-1.png" style="display: block; margin: auto;" />
