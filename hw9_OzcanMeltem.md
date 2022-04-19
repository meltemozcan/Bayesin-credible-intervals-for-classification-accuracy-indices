Preliminary Analysis
================
Meltem Ozcan
April 12, 2022

The goal of this project is to build Bayesian credible intervals around
selection (classification) accuracy indices computed for the reference
and focal groups in a selection setting as described in Millsap & Kwok,
2004.

### Description of variables

The data I use in the illustrative example for this methodological
project were collected in 1994-1996 as part of efforts to develop the
International personality Item Pool (Donnellan et al., 2006; Goldberg et
al., 2006) and were previously used in studies of measurement invariance
(Lai and Zhang, 2022; Ock et. al, 2020). The dataset contains 564
participants’ item-level scores on the Mini-IPIP, which is a short
personality inventory made up of 5 facets: Agreeableness, Extraversion,
Openness to Experience, Neuroticism, and Conscientiousness. Each
personality dimension is measured with 4 items on a 1-5 Likert scale.

In this preliminary analysis I focus on the Agreeableness factor which
has the following four items: a2, a5, a7, and a9.

-   a2 “Sympathize with others’ feelings.”
-   a5 “Feel others’ emotions.”
-   a7 “Am not really interested in others.” (R)
-   a9 “Am not interested in other people’s problems.” (R)

The fifth variable we have in the ‘agree’ dataset is sex, where males
(N=239) are coded as “1” and females (N=325) as “2”.

### Mathematical expressions of model and priors

The strength of a participant’s endorsement of each item (observed
values) is thoughtto be driven by the participant’s underlying (true)
level of the unobservable (latent) construct of interest. The observed
variables (a2, a5, a7, a9) and the latent variable(s) (A) are connected
via a linear system of equations under the common factor model.

Model specification:

-   *X*: N x J matrix of observed values

-   J: number of observed variables (j=1,…J; here we have J=4)

-   N: number of individuals in dataset (i=1,…,N; here we have N=564)

-   *η*: latent variable values (N x 1 vector)

-   *α*: latent factor mean

-   *ψ*: latent factor variance

-   *λ*: vector of factor loadings (J x 1 vector)

-   *ν*: measurement intercepts (observed) (J x 1 vector)

-   *θ*: diagonal matrix of the unique factor variances

where, plugging in values for hyperparameters as specified by default in
(will be covered later in the document),

### Code for running Bayesian analyses

``` r
data <- read.table("IPIPFFM.dat", header = TRUE)
agree <- data[, c("sex", "a2", "a5", "a7", "a9")]
head(agree)
```

    ##   sex a2 a5 a7 a9
    ## 1   1  5  4  5  5
    ## 2   1  4  4  4  4
    ## 3   1  5  4  5  4
    ## 4   1  3  3  4  3
    ## 5   1  2  4  4  4
    ## 6   1  4  4  4  4

``` r
m_a <- 'A =~  a2 + a5  + a7 + a9
          a2 ~~ a5'
```

``` r
set.seed(7)
agree_fit <- bcfa(m_a, data = agree, group = "sex",
               group.equal = c("loadings", "Intercepts", "residuals"),
               group.partial = c("a2 ~ 1"), 
               std.lv = TRUE, save.lvs = TRUE, n.chains = 3)
```

    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.001676 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 16.76 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1500 [  0%]  (Warmup)
    ## Chain 1: Iteration:  150 / 1500 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  300 / 1500 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  450 / 1500 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  501 / 1500 [ 33%]  (Sampling)
    ## Chain 1: Iteration:  650 / 1500 [ 43%]  (Sampling)
    ## Chain 1: Iteration:  800 / 1500 [ 53%]  (Sampling)
    ## Chain 1: Iteration:  950 / 1500 [ 63%]  (Sampling)
    ## Chain 1: Iteration: 1100 / 1500 [ 73%]  (Sampling)
    ## Chain 1: Iteration: 1250 / 1500 [ 83%]  (Sampling)
    ## Chain 1: Iteration: 1400 / 1500 [ 93%]  (Sampling)
    ## Chain 1: Iteration: 1500 / 1500 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 16.2085 seconds (Warm-up)
    ## Chain 1:                30.1826 seconds (Sampling)
    ## Chain 1:                46.3911 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 2).
    ## Chain 2: 
    ## Chain 2: Gradient evaluation took 0.001385 seconds
    ## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 13.85 seconds.
    ## Chain 2: Adjust your expectations accordingly!
    ## Chain 2: 
    ## Chain 2: 
    ## Chain 2: Iteration:    1 / 1500 [  0%]  (Warmup)
    ## Chain 2: Iteration:  150 / 1500 [ 10%]  (Warmup)
    ## Chain 2: Iteration:  300 / 1500 [ 20%]  (Warmup)
    ## Chain 2: Iteration:  450 / 1500 [ 30%]  (Warmup)
    ## Chain 2: Iteration:  501 / 1500 [ 33%]  (Sampling)
    ## Chain 2: Iteration:  650 / 1500 [ 43%]  (Sampling)
    ## Chain 2: Iteration:  800 / 1500 [ 53%]  (Sampling)
    ## Chain 2: Iteration:  950 / 1500 [ 63%]  (Sampling)
    ## Chain 2: Iteration: 1100 / 1500 [ 73%]  (Sampling)
    ## Chain 2: Iteration: 1250 / 1500 [ 83%]  (Sampling)
    ## Chain 2: Iteration: 1400 / 1500 [ 93%]  (Sampling)
    ## Chain 2: Iteration: 1500 / 1500 [100%]  (Sampling)
    ## Chain 2: 
    ## Chain 2:  Elapsed Time: 16.4984 seconds (Warm-up)
    ## Chain 2:                31.2445 seconds (Sampling)
    ## Chain 2:                47.7429 seconds (Total)
    ## Chain 2: 
    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 3).
    ## Chain 3: 
    ## Chain 3: Gradient evaluation took 0.001343 seconds
    ## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 13.43 seconds.
    ## Chain 3: Adjust your expectations accordingly!
    ## Chain 3: 
    ## Chain 3: 
    ## Chain 3: Iteration:    1 / 1500 [  0%]  (Warmup)
    ## Chain 3: Iteration:  150 / 1500 [ 10%]  (Warmup)
    ## Chain 3: Iteration:  300 / 1500 [ 20%]  (Warmup)
    ## Chain 3: Iteration:  450 / 1500 [ 30%]  (Warmup)
    ## Chain 3: Iteration:  501 / 1500 [ 33%]  (Sampling)
    ## Chain 3: Iteration:  650 / 1500 [ 43%]  (Sampling)
    ## Chain 3: Iteration:  800 / 1500 [ 53%]  (Sampling)
    ## Chain 3: Iteration:  950 / 1500 [ 63%]  (Sampling)
    ## Chain 3: Iteration: 1100 / 1500 [ 73%]  (Sampling)
    ## Chain 3: Iteration: 1250 / 1500 [ 83%]  (Sampling)
    ## Chain 3: Iteration: 1400 / 1500 [ 93%]  (Sampling)
    ## Chain 3: Iteration: 1500 / 1500 [100%]  (Sampling)
    ## Chain 3: 
    ## Chain 3:  Elapsed Time: 17.31 seconds (Warm-up)
    ## Chain 3:                30.6083 seconds (Sampling)
    ## Chain 3:                47.9183 seconds (Total)
    ## Chain 3: 
    ## Computing posterior predictives...

``` r
summary(agree_fit)
```

    ## blavaan (0.4-1) results of 1000 samples after 500 adapt/burnin iterations
    ## 
    ##   Number of observations per group         
    ##   1                                                239
    ##   2                                                325
    ## 
    ##   Number of missing patterns per group     
    ##   1                                                  1
    ##   2                                                  1
    ## 
    ##   Statistic                                 MargLogLik         PPP
    ##   Value                                      -2724.953       0.518
    ## 
    ## 
    ## Group 1 []:
    ## 
    ## Latent Variables:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   A =~                                                                         
    ##     a2      (.p1.)    0.337    0.050    0.242    0.442    1.002    normal(0,10)
    ##     a5      (.p2.)    0.670    0.056    0.566    0.780    1.000    normal(0,10)
    ##     a7      (.p3.)    0.634    0.050    0.538    0.735    0.999    normal(0,10)
    ##     a9      (.p4.)    0.783    0.056    0.678    0.894    1.000    normal(0,10)
    ## 
    ## Covariances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##  .a2 ~~                                                                        
    ##    .a5                0.169    0.036    0.099    0.241    1.001       beta(1,1)
    ## 
    ## Intercepts:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .a2                3.945    0.053    3.839    4.044    1.001    normal(0,32)
    ##    .a5      (.12.)    3.564    0.060    3.447    3.678    1.003    normal(0,32)
    ##    .a7      (.13.)    4.012    0.057    3.899    4.123    1.001    normal(0,32)
    ##    .a9      (.14.)    3.611    0.062    3.483    3.732    1.003    normal(0,32)
    ##     A                 0.000                                                    
    ## 
    ## Variances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .a2      (.p6.)    0.548    0.035    0.483    0.619    1.000 gamma(1,.5)[sd]
    ##    .a5      (.p7.)    0.539    0.043    0.458    0.628    0.999 gamma(1,.5)[sd]
    ##    .a7      (.p8.)    0.436    0.035    0.371    0.509    1.001 gamma(1,.5)[sd]
    ##    .a9      (.p9.)    0.332    0.041    0.254    0.413    1.000 gamma(1,.5)[sd]
    ##     A                 1.000                                                    
    ## 
    ## 
    ## Group 2 []:
    ## 
    ## Latent Variables:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   A =~                                                                         
    ##     a2      (.p1.)    0.337    0.050    0.242    0.442    1.002                
    ##     a5      (.p2.)    0.670    0.056    0.566    0.780    1.000                
    ##     a7      (.p3.)    0.634    0.050    0.538    0.735    0.999                
    ##     a9      (.p4.)    0.783    0.056    0.678    0.894    1.000                
    ## 
    ## Covariances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##  .a2 ~~                                                                        
    ##    .a5                0.099    0.037    0.027    0.173    1.000       beta(1,1)
    ## 
    ## Intercepts:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .a2                4.094    0.058    3.978    4.203    1.002    normal(0,32)
    ##    .a5      (.12.)    3.564    0.060    3.447    3.678    1.003                
    ##    .a7      (.13.)    4.012    0.057    3.899    4.123    1.001                
    ##    .a9      (.14.)    3.611    0.062    3.483    3.732    1.003                
    ##     A                 0.726    0.101    0.518    0.919    1.003    normal(0,10)
    ## 
    ## Variances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .a2      (.p6.)    0.548    0.035    0.483    0.619    1.000                
    ##    .a5      (.p7.)    0.539    0.043    0.458    0.628    0.999                
    ##    .a7      (.p8.)    0.436    0.035    0.371    0.509    1.001                
    ##    .a9      (.p9.)    0.332    0.041    0.254    0.413    1.000                
    ##     A                 0.551    0.094    0.385    0.756    1.000 gamma(1,.5)[sd]

``` r
blavInspect(agree_fit, "est")
```

    ## $`1`
    ## $`1`$lambda
    ##        A
    ## a2 0.337
    ## a5 0.670
    ## a7 0.634
    ## a9 0.783
    ## 
    ## $`1`$theta
    ##    a2    a5    a7    a9   
    ## a2 0.548                  
    ## a5 0.169 0.539            
    ## a7 0.000 0.000 0.436      
    ## a9 0.000 0.000 0.000 0.332
    ## 
    ## $`1`$psi
    ##   A
    ## A 1
    ## 
    ## $`1`$nu
    ##    intrcp
    ## a2  3.945
    ## a5  3.564
    ## a7  4.012
    ## a9  3.611
    ## 
    ## $`1`$alpha
    ##   intrcp
    ## A      0
    ## 
    ## 
    ## $`2`
    ## $`2`$lambda
    ##        A
    ## a2 0.337
    ## a5 0.670
    ## a7 0.634
    ## a9 0.783
    ## 
    ## $`2`$theta
    ##    a2    a5    a7    a9   
    ## a2 0.548                  
    ## a5 0.099 0.539            
    ## a7 0.000 0.000 0.436      
    ## a9 0.000 0.000 0.000 0.332
    ## 
    ## $`2`$psi
    ##   A    
    ## A 0.551
    ## 
    ## $`2`$nu
    ##    intrcp
    ## a2  4.094
    ## a5  3.564
    ## a7  4.012
    ## a9  3.611
    ## 
    ## $`2`$alpha
    ##   intrcp
    ## A  0.726

### A convergence check of MCMC

``` r
# Check for convergence, mixing
a_mcmc <- blavInspect(agree_fit, "mcmc") # mcmc draws for each parameter
a_mcmc3 <- rbind(a_mcmc[[1]], a_mcmc[[2]], a_mcmc[[3]]) # combine chains
a_draws <- a_mcmc3[,c(1:13,18,23,24,28)] # remove duplicates from combined chain
a_mcmc[[1]] <- a_mcmc[[1]][,c(1:13,18,23,24,28)] # remove duplicates in each chain
a_mcmc[[2]] <- a_mcmc[[2]][,c(1:13,18,23,24,28)]
a_mcmc[[3]] <- a_mcmc[[3]][,c(1:13,18,23,24,28)]

colnames(a_draws) <- colnames(a_mcmc[[1]]) <- colnames(a_mcmc[[2]]) <-
  colnames(a_mcmc[[3]]) <- c("lambda1", "lambda2", "lambda3", "lambda4", 
                             "Theta_cov_f", "Theta_var1", "Theta_var2", 
                             "Theta_var3", "Theta_var4", "nu1_f", "nu2", "nu3", 
                             "nu4", "Theta_cov_r", "Psi_r", "nu1_r", "alpha_r" )
```

``` r
mcmc_trace(a_mcmc) 
```

![](hw9_OzcanMeltem_files/figure-gfm/agree-convergence-1.png)<!-- -->

``` r
mcmc_rank_hist(a_mcmc)
```

![](hw9_OzcanMeltem_files/figure-gfm/agree-convergence-2.png)<!-- -->

``` r
a_mcmc %>% summarize_draws() %>% knitr::kable(digits = 2)
```

| variable      | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:--------------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| lambda1       | 0.34 |   0.34 | 0.05 | 0.05 | 0.26 | 0.42 |    1 |   2564.44 |   2371.41 |
| lambda2       | 0.67 |   0.67 | 0.06 | 0.06 | 0.58 | 0.76 |    1 |   2108.66 |   2037.32 |
| lambda3       | 0.63 |   0.63 | 0.05 | 0.05 | 0.55 | 0.72 |    1 |   2597.80 |   2251.53 |
| lambda4       | 0.78 |   0.78 | 0.06 | 0.06 | 0.69 | 0.88 |    1 |   2340.23 |   2174.12 |
| Theta\_cov\_f | 0.17 |   0.17 | 0.04 | 0.04 | 0.11 | 0.23 |    1 |   3400.70 |   1943.05 |
| Theta\_var1   | 0.55 |   0.55 | 0.03 | 0.04 | 0.49 | 0.61 |    1 |   3304.47 |   2541.51 |
| Theta\_var2   | 0.54 |   0.54 | 0.04 | 0.04 | 0.47 | 0.61 |    1 |   3191.95 |   2394.90 |
| Theta\_var3   | 0.44 |   0.44 | 0.04 | 0.04 | 0.38 | 0.50 |    1 |   3600.57 |   2352.12 |
| Theta\_var4   | 0.33 |   0.33 | 0.04 | 0.04 | 0.27 | 0.40 |    1 |   3015.95 |   2096.15 |
| nu1\_f        | 3.95 |   3.95 | 0.05 | 0.05 | 3.86 | 4.03 |    1 |   2648.27 |   2065.63 |
| nu2           | 3.56 |   3.56 | 0.06 | 0.06 | 3.46 | 3.66 |    1 |   1978.23 |   2234.96 |
| nu3           | 4.01 |   4.01 | 0.06 | 0.06 | 3.92 | 4.10 |    1 |   1948.02 |   1551.20 |
| nu4           | 3.61 |   3.61 | 0.06 | 0.06 | 3.51 | 3.71 |    1 |   1858.04 |   2071.19 |
| Theta\_cov\_r | 0.10 |   0.10 | 0.04 | 0.04 | 0.04 | 0.16 |    1 |   3365.72 |   2210.40 |
| Psi\_r        | 0.55 |   0.54 | 0.09 | 0.09 | 0.41 | 0.72 |    1 |   2452.56 |   2454.12 |
| nu1\_r        | 4.09 |   4.09 | 0.06 | 0.06 | 4.00 | 4.19 |    1 |   2630.36 |   2446.91 |
| alpha\_r      | 0.73 |   0.73 | 0.10 | 0.10 | 0.56 | 0.89 |    1 |   1733.31 |   2228.69 |

The trace plots indicate good mixing and suggest that the three chains
sampled the same target distribution as the lines frequently cros and
explore the same region. The rank histograms give support to this
finding as we see that the distributions appear to be approximately
uniform for each model parameter for each chain. The R-hat value
equaling 1 for each parameter suggests little between-chain variability,
indicating convergence. We see that while ESS is slightly lower for the
tail regions, both bulk-ESS and tail-ESS are large enough to allow us to
conclude that the chains converged.

### A table and/or a figure showing the posterior distributions of the key model parameters

``` r
# Now combining the chains
a_draws %>%
    summarize_draws() %>% knitr::kable(digits = 2)
```

| variable      | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:--------------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| lambda1       | 0.34 |   0.34 | 0.05 | 0.05 | 0.26 | 0.42 |    1 |   2553.68 |   2356.09 |
| lambda2       | 0.67 |   0.67 | 0.06 | 0.06 | 0.58 | 0.76 |    1 |   2089.14 |   2025.98 |
| lambda3       | 0.63 |   0.63 | 0.05 | 0.05 | 0.55 | 0.72 |    1 |   2587.94 |   2239.66 |
| lambda4       | 0.78 |   0.78 | 0.06 | 0.06 | 0.69 | 0.88 |    1 |   2321.86 |   2140.95 |
| Theta\_cov\_f | 0.17 |   0.17 | 0.04 | 0.04 | 0.11 | 0.23 |    1 |   3378.26 |   1924.11 |
| Theta\_var1   | 0.55 |   0.55 | 0.03 | 0.04 | 0.49 | 0.61 |    1 |   3243.09 |   2521.33 |
| Theta\_var2   | 0.54 |   0.54 | 0.04 | 0.04 | 0.47 | 0.61 |    1 |   3181.91 |   2385.01 |
| Theta\_var3   | 0.44 |   0.44 | 0.04 | 0.04 | 0.38 | 0.50 |    1 |   3577.29 |   2320.44 |
| Theta\_var4   | 0.33 |   0.33 | 0.04 | 0.04 | 0.27 | 0.40 |    1 |   3007.61 |   2077.01 |
| nu1\_f        | 3.95 |   3.95 | 0.05 | 0.05 | 3.86 | 4.03 |    1 |   2635.14 |   2058.85 |
| nu2           | 3.56 |   3.56 | 0.06 | 0.06 | 3.46 | 3.66 |    1 |   1966.15 |   2221.06 |
| nu3           | 4.01 |   4.01 | 0.06 | 0.06 | 3.92 | 4.10 |    1 |   1923.41 |   1513.25 |
| nu4           | 3.61 |   3.61 | 0.06 | 0.06 | 3.51 | 3.71 |    1 |   1854.22 |   2028.99 |
| Theta\_cov\_r | 0.10 |   0.10 | 0.04 | 0.04 | 0.04 | 0.16 |    1 |   3320.07 |   2190.77 |
| Psi\_r        | 0.55 |   0.54 | 0.09 | 0.09 | 0.41 | 0.72 |    1 |   2437.18 |   2426.67 |
| nu1\_r        | 4.09 |   4.09 | 0.06 | 0.06 | 4.00 | 4.19 |    1 |   2624.90 |   2436.01 |
| alpha\_r      | 0.73 |   0.73 | 0.10 | 0.10 | 0.56 | 0.89 |    1 |   1725.40 |   2209.55 |

``` r
a_draws %>%
    mcmc_areas()
```

![](hw9_OzcanMeltem_files/figure-gfm/areas-1.png)<!-- -->

Determining the priors on model parameters:

``` r
blavInspect(agree_fit, "dp")
```

    ##                nu             alpha            lambda              beta 
    ##    "normal(0,32)"    "normal(0,10)"    "normal(0,10)"    "normal(0,10)" 
    ##             theta               psi               rho             ibpsi 
    ## "gamma(1,.5)[sd]" "gamma(1,.5)[sd]"       "beta(1,1)" "wishart(3,iden)" 
    ##               tau             delta 
    ##   "normal(0,1.5)" "gamma(1,.5)[sd]"

#### Classification accuracy indices, interpretation of results

``` r
# Posterior distribution of the latent variable conditioned on the observed
# variables
agree_fit.lvs <- blavInspect(agree_fit, "lvs") #list of 3, each is 1000 x 564
agree_post_eta <- rbind(agree_fit.lvs[[1]], agree_fit.lvs[[2]],
                        agree_fit.lvs[[3]]) #combine the chains
# The posterior expected value of observed variables conditioned on the sampled
# latent variables
agree_fit.ypred <- blavPredict(agree_fit, type = "ypred")
agree_post_z <- matrix(ncol = 3000, nrow = 564)
rowSums <- vector(mode = "list", length = 3000)
for(i in seq_along(1:length(agree_fit.ypred))){
  agree_post_z[,i] <- t(rowSums(agree_fit.ypred[[i]]))
}

agree_post_CAI <- get_posterior_CAI(agree_post_eta, 
                                    t(agree_post_z), 0.15, agree$sex,
                                    ref = "2", foc = "1")
```

``` r
# Overall classification accuracy indices
agree_post_CAI$overall %>%
    summarize_draws() %>%
    knitr::kable(digits=2)
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| TP       | 0.09 |   0.09 | 0.01 | 0.01 | 0.07 | 0.10 |    1 |   2772.68 |   2852.11 |
| FP       | 0.06 |   0.06 | 0.01 | 0.01 | 0.05 | 0.08 |    1 |   2772.68 |   2876.62 |
| TN       | 0.79 |   0.79 | 0.01 | 0.01 | 0.77 | 0.80 |    1 |   2772.68 |   2852.11 |
| FN       | 0.06 |   0.06 | 0.01 | 0.01 | 0.05 | 0.08 |    1 |   2772.68 |   2876.62 |
| PS       | 0.15 |   0.15 | 0.00 | 0.00 | 0.15 | 0.15 |    1 |   3118.67 |        NA |
| SR       | 0.57 |   0.58 | 0.05 | 0.05 | 0.49 | 0.65 |    1 |   2772.68 |   2852.11 |
| SE       | 0.57 |   0.58 | 0.05 | 0.05 | 0.49 | 0.65 |    1 |   2772.68 |   2852.11 |
| SP       | 0.92 |   0.92 | 0.01 | 0.01 | 0.91 | 0.94 |    1 |   2772.68 |   2852.11 |

``` r
# Reference group classification accuracy indices
agree_post_CAI$reference %>%
    summarize_draws() %>%
    knitr::kable(digits=2)
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| TP       | 0.12 |   0.12 | 0.01 | 0.01 | 0.10 | 0.14 |    1 |   2695.40 |   2601.78 |
| FP       | 0.08 |   0.08 | 0.01 | 0.01 | 0.06 | 0.10 |    1 |   2732.42 |   2841.65 |
| TN       | 0.71 |   0.71 | 0.01 | 0.01 | 0.69 | 0.73 |    1 |   2910.91 |   2785.11 |
| FN       | 0.09 |   0.09 | 0.01 | 0.01 | 0.07 | 0.10 |    1 |   2853.52 |   3000.11 |
| PS       | 0.20 |   0.20 | 0.01 | 0.01 | 0.18 | 0.22 |    1 |   2830.27 |   2816.46 |
| SR       | 0.60 |   0.60 | 0.05 | 0.05 | 0.51 | 0.68 |    1 |   2697.85 |   2804.33 |
| SE       | 0.58 |   0.58 | 0.05 | 0.05 | 0.49 | 0.67 |    1 |   2745.19 |   2976.86 |
| SP       | 0.90 |   0.90 | 0.01 | 0.01 | 0.87 | 0.92 |    1 |   2755.84 |   2927.11 |

``` r
# Focal group classification accuracy indices
agree_post_CAI$focal %>%
    summarize_draws() %>%
    knitr::kable(digits=2)
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| TP       | 0.04 |   0.04 | 0.01 | 0.01 | 0.02 | 0.05 |    1 |   2986.01 |   2879.92 |
| FP       | 0.04 |   0.04 | 0.01 | 0.01 | 0.02 | 0.06 |    1 |   2734.40 |   2717.04 |
| TN       | 0.89 |   0.89 | 0.02 | 0.02 | 0.86 | 0.92 |    1 |   2770.28 |   2780.90 |
| FN       | 0.03 |   0.03 | 0.01 | 0.01 | 0.02 | 0.05 |    1 |   2899.25 |   2894.89 |
| PS       | 0.08 |   0.08 | 0.01 | 0.01 | 0.05 | 0.10 |    1 |   2867.74 |   2836.66 |
| SR       | 0.49 |   0.50 | 0.11 | 0.11 | 0.31 | 0.67 |    1 |   2927.63 |   2886.13 |
| SE       | 0.54 |   0.54 | 0.12 | 0.12 | 0.35 | 0.73 |    1 |   2968.70 |   2690.75 |
| SP       | 0.96 |   0.96 | 0.01 | 0.01 | 0.94 | 0.98 |    1 |   2733.79 |   2711.51 |

We see that PS is higher for the reference group (PS=0.203, 90%
CI=\[0.18, 0.22\]) compared to the focal group (PS=0.09, 90% CI=\[0.05,
0.10\]). Similarly, SR and SE are higher in the reference group than the
focal group, and the CI for these two indices are wider in comparison to
other indices such as PS and SP. Finally, the test appears to be more
successful at excluding unqualified candidates in the focal group
compared to the reference group (as SP\_f = 0.96 with 90% CI\_f =
\[0.94, 0.98\] as opposed to SP\_r = 0.90 with 90% CI\_r = \[0.87,
0.92\]).
