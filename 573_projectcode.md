573project
================
Meltem Ozcan

``` r
library(blavaan) # to fit Bayesian latent variable models
library(dplyr)
library(bayesplot)
library(posterior)
library(here)
source('573_helper.R') #read in helper functions
```

``` r
data <- read.table("IPIPFFM.dat", header = TRUE)
head(data)
```

    ##   sex e1 e2 e3 e4 e5 e6 e7 e8 e9 e10 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 c1 c2 c3 c4
    ## 1   1  4  4  4  3  4  4  4  4  2   4  5  5  1  5  4  4  5  3  5   5  3  2  4  4
    ## 2   1  3  4  2  2  2  2  2  4  2   1  5  4  4  5  4  4  4  4  4   4  4  4  4  4
    ## 3   1  4  5  5  5  5  4  4  4  4   4  4  5  4  4  4  4  5  5  4   4  4  5  5  5
    ## 4   1  2  2  4  4  3  4  3  3  3   3  4  3  5  5  3  3  4  2  3   4  4  2  2  4
    ## 5   1  2  4  3  4  2  4  3  4  3   4  4  2  4  4  4  4  4  5  4   4  3  3  3  3
    ## 6   1  2  4  4  4  4  4  4  4  4   4  5  4  4  4  4  4  4  5  4   4  4  4  4  4
    ##   c5 c6 c7 c8 c9 c10 n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 i1 i2 i3 i4 i5 i6 i7 i8 i9
    ## 1  3  2  3  2  1   4  4  3  3  4  4  3  3  3  1   3  5  4  4  4  5  4  4  5  5
    ## 2  4  4  4  4  4   4  2  4  4  4  2  4  4  4  4   3  4  4  4  4  4  4  4  4  4
    ## 3  5  5  4  2  5   5  1  1  2  2  2  2  2  2  4   1  5  5  5  5  5  5  5  5  5
    ## 4  3  4  2  2  3   4  3  4  3  3  4  4  5  4  3   4  4  4  2  2  4  3  4  2  2
    ## 5  3  3  4  4  3   1  3  2  2  2  2  3  1  2  3   2  3  3  5  3  2  3  3  3  3
    ## 6  4  4  4  5  4   5  2  2  2  2  1  2  2  1  2   1  4  4  4  4  3  4  4  4  3
    ##   i10
    ## 1   5
    ## 2   4
    ## 3   5
    ## 4   5
    ## 5   4
    ## 6   4

``` r
agree_m <- data[data$sex == 1, c("a2", "a5", "a7", "a9")] #239 males
agree_f <- data[data$sex == 2, c("a2", "a5", "a7", "a9") ] #325 females
agree <- data[, c("sex", "a2", "a5", "a7", "a9")]
consc <- data[, c("sex", "c3", "c4", "c8", "c9")]
extr <- data[, c("sex", "e1", "e4", "e6", "e7")]
neur <- data[, c("sex", "n1", "n2", "n6", "n8")]
open <- data[, c("sex", "i2", "i8", "i9", "i10")]  
```

## Fit Bayesian CFA and compute posteriors

``` r
m_a <- 'A =~  a2 + a5  + a7 + a9
          a2 ~~ a5'
m_c <- 'C =~  c3 + c4 + c8 + c9'
m_e <- 'E =~ e1 + e4 + e6 + e7
          e4 ~~ e7  '
m_n <- 'N =~ n1 + n2 + n6 + n8'
m_o <- 'O =~ i2 + i8 + i9 + i10
          i2 ~~ i10
          i8 ~~ i9'
```

### Fit Bayesian CFA and compute posteriors: Agreeableness

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
    ## Chain 1: Gradient evaluation took 0.001635 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 16.35 seconds.
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
    ## Chain 1:  Elapsed Time: 17.0387 seconds (Warm-up)
    ## Chain 1:                31.6003 seconds (Sampling)
    ## Chain 1:                48.639 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 2).
    ## Chain 2: 
    ## Chain 2: Gradient evaluation took 0.001438 seconds
    ## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 14.38 seconds.
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
    ## Chain 2:  Elapsed Time: 16.5184 seconds (Warm-up)
    ## Chain 2:                31.4955 seconds (Sampling)
    ## Chain 2:                48.0139 seconds (Total)
    ## Chain 2: 
    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 3).
    ## Chain 3: 
    ## Chain 3: Gradient evaluation took 0.001389 seconds
    ## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 13.89 seconds.
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
    ## Chain 3:  Elapsed Time: 17.1391 seconds (Warm-up)
    ## Chain 3:                30.331 seconds (Sampling)
    ## Chain 3:                47.4701 seconds (Total)
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

``` r
#Get the posterior distribution of the latent variable conditioned on the
#observed variables
agree_fit.lvs <- blavInspect(agree_fit, "lvs") #list of 3, each is 1000 x 564
agree_post_eta <- rbind(agree_fit.lvs[[1]], agree_fit.lvs[[2]],
                           agree_fit.lvs[[3]]) #combine the chains
#The posterior expected value of observed variables conditioned on the sampled latent variables
agree_fit.ypred <- blavPredict(agree_fit, type = "ypred")
#dim(agree_fit.ypred[[1]]) # agree_fit.ypred = list of 3000, each 564 x 4.

agree_post_z <- matrix(ncol = 3000, nrow = 564)
rowSums <- vector(mode = "list", length = 3000)
for(i in seq_along(1:length(agree_fit.ypred))){
  agree_post_z[,i] <- t(rowSums(agree_fit.ypred[[i]]))
}

agree_post_CAI <- get_posterior_CAI(agree_post_eta, 
                                    t(agree_post_z), 0.15, agree$sex,
                                    ref = "2", foc = "1")
```

### Fit Bayesian CFA and compute posteriors: Conscientiousness

``` r
set.seed(7)

consc_fit <- bcfa(m_c, data = consc, group = "sex",
               group.equal = c("loadings", "Intercepts", "residuals"),
               group.partial = c("c3 ~ 1"), 
               std.lv = TRUE, save.lvs = TRUE, n.chains = 3)
```

    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.001305 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 13.05 seconds.
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
    ## Chain 1:  Elapsed Time: 14.4766 seconds (Warm-up)
    ## Chain 1:                26.6809 seconds (Sampling)
    ## Chain 1:                41.1576 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 2).
    ## Chain 2: 
    ## Chain 2: Gradient evaluation took 0.001243 seconds
    ## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 12.43 seconds.
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
    ## Chain 2:  Elapsed Time: 14.3221 seconds (Warm-up)
    ## Chain 2:                26.5772 seconds (Sampling)
    ## Chain 2:                40.8992 seconds (Total)
    ## Chain 2: 
    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 3).
    ## Chain 3: 
    ## Chain 3: Gradient evaluation took 0.001275 seconds
    ## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 12.75 seconds.
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
    ## Chain 3:  Elapsed Time: 13.5921 seconds (Warm-up)
    ## Chain 3:                24.6641 seconds (Sampling)
    ## Chain 3:                38.2562 seconds (Total)
    ## Chain 3: 
    ## Computing posterior predictives...

``` r
summary(consc_fit)
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
    ##   Value                                      -3140.345       0.350
    ## 
    ## 
    ## Group 1 []:
    ## 
    ## Latent Variables:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   C =~                                                                         
    ##     c3      (.p1.)    0.575    0.065    0.454    0.712    1.000    normal(0,10)
    ##     c4      (.p2.)    0.350    0.040    0.274    0.432    1.000    normal(0,10)
    ##     c8      (.p3.)    0.383    0.046    0.300    0.476    0.999    normal(0,10)
    ##     c9      (.p4.)    0.787    0.081    0.633    0.952    1.000    normal(0,10)
    ## 
    ## Intercepts:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .c3                3.286    0.071    3.141    3.427    1.000    normal(0,32)
    ##    .c4      (.11.)    4.263    0.039    4.187    4.337    1.000    normal(0,32)
    ##    .c8      (.12.)    4.182    0.045    4.095    4.270    1.000    normal(0,32)
    ##    .c9      (.13.)    3.599    0.073    3.459    3.741    1.001    normal(0,32)
    ##     C                 0.000                                                    
    ## 
    ## Variances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .c3      (.p5.)    0.954    0.079    0.805    1.113    0.999 gamma(1,.5)[sd]
    ##    .c4      (.p6.)    0.411    0.031    0.353    0.473    1.000 gamma(1,.5)[sd]
    ##    .c8      (.p7.)    0.642    0.046    0.553    0.731    1.001 gamma(1,.5)[sd]
    ##    .c9      (.p8.)    0.876    0.102    0.671    1.078    1.000 gamma(1,.5)[sd]
    ##     C                 1.000                                                    
    ## 
    ## 
    ## Group 2 []:
    ## 
    ## Latent Variables:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   C =~                                                                         
    ##     c3      (.p1.)    0.575    0.065    0.454    0.712    1.000                
    ##     c4      (.p2.)    0.350    0.040    0.274    0.432    1.000                
    ##     c8      (.p3.)    0.383    0.046    0.300    0.476    0.999                
    ##     c9      (.p4.)    0.787    0.081    0.633    0.952    1.000                
    ## 
    ## Intercepts:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .c3                3.292    0.080    3.134    3.454    1.000    normal(0,32)
    ##    .c4      (.11.)    4.263    0.039    4.187    4.337    1.000                
    ##    .c8      (.12.)    4.182    0.045    4.095    4.270    1.000                
    ##    .c9      (.13.)    3.599    0.073    3.459    3.741    1.001                
    ##     C                 0.135    0.124   -0.114    0.378    1.001    normal(0,10)
    ## 
    ## Variances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .c3      (.p5.)    0.954    0.079    0.805    1.113    0.999                
    ##    .c4      (.p6.)    0.411    0.031    0.353    0.473    1.000                
    ##    .c8      (.p7.)    0.642    0.046    0.553    0.731    1.001                
    ##    .c9      (.p8.)    0.876    0.102    0.671    1.078    1.000                
    ##     C                 1.440    0.280    0.983    2.058    1.000 gamma(1,.5)[sd]

``` r
blavInspect(consc_fit, "est")
```

    ## $`1`
    ## $`1`$lambda
    ##        C
    ## c3 0.575
    ## c4 0.350
    ## c8 0.383
    ## c9 0.787
    ## 
    ## $`1`$theta
    ##    c3    c4    c8    c9   
    ## c3 0.954                  
    ## c4 0.000 0.411            
    ## c8 0.000 0.000 0.642      
    ## c9 0.000 0.000 0.000 0.876
    ## 
    ## $`1`$psi
    ##   C
    ## C 1
    ## 
    ## $`1`$nu
    ##    intrcp
    ## c3  3.286
    ## c4  4.263
    ## c8  4.182
    ## c9  3.599
    ## 
    ## $`1`$alpha
    ##   intrcp
    ## C      0
    ## 
    ## 
    ## $`2`
    ## $`2`$lambda
    ##        C
    ## c3 0.575
    ## c4 0.350
    ## c8 0.383
    ## c9 0.787
    ## 
    ## $`2`$theta
    ##    c3    c4    c8    c9   
    ## c3 0.954                  
    ## c4 0.000 0.411            
    ## c8 0.000 0.000 0.642      
    ## c9 0.000 0.000 0.000 0.876
    ## 
    ## $`2`$psi
    ##   C   
    ## C 1.44
    ## 
    ## $`2`$nu
    ##    intrcp
    ## c3  3.292
    ## c4  4.263
    ## c8  4.182
    ## c9  3.599
    ## 
    ## $`2`$alpha
    ##   intrcp
    ## C  0.135

``` r
#Get the posterior distribution of the latent variable conditioned on the
#observed variables
consc_fit.lvs <- blavInspect(consc_fit, "lvs") #list of 3, each is 1000 x 564
consc_post_eta <- rbind(consc_fit.lvs[[1]], consc_fit.lvs[[2]],
                           consc_fit.lvs[[3]]) #combine the chains
#The posterior expected value of observed variables conditioned on the sampled latent variables
consc_fit.ypred <- blavPredict(consc_fit, type = "ypred")
consc_post_z <- matrix(ncol = 3000, nrow = 564)
rowSums <- vector(mode = "list", length = 3000)
for(i in seq_along(1:length(consc_fit.ypred))){
  consc_post_z[,i] <- t(rowSums(consc_fit.ypred[[i]]))
}

consc_post_CAI <- get_posterior_CAI(consc_post_eta, 
                                    t(consc_post_z), 0.15, consc$sex,
                                    ref = "2", foc = "1")
```

### Fit Bayesian CFA and compute posteriors: Extraversion

``` r
set.seed(7)

extr_fit <- bcfa(m_e, data = extr, group = "sex",
               group.equal = c("loadings", "Intercepts", "residuals"),
               group.partial = c("e1 ~ 1"), 
               std.lv = TRUE, save.lvs = TRUE, n.chains = 3)
```

    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.001657 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 16.57 seconds.
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
    ## Chain 1:  Elapsed Time: 16.1383 seconds (Warm-up)
    ## Chain 1:                27.9671 seconds (Sampling)
    ## Chain 1:                44.1054 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 2).
    ## Chain 2: 
    ## Chain 2: Gradient evaluation took 0.001311 seconds
    ## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 13.11 seconds.
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
    ## Chain 2:  Elapsed Time: 16.3178 seconds (Warm-up)
    ## Chain 2:                26.4237 seconds (Sampling)
    ## Chain 2:                42.7415 seconds (Total)
    ## Chain 2: 
    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 3).
    ## Chain 3: 
    ## Chain 3: Gradient evaluation took 0.001323 seconds
    ## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 13.23 seconds.
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
    ## Chain 3:  Elapsed Time: 15.5135 seconds (Warm-up)
    ## Chain 3:                26.2207 seconds (Sampling)
    ## Chain 3:                41.7341 seconds (Total)
    ## Chain 3: 
    ## Computing posterior predictives...

``` r
summary(extr_fit)
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
    ##   Value                                      -3413.732       0.061
    ## 
    ## 
    ## Group 1 []:
    ## 
    ## Latent Variables:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   E =~                                                                         
    ##     e1      (.p1.)    0.664    0.061    0.547    0.784    1.002    normal(0,10)
    ##     e4      (.p2.)    0.829    0.082    0.674    0.994    1.003    normal(0,10)
    ##     e6      (.p3.)    0.759    0.070    0.624    0.902    1.001    normal(0,10)
    ##     e7      (.p4.)    0.720    0.071    0.584    0.866    1.004    normal(0,10)
    ## 
    ## Covariances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##  .e4 ~~                                                                        
    ##    .e7               -0.150    0.089   -0.329    0.018    1.004       beta(1,1)
    ## 
    ## Intercepts:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .e1                2.314    0.071    2.175    2.449    1.000    normal(0,32)
    ##    .e4      (.12.)    2.852    0.073    2.708    2.991    1.001    normal(0,32)
    ##    .e6      (.13.)    3.183    0.070    3.037    3.320    1.002    normal(0,32)
    ##    .e7      (.14.)    3.028    0.061    2.909    3.146    1.002    normal(0,32)
    ##     E                 0.000                                                    
    ## 
    ## Variances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .e1      (.p6.)    0.785    0.063    0.665    0.912    1.002 gamma(1,.5)[sd]
    ##    .e4      (.p7.)    0.965    0.103    0.762    1.167    1.005 gamma(1,.5)[sd]
    ##    .e6      (.p8.)    0.978    0.083    0.821    1.151    1.002 gamma(1,.5)[sd]
    ##    .e7      (.p9.)    0.583    0.072    0.438    0.722    1.004 gamma(1,.5)[sd]
    ##     E                 1.000                                                    
    ## 
    ## 
    ## Group 2 []:
    ## 
    ## Latent Variables:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   E =~                                                                         
    ##     e1      (.p1.)    0.664    0.061    0.547    0.784    1.002                
    ##     e4      (.p2.)    0.829    0.082    0.674    0.994    1.003                
    ##     e6      (.p3.)    0.759    0.070    0.624    0.902    1.001                
    ##     e7      (.p4.)    0.720    0.071    0.584    0.866    1.004                
    ## 
    ## Covariances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##  .e4 ~~                                                                        
    ##    .e7               -0.113    0.077   -0.274    0.033    1.006       beta(1,1)
    ## 
    ## Intercepts:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .e1                2.133    0.078    1.973    2.284    1.004    normal(0,32)
    ##    .e4      (.12.)    2.852    0.073    2.708    2.991    1.001                
    ##    .e6      (.13.)    3.183    0.070    3.037    3.320    1.002                
    ##    .e7      (.14.)    3.028    0.061    2.909    3.146    1.002                
    ##     E                 0.136    0.106   -0.064    0.347    1.005    normal(0,10)
    ## 
    ## Variances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .e1      (.p6.)    0.785    0.063    0.665    0.912    1.002                
    ##    .e4      (.p7.)    0.965    0.103    0.762    1.167    1.005                
    ##    .e6      (.p8.)    0.978    0.083    0.821    1.151    1.002                
    ##    .e7      (.p9.)    0.583    0.072    0.438    0.722    1.004                
    ##     E                 0.992    0.162    0.719    1.355    1.000 gamma(1,.5)[sd]

``` r
blavInspect(extr_fit, "est")
```

    ## $`1`
    ## $`1`$lambda
    ##        E
    ## e1 0.664
    ## e4 0.829
    ## e6 0.759
    ## e7 0.720
    ## 
    ## $`1`$theta
    ##    e1     e4     e6     e7    
    ## e1  0.785                     
    ## e4  0.000  0.965              
    ## e6  0.000  0.000  0.978       
    ## e7  0.000 -0.150  0.000  0.583
    ## 
    ## $`1`$psi
    ##   E
    ## E 1
    ## 
    ## $`1`$nu
    ##    intrcp
    ## e1  2.314
    ## e4  2.852
    ## e6  3.183
    ## e7  3.028
    ## 
    ## $`1`$alpha
    ##   intrcp
    ## E      0
    ## 
    ## 
    ## $`2`
    ## $`2`$lambda
    ##        E
    ## e1 0.664
    ## e4 0.829
    ## e6 0.759
    ## e7 0.720
    ## 
    ## $`2`$theta
    ##    e1     e4     e6     e7    
    ## e1  0.785                     
    ## e4  0.000  0.965              
    ## e6  0.000  0.000  0.978       
    ## e7  0.000 -0.113  0.000  0.583
    ## 
    ## $`2`$psi
    ##   E    
    ## E 0.992
    ## 
    ## $`2`$nu
    ##    intrcp
    ## e1  2.133
    ## e4  2.852
    ## e6  3.183
    ## e7  3.028
    ## 
    ## $`2`$alpha
    ##   intrcp
    ## E  0.136

``` r
#Get the posterior distribution of the latent variable conditioned on the
#observed variables
extr_fit.lvs <- blavInspect(extr_fit, "lvs") #list of 3, each is 1000 x 564
extr_post_eta <- rbind(extr_fit.lvs[[1]], extr_fit.lvs[[2]],
                           extr_fit.lvs[[3]]) #combine the chains
#The posterior expected value of observed variables conditioned on the sampled latent variables
extr_fit.ypred <- blavPredict(extr_fit, type = "ypred")

extr_post_z <- matrix(ncol = 3000, nrow = 564)
rowSums <- vector(mode = "list", length = 3000)
for(i in seq_along(1:length(extr_fit.ypred))){
  extr_post_z[,i] <- t(rowSums(extr_fit.ypred[[i]]))
}

extr_post_CAI <- get_posterior_CAI(extr_post_eta, 
                                    t(extr_post_z), 0.15, extr$sex,
                                    ref = "2", foc = "1")
```

### Fit Bayesian CFA and compute posteriors: Neuroticism

``` r
set.seed(7)

neur_fit <- bcfa(m_n, data = neur, group = "sex",
               group.equal = c("loadings", "Intercepts", "residuals"),
               group.partial = c("n1 ~ 1"), 
               std.lv = TRUE, save.lvs=TRUE, n.chains = 3)
```

    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.002144 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 21.44 seconds.
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
    ## Chain 1:  Elapsed Time: 13.9345 seconds (Warm-up)
    ## Chain 1:                27.2201 seconds (Sampling)
    ## Chain 1:                41.1546 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 2).
    ## Chain 2: 
    ## Chain 2: Gradient evaluation took 0.001249 seconds
    ## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 12.49 seconds.
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
    ## Chain 2:  Elapsed Time: 13.6994 seconds (Warm-up)
    ## Chain 2:                25.549 seconds (Sampling)
    ## Chain 2:                39.2485 seconds (Total)
    ## Chain 2: 
    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 3).
    ## Chain 3: 
    ## Chain 3: Gradient evaluation took 0.001232 seconds
    ## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 12.32 seconds.
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
    ## Chain 3:  Elapsed Time: 13.3032 seconds (Warm-up)
    ## Chain 3:                26.1963 seconds (Sampling)
    ## Chain 3:                39.4994 seconds (Total)
    ## Chain 3: 
    ## Computing posterior predictives...

``` r
summary(neur_fit)
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
    ##   Value                                      -3335.412       0.025
    ## 
    ## 
    ## Group 1 []:
    ## 
    ## Latent Variables:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   N =~                                                                         
    ##     n1      (.p1.)    0.500    0.053    0.395    0.608    1.000    normal(0,10)
    ##     n2      (.p2.)    0.812    0.064    0.690    0.941    1.001    normal(0,10)
    ##     n6      (.p3.)    0.606    0.055    0.498    0.718    1.000    normal(0,10)
    ##     n8      (.p4.)    0.797    0.066    0.672    0.931    1.001    normal(0,10)
    ## 
    ## Intercepts:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .n1                2.261    0.065    2.131    2.386    0.999    normal(0,32)
    ##    .n2      (.11.)    2.549    0.071    2.407    2.687    1.001    normal(0,32)
    ##    .n6      (.12.)    2.287    0.058    2.173    2.399    1.001    normal(0,32)
    ##    .n8      (.13.)    2.230    0.069    2.095    2.365    1.001    normal(0,32)
    ##     N                 0.000                                                    
    ## 
    ## Variances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .n1      (.p5.)    0.829    0.058    0.721    0.948    1.000 gamma(1,.5)[sd]
    ##    .n2      (.p6.)    0.817    0.076    0.668    0.969    1.000 gamma(1,.5)[sd]
    ##    .n6      (.p7.)    0.732    0.055    0.629    0.843    1.001 gamma(1,.5)[sd]
    ##    .n8      (.p8.)    0.687    0.067    0.551    0.817    1.001 gamma(1,.5)[sd]
    ##     N                 1.000                                                    
    ## 
    ## 
    ## Group 2 []:
    ## 
    ## Latent Variables:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   N =~                                                                         
    ##     n1      (.p1.)    0.500    0.053    0.395    0.608    1.000                
    ##     n2      (.p2.)    0.812    0.064    0.690    0.941    1.001                
    ##     n6      (.p3.)    0.606    0.055    0.498    0.718    1.000                
    ##     n8      (.p4.)    0.797    0.066    0.672    0.931    1.001                
    ## 
    ## Intercepts:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .n1                2.513    0.066    2.383    2.638    1.002    normal(0,32)
    ##    .n2      (.11.)    2.549    0.071    2.407    2.687    1.001                
    ##    .n6      (.12.)    2.287    0.058    2.173    2.399    1.001                
    ##    .n8      (.13.)    2.230    0.069    2.095    2.365    1.001                
    ##     N                 0.056    0.107   -0.150    0.265    1.001    normal(0,10)
    ## 
    ## Variances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .n1      (.p5.)    0.829    0.058    0.721    0.948    1.000                
    ##    .n2      (.p6.)    0.817    0.076    0.668    0.969    1.000                
    ##    .n6      (.p7.)    0.732    0.055    0.629    0.843    1.001                
    ##    .n8      (.p8.)    0.687    0.067    0.551    0.817    1.001                
    ##     N                 1.207    0.195    0.876    1.644    1.001 gamma(1,.5)[sd]

``` r
blavInspect(neur_fit, "est")
```

    ## $`1`
    ## $`1`$lambda
    ##        N
    ## n1 0.500
    ## n2 0.812
    ## n6 0.606
    ## n8 0.797
    ## 
    ## $`1`$theta
    ##    n1    n2    n6    n8   
    ## n1 0.829                  
    ## n2 0.000 0.817            
    ## n6 0.000 0.000 0.732      
    ## n8 0.000 0.000 0.000 0.687
    ## 
    ## $`1`$psi
    ##   N
    ## N 1
    ## 
    ## $`1`$nu
    ##    intrcp
    ## n1  2.261
    ## n2  2.549
    ## n6  2.287
    ## n8  2.230
    ## 
    ## $`1`$alpha
    ##   intrcp
    ## N      0
    ## 
    ## 
    ## $`2`
    ## $`2`$lambda
    ##        N
    ## n1 0.500
    ## n2 0.812
    ## n6 0.606
    ## n8 0.797
    ## 
    ## $`2`$theta
    ##    n1    n2    n6    n8   
    ## n1 0.829                  
    ## n2 0.000 0.817            
    ## n6 0.000 0.000 0.732      
    ## n8 0.000 0.000 0.000 0.687
    ## 
    ## $`2`$psi
    ##   N    
    ## N 1.207
    ## 
    ## $`2`$nu
    ##    intrcp
    ## n1  2.513
    ## n2  2.549
    ## n6  2.287
    ## n8  2.230
    ## 
    ## $`2`$alpha
    ##   intrcp
    ## N  0.056

``` r
#Get the posterior distribution of the latent variable conditioned on the
#observed variables
neur_fit.lvs <- blavInspect(neur_fit, "lvs") #list of 3, each is 1000 x 564
neur_post_eta <- rbind(neur_fit.lvs[[1]], neur_fit.lvs[[2]],
                           neur_fit.lvs[[3]]) #combine the chains
#The posterior expected value of observed variables conditioned on the sampled latent variables
neur_fit.ypred <- blavPredict(neur_fit, type = "ypred")
neur_post_z <- matrix(ncol = 3000, nrow = 564)
rowSums <- vector(mode="list", length = 3000)
for(i in seq_along(1:length(neur_fit.ypred))){
  neur_post_z[,i] <- t(rowSums(neur_fit.ypred[[i]]))
}

neur_post_CAI <- get_posterior_CAI(neur_post_eta, 
                                    t(neur_post_z), 0.15, neur$sex,
                                    ref="2", foc="1")
```

### Fit Bayesian CFA and compute posteriors: Openness to Experience

``` r
set.seed(7)

open_fit <- bcfa(m_o, data = open, group = "sex",
               group.equal = c("loadings", "Intercepts", "residuals"),
               group.partial = c("i2 ~ 1"), 
               std.lv = TRUE, save.lvs=TRUE, n.chains = 3)
```

    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.001625 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 16.25 seconds.
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
    ## Chain 1:  Elapsed Time: 39.8009 seconds (Warm-up)
    ## Chain 1:                72.7411 seconds (Sampling)
    ## Chain 1:                112.542 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 2).
    ## Chain 2: 
    ## Chain 2: Gradient evaluation took 0.001454 seconds
    ## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 14.54 seconds.
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
    ## Chain 2:  Elapsed Time: 34.0881 seconds (Warm-up)
    ## Chain 2:                60.1965 seconds (Sampling)
    ## Chain 2:                94.2846 seconds (Total)
    ## Chain 2: 
    ## 
    ## SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 3).
    ## Chain 3: 
    ## Chain 3: Gradient evaluation took 0.001276 seconds
    ## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 12.76 seconds.
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
    ## Chain 3:  Elapsed Time: 40.5926 seconds (Warm-up)
    ## Chain 3:                80.3871 seconds (Sampling)
    ## Chain 3:                120.98 seconds (Total)
    ## Chain 3:

    ## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
    ## Running the chains for more iterations may help. See
    ## https://mc-stan.org/misc/warnings.html#bulk-ess

    ## Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
    ## Running the chains for more iterations may help. See
    ## https://mc-stan.org/misc/warnings.html#tail-ess

    ## Computing posterior predictives...

    ## Warning: blavaan WARNING: Small effective sample sizes (< 100) for some
    ## parameters.

``` r
summary(open_fit)
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
    ##   Value                                      -3277.044       0.318
    ## 
    ## 
    ## Group 1 []:
    ## 
    ## Latent Variables:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   O =~                                                                         
    ##     i2      (.p1.)    0.587    0.181    0.281    0.906    1.020    normal(0,10)
    ##     i8      (.p2.)    0.457    0.169    0.263    0.914    1.057    normal(0,10)
    ##     i9      (.p3.)    0.514    0.180    0.306    1.006    1.059    normal(0,10)
    ##     i10     (.p4.)    0.639    0.191    0.324    0.986    1.020    normal(0,10)
    ## 
    ## Covariances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##  .i2 ~~                                                                        
    ##    .i10               0.221    0.258   -0.193    0.654    1.023       beta(1,1)
    ##  .i8 ~~                                                                        
    ##    .i9                0.429    0.252   -0.248    0.740    1.054       beta(1,1)
    ## 
    ## Intercepts:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .i2                3.934    0.067    3.807    4.067    1.003    normal(0,32)
    ##    .i8      (.13.)    3.585    0.067    3.459    3.719    1.001    normal(0,32)
    ##    .i9      (.14.)    3.606    0.068    3.477    3.741    1.001    normal(0,32)
    ##    .i10     (.15.)    4.069    0.073    3.929    4.203    1.010    normal(0,32)
    ##     O                 0.000                                                    
    ## 
    ## Variances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .i2      (.p7.)    0.650    0.257    0.180    1.076    1.020 gamma(1,.5)[sd]
    ##    .i8      (.p8.)    1.155    0.244    0.478    1.459    1.051 gamma(1,.5)[sd]
    ##    .i9      (.p9.)    0.989    0.285    0.190    1.315    1.058 gamma(1,.5)[sd]
    ##    .i10     (.10.)    0.700    0.292    0.169    1.179    1.022 gamma(1,.5)[sd]
    ##     O                 1.000                                                    
    ## 
    ## 
    ## Group 2 []:
    ## 
    ## Latent Variables:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   O =~                                                                         
    ##     i2      (.p1.)    0.587    0.181    0.281    0.906    1.020                
    ##     i8      (.p2.)    0.457    0.169    0.263    0.914    1.057                
    ##     i9      (.p3.)    0.514    0.180    0.306    1.006    1.059                
    ##     i10     (.p4.)    0.639    0.191    0.324    0.986    1.020                
    ## 
    ## Covariances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##  .i2 ~~                                                                        
    ##    .i10               0.186    0.267   -0.232    0.628    1.023       beta(1,1)
    ##  .i8 ~~                                                                        
    ##    .i9                0.542    0.259   -0.168    0.840    1.056       beta(1,1)
    ## 
    ## Intercepts:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .i2                4.006    0.095    3.829    4.194    1.014    normal(0,32)
    ##    .i8      (.13.)    3.585    0.067    3.459    3.719    1.001                
    ##    .i9      (.14.)    3.606    0.068    3.477    3.741    1.001                
    ##    .i10     (.15.)    4.069    0.073    3.929    4.203    1.010                
    ##     O                -0.413    0.158   -0.761   -0.140    1.006    normal(0,10)
    ## 
    ## Variances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .i2      (.p7.)    0.650    0.257    0.180    1.076    1.020                
    ##    .i8      (.p8.)    1.155    0.244    0.478    1.459    1.051                
    ##    .i9      (.p9.)    0.989    0.285    0.190    1.315    1.058                
    ##    .i10     (.10.)    0.700    0.292    0.169    1.179    1.022                
    ##     O                 1.647    0.400    1.101    2.629    1.010 gamma(1,.5)[sd]

``` r
blavInspect(open_fit, "est")
```

    ## $`1`
    ## $`1`$lambda
    ##         O
    ## i2  0.587
    ## i8  0.457
    ## i9  0.514
    ## i10 0.639
    ## 
    ## $`1`$theta
    ##     i2    i8    i9    i10  
    ## i2  0.650                  
    ## i8  0.000 1.155            
    ## i9  0.000 0.429 0.989      
    ## i10 0.221 0.000 0.000 0.700
    ## 
    ## $`1`$psi
    ##   O
    ## O 1
    ## 
    ## $`1`$nu
    ##     intrcp
    ## i2   3.934
    ## i8   3.585
    ## i9   3.606
    ## i10  4.069
    ## 
    ## $`1`$alpha
    ##   intrcp
    ## O      0
    ## 
    ## 
    ## $`2`
    ## $`2`$lambda
    ##         O
    ## i2  0.587
    ## i8  0.457
    ## i9  0.514
    ## i10 0.639
    ## 
    ## $`2`$theta
    ##     i2    i8    i9    i10  
    ## i2  0.650                  
    ## i8  0.000 1.155            
    ## i9  0.000 0.542 0.989      
    ## i10 0.186 0.000 0.000 0.700
    ## 
    ## $`2`$psi
    ##   O    
    ## O 1.647
    ## 
    ## $`2`$nu
    ##     intrcp
    ## i2   4.006
    ## i8   3.585
    ## i9   3.606
    ## i10  4.069
    ## 
    ## $`2`$alpha
    ##   intrcp
    ## O -0.413

``` r
#Get the posterior distribution of the latent variable conditioned on the
#observed variables
open_fit.lvs <- blavInspect(open_fit, "lvs") #list of 3, each is 1000 x 564
open_post_eta <- rbind(open_fit.lvs[[1]], open_fit.lvs[[2]],
                           open_fit.lvs[[3]]) #combine the chains
#The posterior expected value of observed variables conditioned on the sampled latent variables
open_fit.ypred <- blavPredict(open_fit, type = "ypred")
open_post_z <- matrix(ncol = 3000, nrow = 564)
rowSums <- vector(mode = "list", length = 3000)
for(i in seq_along(1:length(open_fit.ypred))){
  open_post_z[,i] <- t(rowSums(open_fit.ypred[[i]]))
}

open_post_CAI <- get_posterior_CAI(open_post_eta, 
                                    t(open_post_z), 0.15, open$sex,
                                    ref = "2", foc = "1")
```

## Convergence checks of MCMC

### Convergence checks of MCMC: Agreeableness

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
                             "Theta_cov_r", "Theta_var1", "Theta_var2", 
                             "Theta_var3", "Theta_var4", "nu1_r", "nu2", "nu3", 
                             "nu4", "Theta_cov_f", "Psi_r", "nu1_f", "alpha_f" )
```

``` r
mcmc_trace(a_mcmc, regex_pars=c("lambda", "Psi", "alpha")) 
```

![](573_projectcode_files/figure-gfm/agree-convergence-1.png)<!-- -->

``` r
mcmc_trace(a_mcmc, regex_pars="Theta") 
```

![](573_projectcode_files/figure-gfm/agree-convergence-2.png)<!-- -->

``` r
mcmc_trace(a_mcmc, regex_pars="nu")
```

![](573_projectcode_files/figure-gfm/agree-convergence-3.png)<!-- -->

``` r
mcmc_rank_hist(a_mcmc,regex_pars=c("lambda", "Psi", "alpha"))
```

![](573_projectcode_files/figure-gfm/agree-convergence-4.png)<!-- -->

``` r
mcmc_rank_hist(a_mcmc,regex_pars=c("Theta"))
```

![](573_projectcode_files/figure-gfm/agree-convergence-5.png)<!-- -->

``` r
mcmc_rank_hist(a_mcmc,regex_pars=c("nu"))
```

![](573_projectcode_files/figure-gfm/agree-convergence-6.png)<!-- -->

### Convergence checks of MCMC: Conscientiousness

``` r
# Check for convergence, mixing
c_mcmc <- blavInspect(consc_fit, "mcmc") # mcmc draws for each parameter
c_mcmc3 <- rbind(c_mcmc[[1]], c_mcmc[[2]], c_mcmc[[3]]) # combine chains
c_draws <- c_mcmc3[,c(1:12, 21, 26)] # remove duplicates from combined chain
c_mcmc[[1]] <- c_mcmc[[1]][,c(1:12, 21, 26)] # remove duplicates in each chain
c_mcmc[[2]] <- c_mcmc[[2]][,c(1:12, 21, 26)]
c_mcmc[[3]] <- c_mcmc[[3]][,c(1:12, 21, 26)]

colnames(c_draws) <- colnames(c_mcmc[[1]]) <- colnames(c_mcmc[[2]]) <-
  colnames(c_mcmc[[3]]) <- c("lambda1", "lambda2", "lambda3", "lambda4", 
                            "Theta_var1", "Theta_var2", "Theta_var3", "Theta_var4",
                            "nu1", "nu2", "nu3", "nu4", 
                            "Psi_f", "alpha_f" )
```

``` r
mcmc_trace(c_mcmc, regex_pars=c("lambda", "Psi", "alpha")) 
```

![](573_projectcode_files/figure-gfm/c-convergence-1.png)<!-- -->

``` r
mcmc_trace(c_mcmc, regex_pars="Theta") 
```

![](573_projectcode_files/figure-gfm/c-convergence-2.png)<!-- -->

``` r
mcmc_trace(c_mcmc, regex_pars="nu")
```

![](573_projectcode_files/figure-gfm/c-convergence-3.png)<!-- -->

``` r
mcmc_rank_hist(c_mcmc,regex_pars=c("lambda", "Psi", "alpha"))
```

![](573_projectcode_files/figure-gfm/c-convergence-4.png)<!-- -->

``` r
mcmc_rank_hist(c_mcmc,regex_pars=c("Theta"))
```

![](573_projectcode_files/figure-gfm/c-convergence-5.png)<!-- -->

``` r
mcmc_rank_hist(c_mcmc,regex_pars=c("nu"))
```

![](573_projectcode_files/figure-gfm/c-convergence-6.png)<!-- -->

### Convergence checks of MCMC: Extraversion

``` r
# Check for convergence, mixing
e_mcmc <- blavInspect(extr_fit, "mcmc") # mcmc draws for each parameter
e_mcmc3 <- rbind(e_mcmc[[1]], e_mcmc[[2]], e_mcmc[[3]]) # combine chains
e_draws <- e_mcmc3[,c(1:13,18,23,24,28)] # remove duplicates from combined chain
e_mcmc[[1]] <- e_mcmc[[1]][,c(1:13,18,23,24,28)] # remove duplicates in each chain
e_mcmc[[2]] <- e_mcmc[[2]][,c(1:13,18,23,24,28)]
e_mcmc[[3]] <- e_mcmc[[3]][,c(1:13,18,23,24,28)]

colnames(e_draws) <- colnames(e_mcmc[[1]]) <- colnames(e_mcmc[[2]]) <-
  colnames(e_mcmc[[3]]) <- c("lambda1", "lambda2", "lambda3", "lambda4", 
                             "Theta_cov_r", "Theta_var1", "Theta_var2", 
                             "Theta_var3", "Theta_var4", "nu1_r", "nu2", "nu3", 
                             "nu4", "Theta_cov_f", "Psi_r", "nu1_f", "alpha_f" )
```

``` r
mcmc_trace(e_mcmc, regex_pars=c("lambda", "Psi", "alpha")) 
```

![](573_projectcode_files/figure-gfm/e-convergence-1.png)<!-- -->

``` r
mcmc_trace(e_mcmc, regex_pars="Theta") 
```

![](573_projectcode_files/figure-gfm/e-convergence-2.png)<!-- -->

``` r
mcmc_trace(e_mcmc, regex_pars="nu")
```

![](573_projectcode_files/figure-gfm/e-convergence-3.png)<!-- -->

``` r
mcmc_rank_hist(e_mcmc,regex_pars=c("lambda", "Psi", "alpha"))
```

![](573_projectcode_files/figure-gfm/e-convergence-4.png)<!-- -->

``` r
mcmc_rank_hist(e_mcmc,regex_pars=c("Theta"))
```

![](573_projectcode_files/figure-gfm/e-convergence-5.png)<!-- -->

``` r
mcmc_rank_hist(e_mcmc,regex_pars=c("nu"))
```

![](573_projectcode_files/figure-gfm/e-convergence-6.png)<!-- -->

### Convergence checks of MCMC: Neuroticsm

``` r
# Check for convergence, mixing
n_mcmc <- blavInspect(neur_fit, "mcmc") # mcmc draws for each parameter
n_mcmc3 <- rbind(n_mcmc[[1]], n_mcmc[[2]], n_mcmc[[3]]) # combine chains
n_draws <- n_mcmc3[,c(1:12, 21, 26)] # remove duplicates from combined chain
n_mcmc[[1]] <- n_mcmc[[1]][,c(1:12, 21, 26)] # remove duplicates in each chain
n_mcmc[[2]] <- n_mcmc[[2]][,c(1:12, 21, 26)]
n_mcmc[[3]] <- n_mcmc[[3]][,c(1:12, 21, 26)]

colnames(n_draws) <- colnames(n_mcmc[[1]]) <- colnames(n_mcmc[[2]]) <-
  colnames(n_mcmc[[3]]) <- c("lambda1", "lambda2", "lambda3", "lambda4", 
                            "Theta_var1", "Theta_var2", "Theta_var3", "Theta_var4",
                            "nu1", "nu2", "nu3", "nu4", 
                            "Psi_f", "alpha_f" )
```

``` r
mcmc_trace(n_mcmc, regex_pars=c("lambda", "Psi", "alpha")) 
```

![](573_projectcode_files/figure-gfm/n-convergence-1.png)<!-- -->

``` r
mcmc_trace(n_mcmc, regex_pars="Theta") 
```

![](573_projectcode_files/figure-gfm/n-convergence-2.png)<!-- -->

``` r
mcmc_trace(n_mcmc, regex_pars="nu")
```

![](573_projectcode_files/figure-gfm/n-convergence-3.png)<!-- -->

``` r
mcmc_rank_hist(n_mcmc,regex_pars=c("lambda", "Psi", "alpha"))
```

![](573_projectcode_files/figure-gfm/n-convergence-4.png)<!-- -->

``` r
mcmc_rank_hist(n_mcmc,regex_pars=c("Theta"))
```

![](573_projectcode_files/figure-gfm/n-convergence-5.png)<!-- -->

``` r
mcmc_rank_hist(n_mcmc,regex_pars=c("nu"))
```

![](573_projectcode_files/figure-gfm/n-convergence-6.png)<!-- -->

### Convergence checks of MCMC: Openness to Experience

``` r
# Check for convergence, mixing
o_mcmc <- blavInspect(open_fit, "mcmc") # mcmc draws for each parameter
o_mcmc3 <- rbind(o_mcmc[[1]], o_mcmc[[2]], o_mcmc[[3]]) # combine chains
o_draws <- o_mcmc3[,c(1:14, 19, 20,  25, 26, 30)] # remove duplicates from combined chain
o_mcmc[[1]] <- o_mcmc[[1]][,c(1:14, 19, 20,  25, 26, 30)] # remove duplicates in each chain
o_mcmc[[2]] <- o_mcmc[[2]][,c(1:14, 19, 20, 25, 26, 30)]
o_mcmc[[3]] <- o_mcmc[[3]][,c(1:14, 19, 20,  25, 26, 30)]

colnames(o_draws) <- colnames(o_mcmc[[1]]) <- colnames(o_mcmc[[2]]) <-
  colnames(o_mcmc[[3]]) <- c("lambda1", "lambda2", "lambda3", "lambda4", 
                             "Theta_cov_r1", "Theta_cov_r2",
                             "Theta_var1", "Theta_var2","Theta_var3", "Theta_var4",
                             "nu1_r", "nu2", "nu3", "nu4", 
                             "Theta_cov_f1","Theta_cov_f2", 
                             "Psi_r", "nu1_f", "alpha_f" )
```

``` r
mcmc_trace(o_mcmc, regex_pars=c("lambda", "Psi", "alpha")) 
```

![](573_projectcode_files/figure-gfm/o-convergence-1.png)<!-- -->

``` r
mcmc_trace(o_mcmc, regex_pars="Theta") 
```

![](573_projectcode_files/figure-gfm/o-convergence-2.png)<!-- -->

``` r
mcmc_trace(o_mcmc, regex_pars="nu")
```

![](573_projectcode_files/figure-gfm/o-convergence-3.png)<!-- -->

``` r
mcmc_rank_hist(o_mcmc,regex_pars=c("lambda", "Psi", "alpha"))
```

![](573_projectcode_files/figure-gfm/o-convergence-4.png)<!-- -->

``` r
mcmc_rank_hist(o_mcmc,regex_pars=c("Theta"))
```

![](573_projectcode_files/figure-gfm/o-convergence-5.png)<!-- -->

``` r
mcmc_rank_hist(o_mcmc,regex_pars=c("nu"))
```

![](573_projectcode_files/figure-gfm/o-convergence-6.png)<!-- -->

``` r
# default priors, same for all
blavInspect(agree_fit, "dp")
```

    ##                nu             alpha            lambda              beta 
    ##    "normal(0,32)"    "normal(0,10)"    "normal(0,10)"    "normal(0,10)" 
    ##             theta               psi               rho             ibpsi 
    ## "gamma(1,.5)[sd]" "gamma(1,.5)[sd]"       "beta(1,1)" "wishart(3,iden)" 
    ##               tau             delta 
    ##   "normal(0,1.5)" "gamma(1,.5)[sd]"

``` r
blavInspect(consc_fit, "dp")
```

    ##                nu             alpha            lambda              beta 
    ##    "normal(0,32)"    "normal(0,10)"    "normal(0,10)"    "normal(0,10)" 
    ##             theta               psi               rho             ibpsi 
    ## "gamma(1,.5)[sd]" "gamma(1,.5)[sd]"       "beta(1,1)" "wishart(3,iden)" 
    ##               tau             delta 
    ##   "normal(0,1.5)" "gamma(1,.5)[sd]"

``` r
blavInspect(extr_fit, "dp")
```

    ##                nu             alpha            lambda              beta 
    ##    "normal(0,32)"    "normal(0,10)"    "normal(0,10)"    "normal(0,10)" 
    ##             theta               psi               rho             ibpsi 
    ## "gamma(1,.5)[sd]" "gamma(1,.5)[sd]"       "beta(1,1)" "wishart(3,iden)" 
    ##               tau             delta 
    ##   "normal(0,1.5)" "gamma(1,.5)[sd]"

``` r
blavInspect(neur_fit, "dp")
```

    ##                nu             alpha            lambda              beta 
    ##    "normal(0,32)"    "normal(0,10)"    "normal(0,10)"    "normal(0,10)" 
    ##             theta               psi               rho             ibpsi 
    ## "gamma(1,.5)[sd]" "gamma(1,.5)[sd]"       "beta(1,1)" "wishart(3,iden)" 
    ##               tau             delta 
    ##   "normal(0,1.5)" "gamma(1,.5)[sd]"

``` r
blavInspect(open_fit, "dp")
```

    ##                nu             alpha            lambda              beta 
    ##    "normal(0,32)"    "normal(0,10)"    "normal(0,10)"    "normal(0,10)" 
    ##             theta               psi               rho             ibpsi 
    ## "gamma(1,.5)[sd]" "gamma(1,.5)[sd]"       "beta(1,1)" "wishart(3,iden)" 
    ##               tau             delta 
    ##   "normal(0,1.5)" "gamma(1,.5)[sd]"

## Posterior summaries of MCMC draws for model parameters

### Posterior summaries: Agreeableness model parameters

``` r
a_mcmc %>% summarize_draws() %>% knitr::kable(digits = 2)
```

| variable      | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:--------------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| lambda1       | 0.34 |   0.34 | 0.05 | 0.05 | 0.26 | 0.42 |    1 |   2564.44 |   2371.41 |
| lambda2       | 0.67 |   0.67 | 0.06 | 0.06 | 0.58 | 0.76 |    1 |   2108.66 |   2037.32 |
| lambda3       | 0.63 |   0.63 | 0.05 | 0.05 | 0.55 | 0.72 |    1 |   2597.80 |   2251.53 |
| lambda4       | 0.78 |   0.78 | 0.06 | 0.06 | 0.69 | 0.88 |    1 |   2340.23 |   2174.12 |
| Theta\_cov\_r | 0.17 |   0.17 | 0.04 | 0.04 | 0.11 | 0.23 |    1 |   3400.70 |   1943.05 |
| Theta\_var1   | 0.55 |   0.55 | 0.03 | 0.04 | 0.49 | 0.61 |    1 |   3304.47 |   2541.51 |
| Theta\_var2   | 0.54 |   0.54 | 0.04 | 0.04 | 0.47 | 0.61 |    1 |   3191.95 |   2394.90 |
| Theta\_var3   | 0.44 |   0.44 | 0.04 | 0.04 | 0.38 | 0.50 |    1 |   3600.57 |   2352.12 |
| Theta\_var4   | 0.33 |   0.33 | 0.04 | 0.04 | 0.27 | 0.40 |    1 |   3015.95 |   2096.15 |
| nu1\_r        | 3.95 |   3.95 | 0.05 | 0.05 | 3.86 | 4.03 |    1 |   2648.27 |   2065.63 |
| nu2           | 3.56 |   3.56 | 0.06 | 0.06 | 3.46 | 3.66 |    1 |   1978.23 |   2234.96 |
| nu3           | 4.01 |   4.01 | 0.06 | 0.06 | 3.92 | 4.10 |    1 |   1948.02 |   1551.20 |
| nu4           | 3.61 |   3.61 | 0.06 | 0.06 | 3.51 | 3.71 |    1 |   1858.04 |   2071.19 |
| Theta\_cov\_f | 0.10 |   0.10 | 0.04 | 0.04 | 0.04 | 0.16 |    1 |   3365.72 |   2210.40 |
| Psi\_r        | 0.55 |   0.54 | 0.09 | 0.09 | 0.41 | 0.72 |    1 |   2452.56 |   2454.12 |
| nu1\_f        | 4.09 |   4.09 | 0.06 | 0.06 | 4.00 | 4.19 |    1 |   2630.36 |   2446.91 |
| alpha\_f      | 0.73 |   0.73 | 0.10 | 0.10 | 0.56 | 0.89 |    1 |   1733.31 |   2228.69 |

### Posterior summaries: Conscientiousness model parameters

``` r
c_mcmc %>% summarize_draws() %>% knitr::kable(digits = 2)
```

| variable    | mean | median |   sd |  mad |    q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:------------|-----:|-------:|-----:|-----:|------:|-----:|-----:|----------:|----------:|
| lambda1     | 0.58 |   0.57 | 0.06 | 0.06 |  0.47 | 0.69 |    1 |   1804.25 |   1989.16 |
| lambda2     | 0.35 |   0.35 | 0.04 | 0.04 |  0.28 | 0.42 |    1 |   2197.37 |   2093.40 |
| lambda3     | 0.38 |   0.38 | 0.05 | 0.04 |  0.31 | 0.46 |    1 |   2203.70 |   2312.89 |
| lambda4     | 0.79 |   0.78 | 0.08 | 0.08 |  0.66 | 0.92 |    1 |   1643.96 |   1610.98 |
| Theta\_var1 | 0.95 |   0.95 | 0.08 | 0.08 |  0.83 | 1.09 |    1 |   2458.12 |   1696.94 |
| Theta\_var2 | 0.41 |   0.41 | 0.03 | 0.03 |  0.36 | 0.46 |    1 |   3416.64 |   2255.82 |
| Theta\_var3 | 0.64 |   0.64 | 0.05 | 0.04 |  0.57 | 0.72 |    1 |   2882.98 |   2221.92 |
| Theta\_var4 | 0.88 |   0.88 | 0.10 | 0.10 |  0.70 | 1.05 |    1 |   2149.70 |   1812.53 |
| nu1         | 3.29 |   3.29 | 0.07 | 0.07 |  3.17 | 3.41 |    1 |   2799.61 |   2194.97 |
| nu2         | 4.26 |   4.26 | 0.04 | 0.04 |  4.20 | 4.33 |    1 |   2242.52 |   2386.79 |
| nu3         | 4.18 |   4.18 | 0.04 | 0.04 |  4.11 | 4.26 |    1 |   2428.48 |   2177.27 |
| nu4         | 3.60 |   3.60 | 0.07 | 0.07 |  3.48 | 3.72 |    1 |   1850.86 |   2236.79 |
| Psi\_f      | 1.44 |   1.41 | 0.28 | 0.27 |  1.04 | 1.94 |    1 |   1569.85 |   1677.73 |
| alpha\_f    | 0.13 |   0.13 | 0.12 | 0.12 | -0.07 | 0.34 |    1 |   1979.39 |   2370.30 |

### Posterior summaries: Extraversion model parameters

``` r
e_mcmc %>% summarize_draws() %>% knitr::kable(digits = 2)
```

| variable      |  mean | median |   sd |  mad |    q5 |   q95 | rhat | ess\_bulk | ess\_tail |
|:--------------|------:|-------:|-----:|-----:|------:|------:|-----:|----------:|----------:|
| lambda1       |  0.66 |   0.66 | 0.06 | 0.06 |  0.57 |  0.77 | 1.00 |   1936.44 |   2169.25 |
| lambda2       |  0.83 |   0.83 | 0.08 | 0.08 |  0.70 |  0.97 | 1.00 |   1165.14 |   1770.37 |
| lambda3       |  0.76 |   0.76 | 0.07 | 0.07 |  0.64 |  0.88 | 1.00 |   1677.85 |   2126.18 |
| lambda4       |  0.72 |   0.72 | 0.07 | 0.07 |  0.61 |  0.84 | 1.00 |   1069.11 |   1773.44 |
| Theta\_cov\_r | -0.15 |  -0.15 | 0.09 | 0.09 | -0.30 | -0.01 | 1.00 |   1111.32 |   1504.79 |
| Theta\_var1   |  0.79 |   0.78 | 0.06 | 0.06 |  0.68 |  0.89 | 1.00 |   1960.57 |   2334.87 |
| Theta\_var2   |  0.97 |   0.96 | 0.10 | 0.10 |  0.79 |  1.13 | 1.01 |    948.69 |   1225.74 |
| Theta\_var3   |  0.98 |   0.97 | 0.08 | 0.08 |  0.85 |  1.12 | 1.00 |   1568.20 |   1912.22 |
| Theta\_var4   |  0.58 |   0.59 | 0.07 | 0.07 |  0.46 |  0.70 | 1.00 |   1081.15 |   1403.30 |
| nu1\_r        |  2.31 |   2.31 | 0.07 | 0.07 |  2.20 |  2.43 | 1.00 |   2054.33 |   2198.00 |
| nu2           |  2.85 |   2.85 | 0.07 | 0.07 |  2.73 |  2.97 | 1.00 |   1256.22 |   1832.09 |
| nu3           |  3.18 |   3.18 | 0.07 | 0.07 |  3.06 |  3.29 | 1.00 |   1249.77 |   1936.65 |
| nu4           |  3.03 |   3.03 | 0.06 | 0.06 |  2.93 |  3.13 | 1.00 |   1248.97 |   1911.06 |
| Theta\_cov\_f | -0.11 |  -0.11 | 0.08 | 0.08 | -0.24 |  0.01 | 1.01 |    930.96 |    943.84 |
| Psi\_r        |  0.99 |   0.98 | 0.16 | 0.15 |  0.76 |  1.28 | 1.00 |   1782.45 |   1808.38 |
| nu1\_f        |  2.13 |   2.13 | 0.08 | 0.08 |  2.00 |  2.26 | 1.00 |    971.48 |   1348.81 |
| alpha\_f      |  0.14 |   0.13 | 0.11 | 0.11 | -0.03 |  0.31 | 1.01 |    960.58 |   1664.32 |

### Posterior summaries: Neuroticism model parameters

``` r
n_mcmc %>% summarize_draws() %>% knitr::kable(digits = 2)
```

| variable    | mean | median |   sd |  mad |    q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:------------|-----:|-------:|-----:|-----:|------:|-----:|-----:|----------:|----------:|
| lambda1     | 0.50 |   0.50 | 0.05 | 0.05 |  0.42 | 0.59 |    1 |   2707.92 |   2243.93 |
| lambda2     | 0.81 |   0.81 | 0.06 | 0.06 |  0.71 | 0.92 |    1 |   2124.16 |   2073.25 |
| lambda3     | 0.61 |   0.60 | 0.06 | 0.06 |  0.52 | 0.70 |    1 |   2140.82 |   2366.43 |
| lambda4     | 0.80 |   0.79 | 0.07 | 0.07 |  0.69 | 0.91 |    1 |   1720.86 |   2320.60 |
| Theta\_var1 | 0.83 |   0.83 | 0.06 | 0.06 |  0.74 | 0.93 |    1 |   3230.44 |   2180.44 |
| Theta\_var2 | 0.82 |   0.82 | 0.08 | 0.07 |  0.70 | 0.95 |    1 |   2743.43 |   2007.83 |
| Theta\_var3 | 0.73 |   0.73 | 0.05 | 0.05 |  0.65 | 0.83 |    1 |   2468.54 |   2135.84 |
| Theta\_var4 | 0.69 |   0.69 | 0.07 | 0.07 |  0.58 | 0.80 |    1 |   2293.02 |   1977.10 |
| nu1         | 2.26 |   2.26 | 0.07 | 0.07 |  2.15 | 2.37 |    1 |   2744.40 |   2457.59 |
| nu2         | 2.55 |   2.55 | 0.07 | 0.07 |  2.43 | 2.66 |    1 |   1412.72 |   1806.19 |
| nu3         | 2.29 |   2.29 | 0.06 | 0.06 |  2.19 | 2.38 |    1 |   1546.06 |   1460.41 |
| nu4         | 2.23 |   2.23 | 0.07 | 0.07 |  2.12 | 2.35 |    1 |   1498.60 |   1856.53 |
| Psi\_f      | 1.21 |   1.19 | 0.20 | 0.19 |  0.92 | 1.56 |    1 |   1912.83 |   2182.63 |
| alpha\_f    | 0.06 |   0.05 | 0.11 | 0.11 | -0.12 | 0.23 |    1 |   1530.98 |   1807.22 |

### Posterior summaries: Openness to Experience model parameters

``` r
o_mcmc %>% summarize_draws() %>% knitr::kable(digits = 2)
```

| variable       |  mean | median |   sd |  mad |    q5 |   q95 | rhat | ess\_bulk | ess\_tail |
|:---------------|------:|-------:|-----:|-----:|------:|------:|-----:|----------:|----------:|
| lambda1        |  0.59 |   0.58 | 0.18 | 0.22 |  0.31 |  0.88 | 1.02 |    117.62 |    159.56 |
| lambda2        |  0.46 |   0.40 | 0.17 | 0.11 |  0.28 |  0.86 | 1.03 |    102.49 |     79.12 |
| lambda3        |  0.51 |   0.46 | 0.18 | 0.11 |  0.33 |  0.95 | 1.04 |     97.56 |     74.96 |
| lambda4        |  0.64 |   0.62 | 0.19 | 0.23 |  0.35 |  0.95 | 1.02 |    116.79 |    178.82 |
| Theta\_cov\_r1 |  0.22 |   0.24 | 0.26 | 0.33 | -0.17 |  0.62 | 1.03 |    108.95 |    121.49 |
| Theta\_cov\_r2 |  0.43 |   0.50 | 0.25 | 0.17 | -0.20 |  0.71 | 1.04 |    105.76 |     69.32 |
| Theta\_var1    |  0.65 |   0.66 | 0.26 | 0.31 |  0.23 |  1.03 | 1.02 |    116.97 |    167.71 |
| Theta\_var2    |  1.16 |   1.22 | 0.24 | 0.16 |  0.58 |  1.43 | 1.03 |    113.78 |     83.10 |
| Theta\_var3    |  0.99 |   1.07 | 0.28 | 0.19 |  0.27 |  1.29 | 1.03 |    102.01 |     69.33 |
| Theta\_var4    |  0.70 |   0.72 | 0.29 | 0.35 |  0.22 |  1.15 | 1.03 |    110.37 |    131.63 |
| nu1\_r         |  3.93 |   3.93 | 0.07 | 0.07 |  3.83 |  4.04 | 1.00 |    949.35 |   1653.05 |
| nu2            |  3.59 |   3.58 | 0.07 | 0.07 |  3.48 |  3.70 | 1.00 |    525.44 |   1193.05 |
| nu3            |  3.61 |   3.61 | 0.07 | 0.07 |  3.50 |  3.72 | 1.00 |    489.81 |   1644.76 |
| nu4            |  4.07 |   4.07 | 0.07 | 0.08 |  3.95 |  4.19 | 1.01 |    292.95 |    643.83 |
| Theta\_cov\_f1 |  0.19 |   0.20 | 0.27 | 0.35 | -0.21 |  0.60 | 1.03 |    104.98 |     94.25 |
| Theta\_cov\_f2 |  0.54 |   0.62 | 0.26 | 0.17 | -0.11 |  0.81 | 1.04 |    104.88 |     69.06 |
| Psi\_r         |  1.65 |   1.56 | 0.40 | 0.34 |  1.15 |  2.40 | 1.01 |    540.23 |   1247.25 |
| nu1\_f         |  4.01 |   4.01 | 0.10 | 0.10 |  3.85 |  4.17 | 1.01 |    225.02 |    347.85 |
| alpha\_f       | -0.41 |  -0.40 | 0.16 | 0.15 | -0.69 | -0.18 | 1.01 |    409.66 |   1160.11 |

``` r
a_cai_summary <- summarize_draws(agree_post_CAI$overall)
a_cai_summary[,2:8] <- round(a_cai_summary[,2:8],2)
a_cai_summary[,9:10] <- round(a_cai_summary[,9:10])
a_cai_summary%>% knitr::kable()
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| PS       | 0.15 |   0.15 | 0.00 | 0.00 | 0.15 | 0.15 |    1 |      3119 |        NA |
| SR       | 0.57 |   0.58 | 0.05 | 0.05 | 0.49 | 0.65 |    1 |      2773 |      2852 |
| SE       | 0.57 |   0.58 | 0.05 | 0.05 | 0.49 | 0.65 |    1 |      2773 |      2852 |
| SP       | 0.92 |   0.92 | 0.01 | 0.01 | 0.91 | 0.94 |    1 |      2773 |      2852 |

``` r
c_cai_summary <- summarize_draws(consc_post_CAI$overall)
c_cai_summary[,2:8] <- round(c_cai_summary[,2:8],2)
c_cai_summary[,9:10] <- round(c_cai_summary[,9:10])
c_cai_summary%>% knitr::kable()
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| PS       | 0.15 |   0.15 | 0.00 | 0.00 | 0.15 | 0.15 |    1 |      2874 |        NA |
| SR       | 0.58 |   0.59 | 0.05 | 0.05 | 0.51 | 0.66 |    1 |      2423 |      2764 |
| SE       | 0.58 |   0.59 | 0.05 | 0.05 | 0.51 | 0.66 |    1 |      2423 |      2764 |
| SP       | 0.93 |   0.93 | 0.01 | 0.01 | 0.91 | 0.94 |    1 |      2423 |      2764 |

``` r
e_cai_summary <- summarize_draws(extr_post_CAI$overall)
e_cai_summary[,2:8] <- round(e_cai_summary[,2:8],2)
e_cai_summary[,9:10] <- round(e_cai_summary[,9:10])
e_cai_summary%>% knitr::kable()
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| PS       | 0.15 |   0.15 | 0.00 | 0.00 | 0.15 | 0.15 |    1 |      2882 |        NA |
| SR       | 0.66 |   0.66 | 0.04 | 0.03 | 0.60 | 0.73 |    1 |      2399 |      2849 |
| SE       | 0.66 |   0.66 | 0.04 | 0.03 | 0.60 | 0.73 |    1 |      2399 |      2849 |
| SP       | 0.94 |   0.94 | 0.01 | 0.01 | 0.93 | 0.95 |    1 |      2399 |      2849 |

``` r
n_cai_summary <- summarize_draws(neur_post_CAI$overall)
n_cai_summary[,2:8] <- round(n_cai_summary[,2:8],2)
n_cai_summary[,9:10] <- round(n_cai_summary[,9:10])
n_cai_summary%>% knitr::kable()
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| PS       | 0.15 |   0.15 | 0.00 | 0.00 | 0.15 | 0.15 |    1 |      2595 |        NA |
| SR       | 0.71 |   0.71 | 0.04 | 0.03 | 0.64 | 0.76 |    1 |      2526 |      1908 |
| SE       | 0.71 |   0.71 | 0.04 | 0.03 | 0.64 | 0.76 |    1 |      2526 |      1908 |
| SP       | 0.95 |   0.95 | 0.01 | 0.01 | 0.94 | 0.96 |    1 |      2526 |      1908 |

``` r
o_cai_summary <- summarize_draws(open_post_CAI$overall)
o_cai_summary[,2:8] <- round(o_cai_summary[,2:8],2)
o_cai_summary[,9:10] <- round(o_cai_summary[,9:10])
o_cai_summary%>% knitr::kable()
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| PS       | 0.15 |   0.15 | 0.00 | 0.00 | 0.15 | 0.15 |    1 |      2608 |        NA |
| SR       | 0.51 |   0.52 | 0.06 | 0.07 | 0.40 | 0.61 |    1 |       427 |       780 |
| SE       | 0.51 |   0.52 | 0.06 | 0.07 | 0.40 | 0.61 |    1 |       427 |       780 |
| SP       | 0.91 |   0.91 | 0.01 | 0.01 | 0.89 | 0.93 |    1 |       427 |       780 |

``` r
a_cai_summaryr <- summarize_draws(agree_post_CAI$reference)
a_cai_summaryr[,2:8] <- round(a_cai_summaryr[,2:8],2)
a_cai_summaryr[,9:10] <- round(a_cai_summaryr[,9:10])
a_cai_summaryr%>% knitr::kable()
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| PS       | 0.20 |   0.20 | 0.01 | 0.01 | 0.18 | 0.22 |    1 |      2830 |      2816 |
| SR       | 0.60 |   0.60 | 0.05 | 0.05 | 0.51 | 0.68 |    1 |      2698 |      2804 |
| SE       | 0.58 |   0.58 | 0.05 | 0.05 | 0.49 | 0.67 |    1 |      2745 |      2977 |
| SP       | 0.90 |   0.90 | 0.01 | 0.01 | 0.87 | 0.92 |    1 |      2756 |      2927 |

``` r
c_cai_summaryr <- summarize_draws(consc_post_CAI$reference)
c_cai_summaryr[,2:8] <- round(c_cai_summaryr[,2:8],2)
c_cai_summaryr[,9:10] <- round(c_cai_summaryr[,9:10])
c_cai_summaryr%>% knitr::kable()
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| PS       | 0.17 |   0.17 | 0.01 | 0.01 | 0.15 | 0.19 |    1 |      2406 |      2782 |
| SR       | 0.63 |   0.63 | 0.06 | 0.06 | 0.53 | 0.73 |    1 |      2481 |      2373 |
| SE       | 0.60 |   0.60 | 0.06 | 0.06 | 0.50 | 0.69 |    1 |      2483 |      2798 |
| SP       | 0.92 |   0.92 | 0.01 | 0.01 | 0.90 | 0.94 |    1 |      2908 |      2848 |

``` r
e_cai_summaryr <- summarize_draws(extr_post_CAI$reference)
e_cai_summaryr[,2:8] <- round(e_cai_summaryr[,2:8],2)
e_cai_summaryr[,9:10] <- round(e_cai_summaryr[,9:10])
e_cai_summaryr%>% knitr::kable()
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| PS       | 0.16 |   0.16 | 0.01 | 0.01 | 0.14 | 0.18 |    1 |      1899 |      2109 |
| SR       | 0.68 |   0.69 | 0.06 | 0.06 | 0.58 | 0.78 |    1 |      2524 |      2848 |
| SE       | 0.66 |   0.65 | 0.06 | 0.06 | 0.56 | 0.75 |    1 |      2273 |      2802 |
| SP       | 0.94 |   0.94 | 0.01 | 0.01 | 0.92 | 0.96 |    1 |      2818 |      2943 |

``` r
n_cai_summaryr <- summarize_draws(neur_post_CAI$reference)
n_cai_summaryr[,2:8] <- round(n_cai_summaryr[,2:8],2)
n_cai_summaryr[,9:10] <- round(n_cai_summaryr[,9:10])
n_cai_summaryr%>% knitr::kable()
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| PS       | 0.17 |   0.17 | 0.01 | 0.01 | 0.15 | 0.19 |    1 |      2675 |      2776 |
| SR       | 0.72 |   0.72 | 0.05 | 0.05 | 0.63 | 0.81 |    1 |      2598 |      2319 |
| SE       | 0.74 |   0.74 | 0.05 | 0.05 | 0.65 | 0.82 |    1 |      2755 |      3122 |
| SP       | 0.94 |   0.94 | 0.01 | 0.01 | 0.92 | 0.96 |    1 |      2708 |      2562 |

``` r
o_cai_summaryr <- summarize_draws(open_post_CAI$reference)
o_cai_summaryr[,2:8] <- round(o_cai_summaryr[,2:8],2)
o_cai_summaryr[,9:10] <- round(o_cai_summaryr[,9:10])
o_cai_summaryr%>% knitr::kable()
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| PS       | 0.14 |   0.14 | 0.01 | 0.01 | 0.12 | 0.16 |    1 |      2253 |      2532 |
| SR       | 0.52 |   0.52 | 0.09 | 0.10 | 0.37 | 0.67 |    1 |       476 |      1303 |
| SE       | 0.52 |   0.52 | 0.09 | 0.09 | 0.37 | 0.66 |    1 |       518 |      1156 |
| SP       | 0.92 |   0.92 | 0.02 | 0.02 | 0.89 | 0.95 |    1 |       533 |      1415 |

``` r
a_cai_summaryf <- summarize_draws(agree_post_CAI$focal)
a_cai_summaryf[,2:8] <- round(a_cai_summaryf[,2:8],2)
a_cai_summaryf[,9:10] <- round(a_cai_summaryf[,9:10])
a_cai_summaryf%>% knitr::kable()
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| PS       | 0.08 |   0.08 | 0.01 | 0.01 | 0.05 | 0.10 |    1 |      2868 |      2837 |
| SR       | 0.49 |   0.50 | 0.11 | 0.11 | 0.31 | 0.67 |    1 |      2928 |      2886 |
| SE       | 0.54 |   0.54 | 0.12 | 0.12 | 0.35 | 0.73 |    1 |      2969 |      2691 |
| SP       | 0.96 |   0.96 | 0.01 | 0.01 | 0.94 | 0.98 |    1 |      2734 |      2712 |

``` r
c_cai_summaryf <- summarize_draws(consc_post_CAI$focal)
c_cai_summaryf[,2:8] <- round(c_cai_summaryf[,2:8],2)
c_cai_summaryf[,9:10] <- round(c_cai_summaryf[,9:10])
c_cai_summaryf%>% knitr::kable()
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| PS       | 0.12 |   0.12 | 0.02 | 0.02 | 0.09 | 0.15 |    1 |      2394 |      2456 |
| SR       | 0.49 |   0.48 | 0.09 | 0.09 | 0.34 | 0.63 |    1 |      2869 |      2522 |
| SE       | 0.55 |   0.56 | 0.10 | 0.09 | 0.39 | 0.71 |    1 |      2442 |      2974 |
| SP       | 0.93 |   0.93 | 0.02 | 0.02 | 0.91 | 0.95 |    1 |      2932 |      2805 |

``` r
e_cai_summaryf <- summarize_draws(extr_post_CAI$focal)
e_cai_summaryf[,2:8] <- round(e_cai_summaryf[,2:8],2)
e_cai_summaryf[,9:10] <- round(e_cai_summaryf[,9:10])
e_cai_summaryf%>% knitr::kable()
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| PS       | 0.14 |   0.14 | 0.02 | 0.02 | 0.11 | 0.17 |    1 |      1934 |      2300 |
| SR       | 0.63 |   0.63 | 0.07 | 0.08 | 0.51 | 0.75 |    1 |      2879 |      2790 |
| SE       | 0.67 |   0.68 | 0.08 | 0.08 | 0.55 | 0.80 |    1 |      3074 |      2912 |
| SP       | 0.94 |   0.94 | 0.01 | 0.01 | 0.92 | 0.96 |    1 |      2559 |      2859 |

``` r
n_cai_summaryf <- summarize_draws(neur_post_CAI$focal)
n_cai_summaryf[,2:8] <- round(n_cai_summaryf[,2:8],2)
n_cai_summaryf[,9:10] <- round(n_cai_summaryf[,9:10])
n_cai_summaryf%>% knitr::kable()
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| PS       | 0.13 |   0.13 | 0.01 | 0.01 | 0.10 | 0.15 | 1.00 |      2626 |      2868 |
| SR       | 0.68 |   0.68 | 0.08 | 0.08 | 0.55 | 0.80 | 1.00 |      3059 |      2816 |
| SE       | 0.65 |   0.66 | 0.08 | 0.08 | 0.52 | 0.77 | 1.00 |      2863 |      3015 |
| SP       | 0.95 |   0.95 | 0.01 | 0.01 | 0.93 | 0.97 | 1.01 |      2762 |      2873 |

``` r
o_cai_summaryf <- summarize_draws(open_post_CAI$focal)
o_cai_summaryf[,2:8] <- round(o_cai_summaryf[,2:8],2)
o_cai_summaryf[,9:10] <- round(o_cai_summaryf[,9:10])
o_cai_summaryf%>% knitr::kable()
```

| variable | mean | median |   sd |  mad |   q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|-----:|-------:|-----:|-----:|-----:|-----:|-----:|----------:|----------:|
| PS       | 0.16 |   0.16 | 0.02 | 0.02 | 0.13 | 0.19 |    1 |      2267 |      2783 |
| SR       | 0.50 |   0.50 | 0.08 | 0.08 | 0.36 | 0.63 |    1 |       842 |      1241 |
| SE       | 0.50 |   0.50 | 0.08 | 0.08 | 0.36 | 0.63 |    1 |       922 |      1363 |
| SP       | 0.90 |   0.90 | 0.02 | 0.02 | 0.87 | 0.93 |    1 |      1088 |      2676 |

## Posterior distributions of classification accuracy indices

### Agreeableness, group: F

``` r
mcmc_areas(agree_post_CAI$reference, pars = c( "SR", "SE", "SP"))
```

![](573_projectcode_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

### Agreeableness, group: M

``` r
mcmc_areas(agree_post_CAI$focal, pars = c( "SR", "SE", "SP"))
```

![](573_projectcode_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

### Conscientiousness, group: F

``` r
mcmc_areas(consc_post_CAI$reference, pars = c( "SR", "SE", "SP"))
```

![](573_projectcode_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

### Conscientiousness, group: M

``` r
mcmc_areas(consc_post_CAI$focal, pars = c( "SR", "SE", "SP"))
```

![](573_projectcode_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Extraversion, group: F

``` r
mcmc_areas(extr_post_CAI$reference, pars = c( "SR", "SE", "SP"))
```

![](573_projectcode_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

### Extraversion, group: M

``` r
mcmc_areas(extr_post_CAI$focal, pars = c( "SR", "SE", "SP"))
```

![](573_projectcode_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

### Neuroticism, group: F

``` r
mcmc_areas(neur_post_CAI$reference, pars = c( "SR", "SE", "SP"))
```

![](573_projectcode_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### Neuroticism, group: M

``` r
mcmc_areas(neur_post_CAI$focal, pars = c( "SR", "SE", "SP"))
```

![](573_projectcode_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### Openness to Experience, group: F

``` r
mcmc_areas(open_post_CAI$reference, pars = c( "SR", "SE", "SP"))
```

![](573_projectcode_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

### Openness to Experience, group: M

``` r
mcmc_areas(open_post_CAI$focal, pars = c( "SR", "SE", "SP"))
```

![](573_projectcode_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## Posterior distribution of Cohen’s h computed for the difference in classification accuracy indices

``` r
# Agreeableness
agree_h <- summarize_draws(agree_post_CAI$h)
agree_h[,2:8] <- round(agree_h[,2:8],2)
agree_h[,9:10] <- round(agree_h[,9:10])
agree_h %>% knitr::kable()
```

| variable |  mean | median |   sd |  mad |    q5 |   q95 | rhat | ess\_bulk | ess\_tail |
|:---------|------:|-------:|-----:|-----:|------:|------:|-----:|----------:|----------:|
| h(PS)    |  0.37 |   0.36 | 0.08 | 0.07 |  0.23 |  0.51 |    1 |      2858 |      2816 |
| h(SR)    |  0.22 |   0.22 | 0.27 | 0.25 | -0.21 |  0.65 |    1 |      2693 |      2768 |
| h(SE)    |  0.08 |   0.08 | 0.28 | 0.28 | -0.39 |  0.55 |    1 |      3010 |      2742 |
| h(SP)    | -0.24 |  -0.24 | 0.09 | 0.09 | -0.39 | -0.09 |    1 |      2714 |      2841 |

``` r
# Conscientiousness
consc_h <- summarize_draws(consc_post_CAI$h)
consc_h[,2:8] <- round(consc_h[,2:8],2)
consc_h[,9:10] <- round(consc_h[,9:10])
consc_h %>% knitr::kable()
```

| variable |  mean | median |   sd |  mad |    q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|------:|-------:|-----:|-----:|------:|-----:|-----:|----------:|----------:|
| h(PS)    |  0.15 |   0.14 | 0.09 | 0.09 |  0.00 | 0.30 |    1 |      2406 |      2782 |
| h(SR)    |  0.30 |   0.30 | 0.25 | 0.25 | -0.11 | 0.70 |    1 |      2631 |      2851 |
| h(SE)    |  0.09 |   0.08 | 0.25 | 0.25 | -0.31 | 0.48 |    1 |      2428 |      2979 |
| h(SP)    | -0.03 |  -0.03 | 0.09 | 0.09 | -0.19 | 0.12 |    1 |      3169 |      3111 |

``` r
# Extraversion
extra_h <- summarize_draws(extr_post_CAI$h)
extra_h[,2:8] <- round(extra_h[,2:8],2)
extra_h[,9:10] <- round(extra_h[,9:10])
extra_h %>% knitr::kable()
```

| variable |  mean | median |   sd |  mad |    q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|------:|-------:|-----:|-----:|------:|-----:|-----:|----------:|----------:|
| h(PS)    |  0.04 |   0.04 | 0.08 | 0.09 | -0.08 | 0.19 |    1 |      1899 |      2349 |
| h(SR)    |  0.10 |   0.10 | 0.23 | 0.23 | -0.27 | 0.47 |    1 |      2749 |      2980 |
| h(SE)    | -0.04 |  -0.05 | 0.23 | 0.23 | -0.42 | 0.33 |    1 |      2720 |      2899 |
| h(SP)    |  0.00 |   0.00 | 0.09 | 0.09 | -0.15 | 0.15 |    1 |      2885 |      2740 |

``` r
# Neuroticism
neur_h <- summarize_draws(neur_post_CAI$h)
neur_h[,2:8] <- round(neur_h[,2:8],2)
neur_h[,9:10] <- round(neur_h[,9:10])
neur_h %>% knitr::kable()
```

| variable |  mean | median |   sd |  mad |    q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|------:|-------:|-----:|-----:|------:|-----:|-----:|----------:|----------:|
| h(PS)    |  0.12 |   0.12 | 0.07 | 0.06 |  0.00 | 0.25 |    1 |      2694 |      2776 |
| h(SR)    |  0.09 |   0.09 | 0.22 | 0.22 | -0.27 | 0.44 |    1 |      3212 |      2742 |
| h(SE)    |  0.19 |   0.19 | 0.22 | 0.21 | -0.18 | 0.56 |    1 |      2979 |      3110 |
| h(SP)    | -0.05 |  -0.05 | 0.09 | 0.09 | -0.20 | 0.10 |    1 |      2987 |      2784 |

``` r
# Openness to Experience
open_h <- summarize_draws(open_post_CAI$h)
open_h[,2:8] <- round(open_h[,2:8],2)
open_h[,9:10] <- round(open_h[,9:10])
open_h %>% knitr::kable()
```

| variable |  mean | median |   sd |  mad |    q5 |  q95 | rhat | ess\_bulk | ess\_tail |
|:---------|------:|-------:|-----:|-----:|------:|-----:|-----:|----------:|----------:|
| h(PS)    | -0.07 |  -0.06 | 0.09 | 0.09 | -0.20 | 0.08 |    1 |      2259 |      2419 |
| h(SR)    |  0.03 |   0.04 | 0.25 | 0.25 | -0.38 | 0.42 |    1 |      1090 |      1972 |
| h(SE)    |  0.05 |   0.04 | 0.22 | 0.23 | -0.31 | 0.41 |    1 |      1972 |      2290 |
| h(SP)    |  0.07 |   0.07 | 0.10 | 0.09 | -0.09 | 0.22 |    1 |      1698 |      2487 |

### Posterior distribution of Cohen’s h

``` r
par(mfrow=c(2,3))
# Agreeableness
mcmc_areas(agree_post_CAI$h)
```

![](573_projectcode_files/figure-gfm/h%20areas-1.png)<!-- -->

``` r
# Conscientiousness
mcmc_areas(consc_post_CAI$h)
```

![](573_projectcode_files/figure-gfm/h%20areas-2.png)<!-- -->

``` r
# Extraversion
mcmc_areas(extr_post_CAI$h)
```

![](573_projectcode_files/figure-gfm/h%20areas-3.png)<!-- -->

``` r
# Neuroticism
mcmc_areas(neur_post_CAI$h)
```

![](573_projectcode_files/figure-gfm/h%20areas-4.png)<!-- -->

``` r
# Openness to Experience
mcmc_areas(open_post_CAI$h)
```

![](573_projectcode_files/figure-gfm/h%20areas-5.png)<!-- -->
