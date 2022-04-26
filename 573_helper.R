# Helper functions for 573 project

# Function to compute classification accuracy indices
compute_cai <- function(eta_vector, z_vector, cut_eta, cut_z) {
  tab <- table(eta_vector > cut_eta, z_vector > cut_z)
  total <- tab[1, 1] + tab[1, 2] + tab[2, 1] + tab[2, 2]
  TP <- tab[2, 2] / total # A
  FP <- tab[1, 2] / total # B
  TN <- tab[1, 1] / total # C
  FN <- tab[2, 1] / total # D
  
  PS <- (TP + FP)/(TP + FP + TN + FN)
  SE <- TP / (TP + FN)
  SR <- TP / (TP + FP)
  SP <- TN / (TN + FP)
#  return(c("TP" = TP, "FP" = FP, "TN" = TN, "FN" = FN,
#           "PS" = PS, "SR" = SR, "SE" = SE, "SP" = SP))
  return(c("PS" = PS, "SR" = SR, "SE" = SE, "SP" = SP))
}

# Function that returns posterior distributions of the overall classification
# accuracy indices, the reference group classification indices, and the focal 
# group classification indices.
get_posterior_CAI <- function(latent, obs, propsel, groupingvector, ref, foc) {
  # df_overall <- df_foc <- df_ref <- data.frame(matrix(ncol= 8,
  #                                                     nrow = dim(latent)[1]))
  df_overall <- df_foc <- df_ref <- h <- data.frame(matrix(ncol= 4, 
                                                           nrow = dim(latent)[1]))
  
  propsel <- 1 - propsel
  
  for(i in seq_along(1:dim(latent)[1])) {
    eta_c <- quantile(latent[i,], propsel)
    z_c <- quantile(obs[i,], propsel)
    df_overall[i, ] <- compute_cai(latent[i,], obs[i,], eta_c, z_c)
    df_ref[i, ] <- compute_cai(latent[i,][groupingvector==ref],
                               obs[i,][groupingvector==ref], eta_c, z_c)
    df_foc[i, ] <- compute_cai(latent[i,][groupingvector==foc],
                               obs[i,][groupingvector==foc], eta_c, z_c)
    #h[i,] <- cohens_h(df_ref[i, 5:8], df_foc[i, 5:8])
    h[i,] <- cohens_h(df_ref[i,], df_foc[i,])
    
  }
  
 # cols <- c("TP", "FP", "TN", "FN",  "PS", "SR", "SE", "SP")
  cols <- c("PS", "SR", "SE", "SP")
  colnames(df_overall) <- cols; colnames(df_ref) <- cols; 
  colnames(df_foc) <- cols
  colnames(h) <- c("h(PS)", "h(SR)", "h(SE)", "h(SP)")
  return(list("overall" = df_overall, "reference" = df_ref, "focal" = df_foc, "h" = h))
}



#Computes Cohen's h (Cohen, 1988) for the difference in two 
# proportions using \eqn{h = 2arcsin(\sqrt{p1}) - 2arcsin(\sqrt{p2})}
cohens_h <- function(p1, p2) {
  h <- 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))
  return(h)
}

