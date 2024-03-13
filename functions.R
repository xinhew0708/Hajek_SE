library(estimatr)

generate_assign <- function(nbs = rep(2,10), nbt = rep(1,10)){
  # generate 0 1 treatment assignment
  # input: nb, nbt (number of treatments)
  z <- rep(0, sum(nbs))
  t <- 0
  for (b in 1:length(nbs)){
    tt <- t + nbs[b]
    z[sample((t + 1) : tt, nbt[b])] <- 1
    t <- tt
  }
  return(z)
}

.get_ms_contrast <- function(a,b){
  t = 0
  for (i in 1:length(a))
    t = t + sum((b - a[i])^2)
  return (t/length(a)/length(b))
}

# Function for calculate the estimated variances from a data frame
estimate_variances <- function(df){
  # calculate two variance estimators
  # input: df, a data, with weight, block, z, and y
  ws <- df$weight
  yobs <- df$y  # observed ys
  zobs <- df$z  # observed zs
  bid <- df$block # block ids
  
  nbk <- sapply(c(0,1), function(t) as.vector(table(bid[zobs == t])))
  # number of units in each block and treatment, B by K
  nbk_all <- nbk[bid, ]  # n by K
  B <- length(unique(bid))  # number of blocks
  K <- ncol(nbk)  # number of treatments
  
  gammas <- rowSums(nbk_all) * ws # gamma for variance estimation, n by K
  gamsbk <- list()  # s^2_b,j, sample variance
  
  nu1 <- matrix(nrow=B, ncol=K-1)  # nu1_b,0k
  nu2 <- matrix(nrow=B, ncol=K-1)  # nu2_b,0k
  
  for (k in 1:K){
    indk <- zobs == (k-1)
    ipwk <- (1 / nbk_all[,k] * rowSums(nbk_all))[indk]
    thetak <- sum(ws[indk] * yobs[indk] * ipwk) / sum(ws[indk] * ipwk)  # ratio estimate rho_z
    gammas[indk] <- gammas[indk] * (yobs[indk] - thetak)
    gamsbk[[k]] <- aggregate(gammas[indk], by = list(bid[indk]), FUN = var)
  }
  gamsbk <- merge(gamsbk[[1]], gamsbk[[2]], by = "Group.1")[,2:(K+1)]  # B by K
  gamsbk[is.na(gamsbk)] <- 0
  
  for (k in 2:K){
    for (b in 1:B){
      indb0 <- (zobs == 0) & (bid == b)
      indbk <- (zobs == k-1) & (bid == b)  # units in block b, treatment k
      
      nu2[b,k-1] <- .get_ms_contrast(gammas[indbk], gammas[indb0]) -
        sum((nbk[b,c(1,k)]-1) / nbk[b,c(1,k)] * gamsbk[b,c(1,k)])
      
      if (gamsbk[b,1] == 0 | gamsbk[b,k] == 0){ # small block
        nu1[b,k-1] <- nu2[b,k-1]
      } else { # large block
        nu1[b,k-1] <- sum(gamsbk[b,c(1,k)] / nbk[b,c(1,k)])
      }
    }
  }
  return(c(sum(nu1[,k-1]), sum(nu2[,k-1])) / sum(ws)^2)
}

# Function used for estimating CI
f_se_simu <- function(tau, df, q, est){
  # est: the estimator that is used for the SE, either 1 or 2
  # return the value of this function of tau:
  # (\hat\tau - tau) / se(tau) - q
  # input: df, a data, with weight, block, z, and y
  ws <- df$weight
  yobs <- df$y # observed ys
  zobs <- df$z # observed zs
  bid <- df$block # block ids
  ipw <- df$ipw # inverse probabilities
  
  nbk <- sapply(c(0,1), function(t) as.vector(table(bid[zobs == t])))
  # number of units in each block and treatment, B by K
  nb <- rowSums(nbk)
  B <- length(unique(bid)) # number of blocks
  
  y1 <- yobs*zobs + (yobs + tau)*(1-zobs)
  y0 <- yobs*(1-zobs) + (yobs - tau)*zobs
  theta <- c(sum(ws * y0) / sum(ws), sum(ws * y1) / sum(ws))
  
  gammas <- rowSums(nbk[bid, ]) * ws # gamma for variance estimation, n by K
  gamsbk <- list() # s^2_b,j, sample variance
  
  for (k in 1:2){
    indk <- zobs == (k-1)
    gammas[indk] <- gammas[indk] * (yobs[indk] - theta[k])
    gamsbk[[k]] <- aggregate(gammas[indk], by = list(bid[indk]), FUN = var)
  }
  gamsbk <- merge(gamsbk[[1]], gamsbk[[2]], by = "Group.1")[,2:3] # B by K
  gamsbk[is.na(gamsbk)] <- 0
  gamsb <- aggregate(gammas, by = list(bid), FUN = var)[,2] # B vector
  
  tauhat <- sum(ws * yobs * ipw * zobs) / sum(ws * ipw * zobs) - 
    sum(ws * yobs * ipw * (1-zobs)) / sum(ws * ipw * (1-zobs))
  if (est == 1 & nbk[B,2] > 1){
    if (nbk[1,2] > 1){
      se_tau <- sqrt(sum(gamsbk / nbk)) / sum(ws)
    } else {
      nu1 <- sum(gamsbk[(B/2+1):B,] / nbk[(B/2+1):B,]) # sum of nu1_b,0k
      ind <- 1:(B/2)
      nu1 <- nu1 + sum(2/nbk[ind,1]/nbk[ind,2] * choose(nb[ind],2) * gamsb[ind]) - 
        sum((1/nbk[ind,1] + 1/nbk[ind,2]) * 
              ((nbk[ind,1]-1) * gamsbk[ind,1] + (nbk[ind,2]-1) * gamsbk[ind,2]))
      se_tau <- sqrt(nu1) / sum(ws)
    }
    return((tauhat - tau) / se_tau - q)
  }
  nu2 <- sum(2 / nbk[,1] / nbk[,2] * choose(nb,2) * gamsb) - 
    sum((1/nbk[,1] + 1/nbk[,2]) * 
          ((nbk[,1]-1) * gamsbk[,1] + (nbk[,2]-1) * gamsbk[,2]))
  se_tau <- sqrt(nu2) / sum(ws)
  return((tauhat - tau) / se_tau - q)
} 

# Function to generate fixed data based on parameters
generate_fixed_data <- function(nb, B, a = 0, b = 0, ws, ys) {
  # generate weights, block ids, and potential outcomes 
  # input: block size, number of blocks, parameter alpha, beta
  n <- B * nb  # sample size
  #ws <- rgamma(n=n, shape=1/0.5^2, rate=1/30/0.5^2)
  
  alphai <- qnorm(1 - seq(1,n) / (n+1))
  alphab <- rep(qnorm(1 - seq(1,B) / (B+1)), each = nb) # between
  betab <- rep(qnorm(1 - seq(1,nb) / (nb+1)), B)  # within
  ys <- cbind(alphai + ys[,1], 5 + alphai + alphab * a + betab * b + ys[,2])  
  
  bid <- rep(1:B, each = nb)
  data <- data.frame(weight = ws, block = bid, y0 = ys[,1], y1 = ys[,2])
  return(data)
}

# Function to perform a single replication
single_replication <- function(data, nbs, nbt, q=qnorm(0.975)) {
  z <- generate_assign(nbs = nbs, nbt = nbt)
  data <- cbind(data, z = z)
  data$y <- (data$y1 * z + data$y0 * (1-z))
  data$ipw <- rep(nbs / nbt, each=nbs[1]) * z + 
    rep(nbs / (nbs-nbt), each=nbs[1]) * (1-z)
  
  ests <- estimate_variances(data)
  est3 <- lm_robust(y ~ z, data=data, weights=weight*ipw)
  #est4 <- estimate_variances2hat(data)
  
  tauh <- est3$coefficients[2]
  df <- sum(nbs) - 2
  lb1 <- uniroot(f_se_simu, c(tauh-5, tauh), df=data, q=q, est=1)$root
  ub1 <- uniroot(f_se_simu, c(tauh, tauh+5), df=data, q=-q, est=1)$root
  lb2 <- uniroot(f_se_simu, c(tauh-5, tauh), df=data, q=q, est=2)$root
  ub2 <- uniroot(f_se_simu, c(tauh, tauh+5), df=data, q=-q, est=2)$root
  return(list(est1 = ests[1], est2 = ests[2], est3 = est3$std.error[2]^2,
              hajek = tauh, lb1=lb1, ub1=ub1, lb2=lb2, ub2=ub2))
}

# Function to aggregate results from a matrix of estimators to 
# bias, relative bias, SD, CI length, CI coverage
aggregate_est_results <- function(hajek, estimates, tau, dof){
  n = length(hajek)
  vemp = var(hajek)  # empirical variances
  results <- matrix(nrow = 3, ncol = 7)
  q1 = qnorm(0.975)
  q2 = qt(0.975, dof) / q1
  for (j in 1:3){
    results[j,1] = mean(estimates[j,]) - vemp
    results[j,2] = results[j,1] / vemp * 100  # relative bias
    results[j,3] = sd(estimates[j,])  # standard deviation
    
    margin = q1 * sqrt(estimates[j,])
    results[j,4] = 2 * mean(margin)
    results[j,5] = sum(hajek - margin < tau
                       & hajek + margin > tau) / n * 100
    margin = q2 * margin
    results[j,6] = results[j,4] * q2
    results[j,7] = sum(hajek - margin < tau
                       & hajek + margin > tau) / n * 100
  }
  return(results)
}

# Function to aggregate profile CI results from a matrix to
# CI length, CI coverage (%)
aggregate_CI_results <- function(hajek, ci, tau){
  n = length(hajek)
  results <- matrix(nrow = 3, ncol = 2)
  for (j in 1:2){
    results[j,1] <- mean(ci[j*2,] - ci[j*2-1,])
    results[j,2] <- sum(ci[j*2-1,] < tau & tau < ci[j*2,]) / n * 100
  }
  return(results)
}


# Main simulation function
run_simulation <- function(n, nbs, nbt, aa, bb, ws, ys) {
  B <- length(nbs)
  data <- generate_fixed_data(nb = nbs[1], B = B, a = aa, b = bb, ws, ys)
  results <- replicate(n, single_replication(data, nbs, nbt), simplify = FALSE)
  
  hajek <- sapply(results, function(x) x$hajek)
  est1 <- sapply(results, function(x) x$est1)
  est2 <- sapply(results, function(x) x$est2)
  est3 <- sapply(results, function(x) x$est3)
  lb1 <- sapply(results, function(x) x$lb1)
  ub1 <- sapply(results, function(x) x$ub1)
  lb2 <- sapply(results, function(x) x$lb2)
  ub2 <- sapply(results, function(x) x$ub2)
  
  tau <- sum(data$weight * (data$y1 - data$y0)) / sum(data$weight)
  agg_res <- aggregate_est_results(hajek, rbind(est1,est2,est3),
                                   tau, sum(nbs)-2)
  agg_res_ci <- aggregate_CI_results(hajek, rbind(lb1, ub1, lb2, ub2), tau)
  return(cbind(agg_res, agg_res_ci))
}

                
