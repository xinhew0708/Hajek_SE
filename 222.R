library(MASS)

# --------------- generate data ---------------
set.seed(20)
nbs = c(rep(6, 10))  # number of units in each block
n = sum(nbs)  # sample size
B = length(nbs)  # number of blocks

ws = rnorm(n=n, mean=50, sd=10)
#ws = B * ws / sum(ws)  # weights standardized to summation = B
W = sum(ws)
mu = sum(ws) / W

K = 3  # number of treatment arms
r = 0.5  # correlation of potential outcomes
ys = mvrnorm(n=n, mu=rep(0,K), Sigma=(diag(1-r, K) + r)*10)  # potential outcomes

a = 1
b = 1
c = 1
alphab = qnorm(1 - seq(1, B) / (B+1)) * a
betab = 5 + qnorm(1 - seq(1, B) / (B+1)) * b
betab2 = 8 + qnorm(1 - seq(1, B) / (B+1)) * c
ys = ys + cbind(rep(alphab, each=nbs[1]), 
                rep(betab, each=nbs[1]),
                rep(betab2, each=nbs[1]))  # adjust the expectation of y

pi_all = matrix(rep(1/K, n*K), nrow=n)  # propensity scores, n by K
nbk_all = 6 * pi_all
pi = matrix(rep(1/K, B*K), nrow=B)  # propensity scores, B by K
nbk = 6 * pi

rhok = c(sum(ws * ys[,1]), sum(ws * ys[,2]), sum(ws * ys[,3]))  # ratio parameter
rhok = rhok / sum(ws)


# ---------- repeat sampling and estimation ----------
iter = 5000

theta = matrix(nrow=iter, ncol=K)
varest1 = matrix(nrow=iter, ncol=K-1)
varest3 = matrix(nrow=iter, ncol=K-1)

covest1 = matrix(nrow=iter, ncol=(K-1)/2)
covest22 = matrix(nrow=iter, ncol=(K-1)/2)  # the new cov est for small strata 
covest2 = matrix(nrow=iter, ncol=(K-1))  
covest3 = matrix(nrow=iter, ncol=(K-1)) 
#covest4 = matrix(nrow=iter, ncol=(K-1))
covest5 = matrix(nrow=iter, ncol=(K-1))  # the new cov est for large strata


mean_ss <- function(a,b){
  # calculate the cross term difference squared
  t = 0
  for (i in 1:length(a)){
    for (j in 1:length(b)){
      t = t + (a[i] - b[j])^2
    }
  }
  return (t/length(a)/length(b))
}


scd_term <- function(a,nb){
  # calculate the second term in the covariance est 3
  t = 0
  nbl = length(a)
  for (i in 1:nbl){
    t = t + a[i]^2
    for (j in 1:nbl){
      if (j != i){
        t = t + (nb-1) / (nbl-1) * a[i]*a[j]
      }
    }
  }
  return (t * nb / nbl)
}

thd_term_1 <- function(a,b,pil,pim,nb){
  # calculate the third term in the covariance est 3
  t = 0
  nbl = length(a)
  nbm = length(b)
  for (i in 1:nbl){
    for (j in 1:nbm)
      t = t + a[i] * b[j]
  }
  t = nb * (nb-1) / nbl / nbm * (nbl / nb) * nbm / (nb-1) / pil / pim * t
  return (t)
}

thd_term <- function(a,b,c,pi1,pi2,pi3,nb){
  # calculate the third term in the covariance est 3
  t = - thd_term_1(a,b,pi1,pi2,nb) - thd_term_1(a,c,pi1,pi3,nb) + 
    thd_term_1(b,c,pi2,pi3,nb)
  return (t)
}


for (i in 1:iter){
  yobs = c()  # observed y's
  zobs = c()  # treatment assignment, ctl 1, trt 2,3
  for (b in 1:B){
    assig = sample(c(1,1,2,2,3,3))
    zobs = c(zobs, assig)
    yobs = c(yobs, diag(ys[(sum(nbs[1:b])-nbs[b]+1):(sum(nbs[1:b])), assig]))
    
  }
  # Hajek estimators
  thetak = c()
  
  gammas = cbind(nbk_all[,1] * ws, nbk_all[,2] * ws, nbk_all[,3] * ws) / mu
  # gamma (or xps and yps) n by K
  gamsbk = matrix(nrow=B, ncol=K)  # s^2 b,j
  sigmab = matrix(nrow=B, ncol=K)  # hat sigma b,j
  
  nu1 = matrix(nrow=B, ncol=K-1)  # nu1_b,0k
  nu3 = matrix(nrow=B, ncol=K-1)  # nu3_b,0k
  
  gammas_cov = cbind(ws, ws, ws) / mu  # gamma' for covariance estimation
  bd3 = matrix(nrow=B, ncol=K-2)  # bound 2 0,k,j
  bd2 = matrix(nrow=B, ncol=2*(K-2))  # bound 2.2, 0,k,j
  #bd4 = matrix(nrow=B, ncol=2*(K-2))  # bound 3 0,k,j
  bd5 = matrix(nrow=B, ncol=2*(K-2))  # bound for large strata
  
  for (k in 1:K){
    indk = (zobs == k)
    tk = sum(ws[indk] * yobs[indk] / pi_all[indk, k]) /
      sum(ws[indk] / pi_all[indk, k])
    thetak = c(thetak, tk)
    
    gammas[indk] = gammas[indk] / pi_all[indk, k] * (yobs[indk] - thetak[k])
    gammas_cov[indk] = gammas_cov[indk] * (yobs[indk] - thetak[k])
    
    for (b in 1:B){
      in_b = (sum(nbs[1:b])-nbs[b]+1):(sum(nbs[1:b]))  # indices of units in block b
      indbk = (zobs[(nbs[1]*b-5):(nbs[1]*b)] == k)
      gamsbk[b,k] = sd(gammas[in_b,k][indbk])^2
      sigmab[b,k] = gamsbk[b,k] * (nbs[b] - nbs[b]*pi[b,k]) / nbs[b]^2 / pi[b,k]
      
      # variance estimators by block
      if (k > 1){
        indb0 = (zobs[(nbs[1]*b-5):(nbs[1]*b)] == 1)
        nu1[b,k-1] = nbs[b] / (nbs[b] - nbk[b,1]) * sigmab[b,1] +
          nbs[b] / (nbs[b] - nbk[b,k]) * sigmab[b,k]
        
        nu3[b,k-1] = 
          mean_ss(gammas[in_b,k][indbk], gammas[in_b,1][indb0]) - 
          (nbk[b,1]-1) / nbk[b,1] * gamsbk[b,1] - 
          (nbk[b,k]-1) / nbk[b,k] * gamsbk[b,k] 
      }
      
      # covariance estimators by block
      if (k > 2){
        indbj = (zobs[(nbs[1]*b-5):(nbs[1]*b)] == (k-1))
        
        bd3[b,k-2] = sqrt(sigmab[b,1]*sigmab[b,k-1]) + 
          sqrt(sigmab[b,1]*sigmab[b,k]) + sqrt(sigmab[b,k-1]*sigmab[b,k])
        
        bd2[b,k-2] = - 
          mean_ss(gammas_cov[in_b,k][indbk]/pi[b,k],
                  gammas_cov[in_b,1][indb0]/pi[b,1]) /
          (2 + 2 / (nbs[b] - 1)) -
          mean_ss(gammas_cov[in_b,k-1][indbj]/pi[b,k-1],
                  gammas_cov[in_b,1][indb0]/pi[b,1]) /
          (2 + 2 / (nbs[b] - 1)) - 
          mean_ss(gammas_cov[in_b,k][indbk]/pi[b,k],
                  gammas_cov[in_b,k-1][indbj]/pi[b,k-1]) /
          (2 - 2 / (nbs[b] - 1))
        
        bd2[b,k-1] =  
          mean_ss(gammas_cov[in_b,k][indbk]/pi[b,k],
                  gammas_cov[in_b,1][indb0]/pi[b,1]) /
          (2 - 2 / (nbs[b] - 1)) +
          mean_ss(gammas_cov[in_b,k-1][indbj]/pi[b,k-1],
                  gammas_cov[in_b,1][indb0]/pi[b,1]) /
          (2 - 2 / (nbs[b] - 1)) + 
          mean_ss(gammas_cov[in_b,k][indbk]/pi[b,k],
                  gammas_cov[in_b,k-1][indbj]/pi[b,k-1]) /
          (2 + 2 / (nbs[b] - 1))
        
        #bd4[b,k-2] = 
        #  scd_term(gammas_cov[in_b,1][indb0], nbk[b,1] + nbk[b,k]) +
        #  scd_term(gammas_cov[in_b,k-1][indbj], nbk[b,1] + nbk[b,k]) +
        #  scd_term(gammas_cov[in_b,k][indbk], nbk[b,1] + nbk[b,k])
        #bd4[b,k-1] = thd_term(
        #  gammas_cov[in_b,1][indb0],
        #  gammas_cov[in_b,k-1][indbj],
        #  gammas_cov[in_b,k][indbk],
        #  pi[b,1],
        #  pi[b,k-1],
        #  pi[b,k],
        #  nbk[b,1] + nbk[b,k]
        #)
        
        bd5[b,k-2] = sigmab[b,1] - 
          2/nbs[b] * sum(sqrt(gamsbk[b,] * gamsbk[b,c(2,3,1)]))
        bd5[b,k-1] = sigmab[b,1] + 
          2/nbs[b] * sum(sqrt(gamsbk[b,] * gamsbk[b,c(2,3,1)]))
      }
    }
    
    # variance estimators
    if (k > 1){
      varest1[i,k-1] = sum(nu1[,k-1]) / W^2
      varest3[i,k-1] = sum(nu3[,k-1]) / W^2
    }
    # covariance estimators
    if (k > 2){
      covest1[i,k-2] = sqrt(varest1[i,k-2] * varest1[i,k-1])
      
      covest22[i,k-2] = sum(sqrt(nu3[,1] * nu3[,2])) / W^2
      
      covest3[i,k-2] = (sum(sigmab[,1]) - sum(bd3[,k-2])) / W^2
      covest3[i,k-1] = (sum(sigmab[,1]) + sum(bd3[,k-2])) / W^2
      
      covest2[i,k-2] = (sum(sigmab[,1]) + sum(bd2[,k-2])) / W^2
      covest2[i,k-1] = (sum(sigmab[,1]) + sum(bd2[,k-1])) / W^2
      
      #covest4[i,k-2] = mean(sigmab[,1]) - mean(bd3[,k-2]) + mean(bd4[,k-1])
      #covest4[i,k-1] = mean(sigmab[,1]) + mean(bd3[,k-2]) + mean(bd4[,k-1])
      
      covest5[i,k-2] = sum(bd5[,k-2]) / W^2
      covest5[i,k-1] = sum(bd5[,k-1]) / W^2
    }
    
  }
  
  theta[i,] = thetak
  
}



# --------------- simulation results ---------------
library(xtable)
# point estimate (Hajek estimator)
rhok  # true ratios
colMeans(theta)
bias = c(mean(theta[,2] - theta[,1]) - (rhok[2] - rhok[1]),
         mean(theta[,3] - theta[,1]) - (rhok[3] - rhok[1]))
std = c(sd(theta[,2] - theta[,1]), sd(theta[,3] - theta[,1]))
rmse = sqrt(bias^2 + std^2)

hajek_res = cbind(bias, std, rmse)  # results, 3 x 2

# confidence interval lengths
CI_length = c(mean(qnorm(0.975) * sqrt(varest1[,1])),
              mean(qnorm(0.975) * sqrt(varest1[,2])))
CI_length2 = c(mean(qnorm(0.975) * sqrt(varest3[,1])),
               mean(qnorm(0.975) * sqrt(varest3[,2])))

# confidence intervals coverage
Coverage = c(
  sum(theta[,2] - theta[,1] - qnorm(0.975) * sqrt(varest1[,1]) < rhok[2] - rhok[1]
      &
        theta[,2] - theta[,1] + qnorm(0.975) * sqrt(varest1[,1]) > rhok[2] - rhok[1]
  ) / iter * 100,
  sum(theta[,3] - theta[,1] - qnorm(0.975) * sqrt(varest1[,2]) < rhok[3] - rhok[1]
      & 
        theta[,3] - theta[,1] + qnorm(0.975) * sqrt(varest1[,2]) > rhok[3] - rhok[1]
  ) / iter * 100)

Coverage2 = c(
  sum(theta[,2] - theta[,1] - qnorm(0.975) * sqrt(varest3[,1]) < rhok[2] - rhok[1]
      &
        theta[,2] - theta[,1] + qnorm(0.975) * sqrt(varest3[,1]) > rhok[2] - rhok[1]
  ) / iter * 100,
  sum(theta[,3] - theta[,1] - qnorm(0.975) * sqrt(varest3[,2]) < rhok[3] - rhok[1]
      & 
        theta[,3] - theta[,1] + qnorm(0.975) * sqrt(varest3[,2]) > rhok[3] - rhok[1]
  ) / iter * 100)

xtable(x = rbind(
  cbind(matrix(rep(hajek_res, each=2), nrow=4)[1:2,], CI_length, Coverage), 
  cbind(matrix(rep(hajek_res, each=2), nrow=4)[3:4,], CI_length2, Coverage2)),
  caption = "Simulation results",
  label = "tab:",
  align = rep("c", 6),
  digits = c(rep(4,5), 2)
)


# variance estimates
var1true = sd(theta[,2] - theta[,1])^2  # empirical variances
var2true = sd(theta[,3] - theta[,1])^2

est_infl = c(abs(colMeans(varest1)[1] - var1true) / var1true * 100,
             abs(colMeans(varest1)[2] - var2true) / var2true * 100,
             abs(colMeans(varest3)[1] - var1true) / var1true * 100,
             abs(colMeans(varest3)[2] - var2true) / var2true * 100)  # estimation inflation

var_std = c(sd(varest1[,1]),  # standard deviations of variance estimators
            sd(varest1[,2]),
            sd(varest3[,1]),
            sd(varest3[,2]))

xtable(x=cbind(est_infl, var_std, 
               c(CI_length, CI_length2), c(Coverage, Coverage2)),
       caption = "Simulation results",
       label = "tab:",
       align = rep("c", 5),
       digits = rep(2, 5))


# covariance estimates
library(ggplot2)
covtrue = cov(theta[,2] - theta[,1], theta[,3] - theta[,1])
# empirical covariance

colMeans(covest1)
colMeans(covest22)
colMeans(covest2)
colMeans(covest3)
colMeans(covest5)
#colMeans(covest4)

ggplot(data.frame(cov_estimator=1:5,
                  covariance=rep(covtrue,5),
                  lower=c(-mean(covest1[,1]), -mean(covest22[,1]), mean(covest2[,1]), 
                          mean(covest3[,1]), mean(covest5[,1])),
                  upper=c(mean(covest1[,1]), mean(covest22[,1]), mean(covest2[,2]), 
                          mean(covest3[,2]), mean(covest5[,2]))
), aes(cov_estimator, covariance)) +   # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  ylim(-2,2)

sd(covest1)
sd(covest22[,1])  # standard deviation of covariance estimators
sd(covest3[,1])
sd(covest3[,2])
sd(covest5[,1])
sd(covest5[,2])

# coverage rate
sum(-covest1[,1] < covtrue & covest1[,1] > covtrue) / iter * 100
sum(-covest22[,1] < covtrue & covest22[,1] > covtrue) / iter * 100
sum(covest2[,1] < covtrue & covest2[,2] > covtrue) / iter * 100
sum(covest3[,1] < covtrue & covest3[,2] > covtrue) / iter * 100
#sum(covest4[,1] < covtrue & covest4[,2] > covtrue) / iter * 100
sum(covest5[,1] < covtrue & covest5[,2] > covtrue) / iter * 100

