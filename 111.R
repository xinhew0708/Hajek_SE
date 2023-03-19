library(MASS)

# --------------- generate data ---------------
set.seed(20)
nbs = c(rep(3, 10))  # number of units in each block
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

pi_all = matrix(c(rep(1/3,n), rep(1/3,n*(K-1))), nrow=n)
# propensity scores for all units, n by K
nbk_all = nbs[1] * pi_all  # nbk for all units, n by K
pi = matrix(c(rep(1/3,B), rep(1/3,B*(K-1))), nrow=B)
# propensity scores for each block and treatment, B by K
nbk = nbs[1] * pi  # nbk for each block and treatment, B by K

rhok = c(sum(ws * ys[,1]), sum(ws * ys[,2]), sum(ws * ys[,3]))  # ratio parameter
rhok = rhok / sum(ws)


# ---------- repeat sampling and estimation ----------
iter = 5000

theta = matrix(nrow=iter, ncol=K)
varest3 = matrix(nrow=iter, ncol=K-1)  # var est for small strata

covest1 = matrix(nrow=iter, ncol=(K-1)/2)  # baseline comparison
covest2 = matrix(nrow=iter, ncol=(K-1)/2)  # cov est for small strata

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
    assig = sample(c(1,2,3))
    zobs = c(zobs, assig)
    yobs = c(yobs, diag(ys[(sum(nbs[1:b])-nbs[b]+1):(sum(nbs[1:b])), assig]))
  }
  # Hajek estimators
  thetak = c()
  
  gammas = cbind(nbk_all[,1] * ws, nbk_all[,2] * ws, nbk_all[,3] * ws)
  # gamma (or xps and yps) n by K
  gamsbk = matrix(nrow=B, ncol=K)  # s^2 b,j
  #sigmab = matrix(nrow=B, ncol=K)  # hat sigma b,j
  
  nu3 = matrix(nrow=B, ncol=K-1)  # nu3_b,0k
  
  for (k in 1:K){
    indk = (zobs == k)
    tk = sum(ws[indk] * yobs[indk] / pi_all[indk, k]) /
      sum(ws[indk] / pi_all[indk, k])  # ratio estimates
    thetak = c(thetak, tk)
    
    gammas[indk,k] = gammas[indk,k] / pi_all[indk, k] * (yobs[indk] - thetak[k])
    
    for (b in 1:B){
      in_b = (sum(nbs[1:b])-nbs[b]+1):(sum(nbs[1:b]))  # indices of units in block b
      indbk = (zobs[in_b] == k)
      if (sum(indbk) == 1)
        gamsbk[b,k] = 0
      else
        gamsbk[b,k] = sd(gammas[in_b,k][indbk])^2
      
      #sigmab[b,k] = gamsbk[b,k] * (nbs[b] - nbs[b]*pi[b,k]) / nbs[b]^2 / pi[b,k]
      
      # variance estimators by block
      if (k > 1){
        indb0 = (zobs[in_b] == 1)
        
        nu3[b,k-1] = 
          mean_ss(gammas[in_b,k][indbk], gammas[in_b,1][indb0]) -
          (nbk[b,1]-1) / nbk[b,1] * gamsbk[b,1] -
          (nbk[b,k]-1) / nbk[b,k] * gamsbk[b,k]
      }
    }
    
    # variance estimators
    if (k > 1){
      varest3[i,k-1] = sum(nu3[,k-1]) / W^2
    }
    # covariance estimators
    if (k > 2){
      covest1[i,k-2] = sqrt(varest3[i,k-2] * varest3[i,k-1])
      covest2[i,k-2] = sum(sqrt(nu3[,1] * nu3[,2])) / W^2
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
CI_length3 = c(mean(qnorm(0.975) * sqrt(varest3[,1])),
               mean(qnorm(0.975) * sqrt(varest3[,2])))

# confidence intervals coverage
Coverage3 = c(
  sum(
    theta[,2] - theta[,1] - qnorm(0.975) * sqrt(varest3[,1]) < rhok[2] - rhok[1]
    &
      theta[,2] - theta[,1] + qnorm(0.975) * sqrt(varest3[,1]) > rhok[2] - rhok[1]
    ) / iter * 100,
  sum(
    theta[,3] - theta[,1] - qnorm(0.975) * sqrt(varest3[,2]) < rhok[3] - rhok[1]
    &
      theta[,3] - theta[,1] + qnorm(0.975) * sqrt(varest3[,2]) > rhok[3] - rhok[1]
    ) / iter * 100
  )

xtable(x = cbind(hajek_res, CI_length3, Coverage3),
       caption = "Simulation results",
       label = "tab:",
       align = rep("c", 6),
       digits = c(rep(4,5), 2)
)


# variance estimates
var1true = sd(theta[,2] - theta[,1])^2  # empirical variances ("true values")
var2true = sd(theta[,3] - theta[,1])^2 

est_infl = c(abs(colMeans(varest3)[1] - var1true) / var1true * 100,
             abs(colMeans(varest3)[2] - var2true) / var2true * 100)  # estimation inflation

var_std = c(sd(varest3[,1]),  # standard deviations of variance estimators
            sd(varest3[,2]))

xtable(x=cbind(est_infl, var_std, CI_length3, Coverage3),
       caption = "Simulation results",
       label = "tab:",
       align = rep("c", 5),
       digits = rep(2, 5))


# covariance estimates
library(ggplot2)
# empirical covariance ("true values")
covtrue = cov(theta[,2] - theta[,1], theta[,3] - theta[,1])

#covest2 = covest2 / mu^2

colMeans(covest1)
colMeans(covest2)

ggplot(data.frame(cov_estimator=1:2,
                  covariance=rep(covtrue,2),
                  lower=c(-mean(covest1[,1]), -mean(covest2[,1])),
                  upper=c(mean(covest1[,1]), mean(covest2[,1]))
), aes(cov_estimator, covariance)) +   # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper))

sd(covest1)
sd(covest2[,1])  # standard deviation of covariance estimators

# coverage rate
sum(-covest1[,1] < covtrue & covest1[,1] > covtrue) / iter * 100
sum(-covest2[,1] < covtrue & covest2[,1] > covtrue) / iter * 100

#sd(covest2[,1])^2 + abs(colMeans(covest2)[1] - covtrue)^2
#sd(covest2[,2])^2 + abs(colMeans(covest2)[2] - covtrue)^2
