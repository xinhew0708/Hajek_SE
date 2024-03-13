# Run simulation
library(estimatr)

# number of iterations
iter <- 10000
paramlist <- list(c(0,0), c(0.25,0), c(0.5,0), c(0.75,0), c(1,0),
                  c(0,0), c(0,0.25), c(0,0.5), c(0,0.75), c(0,1))
results <- NA
set.seed(20)

# Loop over all parameter combinations
for (B in c(10,50)){
  nblist <- list('1'=rep(2,B), '2'=rep(4,B), '3'=rep(4,B))
  nbtlist <- list('1'=rep(1,B), '2'=rep(2,B), '3'=c(rep(1,B/2), rep(2,B/2)))
  wslist <- list(rgamma(n=2*B, shape=1/0.5^2, rate=1/30/0.5^2), 
                 rgamma(n=4*B, shape=1/0.5^2, rate=1/30/0.5^2))
  yslist <- list(cbind(rnorm(2*B), rnorm(2*B)), 
                 cbind(rnorm(4*B), rnorm(4*B)))
  for (j in 1:3){
    k <- ceiling((j+1)/2)
    ws <- wslist[[k]]
    ys <- yslist[[k]]
    for (ab in paramlist){
      results <- rbind(
        results, 
        cbind(
          run_simulation(iter, nblist[[j]], nbtlist[[j]], ab[1], ab[2], ws, ys),
          rep(ab[1],3), rep(ab[2],3), rep(nblist[[j]][1],3), rep(B,3),
          rep(c("B","B","U")[j],3))
      )
    }
  }
}

write.csv(results[2:nrow(results),], paste0('case1results.csv'))
