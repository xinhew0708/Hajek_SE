# Run simulation for case 2
iter <- 10000
paramlist <- list(c(0,0), c(0.25,0), c(0.5,0), c(0.75,0), c(1,0),
                  c(0,0), c(0,0.25), c(0,0.5), c(0,0.75), c(0,1))
results <- NA
set.seed(20)

# Loop over all parameter combinations
for (B in c(2)){
  nblist <- list('1'=rep(10,B), '2'=rep(20,B), '3'=rep(50,B))
  nbtlist <- list('1'=rep(2,B), '2'=rep(5,B), '3'=rep(4,B),
                  '4'=rep(10,B), '5'=rep(10,B), '6'=rep(25,B))
  wslist <- list(rgamma(n=10*B, shape=1/0.5^2, rate=1/30/0.5^2), 
                 rgamma(n=20*B, shape=1/0.5^2, rate=1/30/0.5^2),
                 rgamma(n=50*B, shape=1/0.5^2, rate=1/30/0.5^2))
  yslist <- list(cbind(rnorm(10*B), rnorm(10*B)), 
                 cbind(rnorm(20*B), rnorm(20*B)), 
                 cbind(rnorm(50*B), rnorm(50*B)))
  for (j in 1:6){
    k <- ceiling(j/2)
    ws <- wslist[[k]]
    ys <- yslist[[k]]
    for (ab in paramlist){
      results <- rbind(
        results, 
        cbind(
          run_simulation(iter, nblist[[k]], nbtlist[[j]], ab[1], ab[2], ws, ys),
          rep(ab[1],2), rep(ab[2],2), rep(nblist[[k]][1],2), rep(B,2),
          rep(c("U","B","U","B","U","B")[j],2))
      )
    }
  }
}

write.csv(results[2:nrow(results),], paste0('case2results10000.csv'))
