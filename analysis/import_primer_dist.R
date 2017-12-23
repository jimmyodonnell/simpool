################################################################################
# get primer distances from real dataset

primer_opts <- c("Ford2016", "Miya2015")

primer_set <- primer_opts[2]

if(primer_set == "Ford2016"){
  primername_f <- "16s-prey-f-dg"
  primername_r <- "_R_16s-prey-r"
}else if(primer_set == "Miya2015"){
  primername_f <- "mifish-u-f"
  primername_r <- "_R_mifish-u-r"
}


# forward:
distfile_f <- paste0(
  "/Users/jimmy.odonnell/Desktop/mitofish/v3.28/alignments/primer_mismatch/", 
  primername_f, "/primer_dist_", primername_f, ".csv"
)

dist_f <- fread(distfile_f)

# reverse:
distfile_r <- paste0(
  "/Users/jimmy.odonnell/Desktop/mitofish/v3.28/alignments/primer_mismatch/", 
  primername_r, "/primer_dist_", primername_r, ".csv"
)

dist_r <- fread(distfile_r)

dist_b <- merge(dist_f[,.(seq_id, label, taxon, distf = dist)], dist_r[,.(seq_id, distr = dist)], by = 'seq_id')

# calculate average of F and R distances
dist_b[ , eff := 1-((distf + distr)/2) ]

plot(dist_b$eff)
plot(density(dist_b$eff, na.rm = TRUE, bw = 0.02), main = "")

hist(dist_b$eff, 
  xlab = 'primer efficiency', 
  main = paste0('primer set: ', primer_set))

plot(scale(dist_b$eff))
plot(density(scale(dist_b$eff), na.rm = TRUE, bw = 0.5), main = "")

hist(scale(dist_b$eff), prob = TRUE)
lines(density(scale(dist_b$eff), na.rm = TRUE, bw = 0.5), col = 4, lwd = 2)

mymu <- mean(dist_b$eff, na.rm = TRUE)
myvar <- var(dist_b$eff, na.rm = TRUE)
estBetaParams <- function(mu, var){
  # estimate parameters of beta distribution given mean and variance
  # https://stats.stackexchange.com/a/12239
  alpha <- ((1 - mu) / var - 1/mu) * mu^2
  beta  <- alpha * (1 / mu - 1)
  return(c(alpha = alpha, beta = beta))
}
myparams <- estBetaParams(mymu, myvar)

(param.scale <- sum(myparams))
myparams/param.scale

hist(rbeta(2500, myparams[1], myparams[2]))