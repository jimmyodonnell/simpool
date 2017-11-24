evenness.lev <- c('even', 'skew.lin', 'skew.low', 'skew.med', 'skew.hi')
even.col <- c(grey(0.7), 'chartreuse', 'gold', 'coral', 'mediumorchid')

p <- function(Y){ Y/sum(Y) }

# richness.lev <- c(10, 100, 1000)
N_sp <- 25

x <- 1:N_sp

par(mar = c(4,5,1,1))
plot(x, x, ylim = c(0, max(p(1/x))), 
     type = 'n', xlab = 'rank', ylab = 'p(S)', las = 1)

# even
y <- rep(1/N_sp, length(x))
points(x, p(y), type = 'b', col = even.col[1], lwd = 2)

# Here is a linearly decreasing option
B0 <- (N_sp+1)/N_sp
B1 <- -1/(2*N_sp)
y <- B0 + B1*x
points(x, p(y), type = 'b', col = even.col[2], lwd = 2)

# skew.low
Beta <- 1/(N_sp^0.5) # 0.8 = 0.125*10
y <- 1/((N_sp*Beta) + x)
points(x, p(y), type = 'b', col = even.col[3], lwd = 2)

# skew.med
Beta <- 1/(N_sp^1)# 0.01*length(x) # 0.1 = 0.01 * 10
y <- 1/((N_sp*Beta) + x)
points(x, p(y), type = 'b', col = even.col[4], lwd = 2)

# skew.hi
Beta <- 0 # same as: y <- 1/(x)
y <- 1/((N_sp*Beta) + x) # same as: y <- 1/(x)
points(x, p(y), type = 'b', col = even.col[5], lwd = 2)

legend('topright', 
       legend = evenness.lev, lwd = 2, col = even.col, pch = 1,
       bty = 'n')




################################################################################
# earlier mucking about
library(rmutil) # rmulti::rpareto

# columns in output: sample species templates

N_samples <- 8

richness.lev <- c(10, 100, 1000)

evenness.lev <- c("constant", "unif", "skew.low", "skew.med", "skew.hi")

x <- 1:100
y <- 1/(x)
plot(x, y, type = 'b', log = 'xy')


# even but varying only as draws from multinomial
# even but stochastic: 1, 1
hist(rbeta(100000, 1, 1))

# perfectly diagonal
hist(rbeta(100000, 1, 2))

hist(rbeta(100000, 1, 5), xlim = c(0,1))
hist(rbeta(100000, 0.3, 1), xlim = c(0,1))

hist(rbeta(100000, 1, 10), xlim = c(0,1))
rmultinom(1, 2, prob = rep(1, 6)/6)
probs <- rexp(20, 1)
probs <- rbeta(20, 2, 1)
hist(probs)
plot(probs/sum(probs))
hist(rmultinom(1, size = 1000, prob = probs), breaks = 20)
hist(rpareto(1000, 2, 2))

hist(rpareto(1000, 100, 2))


plotpareto <- function(mean, shape, type = 'comm'){
  dat <- replicate(100, sort(rpareto(20, mean, shape)))
  if(type == 'comm'){
    pldat <- dat
  }else if(type == 'sp'){
    pldat <- t(dat)
  }
  boxplot(pldat)
  # points(20, mean, col = "red", lwd = 2)
}
plotpareto(20, 6, type = 'sp')

plotmulti <- function(Nsp = 20, beta1 = 1, beta2 = 1, Nind = 1000, type = 'comm', ...){
  probs <- rbeta(Nsp, beta1, beta2)
  dat <- replicate(100, 
                   sort(rmultinom(n = 1, size = Nind, prob = 1/(1*(1:Nsp))))) # rbeta(Nsp, beta1, beta2)
  if(type == 'comm'){
    pldat <- dat
  }else if(type == 'sp'){
    pldat <- t(dat)
  }
  print(min(dat))
  par(mar = c(4,5,1,1))
  boxplot(pldat, ...)
  # points(20, mean, col = "red", lwd = 2)
}
par(mfrow = c(1,2))
plotmulti(beta2 = 100, type = 'sp', ylim = c(0, 500))
plotmulti(beta2 = 100, type = 'comm', ylim = c(0, 500))

N_species <- 10

communities <- replicate(100, simplify = 'matrix', 
                         sort(rmultinom(n = 1, size = 1000, prob = 1/(1:N_species) ))) # prob = rbeta(10, 1, 0)

boxplot(communities, ylim = c(0,500), outpch = 21, outcol=1:10, outbg=1:10)

# total templates per sample is constant
N_templates <- 1000



div_opt <- c("even", "exp")

div_scenario <- div_opt[1]

templates <- c()
if(div_scenario == "exp"){
  probs <- rexp(N_species, rate = 1)
  for(i in 1:N_samples){
    templates <- c(templates, rmultinom(1, N_templates, probs)[,1])
  }
}else if(div_scenario == "even"){
  # probs <- rep(1, N_species)
  for(i in 1:N_samples){
    templates <- c(templates, rep(N_templates/N_species, N_species))
    # templates <- c(templates, rmultinom(1, N_templates, probs)[,1])
  }
}

temp <- expand.grid(1:N_species, 1:N_samples)[,c(2,1)]
colnames(temp) <- c("sample", "species")
template_dat <- data.table(temp, templates = templates)
rm(temp)

template_dat

boxplot(templates ~ species, data = template_dat)
