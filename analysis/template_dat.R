################################################################################
# generate a bunch of 'purified dna' samples
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
