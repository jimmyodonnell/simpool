# generate a bunch of 'purified dna' samples
library(data.table)

functions.f <- list.files("functions", full.names = TRUE)
sapply(functions.f, source, echo = FALSE)

even.levs <- c('even', 'skew.lin', 'skew.low', 'skew.med', 'skew.hi')
N_even <- length(even.levs)

rich.levs <- c(10, 100, 1000)

reps.each <- 10
reps.levs <- 1:reps.each

N_comm <- length(even.levs) * length(rich.levs)


for(i in 1:length(rich.levs)){
  
}

template_dat <- data.table(
  # comm = rep(1:N_even, each = rich.levs[1]*reps.each), 
  even = rep(even.levs, each = rich.levs[1]), # rich.levs[1]*reps.each
  # rich = rep(rich.levs[1], rich.levs[1]*reps.each*N_even), 
  sample = rep(1:(N_even*rich.levs[1]), each = reps.each), 
  species = rep(1:rich.levs[1], N_even*reps.each), 
  count = NA_integer_
)[1:100,]

N_even*rich.levs[1]

for(i in template_dat$sample){
  
}
template_dat[sample == 1, count] <- 
sim_templates(N_sp = rich.levs[1], N_out = 1e6, evenness = even.levs[1], stochastic = TRUE)

temp <- expand.grid(reps.levs, even.levs, rich.levs)
colnames(temp) <- c('rep', 'even', 'rich')
template_dat <- data.table(comm = rep(1:N_comm, each = reps.each), sample = 1:nrow(temp), temp, templates = NA_integer_)
rm(temp)
# Dt[sample == I, unique(evenness)]

template_dat

