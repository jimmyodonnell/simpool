# generate a bunch of 'purified dna' samples
library(data.table)

functions.f <- list.files("functions", full.names = TRUE)
sapply(functions.f, source, echo = FALSE)

even.levs <- c('even', 'skew.lin', 'skew.low', 'skew.med', 'skew.hi')
N_even <- length(even.levs)

rich.levs <- c(10, 100, 1000)
N_rich <- length(rich.levs)

reps.each <- 10
reps.levs <- 1:reps.each

N_comm <- length(even.levs) * length(rich.levs)

N_templates <- 1e5

temp <- expand.grid(reps.levs, even.levs, rich.levs)
colnames(temp) <- c('rep', 'even', 'rich')
sample_dat <- data.table(sample = 1:nrow(temp), comm = rep(1:N_comm, each = reps.each), temp)
rm(temp)

# simulate template counts
template_dat <- sample_dat[ , 
  list(templates = sim_templates(N_sp = rich, N_out = N_templates, 
    evenness = even, stochastic = TRUE, sort = TRUE)), 
  by = sample]

template_dat[, species := seq_along(templates), by = sample]

template_dat <- merge(x = sample_dat, y = template_dat, by = 'sample')

template_dat
