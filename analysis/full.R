library(data.table)

# 1. We have a pool of different types of organisms. We'll classify them into 
# groups based on the identity at a given region of the mitochondrial genome. 
# Thus, each 'variant' of that sequence corresponds to a type (aka OTU).
N_seq_var <- 50
seq_var <- 1:N_seq_var
dim_pool <- c(50, 25, 3)

# Each type of organism occurs in different abundances. 
# We could just simulate some random numbers of abundances, 
# or be more realistic and incorporate information about the mean body size of 
# each type (e.g. smaller organisms are more abundant). 
# To do that, we need to pick some mean body sizes for each of our 
# organism types. 
set.seed(2017)
# choose between "realistic", "uniform", or give a single number for all. 
BODYSIZE_DIST <- "realistic"
if(BODYSIZE_DIST == "realistic"){
  bodysize_mean <- rlnorm(N_seq_var)
} else if(BODYSIZE_DIST == "uniform"){
  bodysize_mean <- runif(n = N_seq_var, min = 0, max = 20)
} else if(is.numeric(BODYSIZE_DIST)){
  bodysize_mean <- rep(BODYSIZE_DIST, N_seq_var)
}
par(mar = c(4,4,1,1))
plot(bodysize_mean, xlab = 'organism type', ylab = 'mean body mass (g)', las = 1)
# OK. Lots of small things, a few big ones.
# Now, let's simulate counts of those things based on their body size. 
# First we'll establish a relationship between body mass and expected counts. 
COEFF <- 1000
counts_mean <- COEFF/bodysize_mean
plot(bodysize_mean, counts_mean, 
     col = grey(0.5, alpha = 0.5), lwd = 2, 
     log = 'y', 
     xlab = 'body mass (g)', ylab = 'expected mean number of individuals', 
     las = 1)
pts.df <- data.frame(x = range(bodysize_mean), y = rev(range(counts_mean)))
text(pts.df, labels = round(pts.df$y), col = 'darkorchid', pos = c(4,3), font = 2)
pts.df

# The smallest 'species' has a mean body mass of 0.14 g and around 7000 individuals.
# The largest has a mean body mass of 23 g and around 43 individuals.
# From this information we can generate some simulated but realistic count data. 
ind_per_type <- rpois(n = length(counts_mean), lambda = counts_mean)

# And to bring it back around, we can simulate a mass for each individual
realsd <- 0.00001/bodysize_mean
logmean <- log(bodysize_mean/sqrt(1+(sqrt(realsd)/(bodysize_mean ^2))))
logsd <- log(1+((sqrt(realsd)/(bodysize_mean ^2))))
mass_ind <- rlnorm(
  n = sum(ind_per_type), 
  meanlog = rep(logmean, times = ind_per_type), 
  sdlog = rep(logsd, times = ind_per_type)
  )
# SCALE <- 1
# SHAPE <- bodysize_mean/SCALE
# SHAPES <- rep(SHAPE, times = ind_per_type)
# mass_ind <- rgamma(n = sum(ind_per_type), shape = SHAPES, scale = SCALE)
inds <- data.table(
  # ID = 1:length(mass_ind), 
  type = rep(seq_var, times = ind_per_type), 
  mass = mass_ind
)
boxplot(split(mass_ind, f = inds$type), log = 'y', cex = 0.5)
points(x = 1:N_seq_var, y = bodysize_mean, col = 'red')

# I'm going to call place_inds() multiple times so each type has a different 
# distribution in space. 
inds[,c('x', 'y', 'z') := data.table(place_inds(.N, dim_pool)), by = type]

ncols <- 3
mycols <- hsv(h = seq(from = 1/ncols, to = 1, length.out = ncols), s = 0.7, alpha = 0.7)
# mycols <- colorRampPalette(c("gold", "slateblue", 'green'))(ncols)
# mycols <- viridis(3, end = 1, option = 'magma')
with(inds[type == 2, ], 
     plot(x,z, col = mycols[cut(z, breaks = ncols)], lwd = 2, las = 1)
     )
legend('bottomright', legend = c('upper', 'middle', 'lower'), 
       pch = 1, pt.lwd = 2, col = rev(mycols), bty = 'n')
with(inds[type == 2, ], 
     plot(x,y, col = mycols[cut(z, breaks = ncols)], lwd = 2, las = 1)
)
legend('bottomright', legend = c('upper', 'middle', 'lower'), pch = 1, pt.lwd = 2, col = rev(mycols), bty = 'n')

# NOW we can sample from the environment.
# Let's start by deciding where to take a sample. 
corner_min <- c(0,0,0) # in the corner
pool_range <- cbind(corner_min, dim_pool)
sample_points <- as.matrix(expand.grid(pool_range[1,], pool_range[2,], pool_range[3,]))

# In the following lines I'll be doing things in an inefficient way, but in the 
# hopes it will illustrate each step individually.
# First, calculate the distance that each sample is from each individual.
for(i in 1:nrow(sample_points)){
  inds[ , paste0('dist.', i)] <- dist_from_samp(sample_points[i,], as.matrix(inds[,.(x,y,z)]))
}

# Then, sample cells from the environment based on density calculated based on 
# mass and distance. 
for(i in 1:nrow(sample_points)){
  distvar <- paste0('dist.', i)
  cellvar <- paste0('cells.', i)
  inds[ , (cellvar)] <- cells_by_dist(inds[[distvar]], inds[,mass])
}

temp <- melt.data.table(inds, id.vars = 'type', 
                measure.vars = paste0('cells.', 1:nrow(sample_points)), 
                variable.name = 'sample', value.name = 'count')
temp[,sample := gsub('cells.', replacement = '', x = sample, fixed = TRUE)]

cells_in_samples <- temp[, sum(count), by = .(type, sample)]
rm(temp)

with(cells_in_samples[sample == 1], plot(V1))
with(cells_in_samples[sample == 2], points(V1, col = 'red'))
abline(v = 1:50, col = grey(0.9))

cells_in_samples[, count := copies_per_cell(V1)]
template_abun <- cells_in_samples
template_abun[,V1:= NULL]

## TODO PICK UP HERE

N_samples <- 10
samples <- 1:N_samples

# the templates are present in a certain abundance in the extraction
template_abun <- matrix(NA, nrow = N_samples, ncol = N_seq_var)
for(i in 1:N_samples){
  template_abun[i,] <- rnbinom(n = N_seq_var, mu = 10, size = 0.6)
}


# each template has some primer efficiency / bias
primer_efficiencies <- rbeta(n = N_seq_var, 1, 1)

true_pcr_products <- expand.grid(seq_var, samples, NA)[,c(2,1,3)]
names(true_pcr_products) <- c('sample', 'seq.var', 'count')
for(i in 1:N_samples){
  current_max <- i * N_seq_var
  current_min <- current_max - N_seq_var + 1
  true_pcr_products[current_min:current_max,'count'] <- do_pcr(
    template_copies = template_abun[i,], template_effs = primer_efficiencies, 
    ncycles = 40, inflection = 10, slope = 0.2
  )
}

# when we pool our samples, the number of molecules will be influenced by:
# concentration measurement error:
DT <- data.table(true_pcr_products)
# samples contain different numbers of fragments:
molecules_per_sample_product <- DT[,sum(count), by = sample]$V1
# and our measurement of them is not accurate, both because of pipetting:
qubit_vol_exp <- rep(1, N_samples)
qubit_vol_act <- pipette(qubit_vol_exp)
pcr_vol   <- 50
lambderp <- molecules_per_sample_product * qubit_vol_act/pcr_vol
molecules_into_qubit <- rpois(N_samples, lambda = lambderp)
# and because of the measurement error of the qubit
qubit_mass <- rpois(N_samples, lambda = molecules_into_qubit)
min_vol_to_pool <- 10
vol_to_pool <- 10 * (min(qubit_mass)/qubit_mass)

# and again, when we pool them, there is error inherent to pipetting
mol_pipette_exp <- molecules_per_sample_product * (vol_to_pool/(pcr_vol))
mol_pipette_actual <- rpois(n = N_samples, lambda = mol_pipette_exp)

mols_to_seq <- expand.grid(seq_var, samples, NA)[,c(2,1,3)]
names(mols_to_seq) <- c('sample', 'seq.var', 'count')
for(i in 1:N_samples){
  current_max <- i * N_seq_var
  current_min <- current_max - N_seq_var + 1
  mols_to_seq[current_min:current_max,'count'] <- rmultinom(
    n = 1, size = mol_pipette_actual[i], 
    prob = true_pcr_products[true_pcr_products$sample == i, 'count'])
}

# The sequencer has some expected output ('depth')
seq_depth <- rpois(n = 1, lambda = 1e6)

# molecules (sample and variant) are sequenced according to their abundance
seq_prob <- mols_to_seq$count/sum(mols_to_seq$count)

# and the number of reads of each taxon from each sample might be:
seq_dat <- expand.grid(seq_var, samples, NA)[,c(2,1,3)]
names(seq_dat) <- c('sample', 'seq.var', 'count')
seq_dat[,'count'] <- rmultinom(n = 1, size = seq_depth, prob = seq_prob)