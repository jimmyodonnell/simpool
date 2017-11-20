library(data.table)

functions.f <- list.files("functions", full.names = TRUE)
sapply(functions.f, source, echo = FALSE)

# 1. We have a pool of different types of organisms. We'll classify them into 
# groups based on the identity at a given region of the mitochondrial genome. 
# Thus, each 'variant' of that sequence corresponds to a type (aka species/OTU).
N_species <- 10 # was 50
species <- 1:N_species
dim_pool <- c(50, 25, 3)

# Each type of organism occurs in different abundances. 
# We could just simulate some random numbers of abundances, 
# or be more realistic and incorporate information about the mean body size of 
# each type (e.g. smaller organisms are more abundant). 
# To do that, we need to pick some mean body sizes for each of our 
# organism types. 
set.seed(1) # was 2017
# choose between "realistic", "uniform", or give a single number for all. 
BODYSIZE_DIST <- "realistic"
if(BODYSIZE_DIST == "realistic"){
  bodysize_mean <- rlnorm(N_species)
} else if(BODYSIZE_DIST == "uniform"){
  bodysize_mean <- runif(n = N_species, min = 0, max = 20)
} else if(is.numeric(BODYSIZE_DIST)){
  bodysize_mean <- rep(BODYSIZE_DIST, N_species)
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

# The smallest 'species' has a mean body mass of 0.43 g and around 2306 individuals.
# The largest has a mean body mass of 4.9 g and around 203 individuals.
# From this information we can generate some simulated but realistic count data. 
ind_per_species <- rpois(n = length(counts_mean), lambda = counts_mean)

# And to bring it back around, we can simulate a mass for each individual
set.seed(1)
mass_ind <- rgamma(
  n = sum(ind_per_species), 
  shape = rep(bodysize_mean, times = ind_per_species), 
  rate = 1
  )
inds <- data.table(
  species = rep(species, times = ind_per_species), 
  mass = mass_ind
)
boxplot(split(mass_ind, f = inds$species), 
  xlab = "species", ylab = "individual body mass (g)", 
  cex = 0.5, pch = '.', las = 1)
points(x = 1:N_species, y = bodysize_mean, col = 'red')

# I'm going to call place_inds() multiple times so each species has a different 
# distribution in space. 
inds[,c('x', 'y', 'z') := data.table(place_inds(.N, dim_pool)), by = species]

# plot distribution of a single species in the pool
sp.plot <- 4 # 4
par(mfrow = c(2,1))
ncols <- 3
mycols <- hsv(h = seq(from = 1/ncols, to = 1, length.out = ncols), s = 0.7, alpha = 0.7)
# mycols <- colorRampPalette(c("gold", "slateblue", 'green'))(ncols)
# mycols <- viridis(ncols, end = 1, option = 'magma')
with(inds[species == sp.plot, ], 
     plot(x, z, col = mycols[cut(z, breaks = ncols)], lwd = 2, 
          xlab = "length", ylab = "depth", las = 1))
legend('bottomright', legend = c('upper', 'middle', 'lower'), 
       pch = 1, pt.lwd = 2, col = rev(mycols), bty = 'n')

with(inds[species == sp.plot, ], 
     plot(x,y, col = mycols[cut(z, breaks = ncols)], lwd = 2, 
          xlab = "length", ylab = "width", las = 1))
legend('bottomright', legend = c('upper', 'middle', 'lower'), 
       pch = 1, pt.lwd = 2, col = rev(mycols), bty = 'n')

# Now we can sample from the environment.
# Let's start by deciding where to take our samples. 
corner_min <- c(0,0,0) # in the corner
pool_range <- cbind(corner_min, dim_pool)
sample_points <- as.matrix(expand.grid(pool_range[1,], pool_range[2,], pool_range[3,]))

samples <- nrow(sample_points)
N_samples <- seq_along(samples)

# In the following lines I'll be doing things in an inefficient way, but in the 
# hopes it will illustrate each step individually.

# First, calculate the distance that each sample is from each individual.
for(i in 1:nrow(sample_points)){
  # add columns for the distance between each individual and each sample location
  inds[ , paste0('dist.', i)] <- dist_from_samp(sample_points[i,], as.matrix(inds[,.(x,y,z)]))
}

# Then, sample cells from the environment based on density calculated based on 
# mass and distance. 
for(i in 1:nrow(sample_points)){
  distvar <- paste0('dist.', i)
  cellvar <- paste0('cells.', i)
  # for each individual, calculate how many cells end up in each sample
  inds[ , (cellvar)] <- cells_by_dist(inds[[distvar]], inds[,mass])
}

# temp data frame; very long: N cells in each sample from each *individual*
temp <- melt.data.table(inds, id.vars = 'species', 
                measure.vars = paste0('cells.', 1:nrow(sample_points)), 
                variable.name = 'sample', value.name = 'count')
temp[,sample := gsub('cells.', replacement = '', x = sample, fixed = TRUE)]

# now collapse that into N cells from each *type*
samples.env <- temp[, sum(count), by = .(sample, species)]
rm(temp)
colnames(samples.env)[3] <- "cells"

# plot cells per sample by species
par(mfrow = c(1,1))
pldat <- with(samples.env, split(cells, species))
boxplot(pldat, 
  xlab = 'number of cells in samples', ylab = 'species', 
  outpch = 21, outcol="slateblue", outbg="red", 
  horizontal = TRUE, las = 1)

# Compute number of templates per sample, 
# given number of cells and range of templates/cell
set.seed(1)
samples.env[, templates := copies_per_cell(cells)]
pldat <- with(samples.env, split(templates, species))
boxplot(pldat, 
  xlab = 'number of templates in samples', ylab = 'species', 
  outpch = 21, outcol="slateblue", outbg="red", 
  horizontal = TRUE, las = 1)

EXPORT <- FALSE
if(EXPORT){
  fwrite(samples.env, "../data/samples_env.csv")
}

# This is the end of environmental sampling
################################################################################



################################################################################
# This is where lab work begins.
# Note that as of now, there is no DNA extraction.

# pcr replication:
pcr_per_sample <- 3
N_pcr  <- N_samples * pcr_per_sample
pcr_id <- 1:N_pcr


# each template has some primer efficiency / bias
set.seed(3)
primer_efficiencies <- rbeta(n = N_species, 10, 1)
plot(primer_efficiencies, ylim = c(0,1), 
  xlab = "Species", ylab = "Primer Efficiency", las = 1)

# set up PCR output
temp <- expand.grid(species, 1:pcr_per_sample, samples)[,c(3,2,1)]
colnames(temp) <- c("sample.env", "pcr.rep", "species")
pcr_output <- data.table(pcr.id = rep(pcr_id, each = N_species), temp, count = NA_real_)
rm(temp)

# Run PCRs
set.seed(1)
for(i in 1:N_pcr){
  current_rows <- which(pcr_output[,pcr.id] == i)
  current_samp <- pcr_output[current_rows,unique(sample.env)]
  pcr_output[current_rows, count := do_pcr(
    template_copies = samples.env[sample == current_samp,templates], 
    template_effs = primer_efficiencies, 
    ncycles = 40, inflection = 10, slope = 0.2)
  ]
}

# when we pool our samples, the number of molecules will be influenced by:
# concentration measurement error:
DT <- data.table(true_pcr_products)
# samples contain different numbers of fragments:
molecules_per_pcr_product <- pcr_output[,sum(count), by = pcr.id]$V1
# and our measurement of them is not accurate, both because of pipetting:
qubit_vol_exp <- rep(1, N_pcr)
qubit_vol_act <- pipette(qubit_vol_exp)
pcr_vol   <- 50
lambderp <- molecules_per_pcr_product * qubit_vol_act/pcr_vol
molecules_into_qubit <- rpois(N_pcr, lambda = lambderp)
# and because of the measurement error of the qubit
qubit_mass <- rpois(N_pcr, lambda = molecules_into_qubit)
min_vol_to_pool <- 10
vol_to_pool <- 10 * (min(qubit_mass)/qubit_mass)

# and again, when we pool them, there is error inherent to pipetting
mol_pipette_exp <- molecules_per_pcr_product * (vol_to_pool/(pcr_vol))
mol_pipette_actual <- rpois(n = N_pcr, lambda = mol_pipette_exp)

# The molecules go onto the sequencer
mols_to_seq <- copy(pcr_output)
mols_to_seq[,count:= NA_real_]
for(i in 1:N_pcr){
  current_rows <- which(mols_to_seq[,pcr.id] == i)
  mols_to_seq[current_rows, count := 
    as.numeric(rmultinom(n = 1, size = mol_pipette_actual[i], 
      prob = pcr_output[current_rows, count]))
  ][]
}

# The sequencer has some expected output ('depth')
seq_depth <- rpois(n = 1, lambda = 1e6)

# molecules (sample and variant) are sequenced according to their abundance
seq_prob <- mols_to_seq[,count/sum(count)]

# and the number of reads of each taxon from each sample might be:
seq_dat <- copy(mols_to_seq)[,count:= NA][]
set.seed(1)
seq_dat[,count := rmultinom(n = 1, size = seq_depth, prob = seq_prob)[,1]]

EXPORT <- FALSE
if(EXPORT){
  fwrite(seq_dat, file = "../data/seq_dat.csv")
}
