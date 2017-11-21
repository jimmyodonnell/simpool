################################################################################
# generate a bunch of 'purified dna' samples

# columns: sample species cells templates

N_species <- 10

N_samples <- 8

# total templates per sample
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
