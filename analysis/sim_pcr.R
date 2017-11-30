################################################################################
# simulate PCR
################################################################################

library(data.table)

# Primer efficiency will be modified, though original sample will not, 
# so add primer efficiency now rather than generate fixed efficiences upfront
template_dat[, eff := primer_eff(N = templates, mmv = 'low'), by = sample]

pcr_reps <- 10
temp <- list()
for(i in 1:pcr_reps){
  temp[[i]] <- template_dat[,list(
    rep.pcr = i, 
    species, 
    templates = do_pcr(template_copies = templates, 
                       template_effs = eff, 
                       ncycles = 30, inflection = 15, slope = 0.5)
    ), by = sample]
}
pcr_dat <- rbindlist(temp)
rm(temp)
