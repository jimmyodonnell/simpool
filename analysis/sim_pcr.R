################################################################################
# simulate PCR
################################################################################

library(data.table)

# Primer efficiency will be modified, though original sample will not, 
# so add primer efficiency now rather than generate fixed efficiences upfront
template_dat[, eff := primer_eff(N = templates, mmv = 'low'), by = sample]

pcr_dat <- template_dat[,list(
  species, 
  templates = do_pcr(template_copies = templates, template_effs = eff, ncycles = 30)
  ), by = sample]
