library(data.table)

# Primer efficiency will be modified, though original sample will not
# so, add primer efficiency now rather than generate fixed efficiences upfront
template_dat[, eff := primer_eff(N = templates, mmv = 'low'), by = sample]

template_dat[,
  do_pcr(template_copies = templates, template_effs = eff, ncycles = 40), by = sample
]
with(template_dat[sample == 1], do_pcr(templates, eff, ncycles = 30))
tmp <- template_dat[sample == 1, ]
do_pcr(tmp$templates, tmp$eff, ncycles = 30)
# do_pcr(template_copies = , template_effs = , ncycles = )
template_dat[ , rbeta(n = length(templates), shape1 = , shape2 = ), by = sample]


hist(rbeta(1000, 100, 1))

