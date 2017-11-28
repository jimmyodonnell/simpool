library(data.table)

# Primer efficiency will be modified, though original sample will not
# so, add primer efficiency now rather than generate fixed efficiences upfront
template_dat[, eff := primer_eff(N = templates, mmv = 'low'), by = sample]

# running into a very weird error; indicating there might be a maximum to the 
# 'size' argument of rmultinom...
rmultinom(1, size = 2.147e9, prob = 1:3)
rmultinom(1, size = 2.148e9, prob = 1:3) # 'invalid second argument 'size'

template_dat[,
  do_pcr(template_copies = templates, template_effs = eff, ncycles = 30), by = sample
]
with(template_dat[sample == 1], do_pcr(templates, eff, ncycles = 30))
tmp <- template_dat[sample == 1, ]
do_pcr(tmp$templates, tmp$eff, ncycles = 30)
