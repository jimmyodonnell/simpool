################################################################################
# compare estimates from two samples

library(vegan)


div.in <- samples.env[ , .(simp = diversity(templates, index = "simpson")), by = sample]
div.in <- template_dat[ , .(simp = diversity(templates, index = "simpson")), by = sample]

################################################################################
div.out <- seq_dat[ , .(simp = diversity(count, index = "simpson")), 
  by = .(sample.env, pcr.id)]

div.in
div.out



boxplot(simp ~ sample.env, data = div.out, ylim = range(c(div.in[,simp], div.out[,simp])))
points(div.in[,simp], col = 2, lwd = 2)

# calculate differences between input and output
div.out[,simp]