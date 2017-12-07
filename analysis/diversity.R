################################################################################
# compare estimates from two samples

library(vegan)


# div.in <- samples.env[ , .(simp = diversity(templates, index = "simpson")), by = sample]
div.in <- template_dat[ , .(
  simp.in = diversity(templates, index = "simpson"), 
  shan.in = diversity(templates, index = "shannon")
  ), by = sample]

################################################################################
div.out <- pcr_dat[ , .(
  simp.out = diversity(amplicons, index = "simpson"), 
  shan.out = diversity(amplicons, index = "shannon")
  ), by = .(sample, rep.pcr)]

div.in
div.out

div.full <- merge(div.in, div.out, by = 'sample')

# calculate differences between input and output
div.full[, simp.diff := simp.in - simp.out]
div.full[, shan.diff := shan.in - shan.out]
div.full <- merge(div.full, unique(template_dat[,.(sample, even, rich)]), by = 'sample')

div.full[ , div.scen := paste(even, rich, sep = '.')]

divplot <- function(mmv, metric){
  if(metric == "simpson"){
  	divcol <- quote(simp.diff)
    XLAB <- "Delta Simpson"
  }
  if(metric == "shannon"){
  	divcol <- quote(shan.diff)
  	XLAB <- "Delta Shannon"
  }
  if(metric == "richness"){
  	divcol <- quote(rich.diff)
  	XLAB <- "Delta Richness"
  }
  maintext <- paste0('Mismatch Variation = ', mmv)
  par(mar = c(4, 7, 3, 1))
  boxplot(div.full[,eval(divcol)] ~ div.full[,div.scen], 
    # data = div.full, 
    # ylim = range(c(div.full[,simp.diff], div.full[,simp.diff])), 
    horizontal = TRUE, las = 1
  )
  title(main = maintext, xlab = XLAB)
  # points(div.full[, simp.in], col = 2, lwd = 2)
  
}
divplot(mmv = 'Low', metric = "simpson")

