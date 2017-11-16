counts_in <- c(0, 1, 2, 3, 4, 5, 10, 20, 50, 100)
pcr_effs  <- c(1, 1, 1, 1, 1, 0.5, 0.01, 0.01, 0, 0.01)
CYCLES <- 30
pcr_out <- do_pcr(
  template_copies = counts_in, template_effs = pcr_effs, 
  ncycles = CYCLES, inflection = CYCLES/3, slope = 0.2
)

pcr_1 <- do_pcr(
  template_copies = 2, template_effs = 0.5, 
  ncycles = 40, inflection = 1, slope = 0.9)

PLOT <- FALSE
if(PLOT){
  par(mar = c(4,5,1,1))
  plot_dat <- pcr_out
  plot(x = 0:CYCLES, y = seq(1e-2, max(plot_dat), length.out = CYCLES+1), 
       xlab = 'Cycle', ylab = '', las = 1, type = 'n')
  title(ylab = 'Copies', line = 4)
  for(i in 1:ncol(plot_dat)){
    hue_inc <- 1/ncol(plot_dat)
    hue <- i * hue_inc
    points(plot_dat[,i], type = 'l', lwd = 3, col = hsv(hue, 0.5, 1))
    text(x = CYCLES+2, y = plot_dat[CYCLES,i], 
         labels = LETTERS[i], col = hsv(hue, 0.8, 1), font = 2, cex = 1.2, xpd = TRUE)
  }
}
