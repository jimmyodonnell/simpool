evenness.lev <- c('even', 'skew.lin', 'skew.low', 'skew.med', 'skew.hi')
even.col <- c(grey(0.7), 'chartreuse', 'gold', 'coral', 'mediumorchid')

p <- function(Y){ Y/sum(Y) }

# richness.lev <- c(10, 100, 1000)
N_sp <- 25

x <- 1:N_sp

par(mar = c(4,5,1,1))
plot(x, x, ylim = c(0, max(p(1/x))), 
     type = 'n', xlab = 'rank', ylab = 'p(S)', las = 1)

# even
y <- rep(1/N_sp, length(x))
points(x, p(y), type = 'b', col = even.col[1], lwd = 2)

# Here is a linearly decreasing option
B0 <- (N_sp+1)/N_sp
B1 <- -1/(2*N_sp)
y <- B0 + B1*x
points(x, p(y), type = 'b', col = even.col[2], lwd = 2)

# skew.low
Beta <- 1/(N_sp^0.5) # 0.8 = 0.125*10
y <- 1/((N_sp*Beta) + x)
points(x, p(y), type = 'b', col = even.col[3], lwd = 2)

# skew.med
Beta <- 1/(N_sp^1)# 0.01*length(x) # 0.1 = 0.01 * 10
y <- 1/((N_sp*Beta) + x)
points(x, p(y), type = 'b', col = even.col[4], lwd = 2)

# skew.hi
Beta <- 0 # same as: y <- 1/(x)
y <- 1/((N_sp*Beta) + x) # same as: y <- 1/(x)
points(x, p(y), type = 'b', col = even.col[5], lwd = 2)

legend('topright', 
       legend = evenness.lev, lwd = 2, col = even.col, pch = 1,
       bty = 'n')

