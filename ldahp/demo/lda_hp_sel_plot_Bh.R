################################################################################# 
## This script is to create EPS plots of the likelihood ratios' (B(h)) surface 
## from loaded execution data    
## 
## Required: lattice package 
################################################################################# 


eps.bh.surface <- function(ratios, start, end, interval, xlabel, ylabel, zlabel, plot.title, plot.file){
  
  library(lattice)
  
  x <- seq(start, end, by=interval); # e.g. (0.4, 12, 0.4)
  y <- x;
  L <- length(x);
  
  z <- array(0, dim=c(L, L));
  count <- 1;
  for (i in 1:L){
    for (j in 1:L){
      z[i,j] <- ratios[count];
      count <- count + 1;
    }
  }
  
  op <- par(bg = "white")
  
  trellis.device(postscript, file=plot.file,
                 height=5.5, width=6.5, horiz=F,
                 title=plot.title, onefile=T)
  par(mar=c(4-3, 4, 4-3, 2-2) + .1) # c(bottom, left, top, right)
  persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
        ltheta = 120, shade = 0.3, ticktype = "detailed",
        xlab = xlabel, ylab = ylabel, zlab = zlabel, 
        cex = 1.5, cex.lab = 1.35) -> res
  
  par(op)
  dev.off()
  
  round(res, 3);
  
}

## Plots surface for different experiments 


xlabel     <- "alpha"
ylabel     <- "eta"
zlabel     <- "Estimate of B(h)"


# 1 
rdata.file     <- "/home/clintpg/results/fg_ae33.RData"
plot.file      <- "/home/clintpg/results/fg_ae33.eps"
plot.title     <- "true h = (3, 3)"
load(rdata.file)
eps.bh.surface(ratios, start, end, interval, xlabel, ylabel, zlabel, plot.title, plot.file)

# 2 
rdata.file     <- "/home/clintpg/results/fg_ae37.RData"
plot.file      <- "/home/clintpg/results/fg_ae37.eps"
plot.title     <- "true h = (3, 7)"
load(rdata.file)
eps.bh.surface(ratios, start, end, interval, xlabel, ylabel, zlabel, plot.title, plot.file)

# 3
rdata.file     <- "/home/clintpg/results/fg_ae77.RData"
plot.file      <- "/home/clintpg/results/fg_ae77.eps"
plot.title     <- "true h = (7, 7)"
load(rdata.file)
eps.bh.surface(ratios, start, end, interval, xlabel, ylabel, zlabel, plot.title, plot.file)

# 4 
rdata.file     <- "/home/clintpg/results/fg_ae108.RData"
plot.file      <- "/home/clintpg/results/fg_ae108.eps"
plot.title     <- "true h = (10, 8)"
load(rdata.file)
eps.bh.surface(ratios, start, end, interval, xlabel, ylabel, zlabel, plot.title, plot.file)







# rdata.file     <- "~/workspace/tm/data/fg_synth_cfg11.RData"
# plot.file      <- "~/workspace/tm/data/fg_synth_cfg11.eps"
# plot.title     <- "true h = (3, 3)"
# load(rdata.file)
# eps.bh.surface(ratios, start, end, interval, xlabel, ylabel, zlabel, plot.title, plot.file)


# postscript(file=plot.file, title=plot.title)  
# persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
#       ltheta = 120, shade = 0.75, ticktype = "detailed",
#       xlab = xlabel, ylab = ylabel, zlab = zlabel, 
#       cex = 1.5, cex.lab = 1.5) -> res
