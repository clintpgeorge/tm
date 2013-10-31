# ##############################################################################################
# Generates Q-Q plots and histograms for the p-values over all the words in all
# the documents, for four hyperparameter settings. The plots compare the empirical 
# quantiles of the p-values with the quantiles of the uniform distribution on (0, 1).
# 
# Dependency: 
#     It uses the outputs of the R script lda_fgs_acgs_analyze_z 
# ##############################################################################################

plot.p.values <- function(rdata.file, 
                          pvalues.hist.file, 
                          pvalues.hist.title, 
                          pvalues.qqplot.file,
                          pvalues.qqplot.title){
  
  load(rdata.file)
  
  ## Generates the histograms of p-values 
  
  postscript(file=pvalues.hist.file, title=pvalues.hist.title, horiz=F, height=5, width=7.5) 
  par(mar=c(5-1, 4+1, 4, 2) +.1) # c(bottom,left,top,right)
  hist(p.values, cex = 2, cex.lab = 1.8, cex.axis = 1.4, main="", xlab="p-values", breaks=40)
  dev.off()
  
  ## Generating Q-Q plots  
  
  q <- ((1:total.num.words) - 0.5) / total.num.words; ## the theoretical quantile values 
  postscript(file=pvalues.qqplot.file, title=pvalues.qqplot.title, horiz=F, height=5.5, width=5.5) 
  par(mar=c(5, 4, 4, 2) +.1) # c(bottom,left,top,right)
  qqplot(q, p.values, cex = 0.7, cex.lab = 1.3, cex.axis = 1.2, main="", 
         ylab = "Empirical quantiles of the p-values", 
         xlab = "Quantiles of the uniform distribution"); 
  abline(0,1,lwd=2); ## a 45-degree reference line is plotted
  dev.off()
  
}



rdata.file <- "fg_cg_ea108_pvalues.RData"
pvalues.hist.file <- "fg_cg_ea108_pvalues.eps"
pvalues.hist.title <- "h=(10,8): p-values"
pvalues.qqplot.file <- "fg_cg_ea108_pvalues_qqplot.eps"
pvalues.qqplot.title <- "h=(10,8): QQ plot of the p-values"
plot.p.values(rdata.file, 
              pvalues.hist.file, 
              pvalues.hist.title, 
              pvalues.qqplot.file,
              pvalues.qqplot.title);

rdata.file <- "fg_cg_ea33_pvalues.RData"
pvalues.hist.file <- "fg_cg_ea33_pvalues.eps"
pvalues.hist.title <- "h=(3,3): p-values"
pvalues.qqplot.file <- "fg_cg_ea33_pvalues_qqplot.eps"
pvalues.qqplot.title <- "h=(3,3): QQ plot of the p-values"
plot.p.values(rdata.file, 
              pvalues.hist.file, 
              pvalues.hist.title, 
              pvalues.qqplot.file,
              pvalues.qqplot.title);

rdata.file <- "fg_cg_ea37_pvalues.RData"
pvalues.hist.file <- "fg_cg_ea37_pvalues.eps"
pvalues.hist.title <- "h=(3,7): p-values"
pvalues.qqplot.file <- "fg_cg_ea37_pvalues_qqplot.eps"
pvalues.qqplot.title <- "h=(3,7): QQ plot of the p-values"
plot.p.values(rdata.file, 
              pvalues.hist.file, 
              pvalues.hist.title, 
              pvalues.qqplot.file,
              pvalues.qqplot.title);

rdata.file <- "fg_cg_ea77_pvalues.RData"
pvalues.hist.file <- "fg_cg_ea77_pvalues.eps"
pvalues.hist.title <- "h=(7,7): p-values"
pvalues.qqplot.file <- "fg_cg_ea77_pvalues_qqplot.eps"
pvalues.qqplot.title <- "h=(7,7): QQ plot of the p-values"
plot.p.values(rdata.file, 
              pvalues.hist.file, 
              pvalues.hist.title, 
              pvalues.qqplot.file,
              pvalues.qqplot.title);


rdata.file <- "fg_cg_ea73_pvalues.RData"
pvalues.hist.file <- "fg_cg_ea73_pvalues.eps"
pvalues.hist.title <- "h=(7,3): p-values"
pvalues.qqplot.file <- "fg_cg_ea73_pvalues_qqplot.eps"
pvalues.qqplot.title <- "h=(7,3): QQ plot of the p-values"
plot.p.values(rdata.file, 
              pvalues.hist.file, 
              pvalues.hist.title, 
              pvalues.qqplot.file,
              pvalues.qqplot.title);


# ##############################################################################################
