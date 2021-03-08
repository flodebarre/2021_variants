rm(region_index)
figures <- c(
             "~/ownCloud/coronavirus/variant_France/inferred_transmission_advantage",
             "~/ownCloud/coronavirus/variant_France/cases_per_region_France",
             "~/ownCloud/coronavirus/variant_France/all_figures_France",
             "~/ownCloud/coronavirus/variant_France/cases_per_region_France_extrapolated",
             "~/ownCloud/coronavirus/variant_France/frequency_per_region_France",
             "~/ownCloud/coronavirus/variant_France/sigmoid_per_region_France"
)
filedates <- "~/ownCloud/coronavirus/variant_France/date_50pc"

figures <- paste0(figures, "_", Sys.Date(), "_", upper_limit_rho, "_", meangt, ".pdf")
filedates <- paste0(filedates, "_", Sys.Date(), "_", upper_limit_rho, "_", meangt, ".csv")

########################################################################### PLOT RESULTS ###########################################################################


# check the mcmc chain:
par(mfrow = c(1,1))
plot(mcmc$all.lik, type = "l")
for(i in 1:nregions) plot(mcmc$chain[, 5*(i-1) + 2], type = "l", ylab = "initial R wild type", main = allregions[i])
for(i in 1:nregions) plot(mcmc$chain[, 5*(i-1) + 3], type = "l", ylab = "final R wild type", main = allregions[i])
plot(mcmc$chain[, 5*nregions+1], type = "l", ylab = "transmissibility advantage")

pdf(figures[1], width = 4, height = 3)
par(mfrow = c(1,1), mar = c(4, 4, 1, 1))
hist(mcmc$chain[(burnin_fraction*niter):niter, 5*nregions+1], freq = F, main = "inferred transmission advantage", xlim = c(0, 1), las = 1, ylab = "density", xlab = "transmission advantage") # now prior is unifomr
#points(x = seq(0, 2, 0.01), y = dnorm(x = seq(0, 2, 0.01), mean = 0.5, sd = 0.2), col = "black", type = "l", lty = 1, lwd = 1)
#legend("topright", col = "black", lty = 1, legend = "prior", bty = "n")
dev.off()


datelabels <- c("04/01", "11/01", "18/01", "25/01", "01/02", "08/02", "15/02", "22/02")
datenumbers <- seq(370, 419, 7)

# get predictions across mcmc samples
pdf(figures[2], width = 3*4*0.7, height = 4*3*0.7)
par(mfrow = c(4,3), mar = c(4,4,1,1))
all_date_50 <- c() # date to 50%
for(rr in allregions){
  
  #region_index <- which(allregions==rr)
  mymain <- rr
  #if(rr=="PROVENCE ALPES COTE D AZUR") mymain <- "PACA" else mymain <- rr
  plot(a2$dateday[a2$reg2==rr & a2$dateday>=minday], (a2$Psmoothed[a2$reg2==rr & a2$dateday>=minday]), type = "o", pch = 20, xaxs = "i", yaxs = "i",
       main = mymain, cex = 1, xlab = "date", ylab = "daily number of cases", ylim = c(-100, 4500), axes = F, xlim = c(370, maxday+7))
  axis(side = 1, at = datenumbers, labels = datelabels)
  axis(side = 2, at = seq(0, 4500, 1000), las = 1)
  
  myattach(create_enveloppe(rr = rr))
  
  env_TOT <- apply(env_TOT, 2, quantile, c(0.025, 0.5, 0.975))
  env_VOC <- apply(env_VOC, 2, quantile, c(0.025, 0.5, 0.975))
  all_date_50 <- rbind(all_date_50, quantile(date_50, c(0.025, 0.5, 0.975), na.rm = T))
  polygon(x = c(truc$dateday, rev(truc$dateday)), y = c(env_TOT[1,], rev(env_TOT[3,])), col = "gray", border = NA)
  points(truc$dateday, y = env_TOT[2,], col = "black", type = "l")
  points(a2$dateday[a2$reg2==rr & a2$dateday>=minday], (a2$Psmoothed[a2$reg2==rr & a2$dateday>=minday]), type = "o", pch = 20)
  cols <- RColorBrewer::brewer.pal(12, "Paired")[5:6]
  polygon(x = c(truc$dateday, rev(truc$dateday)), y = c(env_VOC[1,], rev(env_VOC[3,])), col = cols[1], border = NA)
  points(truc$dateday, y = env_VOC[2,], col = cols[2], type = "l")
  
}
dev.off()

# figure all of France
pdf(figures[3], width = 4*1., height = 3*3*1.)
par(mfrow = c(3,1), mar = c(4,5,1,1))

mymain <- ""
plot(epitable0$dateday, (epitable0$Psmoothed), type = "o", pch = 20, xaxs = "i", yaxs = "i",
     main = mymain, cex = 1, xlab = "date", ylab = "", ylim = c(-100, 25000), axes = F, xlim = c(370, maxday+7))
title(ylab="daily number of cases", mgp=c(3.5,1,0))
axis(side = 1, at = datenumbers, labels = datelabels)
axis(side = 2, at = seq(0, 25000, 5000), las = 1)
samp <- sample((burnin_fraction*niter):niter, 1000)

# add enveloppes:
polygon(x = c(truc$dateday, rev(truc$dateday)), y = c(env_TOT0[1,], rev(env_TOT0[3,])), col = "gray", border = NA)
points(truc$dateday, y = env_TOT0[2,], col = "black", type = "l")
points(epitable0$dateday, epitable0$Psmoothed, type = "o", pch = 20)
cols <- RColorBrewer::brewer.pal(12, "Paired")[5:6]
polygon(x = c(truc$dateday, rev(truc$dateday)), y = c(env_VOC0[1,], rev(env_VOC0[3,])), col = cols[1], border = NA)
points(truc$dateday, y = env_VOC0[2,], col = cols[2], type = "l")

# frequency:
plot(NULL, type = "o", pch = 20, xaxs = "i", yaxs = "i",
     main = mymain, xlab = "date", ylab = "frequency of the variant", ylim = c(0, 1), axes = F, xlim = c(370, maxday+7))
axis(side = 1, at = datenumbers, labels = datelabels)
axis(side = 2, at = seq(0, 1, 0.2), las = 1)
polygon(x = c(truc$dateday, rev(truc$dateday)), y = c(env_VOC_f_0[1,], rev(env_VOC_f_0[3,])), col = cols[1], border = NA)
points(truc$dateday, y = env_VOC_f_0[2,], col = cols[2], pch = 20, type = "l", lwd = 3)
freqdatatable0$daymid <- 0.5 * (as.numeric(freqdatatable0$daymin) + as.numeric(freqdatatable0$daymax))
freqdatatable0$fVOC <- freqdatatable0$NVOC/freqdatatable0$N
freqdatatable0$lower <- freqdatatable0$fVOC - 1.96 * sqrt(freqdatatable0$fVOC * (1-freqdatatable0$fVOC) / freqdatatable0$N)
freqdatatable0$upper <- freqdatatable0$fVOC + 1.96 * sqrt(freqdatatable0$fVOC * (1-freqdatatable0$fVOC) / freqdatatable0$N)
for(k in 1:3){
  points(freqdatatable0[k, "daymid"], freqdatatable0[k, "fVOC"], pch = 0, col = cols[2], cex = 2)
  segments(x0 = freqdatatable0[k, "daymid"], x1 = freqdatatable0[k, "daymid"], y0 = freqdatatable0[k, "lower"], y1 = freqdatatable0[k, "upper"], col = cols[2], lwd = 3)
}
segments(x0 = date_50_0["50%"], x1 = date_50_0["50%"], y0 = 0, y1 = 0.5, col = "red", lwd =2, lty = 2)
segments(x0 = minday, x1 = date_50_0["50%"], y0 = 0.5, y1 = 0.5, col = "red", lwd =2, lty = 2)

# sigmoid:
plot(NULL, type = "o", pch = 20, xaxs = "i", yaxs = "i",
     main = "", xlab = "date", ylab = "effective reproduction number", ylim = c(0.6, 2), axes = F, xlim = c(370, maxday+7))
axis(side = 1, at = datenumbers, labels = datelabels)
axis(side = 2, at = seq(0, 2, 0.2), las = 1)

cols <- RColorBrewer::brewer.pal(12, "Paired")[1:2]
polygon(x = c(minday:(maxday+90-1), (maxday+90-1):minday), y = c(env_RWT0[1,], rev(env_RWT0[3,])), col = cols[1], border = NA)
points(minday:(maxday+90-1), y = env_RWT0[2,], col = cols[2], pch = 20, type = "l", lwd = 3)

cols <- RColorBrewer::brewer.pal(12, "Paired")[5:6]
polygon(x = c(minday:(maxday+90-1), (maxday+90-1):minday), y = c(env_RVOC0[1,], rev(env_RVOC0[3,])), col = cols[1], border = NA)
points(minday:(maxday+90-1), y = env_RVOC0[2,], col = cols[2], pch = 20, type = "l", lwd = 3)

abline(h = 1, lty = 3, lwd = 1)
abline(v = 366+16, lty = 3, lwd = 1)
dev.off()

datelabels2 <- c("04/01", "18/01", "01/02", "15/02", "01/03", "15/03", "29/03")
datenumbers2 <- seq(370, 454, 14)

pdf(figures[4], width = 3*4*0.7, height = 4*3*0.7)
par(mfrow = c(4,3), mar = c(4,4,1,1))
all_date_50 <- c() # date to 50%
for(rr in allregions){
  #region_index <- which(allregions==rr)
  mymain <- rr
  #if(rr=="PROVENCE ALPES COTE D AZUR") mymain <- "PACA" else mymain <- rr
  plot(a2$dateday[a2$reg2==rr & a2$dateday>=minday], (a2$Psmoothed[a2$reg2==rr & a2$dateday>=minday]), type = "o", pch = 20, xaxs = "i", yaxs = "i",
       main = mymain, cex = 1, xlab = "date", ylab = "daily number of cases", ylim = c(-100, 10000), axes = F, xlim = c(370, maxday+56))
  axis(side = 1, at = datenumbers2, labels = datelabels2)
  axis(side = 2, at = seq(0, 10000, 2000), las = 1)
  
  myattach(create_enveloppe(rr))
  env_TOT <- apply(env_TOT, 2, quantile, c(0.025, 0.5, 0.975))
  env_VOC <- apply(env_VOC, 2, quantile, c(0.025, 0.5, 0.975))
  env_WT <- apply(env_WT, 2, quantile,  c(0.025, 0.5, 0.975))
  all_date_50 <- rbind(all_date_50, quantile(date_50, c(0.025, 0.5, 0.975), na.rm = T))
  
  polygon(x = c(truc$dateday, rev(truc$dateday)), y = c(env_TOT[1,], rev(env_TOT[3,])), col = "gray", border = NA)
  points(truc$dateday, y = env_TOT[2,], col = "black", type = "l")
  cols <- RColorBrewer::brewer.pal(12, "Paired")[1:2]
  polygon(x = c(truc$dateday, rev(truc$dateday)), y = c(env_WT[1,], rev(env_WT[3,])), col = cols[1], border = NA)
  points(truc$dateday, y = env_WT[2,], col = cols[2], type = "l")
  points(a2$dateday[a2$reg2==rr & a2$dateday>=minday], (a2$Psmoothed[a2$reg2==rr & a2$dateday>=minday]), type = "o", pch = 20)
  cols <- RColorBrewer::brewer.pal(12, "Paired")[5:6]
  polygon(x = c(truc$dateday, rev(truc$dateday)), y = c(env_VOC[1,], rev(env_VOC[3,])), col = cols[1], border = NA)
  points(truc$dateday, y = env_VOC[2,], col = cols[2], type = "l")
  
}
dev.off()

# get dates:
all_date_50 <- data.frame(all_date_50)
all_date_50$region <- allregions
names(all_date_50)[1:3] <- c("date_lower", "date_median", "date_upper")
for(mycol in names(all_date_50)[1:3]){
  all_date_50[, mycol] <- as.Date(floor(all_date_50[, mycol]), origin = "2019-12-31")
}
all_date_50 <- all_date_50[, c("region", "date_median", "date_lower", "date_upper")]
write.csv(x = all_date_50, file = filedates, row.names = F)

pdf(figures[5], width = 3*4*0.7, height = 4*3*0.7)
par(mfrow = c(4,3), mar = c(4,4,1,1))

for(rr in allregions){
  region_index <- which(allregions==rr)
  plot(NULL, type = "o", pch = 20, xaxs = "i", yaxs = "i",
       main = rr, xlab = "date", ylab = "frequency of the variant", ylim = c(0, 0.5), axes = F, xlim = c(370, maxday+7))
  axis(side = 1, at = datenumbers, labels = datelabels)
  axis(side = 2, at = seq(0, 1, 0.2), las = 1)
  
  myattach(create_enveloppe(rr = rr))
  
  # use the enveloppe of frequency:
  env_VOC_f <- apply(env_VOC_f, 2, quantile, c(0.025, 0.5, 0.975))
  polygon(x = c(truc$dateday, rev(truc$dateday)), y = c(env_VOC_f[1,], rev(env_VOC_f[3,])), col = cols[1], border = NA)
  points(truc$dateday, y = env_VOC_f[2,], col = cols[2], pch = 20, type = "l", lwd = 3)
  
  sub <- which(freqs4$region==rr)
  if(length(sub) > 0){
    for(k in sub){
      points(freqs4[k, "daymid"], freqs4[k, "fVOC"], pch = 0, col = cols[2], cex = 2)
      segments(x0 = freqs4[k, "daymid"], x1 = freqs4[k, "daymid"], y0 = freqs4[k, "lower"], y1 = freqs4[k, "upper"], col = cols[2], lwd = 3)
    }
  }
  
}
dev.off()

#####   PLOT THE SIGMOID PER REGION  #####

pdf(figures[6], width = 3*4*0.7, height = 4*3*0.7)
par(mfrow = c(4,3), mar = c(4,4,1,1))

for(rr in allregions){
  #region_index <- which(allregions==rr)
  
  plot(NULL, type = "o", pch = 20, xaxs = "i", yaxs = "i",
       main = rr, xlab = "date", ylab = "effective reproduction number", ylim = c(0.6, 2), axes = F, xlim = c(370, maxday+7))
  axis(side = 1, at = datenumbers, labels = datelabels)
  axis(side = 2, at = seq(0, 2, 0.2), las = 1)
  
  myattach(create_enveloppe(rr = rr))
  # create enveloppe for R
  env_RWT <- apply(env_RWT, 2, quantile, c(0.025, 0.5, 0.975))
  env_RVOC <- apply(env_RVOC, 2, quantile, c(0.025, 0.5, 0.975))
  
  cols <- RColorBrewer::brewer.pal(12, "Paired")[1:2]
  polygon(x = c(minday:(maxday-1), (maxday-1):minday), y = c(env_RWT[1,], rev(env_RWT[3,])), col = cols[1], border = NA)
  points(minday:(maxday-1), y = env_RWT[2,], col = cols[2], pch = 20, type = "l", lwd = 3)
  
  cols <- RColorBrewer::brewer.pal(12, "Paired")[5:6]
  polygon(x = c(minday:(maxday-1), (maxday-1):minday), y = c(env_RVOC[1,], rev(env_RVOC[3,])), col = cols[1], border = NA)
  points(minday:(maxday-1), y = env_RVOC[2,], col = cols[2], pch = 20, type = "l", lwd = 3)
  
  abline(h = 1, lty = 2, lwd = 1)
  abline(v = 366+16, lty = 3, lwd = 1)
  
}
dev.off()







