plotSmoothersIk <- function(m, nPoints = 100, ...){
  data <- m$model
  y <- data$y
  #construct time variable based on cell assignments.
  nCurves <- length(m$smooth)
  timeAll <- c()
  col <- c()
  for (jj in seq_len(nCurves)) {
    for (ii in 1:nrow(data)) {
      if (data[ii, paste0("l", jj)] == 1) {
        timeAll[ii] <- data[ii, paste0("t", jj)]
        col[ii] <- jj
      } else {
        next
      }
    }
  }

  # plot raw data
    #cols <- c("#E7298A", "#FF7F00", "#1F78B4")
    cols <- c("#FF7F00", "#1F78B4", "#E7298A")
  plot(x = timeAll, y = log(y + 1), col = alpha(cols[col],2/3), pch = 16, cex = 2 / 3,
       ylab = "log(count + 1)", xlab = "Pseudotime", ...)

  #predict and plot smoothers across the range
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(m, jj, nPoints = nPoints)
  yhat <- predict(m, newdata = df, type = "response")
  lines(x = df[, paste0("t", jj)], y = log(yhat + 1), col = cols[jj], lwd = 2)
  }
  # knots
  #  abline(v=gamList[[1]]$smooth[[1]]$xp[2], lty=2, col="black", lwd=1.5)
  #  abline(v=gamList[[1]]$smooth[[1]]$xp[4], lty=2, col="black", lwd=1.5)
  legend("bottomleft", c("Neuronal", "Microvillous", "Sustentacular"),col = cols,
         lty = 1, lwd = 2, bty = "n", cex = 4 / 5)
}
