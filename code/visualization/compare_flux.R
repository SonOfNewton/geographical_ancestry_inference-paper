##### parameters #####
WORLD = my_world
map = c("friction", "naive")

for (MAP in map){
  x = read.csv(sprintf("data/flux/flux_strait_%s_%s_%s_%s.csv", WORLD, source_pop, end_gen, MAP))
  n = nrow(x)
  
  # true flux
  true_sinai = numeric(n)
  true_mandeb = numeric(n)
  true_gibraltar = numeric(n)
  
  for (i in 1:n) {
    rep_id = x$rep[i]
    f_data = read.csv(sprintf("data/flux/flux_%s_%s_%s_%s_%s.csv", WORLD, source_pop, end_gen, MAP, rep_id), row.names=1)
    true_mat = matrix(f_data$true_flux, nrow=177, ncol=177)
    true_sinai[i] = sum(true_mat[38, 69:70])
    true_mandeb[i] = true_mat[67, 68]
    true_gibraltar[i] = true_mat[155, 7]
  }
  
  # sort
  flux_features <- data.frame(sinai = true_sinai, mandeb = true_mandeb, gibraltar = true_gibraltar)
  dist_matrix <- dist(flux_features)
  mds_1d <- cmdscale(dist_matrix, k = 1)
  smooth_order <- order(mds_1d[, 1])
  
  true_sinai <- true_sinai[smooth_order]
  true_mandeb <- true_mandeb[smooth_order]
  true_gibraltar <- true_gibraltar[smooth_order]
  x <- x[smooth_order, ]
  
  # plot
  output_file = sprintf("output/figures/compare_flux_%s_%s_%s_%s.png", WORLD, source_pop, end_gen, MAP)
  png(file=output_file, width=10, height=5, units="in", res=300, bg = "transparent")
  
  par(mfrow=c(1, 1), mar=c(5, 4, 4, 2) + 0.1)
  plot(0, 0, xlim=c(1, n), ylim=c(0, 1), type='n', las=1, xaxt='n', xlab="", ylab="flux")
  mtext("replications", side=1, line=3.5)
  axis(1, at=1:n, labels=1:n, cex.axis=0.8, las=2)
  
  lines(1:n, x$flux_sinai, col="black", lty=1, lwd=1.2)
  points(1:n, x$flux_sinai, pch=19, col="black", cex=0.9)
  
  lines(1:n, x$flux_mandeb, col="brown", lty=1, lwd=1.2)
  points(1:n, x$flux_mandeb, pch=19, col="brown", cex=0.9)
  
  lines(1:n, x$flux_gibraltar, col="orange", lty=1, lwd=1.2)
  points(1:n, x$flux_gibraltar, pch=19, col="orange", cex=0.9)
  
  lines(1:n, true_sinai, col="black", lty=2)
  points(1:n, true_sinai, pch=1, col="black", cex=0.9)
  
  lines(1:n, true_mandeb, col="brown", lty=2)
  points(1:n, true_mandeb, pch=1, col="brown", cex=0.9)
  
  lines(1:n, true_gibraltar, col="orange", lty=2)
  points(1:n, true_gibraltar, pch=1, col="orange", cex=0.9)
  
  legend("topright", 
         col = c("black", "black", "brown", "brown", "orange", "orange"), 
         legend = c("Estimated Sinai", "True Sinai", "Estimated Mandeb", "True Mandeb", "Estimated Gibraltar", "True Gibraltar"), 
         lty = c(1, 2, 1, 2, 1, 2), 
         pch = c(19, 1, 19, 1, 19, 1), 
         bty = 'n', cex = 0.8)
  
  dev.off()
  message(sprintf("saved figure to: %s", output_file))


  ##### plot mean and std of estimation of both maps in one figure #####
  x_friction <- read.csv(sprintf("data/flux/flux_strait_%s_%s_%s_friction.csv", WORLD, source_pop, end_gen))
  x_naive <- read.csv(sprintf("data/flux/flux_strait_%s_%s_%s_naive.csv", WORLD, source_pop, end_gen))
  
  n <- nrow(x_friction)
  true_sinai <- numeric(n)
  true_mandeb <- numeric(n)
  true_gibraltar <- numeric(n)
  
  for (i in 1:n) {
    rep_id <- x_friction$rep[i]
    f_data <- read.csv(sprintf("data/flux/flux_%s_%s_%s_%s_%s.csv",WORLD, source_pop, end_gen, MAP, rep_id), row.names=1)   
    # Note: theoretically, true_flux in flux_WORLD_friction_1.csv should be identical to flux_WORLD_naive_1.
    # However, due to randomness in selection of multiple optimal routes, they vary slightly from each other.
    # Thus we perform both comparisons, using one as result and the other as robustness test.
    true_mat <- matrix(f_data$true_flux, nrow=177, ncol=177)
    true_sinai[i]     <- sum(true_mat[38, 69:70])
    true_mandeb[i]    <- true_mat[67, 68]
    true_gibraltar[i] <- true_mat[155, 7]
  }
  
  # compute estimation error
  summary_stats <- data.frame(
    Route = rep(c("Sinai", "Mandeb", "Gibraltar"), each = 2),
    Method = rep(c("Naive", "Friction"), 3),
    
    Mean_Dev = c(mean(x_naive$flux_sinai - true_sinai),
                 mean(x_friction$flux_sinai - true_sinai),
                 mean(x_naive$flux_mandeb - true_mandeb),
                 mean(x_friction$flux_mandeb - true_mandeb),
                 mean(x_naive$flux_gibraltar - true_gibraltar),
                 mean(x_friction$flux_gibraltar - true_gibraltar)),
    
    SD_Dev = c(sd(x_naive$flux_sinai - true_sinai),
               sd(x_friction$flux_sinai - true_sinai),
               sd(x_naive$flux_mandeb - true_mandeb),
               sd(x_friction$flux_mandeb - true_mandeb),
               sd(x_naive$flux_gibraltar - true_gibraltar),
               sd(x_friction$flux_gibraltar - true_gibraltar))
  )
  
  summary_stats$y <- rep(3:1, each = 2) + rep(c(0.12, -0.12), 3)
  
  # plot
  output_file <- sprintf("output/figures/compare_flux_%s_%s_%s_summary_%s_base.png", WORLD, source_pop, end_gen, MAP)
  png(file = output_file, width = 9.5, height = 9, units = "in", res = 300, bg = "transparent")
  
  par(mar = c(5, 8, 3, 3), las = 1, cex.axis = 1.35, cex.lab = 1.45)
  x_range <- range(c(summary_stats$Mean_Dev - summary_stats$SD_Dev*1.15,
                     summary_stats$Mean_Dev + summary_stats$SD_Dev*1.15))
  plot(0, 0, type = "n",
       xlim = x_range,
       ylim = c(0.1, 3.7),          
       xlab = "deviation from true flux (estimated - true)",
       ylab = "",
       yaxt = "n",
       bty = "n")                   
  
  axis(2, at = 3:1, labels = c("Sinai", "Mandeb", "Gibraltar"), tick = FALSE, line = -0.5, cex.axis = 1.45)
  abline(v = 0, col = "gray40", lty = 2, lwd = 2.5)
  col_naive <- "#66c2a5"
  col_friction <- "#fc8d62"
  y_offset <- -0.25              
  
  for(i in 1:nrow(summary_stats)) {
    y <- summary_stats$y[i] + y_offset   
    col_use <- ifelse(summary_stats$Method[i] == "Naive", col_naive, col_friction)
    segments(summary_stats$Mean_Dev[i] - summary_stats$SD_Dev[i], y,
             summary_stats$Mean_Dev[i] + summary_stats$SD_Dev[i], y,
             col = col_use, lwd = 5)
    points(summary_stats$Mean_Dev[i], y, pch = 19, col = col_use, cex = 2.2)
  }

  max_abs <- max(abs(x_range))
  breaks <- seq(-ceiling(max_abs / 0.2) * 0.2,
                ceiling(max_abs / 0.2) * 0.2,
                by = 0.2)
  axis(1, at = breaks, labels = sprintf("%.1f", breaks), cex.axis = 1.35)
  
  legend("topright", inset = c(0, 0.08),  
         legend = c("true value", 
                    "naive: mean ± SD", 
                    "friction: mean ± SD"),
         col = c("gray40", col_naive, col_friction),
         lty = c(2, 1, 1),
         lwd = c(2.5, 5, 5),
         pch = c(NA, 19, 19),
         bty = "n",
         cex = 1.35)
  
  dev.off()
  message(sprintf("saved figure to: %s", output_file))
  
  
  # MSE
  mse_table <- data.frame(
    Route = rep(c("Sinai", "Mandeb", "Gibraltar"), each = 2),
    Method = rep(c("Naive", "Friction"), 3),
    MSE = c(
      mean((x_naive$flux_sinai - true_sinai)^2),
      mean((x_friction$flux_sinai - true_sinai)^2),
      
      mean((x_naive$flux_mandeb - true_mandeb)^2),
      mean((x_friction$flux_mandeb - true_mandeb)^2),
      
      mean((x_naive$flux_gibraltar - true_gibraltar)^2),
      mean((x_friction$flux_gibraltar - true_gibraltar)^2)
    )
  )
  output_file <- sprintf("output/tables/flux_mse_%s_%s_%s_%s_base.csv",WORLD, source_pop, end_gen, MAP)
  write.csv(mse_table, file = output_file, row.names = FALSE)
  message(sprintf("saved MSE table to: %s", output_file))
}
