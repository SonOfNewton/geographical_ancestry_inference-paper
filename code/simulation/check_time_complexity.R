# simple test for time complexity of MPR estimation

WORLD <- "afro-eurasia"
SOURCE_POP <- "58"
END_GEN <- "10000"
REP <- "1"
MAP <- "naive"

ts <- treeseq_load(sprintf("data/trees/empirical_tree_%s_%s_%s_%s.trees",WORLD, SOURCE_POP, END_GEN, REP))
cost.mat <- data.matrix(read.csv(sprintf("data/geo/landgrid_costmat_%s_%s.csv", MAP, WORLD)))
neighbor.mat <- data.matrix(read.csv(sprintf("data/geo/landgrid_adjmat_%s_%s.csv", MAP, WORLD),row.names = 1))

dimnames(cost.mat) <- NULL
dimnames(neighbor.mat) <- NULL

nodes <- treeseq_nodes(ts)

# sample sizes to test
sample_sizes <- c(10, 50, 100, 200, 300, 400, 500)

results <- data.frame(
  sample_size = sample_sizes,
  runtime_sec = NA_real_
)

set.seed(123)
for (i in seq_along(sample_sizes)) {
  k <- sample_sizes[i]
  cat("\n============================\n")
  cat("Sample size:", k, "\n")

  sample_nodes <- nodes[nodes$is_sample == 1, ]
  sample_nodes_to_sample <- sample_nodes[sample(seq_len(nrow(sample_nodes)), k),]
  to_keep <- sample_nodes_to_sample$node_id
  
  ts2 <- treeseq_simplify(
    ts,
    to_keep,
    filter.populations = FALSE
  )
  nodes2 <- treeseq_nodes(ts2)
  sample_nodes2 <- nodes2[nodes2$is_sample == 1, ]
  sample_locations <- cbind(
    node_id = sample_nodes2[, 1],
    state_id = sample_nodes2[, 4] + 1L
  )
  
  # time cost of MPR
  t0 <- Sys.time()
  mpr <- treeseq_discrete_mpr(ts2,sample_locations,cost.mat)
  t1 <- Sys.time()
  
  runtime <- as.numeric(difftime(t1, t0, units = "secs"))
  cat("Time:", runtime, "sec\n")
  results$runtime_sec[i] <- runtime
}

# save results
write.csv(results,file = sprintf("output/tables/mpr_runtime_test_%s_%s_%s.csv",WORLD, MAP, REP),row.names = FALSE)
print(results)

ggplot(results, aes(x = sample_size, y = runtime_sec / 60)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = results$sample_size) +
  labs(
    x = "Sample size",
    y = "Computation time (min)"
  ) +
  theme_minimal(base_size = 14)
