process_one_file <- function(i, mode, model, real_case) {
  
  ts_file <- sprintf("data/trees/math-tree-%s-%s-%s.trees", mode, real_case, i)
  cost_mat_error <- data.matrix(read.csv(sprintf("data/math/costmat_%s_naive.csv", mode)))
  cost.mat <- data.matrix(read.csv(sprintf("data/math/costmat_%s_%s.csv", mode, model)))
  
  if (mode == "square"){
    MAXDEME = 100
  } else if (mode == "line"){
    MAXDEME = 20
  } else if (mode == "cube"){
    MAXDEME = 343
  }
  
  ts <- treeseq_load(ts_file)
  nodes <- treeseq_nodes(ts)
  edges <- treeseq_edges(ts)
  
  sample_nodes <- nodes[nodes$is_sample == 1, ]
  populated_demes <- unique(sample_nodes$population_id)
  
  set.seed(0)
  num_demes_to_sample <- min(20, MAXDEME)
  pops_to_sample <- sample(populated_demes, num_demes_to_sample, replace = FALSE)
  sample_nodes_to_sample <- sample_nodes[sample_nodes$population_id %in% pops_to_sample, ]
  
  safe_sample <- function(x) {
    if (length(x) == 1) return(x)
    return(sample(x, 1, replace = FALSE))
  }
  
  indivs <- unlist(unname(tapply(
    sample_nodes_to_sample$individual_id,
    sample_nodes_to_sample$population_id,
    safe_sample
  )))
  
  to_keep <- sample_nodes_to_sample$node_id[sample_nodes_to_sample$individual_id %in% indivs]
  ts2 <- treeseq_simplify(ts, to_keep, filter.populations = FALSE)
  
  nodes2 <- treeseq_nodes(ts2)
  sample_nodes2 <- nodes2[nodes2$is_sample == 1, ]
  sample_locations <- cbind(node_id = sample_nodes2[, 1], state_id = sample_nodes2[, 4] + 1L)
  
  nodes_df <- treeseq_nodes(ts2)
  edges_tree1 <- treeseq_edges(ts2)
  
  ### GAIA
  mpr <- treeseq_discrete_mpr(ts2, sample_locations, cost.mat)
  estimated_node_states = treeseq_discrete_mpr_minimize(mpr)
  real_node_states = nodes_df$population_id + 1L
  
  comparison_df <- data.frame(
    node_id = nodes_df$node_id,
    time = nodes_df$time,
    real_state = real_node_states,
    est_state = estimated_node_states
  )
  
  comparison_df <- comparison_df[order(comparison_df$time, decreasing = TRUE), ]
  comparison_df$distance <- cost_mat_error[cbind(comparison_df$real_state, comparison_df$est_state)]
  comparison_df$time_bin <- cut(comparison_df$time, breaks = 4, include.lowest = TRUE)
  
  node_counts <- aggregate(distance ~ time_bin, data = comparison_df, FUN = length)
  mean_errors <- aggregate(distance ~ time_bin, data = comparison_df, FUN = function(x) mean(x, na.rm = TRUE))
  
  node_counts <- node_counts[order(node_counts$time_bin, decreasing = TRUE), ]
  mean_errors <- mean_errors[order(mean_errors$time_bin, decreasing = TRUE), ]
  
  row_data <- data.frame(
    file_id = i,
    count_bin1_old = node_counts$distance[1],
    count_bin2 = node_counts$distance[2],
    count_bin3 = node_counts$distance[3],
    count_bin4_new = node_counts$distance[4],
    error_bin1_old = mean_errors$distance[1],
    error_bin2 = mean_errors$distance[2],
    error_bin3 = mean_errors$distance[3],
    error_bin4_new = mean_errors$distance[4]
  )
  
  message(paste("Processed file:", i, " | Mode:", mode, " | Model:", model))
  return(row_data)
}