ts_file <- "data/trees/math-tree-line-naive-1.trees"
ts <- treeseq_load(ts_file)
nodes <- treeseq_nodes(ts)
edges <- treeseq_edges(ts)
#unique(edges$left)   # check recombinations

MAXDEME = 20
cost.mat = data.matrix(read.csv("data/math/costmat_line_naive.csv"))

sample_nodes <- nodes[nodes$is_sample == 1, ]   # all contemporary nodes
pops_to_sample <- unique(sample_nodes$population_id)   # sample from each populated deme

safe_sample <- function(x) {
  if (length(x) == 1) return(x)
  return(sample(x, 1, replace = FALSE))  
}

set.seed(0) 
# sample 1 individual from each deme
indivs <- unlist(unname(tapply(
  sample_nodes$individual_id, 
  sample_nodes$population_id, 
  safe_sample
)))

# extract node id of the sampled individuals and simplify the tree sequence
to_keep <- sample_nodes_to_sample$node_id[sample_nodes_to_sample$individual_id %in% indivs]
ts2 <- treeseq_simplify(ts, to_keep, filter.populations = FALSE)

nodes2 <- treeseq_nodes(ts2)
sample_nodes2 <- nodes2[nodes2$is_sample == 1, ]
sample_locations <- cbind(node_id = sample_nodes2[, 1], state_id = sample_nodes2[, 4] + 1L)

nodes_df <- treeseq_nodes(ts2)
edges_tree1 <- treeseq_edges(ts2)

### GAIA
mpr <- treeseq_discrete_mpr(ts2, sample_locations, cost.mat)  #, use_brlen=TRUE)
estimated_node_states = treeseq_discrete_mpr_minimize(mpr)
real_node_states = nodes_df$population_id + 1L
comparison_df <- data.frame(
  node_id = nodes_df$node_id,
  time = nodes_df$time,
  real_state = real_node_states,
  est_state = estimated_node_states
)

# visualize tree structure
target_pos <- 0    # position on the genome (structure of the tree sequence differs for different genome positions)
local_edges <- edges_tree1[edges_tree1$left <= target_pos & edges_tree1$right > target_pos, ]

g_simple <- graph_from_data_frame(
  data.frame(from = local_edges$parent, to = local_edges$child)
)

node_indices <- match(as.numeric(V(g_simple)$name), nodes_df$node_id)
times <- nodes_df$time[node_indices]
real_demes <- nodes_df$population_id[node_indices] + 1L
est_demes <- estimated_node_states[node_indices]

lay <- layout_as_tree(g_simple)
lay[, 2] <- times  # y axis represent real time in history

# extract and sort contemporary samples
out_degrees <- igraph::degree(g_simple, mode = "out")
leaf_idx <- which(out_degrees == 0)
sorted_leaf_idx <- leaf_idx[order(real_demes[leaf_idx])]
lay[sorted_leaf_idx, 1] <- seq_along(sorted_leaf_idx)

# sort internal nodes
internal_idx <- which(out_degrees > 0)
internal_idx <- internal_idx[order(times[internal_idx])]

# Bottom-up approach: Place the X-coordinate of each ancestor node in the center of the X-coordinates of all its child nodes
for (i in internal_idx) {
  children <- neighbors(g_simple, v = i, mode = "out")
  children_idx <- match(names(children), V(g_simple)$name)
  lay[i, 1] <- mean(lay[children_idx, 1])
}

soft_colors <- c("#4575b4", "#e0f3f8", "#fee090", "#d73027")
pal <- colorRampPalette(soft_colors)(MAXDEME)

V(g_simple)$real_color <- pal[real_demes]
V(g_simple)$est_color <- pal[est_demes]  

is_sample_node <- nodes_df$is_sample[node_indices] == 1
V(g_simple)$frame_color <- ifelse(is_sample_node, "black", "gray80")

# plot
png(file="output/figures/line_tree_structure.png", width=10, height=5, units="in", res=300)
par(mfrow = c(1, 2), mar = c(3, 4, 4, 1))  
# left panel: actual deme
plot(g_simple, layout = lay, vertex.color = V(g_simple)$real_color,
     vertex.frame.color = V(g_simple)$frame_color, vertex.size = 10,
     vertex.label = NA, edge.arrow.size = 0.2, edge.color = "gray60",
     #main = paste("Local Tree at position:", target_pos, "\n(Real Node States)")
     )
# right panel: estimated deme
par(mar = c(3, 1, 4, 4))
plot(g_simple, layout = lay, vertex.color = V(g_simple)$est_color,
     vertex.frame.color = V(g_simple)$frame_color, vertex.size = 10,
     vertex.label = NA, edge.arrow.size = 0.2, edge.color = "gray60",
     #main = paste("Local Tree at position:", target_pos, "\n(Estimated Node States)")
     )
legend_demes <- unique(round(seq(1, MAXDEME, length.out = 5)))
legend("topright", title = "Deme", legend = legend_demes, 
       fill = pal[legend_demes], bty = "n", cex = 1.2, inset = 0.02)
dev.off()



