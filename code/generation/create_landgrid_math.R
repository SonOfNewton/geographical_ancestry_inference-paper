# scenarios: line, square, cube, annulus
# each scenario output 4 files:
# adjmat_scenario.csv (homogeneous, for SLiM)
# costmat_scenario_naive.csv (homogeneous, for GAIA)
# costmat_scenario_friction.csv (heterogeneous, for GAIA)
# scenario_weights.csv (heterogeneous, for SLiM)

########## grid generation for line ########## 
# appart from the heterogeneous line model used in figure 2, this part also prepares for the homogeneous line model used in figure 1.

# naive adjcency matrix and cost matrix
n <- 20  # dimension
adj_mat <- matrix(0, nrow = n, ncol = n)
for(i in 1:(n-1)){
  adj_mat[i, i+1] <- 1
  adj_mat[i+1, i] <- 1
}
write.table(adj_mat, file = "data/math/adjmat_line.csv", sep = ",", row.names = FALSE, col.names = FALSE)
cost.mat <- outer(1:n, 1:n, function(i, j) abs(i - j))
write.table(cost.mat, file = "data/math/costmat_line_naive.csv", sep = ",")

# friction cost matrix scenario 1: cost=1 from deme 1 to deme 10, cost=5 from deme 11 to deme 20
adj_mat <- as.matrix(read.csv("data/math/adjmat_line.csv", header = FALSE))
g_friction <- graph_from_adjacency_matrix(adj_mat, mode = "undirected")
E(g_friction)$weight <- 1
friction_range <- 11:21
edges_in_graph <- as_edgelist(g_friction, names = FALSE)
for(i in 1:nrow(edges_in_graph)) {
  node_a <- edges_in_graph[i, 1]
  node_b <- edges_in_graph[i, 2]
  if (node_a %in% friction_range | node_b %in% friction_range) {
    E(g_friction)$weight[i] <- 5.0
  }
}
cost.mat <- distances(g_friction, weights = E(g_friction)$weight)
write.table(cost.mat, file = "data/math/costmat_line_friction.csv", sep = ",")

# weights (friction and naive)
line_weights <- sapply(1:n, function(i) {
  if (i <= 10) {
    return(1)
  } else {
    return(5)
  }
})
write.table(line_weights, file = "data/math/line_weights.csv", row.names = FALSE, col.names = FALSE)
line_weights <- sapply(1:n, function(i) {return(1)})
write.table(line_weights, file = "data/math/line_weights_naive.csv", row.names = FALSE, col.names = FALSE)

# friction cost matrix scenario 2: cost=1, 2, 5, 10 for deme regions 1-5, 6-10, 11-15, 16-20, respectively
adj_mat <- as.matrix(read.csv("data/math/adjmat_line.csv", header = FALSE))
g_friction <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
E(g_friction)$weight <- 1  
get_multiplier <- function(node_id) {
  if (node_id <= 5) {
    return(1)
  } else if (node_id <= 10) {
    return(2)
  } else if (node_id <= 15) {
    return(5)
  } else {
    return(10)
  }
}
edges_in_graph <- igraph::as_edgelist(g_friction, names = FALSE)
for (i in 1:nrow(edges_in_graph)) {
  node_a <- edges_in_graph[i, 1]
  node_b <- edges_in_graph[i, 2]
  mult_a <- get_multiplier(node_a)
  mult_b <- get_multiplier(node_b)
  edge_mult <- max(mult_a, mult_b)
  E(g_friction)$weight[i] <- edge_mult
}
cost.mat <- igraph::distances(g_friction, weights = E(g_friction)$weight)
write.table(cost.mat, file = "data/math/costmat_line_friction_case2.csv", sep = ",")

line_weights <- sapply(1:n, function(i) {
  if (i <= 5) {
    return(1)
  } else if (i <= 10) {
    return(2)
  } else if (i <= 15) {
    return(5)
  } else {
    return(10)
  }
})
write.table(line_weights, file = "data/math/line_weights_case2.csv", row.names = FALSE, col.names = FALSE)


########## grid generation for square ########## 
# naive adjcency matrix and cost matrix
n_rows <- 10
n_cols <- 10
n <- n_rows * n_cols   # dimension: 10 x 10
grid_rows <- ceiling((1:n) / n_cols)
grid_cols <- ((1:n) - 1) %% n_cols + 1
cost.mat <- abs(outer(grid_rows, grid_rows, "-")) + abs(outer(grid_cols, grid_cols, "-"))
write.table(cost.mat, file = "data/math/costmat_square_naive.csv", sep = ",")
adjmat <- (cost.mat == 1) * 1
write.table(adjmat, file = "data/math/adjmat_square.csv", sep = ",", row.names = FALSE, col.names = FALSE)

# friction cost matrix 
get_weight <- function(r, c) {
  if (r <= 5 & c <= 5) return(1)      # top left: 1
  if (r <= 5 & c > 5)  return(2)      # top right: 2
  if (r > 5  & c <= 5) return(3)      # bottom left: 3
  if (r > 5  & c > 5)  return(4)      # bottom right: 4
}
weights_vec <- sapply(1:n, function(i) get_weight(grid_rows[i], grid_cols[i]))
dist_phys <- abs(outer(grid_rows, grid_rows, "-")) + abs(outer(grid_cols, grid_cols, "-"))
adj_mat <- (dist_phys == 1) * 1
g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected")
edges <- as_edgelist(g, names = FALSE)
edge_weights <- sapply(1:nrow(edges), function(i) {
  target_node <- edges[i, 2]
  return(weights_vec[target_node])
})
cost.mat <- distances(g, weights = edge_weights)
write.table(cost.mat, file = "data/math/costmat_square_friction.csv", sep = ",")
write.table(weights_vec, file = "data/math/square_weights.csv", row.names = FALSE, col.names = FALSE)


########## grid generation for cube ########## 
size <- 7
coords <- expand.grid(x = 1:size, y = 1:size, z = 1:size)
# strict sorting: z is first major order, y is second, x is third, ensuring that the ID corresponds to 0:342 in SLiM.
coords <- coords[order(coords$z, coords$y, coords$x), ]
rownames(coords) <- NULL
N <- nrow(coords) # 343

# divide into 8 blocks and assign weights
get_3d_weight <- function(x, y, z) {
  xi <- if(x <= 4) 0 else 1
  yi <- if(y <= 4) 0 else 2
  zi <- if(z <= 4) 0 else 4
  return(xi + yi + zi + 1)
}
coords$weight <- mapply(get_3d_weight, coords$x, coords$y, coords$z)
write.table(coords$weight, "data/math/cube_weights.csv", row.names=FALSE, col.names=FALSE)

adjmat <- matrix(0, nrow=N, ncol=N)
edges_naive <- list()
edges_friction <- list()

for (i in 1:N) {
  for (j in 1:N) {
    if (i == j) next
    dist <- sum(abs(coords[i,1:3] - coords[j,1:3]))
    if (dist == 1) {
      adjmat[i, j] <- 1
      edges_naive <- append(edges_naive, list(c(i, j, 1)))
      edges_friction <- append(edges_friction, list(c(i, j, coords$weight[i])))
    }
  }
}
write.table(adjmat, "data/math/adjmat_cube.csv", sep=",", row.names=FALSE, col.names=FALSE)   # naive adjcency matrix

calc_cost <- function(edge_list) {
  el <- do.call(rbind, edge_list)
  g <- graph_from_data_frame(d = data.frame(from=el[,1], to=el[,2], weight=el[,3]), directed=TRUE)
  return(distances(g, mode="out"))
}
cost_naive <- calc_cost(edges_naive)
cost_friction <- calc_cost(edges_friction)
write.table(cost_naive, "data/math/costmat_cube_naive.csv", sep=",")   #naive cost matrix
write.table(cost_friction, "data/math/costmat_cube_friction.csv", sep=",")   #friction cost matrix


########## grid generation for annulus ########## 
n_rows <- 10
n_cols <- 10
grid <- expand.grid(row = 1:n_rows,col = 1:n_cols)

# remove 4 x 4 nodes in the middle
grid <- subset(
  grid,
  !(row %in% 4:7 & col %in% 4:7)
)
n <- nrow(grid)   # 84
grid$id <- 1:n

# naive adjmat and costmat
adjmat <- matrix(0, n, n)
for(i in 1:n){
  for(j in 1:n){
    d <- abs(grid$row[i]-grid$row[j]) +
      abs(grid$col[i]-grid$col[j])
    if(d == 1)
      adjmat[i,j] <- 1
  }
}
write.table(adjmat, file = "data/math/adjmat_annulus.csv",sep = ",", row.names = FALSE, col.names = FALSE)
g <- graph_from_adjacency_matrix(adjmat, mode = "undirected")
cost.mat <- distances(g)
write.table(cost.mat, file = "data/math/costmat_annulus_naive.csv", sep = ",")

# friction weights and costmat
get_weight <- function(r, c) {
  if (r <= 5 & c <= 5) return(1)      # top left: 1
  if (r <= 5 & c > 5)  return(2)      # top right: 2
  if (r > 5  & c <= 5) return(3)      # bottom left: 3
  if (r > 5  & c > 5)  return(4)      # bottom right: 4
}
weights_vec <- mapply(get_weight, grid$row, grid$col)
write.table(weights_vec,"data/math/annulus_weights.csv",sep = ",",row.names = FALSE,col.names = FALSE)

edges <- as_edgelist(g, names = FALSE)
edge_weights <- sapply(
  1:nrow(edges),
  function(i){
    target_node <- edges[i,2]
    weights_vec[target_node]
  }
)
cost.mat <- distances(g,weights = edge_weights)
write.table(cost.mat,"data/math/costmat_annulus_friction.csv",sep = ",")

