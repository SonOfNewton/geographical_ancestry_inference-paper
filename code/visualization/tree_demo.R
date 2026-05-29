# parameters
levels <- 5
nodes_per_level <- 2^(0:(levels-1))

node_id <- list()
counter <- 1
for (l in 1:levels) {
  node_id[[l]] <- counter:(counter + nodes_per_level[l] - 1)
  counter <- counter + nodes_per_level[l]
}

# recursive calculation of x coordinate
coords <- list()

# distribution of leaves
leaf_x <- seq(0, 1, length.out = nodes_per_level[levels])
coords[[levels]] <- cbind(leaf_x, rep(0, nodes_per_level[levels]))

# bottom up
for (l in (levels-1):1) {
  n <- nodes_per_level[l]
  xy <- matrix(NA, n, 2)
  for (i in 1:n) {
    left  <- coords[[l+1]][2*i-1, ]
    right <- coords[[l+1]][2*i, ]
    xy[i,1] <- mean(c(left[1], right[1]))
    xy[i,2] <- levels - l
  }
  coords[[l]] <- xy
}

labels <- list(
  c(1),
  c(1,2),
  c(1,1,2,3),
  c(1,1,1,1,2,2,3,4),
  c(1,1,1,1,1,1,1,1,2,2,2,2,3,3,4,5)
)

highlight_leaves <- c(1,9,13,15,16)
get_path <- function(idx) {
  path <- c()
  for (l in levels:1) {
    path <- c(idx, path)
    idx <- ceiling(idx/2)
  }
  path
}

edge_key <- function(l,p,c) paste(l,p,c)
highlight_edges <- c()
for (leaf in highlight_leaves) {
  p <- get_path(leaf)
  for (l in 1:(levels-1)) {
    highlight_edges <- c(highlight_edges, edge_key(l,p[l],p[l+1]))
  }
}

# plot
png("output/figures/tree-demo.png", width=1800, height=1800, res=300)
plot(NULL, xlim=c(0,1), ylim=c(0,levels), axes=FALSE, xlab="", ylab="")

# nodes on true path
true_nodes <- list()
for (l in 1:levels) true_nodes[[l]] <- rep(FALSE, nodes_per_level[l])

for (leaf in highlight_leaves) {
  path <- get_path(leaf)
  for (l in 1:levels) {
    true_nodes[[l]][path[l]] <- TRUE
  }
}

# edges
for (l in 1:(levels-1)) {
  for (i in 1:nodes_per_level[l]) {
    parent <- coords[[l]][i,]
    for (k in 0:1) {
      child_idx <- 2*i - 1 + k
      child <- coords[[l+1]][child_idx,]
      
      key <- edge_key(l, i, child_idx)
      is_highlight <- key %in% highlight_edges
      
      segments(parent[1], parent[2], child[1], child[2],
               col = if (is_highlight) "black" else "grey80",
               lwd = if (is_highlight) 2.5 else 1)
    }
  }
}

# nodes
for (l in 1:levels) {
  for (i in 1:nodes_per_level[l]) {
    points(coords[[l]][i,1], coords[[l]][i,2],
           pch=21,
           bg = if (true_nodes[[l]][i]) "black" else "grey80",
           col = NA,
           cex=1.2)
  }
}

# true labels
for (l in 1:levels) {
  text(coords[[l]],
       labels=labels[[l]],
       pos=3,
       cex=1.1,   # 🔥 放大
       col="black")
}

# estimated labels
x_offset <- 0.05
y_offset <- 0.2

text(coords[[1]][1,1] + x_offset, coords[[1]][1,2] + y_offset,
     labels="1,2", col="red", cex=1.2, font=2)

text(coords[[2]][2,1] + x_offset, coords[[2]][2,2] + y_offset,
     labels="2", col="red", cex=1.2, font=2)

text(coords[[3]][4,1] + x_offset, coords[[3]][4,2] + y_offset,
     labels="3", col="red", cex=1.2, font=2)

text(coords[[4]][8,1] + x_offset, coords[[4]][8,2] + y_offset,
     labels="4", col="red", cex=1.2, font=2)

legend("topright",
       legend=c("True", "Estimated", "Not sampled"),
       lty=c(1,1,1),
       lwd=c(2.5,2.5,2.5),
       col=c("black","red","grey80"),
       bty="n",
       cex=1.2)

dev.off()

