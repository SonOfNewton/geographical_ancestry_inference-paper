# note: this script has not been parallelized yet

##### parameters #####
world <- "afro-eurasia" 
num_reps <- 200   #number of SLiM runs, controlled in slim_empirical.sh

cost.mat = data.matrix(read.csv("data/geo/landgrid_costmat_naive_afro-eurasia.csv", row.names=1))
neighbor.mat = data.matrix(read.csv("data/geo/landgrid_adjmat_naive_afro-eurasia.csv", row.names=1))
dimnames(cost.mat) = NULL
dimnames(neighbor.mat) = NULL
state_sets = 1:nrow(cost.mat)   # each grid represents one state (another choice is setting multiple grid to the same state to account for them as a whole)

# cells to sample from (originate from make-data.R in GAIA paper)
set.seed(1123581321L)
D = data_full(chromosome=18L, world="afro-eurasia", short=TRUE)
data = D$data
write.csv(D$data, file=sprintf("data/genetics/sample_states_%s.csv", world), row.names=FALSE)   # cell ID for all sampled individuals

pops_to_sample = sort(unique(data[,2])-1L)

# extract true flux from all replications
results_list <- list()

for (i in 1:num_reps) {
  ts_file <- sprintf("data/trees/empirical_tree_afro-eurasia-%s.trees", i)
  # skip if missing some files
  if (!file.exists(ts_file)) {
    cat(sprintf("Rep %d file does not exist, skipped\n", i))
    next
  }
  
  ts <- treeseq_load(ts_file)
  nodes <- treeseq_nodes(ts)
  sample_nodes <- nodes[nodes$is_sample==1, ]
  sample_nodes_to_sample <- sample_nodes[sample_nodes$population_id %in% pops_to_sample,]
  
  set.seed(1024)  # make sure the same individuals are selected for all replications so that the tree sequence is comparable
  indivs <- unname(tapply(sample_nodes_to_sample$individual_id, sample_nodes_to_sample$population_id, sample, 1))
  to_keep <- sample_nodes_to_sample$node_id[sample_nodes_to_sample$individual_id %in% indivs]
  
  # prune tree
  ts3 <- treeseq_simplify(ts, to_keep, filter.populations=FALSE, keep.unary=TRUE)
  nodes0 <- treeseq_nodes(ts3)
  sample_sets <- rep(1L, treeseq_num_samples(ts3))
  
  # true flux
  flux0 <- true_flux(nodes0$population_id, ts3, c(0,20000), cost.mat, neighbor.mat, sample_sets, state_sets)
  true_mat <- flux0[,,1,1] 
  
  # true flux across the three corridors
  sinai_val <- sum(true_mat[38, 69:70])
  mandeb_val <- true_mat[67, 68]
  gibraltar_val <- true_mat[155, 7]
  
  results_list[[i]] <- data.frame(rep = i, sinai = sinai_val, mandeb = mandeb_val, gibraltar = gibraltar_val)
  if (i %% 10 == 0) cat(sprintf("processed %d / %d ...\n", i, num_reps))
}

# merge all replications
all_flux <- bind_rows(results_list)
write.csv(all_flux, sprintf("output/tables/all_replication_true_flux_%s.csv", world), row.names=FALSE)

# sample replications
n_each <- 10  # select n_each replications for each corridor with the highest true flux

top_sinai <- all_flux %>%
  arrange(desc(sinai)) %>%
  slice(1:n_each)

top_mandeb <- all_flux %>%
  arrange(desc(mandeb)) %>%
  slice(1:n_each)

top_gibraltar <- all_flux %>%
  arrange(desc(gibraltar)) %>%
  slice(1:n_each)

# merge and remove duplicates
selected_reps <- bind_rows(top_sinai, top_mandeb, top_gibraltar) %>%
  distinct(rep, .keep_all = TRUE)
selected_reps <- selected_reps %>% arrange(rep)   # sort
selection_file <- sprintf("output/tables/selected_reps_%s.csv", world)
write.csv(selected_reps$rep, selection_file, row.names = FALSE)
cat(sprintf("sampled %d replications\n", nrow(selected_reps)))

# visualize
flux_features <- selected_reps[, c("sinai", "mandeb", "gibraltar")]
dist_matrix <- dist(flux_features)
mds_1d <- cmdscale(dist_matrix, k = 1)
smooth_order <- order(mds_1d[, 1])
flux_sorted <- selected_reps[smooth_order, ]
n_reps <- nrow(flux_sorted)

output_file <- sprintf("output/figures/true_flux_sampled_replications_%s.pdf", world)
pdf(file=output_file, width=10, height=6)
par(mfrow=c(1, 1), mar=c(5.5, 4, 4, 2) + 0.1) 
max_y <- max(c(flux_sorted$sinai, flux_sorted$mandeb, flux_sorted$gibraltar)) * 1.1
plot(0, 0, xlim=c(1, n_reps), ylim=c(0, max_y), type='n', las=1, xaxt='n', 
     xlab="", ylab="true ancestry flux")
mtext("replication ID", side=1, line=4)
axis(1, at=1:n_reps, labels=flux_sorted$rep, cex.axis=0.8, las=2)   # original ID

# Sinai
lines(1:n_reps, flux_sorted$sinai, col="black", lty=1, lwd=1.5)
points(1:n_reps, flux_sorted$sinai, pch=19, col="black", cex=1)

# Mandeb
lines(1:n_reps, flux_sorted$mandeb, col="brown", lty=1, lwd=1.5)
points(1:n_reps, flux_sorted$mandeb, pch=19, col="brown", cex=1)

# Gibraltar
lines(1:n_reps, flux_sorted$gibraltar, col="orange", lty=1, lwd=1.5)
points(1:n_reps, flux_sorted$gibraltar, pch=19, col="orange", cex=1)

#title("true flux of sampled replications")
legend("topright", legend=c("Sinai", "Mandeb", "Gibraltar"), 
       col=c("black", "brown", "orange"), lty=1, pch=19, bty='n')

dev.off()
message(sprintf("saved figure to: %s", output_file))
