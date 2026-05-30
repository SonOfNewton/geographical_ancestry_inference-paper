# this script runs on its own (without .sh file)

##### parameters #####
num_subsets = 100L
map = c("naive", "friction")

run_mpr <- function(map, num_subsets) {
cost.mat = data.matrix(read.csv(sprintf("data/geo/landgrid_costmat_%s_asia-americas.csv", map)))

process_mpr_subset <- function(i) {
  ts_file <- sprintf("data/trees/empirical_tree_asia-americas_subset_%d.trees", i)
  data_file <- sprintf("data/genetics/subsets/sample_states_asia-americas_subset_%d.csv", i)
  
  # independently load trees in each process
  ts <- treeseq_load(ts_file)
  data <- data.matrix(read.csv(data_file))
  
  mpr <- treeseq_discrete_mpr(ts, data, cost.mat)   ### GAIA
  # rds_filename <- sprintf("data/mpr/mpr_asia-americas_%s_%d.rds", map, i)
  # saveRDS(mpr, file = rds_filename)

  estimated_node_states <- treeseq_discrete_mpr_minimize(mpr)
  nodes <- treeseq_nodes(ts)

  csv_filename <- sprintf("data/mpr/mpr_asia-americas_%s_%d.csv", map, i)
  write.csv(
    data.frame(node_time = nodes$time, estimated_node_state = estimated_node_states), 
    file = csv_filename, 
    row.names = FALSE
  )
  
  # clear memory
  rm(ts, data, mpr, estimated_node_states, nodes)
  gc(verbose = FALSE, full = TRUE)
  
  return(i)
}

# use mclapply to run in parallel
system.time({
  results <- mclapply(
    1:num_subsets, 
    process_mpr_subset, 
    mc.cores = usecores, 
    mc.preschedule = FALSE   # dynamic memory allocation
  )
})
}

for (mymap in map){
  cat("running mpr estimation for the ", mymap, "case ...")
  run_mpr(mymap, num_subsets)
}
