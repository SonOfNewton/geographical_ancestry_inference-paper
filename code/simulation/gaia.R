#! /usr/local/bin/Rscript --vanilla

# this script runs in a new R session so need to load library of its own
suppressPackageStartupMessages({   
  library(gaia)
  library(here)
})
source(here("code","functions.R"),verbose=FALSE)

args = commandArgs(trailingOnly=TRUE)
REP = args[1]
WORLD = args[2]  # afro-eurasia
MAP = args[3]    # friction, naive

# grid cells to sample
data = read.csv(sprintf("data/genetics/sample_states_%s.csv", WORLD))
pops_to_sample = sort(unique(data[,2])-1L)

if (MAP=="friction"){
cost.mat = data.matrix(read.csv("data/geo/landgrid_costmat_friction_afro-eurasia.csv"))
neighbor.mat = data.matrix(read.csv("data/geo/landgrid_adjmat_friction_afro-eurasia.csv"))
} else if (MAP=="naive"){
cost.mat = data.matrix(read.csv("data/geo/landgrid_costmat_naive_afro-eurasia.csv"))
neighbor.mat = data.matrix(read.csv("data/geo/landgrid_adjmat_afro-eurasia.csv", row.names=1))  #need to specify row names
}
dimnames(cost.mat) = NULL
dimnames(neighbor.mat) = NULL

ts = treeseq_load(sprintf("data/trees/empirical_tree_%s-%s.trees", WORLD, REP))
nodes = treeseq_nodes(ts)
sample_nodes = nodes[nodes$is_sample==1, ]
sample_nodes_to_sample = sample_nodes[sample_nodes$population_id %in% pops_to_sample,]
set.seed(as.integer(REP))  # make sure the same people are sampled for naive and friction in the same run
indivs = unname(tapply(sample_nodes_to_sample$individual_id, sample_nodes_to_sample$population_id, sample, 1))
to_keep = sample_nodes_to_sample$node_id[sample_nodes_to_sample$individual_id %in% indivs]

ts2 = treeseq_simplify(ts, to_keep, filter.populations=FALSE)
nodes = treeseq_nodes(ts2)
sample_nodes = nodes[nodes$is_sample==1, ]
sample_locations = cbind(node_id=sample_nodes[,1], state_id=sample_nodes[,4]+1L)

mpr = treeseq_discrete_mpr(ts2, sample_locations, cost.mat)
estimated_node_states = treeseq_discrete_mpr_minimize(mpr)

write.csv(data.frame(node_time=nodes$time, node_state=nodes$population_id+1L, 
    estimated_node_state=estimated_node_states), 
    file=sprintf("data/mpr/mpr_%s_%s_%s.csv", WORLD, MAP, REP), row.names=FALSE)

sample_sets = rep(1L, treeseq_num_samples(ts2))
state_sets = 1:nrow(cost.mat)

flux = treeseq_discrete_mpr_ancestry_flux(ts2, mpr, cost.mat,
    neighbor.mat, c(0, 20000), state_sets, sample_sets)
attr(flux, "sites") = unique(sample_locations[,2])

ts3 = treeseq_simplify(ts, to_keep, filter.populations=FALSE, keep.unary=TRUE)
nodes0 = treeseq_nodes(ts3)

flux0 = true_flux(nodes0$population_id, ts3, c(0,20000), cost.mat, neighbor.mat,
    sample_sets, state_sets)   

# option 1: consider a->b and b->a as separate flux
f = data.frame(true_flux=c(flux0), estimated_flux=c(flux))

# option 2: consider a->b as a flux with direction ( - b->a )
# true_mat = flux0[,,1,1]
# est_mat = flux[,,1,1]
# net_true_mat = true_mat - t(true_mat)
# net_est_mat = est_mat - t(est_mat)
# upper_idx = upper.tri(net_true_mat)
# f = data.frame(
#   true_flux = net_true_mat[upper_idx],
#   estimated_flux = net_est_mat[upper_idx]
# )

write.csv(f, file=sprintf("data/flux/flux_%s_%s_%s.csv", WORLD, MAP, REP))

if (WORLD=="afro-eurasia"){
  ooa_flux = data.frame(
    rep=REP,
    flux_sinai=sum(flux[38,69:70,1,1]),
    flux_mandeb=flux[67,68,1,1],
    flux_gibraltar=flux[155,7,1,1])
  
  out_file = sprintf("data/flux/flux_strait_%s_%s_%s.csv", WORLD, MAP, REP)
  write.csv(ooa_flux, file=out_file, row.names=FALSE)
}



