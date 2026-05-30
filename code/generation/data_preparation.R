if (my_world == "afro-eurasia"){
  ##### generate geography #####
  # we take the landgrid of afroeurasia from GAIA paper so no need to generate
  plot_landgrid_visualization(world = my_world)
  
  adjmat = read.csv("data/geo/landgrid_adjmat_naive_afro-eurasia.csv", row.names = 1)
  # preparation for SLiM: need column names but delete row names (difficult to remove headers in Eidos)
  write.table(adjmat, "data/geo/landgrid_adjmat_naive_afro-eurasia_norownames.csv",sep = ",", row.names = FALSE)
  # create naive costmat according to naive adjmat
  g <- graph_from_adjacency_matrix(as.matrix(adjmat), mode = "undirected", diag = FALSE)
  costmat <- distances(g, weights = NA)
  write.table(costmat, "data/geo/landgrid_costmat_naive_afro-eurasia.csv",sep = ",")
  
  ##### friction map #####
  create_friction_map(world = my_world, visualize = TRUE)
  
  ##### population capacity #####
  # choice of parameter time: "10000BC", "5000BC", "0AD", "1500AD", "2000AD"
  # choice of parameter scaling: "uniform" "sqrt" "complex" "direct"
  create_population_capacity(world = my_world, time="5000BC", scaling = "sqrt", visualize = TRUE)

} else if (my_world == "asia-americas"){
  ##### generate geography #####
  source(here("code", "generation", "create_landgrid.R"),verbose=FALSE)
  plot_landgrid_visualization(world = my_world)
  
  ##### friction map #####
  create_friction_map(world = my_world, visualize = TRUE)
  
  
} else{
  message(sprintf("no such world setting!"))
}


