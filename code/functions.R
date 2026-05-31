# compute estimation error under a specific model setting
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
  }else if (mode == "annulus"){
    MAXDEME = 84
  }
  
  ts <- treeseq_load(ts_file)
  nodes <- treeseq_nodes(ts)
  edges <- treeseq_edges(ts)
  
  sample_nodes <- nodes[nodes$is_sample == 1, ]
  populated_demes <- unique(sample_nodes$population_id)
  
  set.seed(0)
  num_demes_to_sample <- min(20, MAXDEME)   # sample 20 demes or less
  pops_to_sample <- sample(populated_demes, num_demes_to_sample, replace = FALSE)
  sample_nodes_to_sample <- sample_nodes[sample_nodes$population_id %in% pops_to_sample, ]
  
  safe_sample <- function(x) {
    if (length(x) == 1) return(x)
    return(sample(x, 1, replace = FALSE))
  }
  
  # sample 1 individual from each selected deme
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
  comparison_df$distance <- cost_mat_error[cbind(comparison_df$real_state, comparison_df$est_state)]  # piecewise average error distance
  comparison_df$time_bin <- cut(comparison_df$time, breaks = 4, include.lowest = TRUE)  # bins with equal time span
  
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


# visualize landgrid together with cell ID and adjacency
plot_landgrid_visualization <- function(world = "afro-eurasia", label_cex = 0.7, point_cex = 0.5) {
  landgrid = st_read(sprintf("data/geo/landgrid_%s.gpkg", world), quiet=TRUE)
  adjmat = data.matrix(read.csv(sprintf("data/geo/landgrid_adjmat_naive_%s.csv", world), row.names=1))
  
  centroids = st_centroid(landgrid$geom)
  coords = st_coordinates(centroids)
  output_file <- sprintf("output/figures/landgrid_visualization_%s.pdf", world)
  pdf(output_file, width=16, height=8)
  par(mfrow=c(1,2), mar=c(1,1,3,1))
  
  # left panel: cell id
  plot(landgrid$geom, border="gray50", col="gray95", lwd=0.5, main="landgrid cells with ID numbers")
  text(coords[,1], coords[,2], labels=1:nrow(coords), 
       #col="red", cex=0.3, font=1)
       col="red", cex=label_cex, font=2)
  
  # right panel: adjacency
  plot(landgrid$geom, border="gray80", col="white", lwd=0.5, main="adjacency network")
  num_cells = nrow(adjmat)
  for (i in 1:(num_cells-1)) {
    for (j in (i+1):num_cells) {
      if (adjmat[i, j] == 1) {
        lines(c(coords[i,1], coords[j,1]), 
              c(coords[i,2], coords[j,2]), 
              col="#0000FF55", lwd=1)
      }
    }
  }
  points(coords[,1], coords[,2], pch=20, col="black", cex=point_cex)
  
  dev.off()
  message(sprintf("saved figure to: %s", output_file))
}


# create cost matrix and (weighted) adjacency matrix for a given world setting (plus visualization)
create_friction_map <- function(world = "afro-eurasia", visualize = FALSE) {
  landgrid = st_read(sprintf("data/geo/landgrid_%s.gpkg", world), quiet=TRUE)
  adjmat = data.matrix(read.csv(sprintf("data/geo/landgrid_adjmat_naive_%s.csv", world), row.names=1))
  
  # MAP friction surface
  friction_raster = rast("data/geo/2020_walking_only_friction_surface.geotiff")
  
  # transform the coordinate system of the polygonal mesh to be consistent with that of the raster image
  landgrid_proj = st_transform(landgrid, crs(friction_raster))
  
  # map pixels into grid cells
  friction_extract = terra::extract(friction_raster, vect(landgrid_proj), fun=mean, na.rm=TRUE)
  cell_friction = friction_extract[, 2]   # second column is average
  
  # compute migration cost between adjacent cells
  num_cells = nrow(adjmat)
  raw_adj_cost_matrix = matrix(0, nrow=num_cells, ncol=num_cells)
  
  for (i in 1:num_cells) {
    for (j in 1:num_cells) {
      if (adjmat[i, j] == 1) {
        # distance in meters
        dist_meters = as.numeric(st_distance(st_centroid(landgrid$geom[i]), st_centroid(landgrid$geom[j])))
        # cost = average friction (min/meter) × actual distance (meter)
        raw_adj_cost_matrix[i, j] = ((cell_friction[i] + cell_friction[j]) / 2) * dist_meters
      }
    }
  }
  
  # compute migration cost between all cells
  g = graph_from_adjacency_matrix(raw_adj_cost_matrix, mode="undirected", weighted=TRUE)
  full_cost_matrix = distances(g, mode="all")
  colnames(full_cost_matrix) = 1:num_cells
  rownames(full_cost_matrix) = 1:num_cells
  write.csv(full_cost_matrix, file=sprintf("data/geo/landgrid_costmat_friction_%s.csv", world), row.names=FALSE)
  write.csv(raw_adj_cost_matrix, file=sprintf("data/geo/landgrid_adjmat_friction_%s.csv",world), row.names=FALSE)
  message("matrices saved")
  
  # visualization
  if (visualize == TRUE){
    landgrid$friction = cell_friction
    coords = st_coordinates(st_centroid(landgrid$geom))
    
    p = ggplot(data = landgrid) +
      geom_sf(aes(fill = friction), color = "gray20", linewidth = 0.1) +
      scale_fill_viridis(option = "magma", direction = -1, name = "friction\n(mins/meter)") +
      theme_void() +
      theme(
        text = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 8)),
        plot.subtitle = element_text(size = 16, hjust = 0.5, margin = margin(b = 12)),
        legend.position = "right",
        legend.key.size = unit(1.4, "cm"),          
        legend.key.height = unit(1.8, "cm"),       
        legend.key.width  = unit(1.1, "cm"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17, face = "bold", margin = margin(b = 10)),
        legend.spacing.y = unit(0.6, "cm"),         
        legend.margin = margin(t = 10, r = 15, b = 10, l = 5), 
        plot.margin = margin(10, 20, 10, 10)
      )
    
    output_file = sprintf("output/figures/landgrid_friction_heatmap_%s.pdf", world)
    pdf(output_file, width = 12, height = 8)
    print(p)
    dev.off()
    message(sprintf("saved figure to: %s", output_file))
  }
}


# create population capacity data for a given time (plus visualization)
create_population_capacity <- function(world = "afro-eurasia", time="5000BC", scaling = "sqrt", visualize = FALSE){
  pop_raster = rast(sprintf("data/pop/popc_%s.asc", time))
  landgrid = st_read(sprintf("data/geo/landgrid_%s.gpkg", world), quiet=TRUE)
  
  landgrid_proj = st_transform(landgrid, crs(pop_raster))   # align coordinate system
  pop_extract = terra::extract(pop_raster, vect(landgrid_proj), fun=sum, na.rm=TRUE)   # population within each grid cell
  raw_cell_pop = pop_extract[, 2]
  raw_cell_pop[is.na(raw_cell_pop)] <- 0
  
  # population scaling
  if(scaling == "uniform"){
    TARGET_GLOBAL_K = 100000   # larger population size will result in heavier computation costs in tree sequence generation and simulation
    scale_factor = TARGET_GLOBAL_K / sum(raw_cell_pop)
    scaled_K = raw_cell_pop * scale_factor
  } else if(scaling == "sqrt"){
    scaled_K = sqrt(raw_cell_pop)/1    # magic number used to control computation complexity
  } else if(scaling == "complex"){
    scaled_K = raw_cell_pop/(0.1*sqrt(sum(raw_cell_pop)))    # magic number used to control computation complexity
  } else if (scaling =="direct"){   #only for 5000BC
    scaled_K = raw_cell_pop/100
  }
  
  final_K = round(scaled_K)
  final_K[raw_cell_pop > 0 & final_K == 0] <- 1
  
  write.table(final_K, 
              file=sprintf("data/pop/landgrid_capacity_%s_%s_%s.csv", world, time, scaling), 
              sep=",", 
              row.names=FALSE, 
              col.names=FALSE)
  message("population data saved")
  
  # visualization
  if (visualize == TRUE){
    landgrid$K = final_K
    landgrid$plot_K = ifelse(landgrid$K < 0, NA, landgrid$K)
    
    p = ggplot(data = landgrid) +
      geom_sf(aes(fill = plot_K), color = "gray20", linewidth = 0.1) +
      scale_fill_viridis(
        option = "mako",
        direction = -1,
        name = "carrying capacity\n(individuals)",
        na.value = "grey80"
      ) +
      theme_void() +
      theme(
        text = element_text(size = 16),
        legend.position = "right",
        legend.key.size = unit(1.4, "cm"),           
        legend.key.height = unit(1.8, "cm"),         
        legend.key.width  = unit(1.1, "cm"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17, face = "bold", margin = margin(b = 10)),
        legend.spacing.y = unit(0.6, "cm"),
        legend.margin = margin(t = 10, r = 15, b = 10, l = 5),
        plot.margin = margin(10, 20, 10, 10)
      )
  
    output_file = sprintf("output/figures/capacity_%s_%s_%s.pdf", world, time, scaling)
    pdf(output_file, width = 12, height = 8)
    print(p)
    dev.off()
    message(sprintf("saved figure to: %s", output_file))
  }
}

# extract true flux from SLiM data
true_flux <- function(true_node_states, ts, times, cost.mat, neighbor.mat,
                     sample_sets, state_sets)
{
  storage.mode(cost.mat) = "double"
  adjacency_matrix = Matrix::Matrix(neighbor.mat, sparse = TRUE)
  adjacency_matrix = methods::as(adjacency_matrix, "generalMatrix")
  h = .Call(C_treeseq_discrete_mpr_edge_history, ts@treeseq, 
            true_node_states, cost.mat, adjacency_matrix, FALSE)
  H = data.frame(edge_id = rep.int(0:(length(h[[1]]) - 1L), 
                                   times = h[[3]]), state_id = do.call(c, h[[1]]), time = do.call(c, 
                                                                                                  h[[2]]))
  H = structure(H, node.state = true_node_states, path.offset = c(cumsum(h[[3]]) - 
                                                                    h[[3]], nrow(H)) + as.integer(FALSE))
  num_state_sets = as.numeric(max(state_sets))
  num_sample_sets = as.numeric(max(sample_sets))
  num_time_bins = length(times) - 1
  flux = .Call(C_treeseq_discrete_mpr_ancestry_flux, ts@tree, attr(H, 
                                                                   "path.offset"), H$state_id, H$time, as.integer(num_state_sets), 
               state_sets - 1L, as.integer(num_sample_sets), sample_sets - 
                 1L, times)
  return (flux)
}


# process real tree sequence data for afro-eurasia case study (originated from data_full() in GAIA paper)
source("code/generation/proj.R")
process_genetic_data_afroeurasia = function(chromosome, short=TRUE)
{
  if (short)
  {
    tsfile = sprintf("data/genetics/hgdp_tgp_sgdp_high_cov_ancients_chr%d_p.dated.trees", chromosome)
  }
  else
  {
    tsfile = sprintf("data/genetics/hgdp_tgp_sgdp_high_cov_ancients_chr%d_q.dated.trees", chromosome)
  }
  
  ts = treeseq_load(tsfile)
  samples = treeseq_individuals(ts)
  num_samples = nrow(samples)
  nodes = treeseq_nodes(ts)
  
  # keep georeferenced samples only
  {
    index = logical(num_samples)
    for (i in 1:num_samples)
    {
      s = samples[i, ]
      if (length(s[[2]]))
        index[i] = TRUE
    }
    filtered_samples = samples[index,]
    
    filtered_node_ids = nodes[
      nodes$individual_id %in% unlist(filtered_samples[,1]), "node_id"]
    
    ts_simpl = treeseq_simplify(ts, filtered_node_ids)
    
    samples = treeseq_individuals(ts_simpl)
    nodes = treeseq_nodes(ts_simpl)
    
    sample_ids = unlist(samples[,1])
    
    # lons,lats of georeferenced samples
    xy = t(apply(samples, 1, "[[", 2))[,2:1]
    
    # first pass, limit to african and eurasian samples
    filter1 = xy[,1] > -30 
    xy = xy[filter1, ]
    sample_ids = sample_ids[filter1]
    # finish by dropping australasian samples
    filter2 = !(xy[,1] > 110 & xy[,2] < 18) 
    xy = xy[filter2, ]
    sample_ids = sample_ids[filter2]
    sample_idx = match(sample_ids, nodes$individual_id)
    
    node_ids = nodes[
      c(rbind(
        sample_idx
        , sample_idx+1L))
      , 1L
    ]
    
    ts_simpl2 = treeseq_simplify(ts_simpl, node_ids)
    
    # drop edges to nodes with insane number of parents or children
    edges = treeseq_edges(ts_simpl2)
    etab = table(edges$parent_id)
    ignore_parent = as.integer(names(which(etab > quantile(etab, 0.98))))
    etab = table(
      edges$child_id[edges$child_id >= treeseq_num_samples(ts_simpl2)])
    ignore_child = as.integer(names(which(etab > quantile(etab, 0.98))))
    
    ts_simpl3 = treeseq_drop_edges(ts_simpl2, ignore_parent, ignore_child)
    
    samples = treeseq_individuals(ts_simpl3)
    nodes = treeseq_nodes(ts_simpl3)
    
    sample_ids = unlist(samples[,1])
    
    xy = t(apply(samples, 1, "[[", 2))[,2:1]
    
    sample_idx = match(sample_ids, nodes$individual_id)
    
    node_ids = nodes[
      c(rbind(
        sample_idx
        , sample_idx+1L))
      , 1L
    ]
  }
  
  coords = st_transform(st_as_sf(
    data.frame(
      node_id=node_ids,
      lon=rep(xy[,1],each=2),
      lat=rep(xy[,2],each=2)
    )
    , coords=c("lon", "lat")
    , crs=WGS84
  ), crs=st_crs(GRS80))
  
  landgrid = st_read("data/geo/landgrid_afro-eurasia.gpkg", quiet=TRUE)
  
  sample_locations_hex = cbind(
    node_id=coords$node_id,
    state_id=st_nearest_feature(
      st_transform(coords, st_crs(landgrid))$geometry,
      st_centroid(landgrid$geom)
    )
  )
  
  sample_coords = cbind(node_id=node_ids,
                        st_coordinates(st_transform(coords, crs=st_crs(EEGRS80))))
  
  list(
    ts=ts_simpl3,
    data=sample_locations_hex,
    sample.coords=sample_coords
  )
}


# extract samples from preprocessed tree sequence data (can be used in both case studies)
data_subsample = function(D)
{
  ts = D$ts
  dat = D$data
  nodes = treeseq_nodes(ts)
  indivs = tapply(nodes$individual_id, nodes$population_id, function(i) {
    indiv = unique(i)
    if (length(indiv) > 1)
      return (sample(indiv, 1L))
    else
      return (indiv)
  })
  if (indivs[1L] == -1L)
    indivs = indivs[-1L]
  nodes_to_keep = nodes$node_id[which(nodes$individual_id %in% indivs)]
  idx = which(dat[,1] %in% nodes_to_keep)
  ts2 = treeseq_simplify(ts, nodes_to_keep, node.map=TRUE)
  node_ids = attr(ts2, 'node.map')[dat[idx, 1]+1]
  dat2 = cbind(node_id=node_ids, state_id=dat[idx, 2])
  coords = cbind(node_id=node_ids, D$sample.coords[idx, 2:3])
  list(
    ts=ts2,
    data=dat2,
    sample.coords=coords
  )
}


# extract estimated location of oldest ancestors for all sampled replications
get_y <- function(world="afro-eurasia", map, timedepth){
  y = c()
  selection_file <- sprintf("output/tables/selected_reps_%s.csv", world)
  selected_reps_df <- read.csv(selection_file)
  selected_rep_ids <- selected_reps_df$x
  for (rep in selected_rep_ids) {
    x = read.csv(sprintf("data/mpr/mpr_%s_%s_%s.csv", 
                         world, map, rep))
    idx = which(x$node_time == timedepth-1)
    y = rbind(y, x[idx,2:3])
  }
  return(y)
}


# process real tree sequence data for asia-americas case study 
process_genetic_data_asiaamericas = function(chromosome, short=TRUE)
{
  if (short)
  {
    tsfile = sprintf("data/genetics/hgdp_tgp_sgdp_high_cov_ancients_chr%d_p.dated.trees", chromosome)
  }
  else
  {
    tsfile = sprintf("data/genetics/hgdp_tgp_sgdp_high_cov_ancients_chr%d_q.dated.trees", chromosome)
  }
  ts = treeseq_load(tsfile)
  
  # parse population metadata for ethnic group information
  pops = treeseq_populations(ts)
  pop_metadata_list = pops[, "metadata"]
  
  pop_df <- bind_rows(lapply(seq_along(pop_metadata_list), function(i) {
    raw_bytes <- pop_metadata_list[[i]]
    name <- "Unknown"
    if (length(raw_bytes) > 0) {
      json_str <- rawToChar(as.raw(raw_bytes))
      parsed_json <- tryCatch(fromJSON(json_str), error = function(e) list())
      if (!is.null(parsed_json$name)){
        name <- parsed_json$name
      } 
    }
    data.frame(population_id = i - 1, name = name, stringsAsFactors = FALSE)
  }))
  
  target_pop_ids <- pop_df %>% filter(name %in% target_names) %>% pull(population_id)
  
  # parse individual metadata to obtain lat-long for filtering
  inds = treeseq_individuals(ts)
  loc_list = inds[, "location"]
  
  ind_df <- bind_rows(lapply(seq_along(loc_list), function(i) {
    loc <- loc_list[[i]]
    lat <- NA; lon <- NA
    if (length(loc) >= 2) {
      lat <- as.numeric(loc[1])
      lon <- as.numeric(loc[2])
    }
    data.frame(individual_id = i - 1, lat = lat, lon = lon, stringsAsFactors = FALSE)
  }))
  
  # filter
  nodes = as.data.frame(treeseq_nodes(ts))
  
  valid_nodes <- nodes %>%
    filter(is_sample == 1) %>%
    filter(population_id %in% target_pop_ids) %>%
    left_join(ind_df, by = "individual_id") %>%
    filter(!is.na(lat) & !is.na(lon))
  
  filtered_node_ids <- valid_nodes$node_id
  
  # prune tree sequence
  ts_simpl = treeseq_simplify(ts, filtered_node_ids)
  
  edges = treeseq_edges(ts_simpl)
  etab = table(edges$parent_id)
  ignore_parent = as.integer(names(which(etab > quantile(etab, 0.98))))   # delete nodes with exceptionally high in-degree or out-degree
  etab = table(edges$child_id[edges$child_id >= treeseq_num_samples(ts_simpl)])
  ignore_child = as.integer(names(which(etab > quantile(etab, 0.98))))
  
  ts_simpl3 = treeseq_drop_edges(ts_simpl, ignore_parent, ignore_child)
  final_inds = treeseq_individuals(ts_simpl3)
  final_loc_list = final_inds[, "location"]
  
  new_ind_df <- bind_rows(lapply(seq_along(final_loc_list), function(i) {
    loc <- final_loc_list[[i]]
    lat <- NA; lon <- NA
    if (length(loc) >= 2) {
      lat <- as.numeric(loc[1])
      lon <- as.numeric(loc[2])
    }
    # individual_id start from 0
    data.frame(individual_id = i - 1, lat = lat, lon = lon, stringsAsFactors = FALSE)
  }))
  
  final_nodes = as.data.frame(treeseq_nodes(ts_simpl3))
  
  final_samples = final_nodes %>% 
    filter(is_sample == 1) %>%
    left_join(new_ind_df, by = "individual_id") %>%
    filter(!is.na(lat) & !is.na(lon))
  
  node_ids = final_samples$node_id
  xy = as.matrix(final_samples[, c("lon", "lat")])
  
  # projection
  coords = st_transform(st_as_sf(
    data.frame(node_id = node_ids, lon = xy[,1], lat = xy[,2]), 
    coords = c("lon", "lat"), 
    crs = 4326, agr = "constant"
  ), crs = target_crs)
  
  landgrid = st_read("data/geo/landgrid_asia-americas.gpkg", quiet=TRUE)
  if (st_crs(landgrid)$proj4string != target_crs) {
    landgrid <- st_transform(landgrid, crs = target_crs)
  }
  
  sample_locations_hex = cbind(
    node_id = coords$node_id,
    state_id = st_nearest_feature(coords, st_centroid(landgrid))
  )
  
  sample_coords = cbind(node_id = node_ids, st_coordinates(coords))
  
  list(
    ts = ts_simpl3,
    data = sample_locations_hex,
    sample.coords = sample_coords
  )
}


# summary of sampled tree sequence
print_summary_stats <- function(D) {
  cat("total number of kept nodes: ", treeseq_num_nodes(D$ts), "\n")
  cat("total number of kept edges: ", treeseq_num_edges(D$ts), "\n")
  cat("number of kept samples (chromosomes): ", nrow(D$sample.coords), "\n")
  cat("number of grid cells that the samples cover: ", length(unique(D$data[, "state_id"])), "\n")
  cat("number of local genealogies (recombination sites): ", length(unique(c(treeseq_edges(D$ts)$left, treeseq_edges(D$ts)$right))), "\n")
}


# display information of contemporary samples
inspect_tree_nodes <- function(D, n_show = 10) {
  pops <- treeseq_populations(D$ts)
  pop_metadata_list <- pops[, "metadata"]
  pop_names <- sapply(pop_metadata_list, function(raw_bytes) {
    if (length(raw_bytes) > 0) {
      parsed <- tryCatch(jsonlite::fromJSON(rawToChar(as.raw(raw_bytes))), error = function(e) list())
      if (!is.null(parsed$name)) return(parsed$name)
    }
    return("Unknown")
  })
  pop_df <- data.frame(population_id = 0:(length(pop_names)-1), pop_name = pop_names, stringsAsFactors = FALSE)
  nodes <- as.data.frame(treeseq_nodes(D$ts))
  states <- as.data.frame(D$data)
  coords <- as.data.frame(D$sample.coords)
  
  detail_df <- nodes %>%
    filter(is_sample == 1) %>%
    left_join(pop_df, by = "population_id") %>%
    left_join(states, by = "node_id") %>%
    left_join(coords, by = "node_id") %>%
    select(node_id, individual_id, pop_name, time, state_id, X, Y)
  
  print(head(detail_df, n_show))
  cat(sprintf("... left %d samples ...\n", max(0, nrow(detail_df) - n_show)))
}


# visualization for asia-americas
plot_spatial_distribution <- function(D, title = "", landgrid_path = "data/geo/landgrid_asia-americas.gpkg") {
  landgrid <- st_read(landgrid_path, quiet = TRUE)
  if (st_crs(landgrid)$proj4string != target_crs) {
    landgrid <- st_transform(landgrid, crs = target_crs)
  }
  
  coords_df <- as.data.frame(D$sample.coords)
  coords_sf <- st_as_sf(coords_df, coords = c("X", "Y"), crs = target_crs, agr = "constant")
  
  p <- ggplot() +
    geom_sf(data = landgrid, fill = "gray95", color = "white", linewidth = 0.2) +
    geom_sf(data = coords_sf, color = "darkred", size = 1.8, alpha = 0.7) +
    theme_minimal() +
    labs(
      title = title,
      x = "", y = ""
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      panel.grid.major = element_line(color = "gray90")
    )
  
  return(p)
}