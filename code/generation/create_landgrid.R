if (my_world == "asia-americas"){
  ######## map of the americas and eurasia ########
  sf_use_s2(FALSE)   # disable the S2 spherical geometry engine to avoid topology errors when processing polygons crossing datelines
  
  # build map
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  shift <- 150   # shift the center of the map from 0° to 150°W 
  world_pacific <- world %>%
    st_break_antimeridian(lon_0 = shift)
  
  # select continents, delete islands
  target_continents <- c("North America", "South America", "Europe", "Asia")
  islands_to_remove <- c("Greenland", "United Kingdom", "Ireland", "Japan", "Iceland", 
                         "Madagascar", "Taiwan", "Philippines", "Indonesia", "Sri Lanka")
  
  land_mass <- world_pacific %>%
    filter(continent %in% target_continents) %>%
    filter(!name %in% islands_to_remove)
  
  proj_pacific <- "+proj=laea +lat_0=30 +lon_0=-150 +datum=WGS84 +units=m +no_defs"   # use equal area projection: Lambert Azimuthal Equal Area (LAEA)
  land_proj <- st_transform(land_mass, proj_pacific)
  land_union <- st_make_valid(st_union(land_proj))
  
  # build grid
  a_meters <- 300000   # hexagon edge length 300km 
  cell_size <- a_meters * sqrt(3)
  hex_grid <- st_make_grid(land_union, cellsize = cell_size, square = FALSE)
  hex_sf <- st_as_sf(data.frame(id = 1:length(hex_grid)), geometry = hex_grid)
  
  # crop coastline
  hex_clipped <- st_intersection(st_make_valid(hex_sf), land_union) %>%
    st_collection_extract(type = "POLYGON") %>%
    st_cast("MULTIPOLYGON")
  
  # assign ID
  nb <- poly2nb(hex_clipped)
  comp <- n.comp.nb(nb)
  comp_sizes <- table(comp$comp.id)
  
  # threshold for filter
  main_comp_ids <- as.integer(names(comp_sizes[comp_sizes > 20]))
  hex_mainland <- hex_clipped[comp$comp.id %in% main_comp_ids, ] %>%
    mutate(cell_id = 1:n())
  
  # adjacency matrix
  nb_final <- poly2nb(hex_mainland)
  adj_matrix <- nb2mat(nb_final, style = "B", zero.policy = TRUE)
  rownames(adj_matrix) <- hex_mainland$cell_id
  colnames(adj_matrix) <- hex_mainland$cell_id
  write.csv(adj_matrix, "data/geo/landgrid_adjmat_naive_asia-americas.csv", row.names = TRUE)
  
  # cost matrix
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
  cost_matrix <- distances(g, weights = NA)
  write.csv(cost_matrix, "data/geo/landgrid_costmat_naive_asia-americas.csv", row.names = FALSE)
  
  # gpkg file
  st_write(hex_mainland, "data/geo/landgrid_asia-americas.gpkg", append = FALSE)
  
  # visualization
  p <- ggplot() +
    geom_sf(data = hex_mainland, fill = "gray90", color = "darkred", linewidth = 0.1) +
    theme_minimal() +
    labs(title = "Pacific-centric land grid of Eurasia and the Americas")
  
  output_file = "output/figures/landgrid_3d_asia-americas.pdf"
  pdf(output_file, width = 8, height = 8)
  print(p)
  dev.off()
  message(sprintf("saved figure to: %s", output_file))
  
} else{   # add new landgrid here
  message(sprintf("no such world setting!"))
}

