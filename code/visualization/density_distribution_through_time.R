##### parameters #####
num_subsets = 100L
map = c("naive", "friction")
time_from = c(0, 10000, 81719)   # a pair of values in time_from and time_to corresponds to a time period
time_to = c(1, 20000, 81721)


landgrid = st_read("data/geo/landgrid_asia-americas.gpkg", quiet = TRUE)

for (mymap in map){
  for (mytime in seq_along(time_from)){

    y = c()
    for (i in 1:100) {
      filename = sprintf("data/mpr/mpr_asia-americas_%s_%d.csv", mymap, i) 
      if (file.exists(filename)) {
        x = read.csv(filename)
        state = x[time_from[mytime] <= x$node_time & x$node_time < time_to[mytime], 2]    # max time: 81720   
        y = c(y, state)
      }
    }
    
    # frequency of landing in each grid cell
    state_counts = as.data.frame(table(state_id = y))
    state_counts$state_id = as.integer(as.character(state_counts$state_id))
    
    # noise filtering (optional)
    max_freq <- max(state_counts$Freq, na.rm = TRUE)    
    # various ways to set threshold
    # threshold <- max_freq * 0.001                      
    # mean_freq <- mean(state_counts$Freq, na.rm = TRUE)
    # threshold <- mean_freq * 0.05
    # median_freq <- median(state_counts$Freq, na.rm = TRUE)
    # threshold <- median_freq * 0.1
    threshold=0
    state_counts$Freq[state_counts$Freq < threshold] <- NA 
    
    # plot
    plot_data_poly = landgrid %>%
      left_join(state_counts, by = c("cell_id" = "state_id"))
    
    p = ggplot() +
      geom_sf(data = plot_data_poly, aes(fill = Freq), color = "grey85", linewidth = 0.2) +
      
      scale_fill_gradient(
        name = "estimated\nancestors",
        low = "#FFF59D",
        high = "#B71C1C",
        na.value = "grey95",
        
        breaks = c(min(plot_data_poly$Freq, na.rm = TRUE),
                   max(plot_data_poly$Freq, na.rm = TRUE)),
        labels = c("sparse", "dense")
      ) +
      
      theme_minimal() +
      theme(
        panel.grid = element_blank(), 
        axis.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right", #"none",  
        legend.key.size = unit(1.2, "cm"),      
        legend.key.height = unit(1.5, "cm"),   
        legend.key.width  = unit(1.0, "cm"),      
        legend.text = element_text(size = 14),   
        legend.title = element_text(size = 15, face = "bold", margin = margin(b = 8)),
        legend.spacing.y = unit(0.8, "cm"),   
        legend.margin = margin(t = 10, r = 20, b = 10, l = 10) 
      )
    
    output_file = sprintf("output/figures/ancestry_distribution_asia-americas_%s_from_%s_to_%s.pdf", mymap, time_from[mytime], time_to[mytime])
    pdf(output_file, width = 12, height = 8)
    print(p)
    dev.off()
    message(sprintf("saved figure to: %s", output_file))

  }
}