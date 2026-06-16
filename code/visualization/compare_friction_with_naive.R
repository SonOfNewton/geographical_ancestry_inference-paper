# settings
all_mode = c("line", "square", "cube", "annulus", "annulus2", "annulus3")
real_case = "friction"
all_model = c("friction", "naive")

for (mode in all_mode){
  for (model in all_model){
    cat(sprintf("\n=== Running: mode = %s | model = %s ===\n", mode, model))
    
    if (!exists("usecores")) usecores <- 1
    cl <- makeCluster(usecores)
    clusterExport(cl, varlist = c("mode", "model", "real_case", "process_one_file"), envir = environment())
    clusterEvalQ(cl, {
      library(gaia)
      library(tidyverse)
    })
    
    final_output_list <- parLapply(cl, 1:50, function(i) {
      process_one_file(i, mode, model, real_case)
    })
    
    stopCluster(cl)
    
    summary_table <- do.call(rbind, final_output_list)
    output_path <- sprintf("output/tables/all_results_%s_%s.csv", mode, model)
    write.csv(summary_table, output_path, row.names = FALSE)
    
    cat(sprintf("Finished: %s | Saved to: %s\n", paste(mode, model, sep="_"), output_path))
  }
}


######## test the difference ########
for (mode in all_mode){
  
df_naive <- read.csv(sprintf("output/tables/all_results_%s_naive.csv", mode)) %>% mutate(group = "Naive")
df_friction <- read.csv(sprintf("output/tables/all_results_%s_friction.csv", mode)) %>% mutate(group = "Friction")
df_naive=df_naive[1:50,]
df_friction=df_friction[1:50,]
df_all <- rbind(df_naive, df_friction)

df_long <- df_all %>%
  select(group, error_bin1_old, error_bin2, error_bin3, error_bin4_new) %>%
  pivot_longer(
    cols = starts_with("error"),
    names_to = "Time_Bin",
    values_to = "Error_Distance"
  ) %>%
  mutate(Time_Bin = factor(Time_Bin, 
                           levels = c("error_bin1_old", "error_bin2", "error_bin3", "error_bin4_new"),
                           labels = c("T1 (Oldest)", "T2", "T3", "T4 (Recent)")))

stat_results <- df_long %>%
  group_by(Time_Bin) %>%
  summarise(
    n_Naive     = sum(group == "Naive"),
    n_Friction  = sum(group == "Friction"),
    n_Total     = n(),
    p_value     = wilcox.test(Error_Distance ~ group)$p.value,
    mean_Naive    = mean(Error_Distance[group == "Naive"]),
    mean_Friction = mean(Error_Distance[group == "Friction"])
  ) %>%
  ungroup()
print(stat_results)

# plot
p <- ggplot(df_long, aes(x = Time_Bin, y = Error_Distance, fill = group)) +
  geom_violin(alpha = 0.5, position = position_dodge(width = 0.8), color = NA) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.8), outlier.size = 0.5) +
  scale_fill_manual(values = c("Naive" = "#66c2a5", "Friction" = "#fc8d62"),
                    labels = c("Naive" = "naive", "Friction" = "friction")) +
  
  stat_compare_means(aes(group = group), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     hide.ns = FALSE, 
                     size = 9) +   
  labs(
    x = "temporal segments",
    y = "mean error distance (steps)",
    fill = "model"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),                  
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 18, hjust = 0.5),
    
    axis.title.x = element_text(size = 24, margin = margin(t = 12)),
    axis.title.y = element_text(size = 24, margin = margin(r = 12)),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    
    legend.position = "top",
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 24),
    
    panel.spacing = unit(0.5, "lines"),          
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10), 
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(sprintf("output/figures/error_comparison_%s.png", mode), 
       plot = p, width = 11, height = 7, dpi = 300)

print(p)

}

message(paste0('Error comparison figure (different suffix) saved in:\n',"output/figures/error_comparison_suffix.png",'\n\n'))


