library(tidyverse)
library(ggplot2)
library(ggpubr) 
library(gaia)
library(igraph)
library(parallel)

all_mode = c("line", "square", "cube")
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


# # iterate through all modes and both models
# for (mode in mode){
#   for (model in model){
# 
# final_output_list <- list()
# for (i in 1:50) {
#   ts_file <- sprintf("data/trees/math-tree-%s-%s-%s.trees", mode, real_case, i)
#   cost_mat_error <- data.matrix(read.csv(sprintf("data/math/costmat_%s_naive.csv", mode)))
#   cost.mat <- data.matrix(read.csv(sprintf("data/math/costmat_%s_%s.csv", mode, model)))
#   if (mode == "square"){
#     MAXDEME = 100
#   }else if (mode == "line"){
#     MAXDEME = 20
#   }else if (mode == "cube"){
#     MAXDEME = 343
#   }
#   
#   ts <- treeseq_load(ts_file)
#   nodes <- treeseq_nodes(ts)
#   edges <- treeseq_edges(ts)
#   #unique(edges$left)   # check recombinations
#   
#   # 提取所有存活的现代样本节点
#   sample_nodes <- nodes[nodes$is_sample == 1, ]
#   
#   
#   populated_demes <- unique(sample_nodes$population_id)
#   
#   # 2. 【核心修改】：随机选取 10 个格子（包含安全机制）
#   set.seed(42) # 设定种子，保证每次选中的 10 个格子是固定的
#   num_demes_to_sample <- min(30,MAXDEME)
#   
#   # 随机无放回抽取 
#   pops_to_sample <- sample(populated_demes, num_demes_to_sample, replace = FALSE)
#   
#   # 过滤数据，仅保留这 10 个被选中格子里的个体
#   sample_nodes_to_sample <- sample_nodes[sample_nodes$population_id %in% pops_to_sample, ]
#   
#   # 3. 安全抽样函数（防止单个格子只有1人时的 sample 陷阱）
#   safe_sample <- function(x) {
#     if (length(x) == 1) return(x)
#     return(sample(x, 1, replace = FALSE))  # sample how many from each deme
#   }
#   
#   set.seed(40) # 保证每次每个格子抽到的人一样
#   # 从这 10 个格子中，各随机抽取 1 个个体
#   indivs <- unlist(unname(tapply(
#     sample_nodes_to_sample$individual_id, 
#     sample_nodes_to_sample$population_id, 
#     safe_sample
#   )))
#   
#   # 提取保留个体的 node_id 并简化树序列
#   to_keep <- sample_nodes_to_sample$node_id[sample_nodes_to_sample$individual_id %in% indivs]
#   ts2 <- treeseq_simplify(ts, to_keep, filter.populations = FALSE)
#   
#   nodes2 <- treeseq_nodes(ts2)
#   sample_nodes2 <- nodes2[nodes2$is_sample == 1, ]
#   sample_locations <- cbind(node_id = sample_nodes2[, 1], state_id = sample_nodes2[, 4] + 1L)
#   
#   nodes_df <- treeseq_nodes(ts2)
#   edges_tree1 <- treeseq_edges(ts2)
#   
#   ### GAIA
#   mpr <- treeseq_discrete_mpr(ts2, sample_locations, cost.mat)  #, use_brlen=TRUE)
#   estimated_node_states = treeseq_discrete_mpr_minimize(mpr)
#   real_node_states = nodes_df$population_id + 1L
#   comparison_df <- data.frame(
#     node_id = nodes_df$node_id,
#     time = nodes_df$time,
#     real_state = real_node_states,
#     est_state = estimated_node_states
#   )
#   
#   # 2. 按时间降序排列 (time 值越大，代表历史上越古老)
#   comparison_df <- comparison_df[order(comparison_df$time, decreasing = TRUE), ]
#   
#   #分段平均误差距离
#   comparison_df$distance <- cost_mat_error[cbind(comparison_df$real_state, comparison_df$est_state)]
#   
#   # 强制划分为时间上等长的 4 段（注意：如果某些文件时间范围极小，建议检查 breaks）
#   comparison_df$time_bin <- cut(comparison_df$time, breaks = 4, include.lowest = TRUE)
#   
#   # 分组计算
#   node_counts <- aggregate(distance ~ time_bin, data = comparison_df, FUN = length)
#   mean_errors <- aggregate(distance ~ time_bin, data = comparison_df, FUN = function(x) mean(x, na.rm = TRUE))
#   
#   # 排序：从老到新 (4个分段)
#   node_counts <- node_counts[order(node_counts$time_bin, decreasing = TRUE), ]
#   mean_errors <- mean_errors[order(mean_errors$time_bin, decreasing = TRUE), ]
#   
#   # --- 3. 将结果“拉平”成一行 ---
#   # 我们需要 8 个指标：T1_count, T2_count, T3_count, T4_count, T1_err, T2_err, T3_err, T4_err
#   row_data <- data.frame(
#     file_id = i,
#     # 样本量列
#     count_bin1_old = node_counts$distance[1],
#     count_bin2     = node_counts$distance[2],
#     count_bin3     = node_counts$distance[3],
#     count_bin4_new = node_counts$distance[4],
#     # 误差列
#     error_bin1_old = mean_errors$distance[1],
#     error_bin2     = mean_errors$distance[2],
#     error_bin3     = mean_errors$distance[3],
#     error_bin4_new = mean_errors$distance[4]
#   )
#   
#   final_output_list[[i]] <- row_data
#   message(paste("Processed file:", i))
# }
# 
# # ==========================================
# # 合并并保存最终结果
# # ==========================================
# summary_table <- do.call(rbind, final_output_list)
# 
# # 写入文件
# write.csv(summary_table, sprintf("output/tables/all_results_%s_%s.csv", mode, model), row.names = FALSE)
# 
#   }
# }



######## test the difference ########
for (mode in all_mode){
  
df_naive <- read.csv(sprintf("output/tables/all_results_%s_naive.csv", mode)) %>% mutate(group = "Naive")
df_friction <- read.csv(sprintf("output/tables/all_results_%s_friction.csv", mode)) %>% mutate(group = "Friction")
df_naive=df_naive[1:50,]
df_friction=df_friction[1:50,]
df_all <- rbind(df_naive, df_friction)

# 转换为长格式：提取 Error 列进行比较
df_long <- df_all %>%
  select(group, error_bin1_old, error_bin2, error_bin3, error_bin4_new) %>%
  pivot_longer(
    cols = starts_with("error"),
    names_to = "Time_Bin",
    values_to = "Error_Distance"
  ) %>%
  # 优化标签，使其在图表上更易读
  mutate(Time_Bin = factor(Time_Bin, 
                           levels = c("error_bin1_old", "error_bin2", "error_bin3", "error_bin4_new"),
                           labels = c("T1 (Oldest)", "T2", "T3", "T4 (Recent)")))

# stat_results <- df_long %>%
#   group_by(Time_Bin) %>%
#   summarise(
#     p_value = wilcox.test(Error_Distance ~ group)$p.value,
#     mean_naive = mean(Error_Distance[group == "Naive"]),
#     mean_friction = mean(Error_Distance[group == "Friction"])
#   )
stat_results <- df_long %>%
  group_by(Time_Bin) %>%
  summarise(
    # 样本量
    n_Naive     = sum(group == "Naive"),
    n_Friction  = sum(group == "Friction"),
    n_Total     = n(),
    
    # 统计量
    p_value     = wilcox.test(Error_Distance ~ group)$p.value,
    
    # 均值
    mean_Naive    = mean(Error_Distance[group == "Naive"]),
    mean_Friction = mean(Error_Distance[group == "Friction"])
  ) %>%
  ungroup()
print(stat_results)

# p <- ggplot(df_long, aes(x = Time_Bin, y = Error_Distance, fill = group)) +
#   # 绘制小提琴图展示密度
#   geom_violin(alpha = 0.5, position = position_dodge(width = 0.8), color = NA) +
#   # 叠加箱线图展示中位数和分位数
#   geom_boxplot(width = 0.2, position = position_dodge(width = 0.8), outlier.size = 0.5) +
#   # 使用科学论文常用的配色
#   scale_fill_manual(values = c("Naive" = "#66c2a5", "Friction" = "#fc8d62")) +
#   # 添加显著性标注 (自动比较 Naive vs Friction)
#   stat_compare_means(aes(group = group), label = "p.signif", method = "wilcox.test", hide.ns = FALSE, size = 5.5) +
#   labs(
#     #title = "Ancestry Inference Error: Naive vs Friction Landscape",
#     #subtitle = "Comparing Mean Error Distance across 100 Replications",
#     x = "temporal segments",
#     y = "mean error distance (steps)",
#     fill = "model"
#   ) +
#   theme_minimal() +
#   # theme(
#   #   legend.position = "top",
#   #   panel.grid.minor = element_blank()
#   # )
#   theme(
#     legend.position = "top",
#     legend.title = element_text(size = 18, face = "bold"),
#     legend.text = element_text(size = 16),
#     plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
#     plot.subtitle = element_text(size = 16, hjust = 0.5),
#     axis.title.x = element_text(size = 16),
#     axis.title.y = element_text(size = 16),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank()
#   )
# 
# ggsave(sprintf("output/figures/error_comparison_%s.png", mode), 
#        plot = p, 
#        width = 10, 
#        height = 6.5, 
#        dpi = 600)
p <- ggplot(df_long, aes(x = Time_Bin, y = Error_Distance, fill = group)) +
  geom_violin(alpha = 0.5, position = position_dodge(width = 0.8), color = NA) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.8), outlier.size = 0.5) +
  scale_fill_manual(values = c("Naive" = "#66c2a5", "Friction" = "#fc8d62"),
                    labels = c("Naive" = "naive", "Friction" = "friction")) +
  
  stat_compare_means(aes(group = group), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     hide.ns = FALSE, 
                     size = 9) +   # 显著性标记也加大
  
  labs(
    x = "temporal segments",
    y = "mean error distance (steps)",
    fill = "model"
  ) +
  theme_minimal() +
  theme(
    # ==================== 文字整体加大 ====================
    text = element_text(size = 20),                    # 全局基础字体
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 18, hjust = 0.5),
    
    axis.title.x = element_text(size = 24, margin = margin(t = 12)),
    axis.title.y = element_text(size = 24, margin = margin(r = 12)),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    
    legend.position = "top",
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 24),
    
    # ==================== 缩小图内空白 ====================
    panel.spacing = unit(0.5, "lines"),           # 面板间距
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10),  # 整体外边距
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# 保存图像（建议使用较大尺寸）
ggsave(sprintf("output/figures/error_comparison_%s.png", mode), 
       plot = p, width = 11, height = 7, dpi = 300)

print(p)

}





