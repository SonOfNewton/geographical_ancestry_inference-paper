model = "naive"  # "friction" "naive"


ts_file <- "data/toy/tree-line-test.trees"

MAXDEME = 20
if (model == "naive"){
  cost.mat = outer(1:MAXDEME, 1:MAXDEME, function(i, j) abs(i - j))
} else if (model == "friction"){
  cost.mat = data.matrix(read.csv("data/toy/costmat_line_friction.csv"))
}

ts <- treeseq_load(ts_file)
nodes <- treeseq_nodes(ts)
edges <- treeseq_edges(ts)
#unique(edges$left)   # check recombinations

# 提取所有存活的现代样本节点
sample_nodes <- nodes[nodes$is_sample == 1, ]


populated_demes <- unique(sample_nodes$population_id)

# 2. 【核心修改】：随机选取 10 个格子（包含安全机制）
set.seed(42) # 设定种子，保证每次选中的 10 个格子是固定的
num_demes_to_sample <- min(50,MAXDEME)

# 随机无放回抽取 10 个格子
pops_to_sample <- sample(populated_demes, num_demes_to_sample, replace = FALSE)


# 过滤数据，仅保留这 10 个被选中格子里的个体
sample_nodes_to_sample <- sample_nodes[sample_nodes$population_id %in% pops_to_sample, ]

# 3. 安全抽样函数（防止单个格子只有1人时的 sample 陷阱）
safe_sample <- function(x) {
  if (length(x) == 1) return(x)
  return(sample(x, 1, replace = FALSE))  # sample how many from each deme
}

set.seed(40) # 保证每次每个格子抽到的人一样
# 从这 10 个格子中，各随机抽取 1 个个体
indivs <- unlist(unname(tapply(
  sample_nodes_to_sample$individual_id, 
  sample_nodes_to_sample$population_id, 
  safe_sample
)))

# 提取保留个体的 node_id 并简化树序列
to_keep <- sample_nodes_to_sample$node_id[sample_nodes_to_sample$individual_id %in% indivs]
ts2 <- treeseq_simplify(ts, to_keep, filter.populations = FALSE)

nodes2 <- treeseq_nodes(ts2)
sample_nodes2 <- nodes2[nodes2$is_sample == 1, ]
sample_locations <- cbind(node_id = sample_nodes2[, 1], state_id = sample_nodes2[, 4] + 1L)

nodes_df <- treeseq_nodes(ts2)
edges_tree1 <- treeseq_edges(ts2)

### GAIA
mpr <- treeseq_discrete_mpr(ts2, sample_locations, cost.mat)  #, use_brlen=TRUE)
estimated_node_states = treeseq_discrete_mpr_minimize(mpr)
real_node_states = nodes_df$population_id + 1L
comparison_df <- data.frame(
  node_id = nodes_df$node_id,
  time = nodes_df$time,
  real_state = real_node_states,
  est_state = estimated_node_states
)

# 2. 按时间降序排列 (time 值越大，代表历史上越古老)
comparison_df <- comparison_df[order(comparison_df$time, decreasing = TRUE), ]

# 3. 提取最古老的 100 个祖先节点
oldest_20 <- head(comparison_df, 100)

# 4. 从 cost.mat 中批量提取这 100 个节点的真实与估计位置之间的距离
# 核心技巧：使用 cbind() 将两列合并为坐标矩阵，这样可以直接从 cost.mat 中精准抓取 100 个对应距离，速度极快
oldest_20$distance <- cost.mat[cbind(oldest_20$real_state, oldest_20$est_state)]

# 5. 计算平均误差距离
mean(oldest_20$distance, na.rm = TRUE)

#分段平均误差距离
comparison_df$distance <- cost.mat[cbind(comparison_df$real_state, comparison_df$est_state)]
comparison_df$time_bin <- cut(comparison_df$time, breaks = 4, include.lowest = TRUE)
# 3. 按时间段计算样本量 (节点数)
node_counts <- aggregate(distance ~ time_bin, data = comparison_df, FUN = length)
colnames(node_counts)[2] <- "Node_Count"
# 4. 按时间段计算平均误差距离
mean_errors <- aggregate(distance ~ time_bin, data = comparison_df, FUN = function(x) mean(x, na.rm = TRUE))
colnames(mean_errors)[2] <- "Mean_Error"
result_df <- merge(node_counts, mean_errors, by = "time_bin")
result_df <- result_df[order(result_df$time_bin, decreasing = TRUE), ]
rownames(result_df) <- NULL
#print(result_df)


target_pos <- 0  

local_edges <- edges_tree1[edges_tree1$left <= target_pos & edges_tree1$right > target_pos, ]

# ==========================================
# 2. 构建简易网络图
# ==========================================
g_simple <- graph_from_data_frame(
  data.frame(from = local_edges$parent, to = local_edges$child)
)

# ==========================================
# 3. 提取节点属性
# ==========================================
node_indices <- match(as.numeric(V(g_simple)$name), nodes_df$node_id)
times <- nodes_df$time[node_indices]
real_demes <- nodes_df$population_id[node_indices] + 1L
est_demes <- estimated_node_states[node_indices]

# ==========================================
# 4. 【核心升级 1】自定义美化布局 (Y轴时间，X轴按Deme排序)
# ==========================================
lay <- layout_as_tree(g_simple)
lay[, 2] <- times  # Y轴强制设定为真实世代时间

# 获取所有的叶子节点（出度为 0 的节点，即当代样本）
out_degrees <- igraph::degree(g_simple, mode = "out")
leaf_idx <- which(out_degrees == 0)

# 按照它们所在的真实地理位置 (Deme) 从小到大排序
sorted_leaf_idx <- leaf_idx[order(real_demes[leaf_idx])]

# 给这些排好序的叶子节点分配从左到右的等距 X 坐标
lay[sorted_leaf_idx, 1] <- seq_along(sorted_leaf_idx)

# 对于内部节点（祖先），按时间从新到老排序
internal_idx <- which(out_degrees > 0)
internal_idx <- internal_idx[order(times[internal_idx])]

# 自下而上：把每个祖先节点的 X 坐标，放在它所有子节点 X 坐标的“正中间”，让树枝最优雅
for (i in internal_idx) {
  children <- neighbors(g_simple, v = i, mode = "out")
  # 提取子节点在当前图中的数字索引，用于获取它们已经计算好的 X 坐标
  children_idx <- match(names(children), V(g_simple)$name)
  lay[i, 1] <- mean(lay[children_idx, 1])
}

# ==========================================
# 5. 【核心升级 2】柔和、易读的颜色与形状映射
# ==========================================
# 使用非常经典的科研绘图柔和渐变色 (深蓝 -> 浅蓝 -> 柔黄 -> 砖红)
soft_colors <- c("#4575b4", "#e0f3f8", "#fee090", "#d73027")
pal <- colorRampPalette(soft_colors)(MAXDEME)

V(g_simple)$real_color <- pal[real_demes]
V(g_simple)$est_color <- pal[est_demes]  

is_sample_node <- nodes_df$is_sample[node_indices] == 1
V(g_simple)$frame_color <- ifelse(is_sample_node, "black", "gray80")

# ==========================================
# 6. 并排画图与图例
# ==========================================
#png(file="output/figures/line_tree-demo-20.png", width=10, height=5, units="in", res=300)
# 2. 这里的绘图逻辑和你之前的一模一样
par(mfrow = c(1, 2), mar = c(3, 4, 4, 1))  

# --- 左图 ---
plot(g_simple, layout = lay, vertex.color = V(g_simple)$real_color,
     vertex.frame.color = V(g_simple)$frame_color, vertex.size = 10,
     vertex.label = NA, edge.arrow.size = 0.2, edge.color = "gray60",
     #main = paste("Local Tree at position:", target_pos, "\n(Real Node States)")
     )

# --- 右图 ---
par(mar = c(3, 1, 4, 4))
plot(g_simple, layout = lay, vertex.color = V(g_simple)$est_color,
     vertex.frame.color = V(g_simple)$frame_color, vertex.size = 10,
     vertex.label = NA, edge.arrow.size = 0.2, edge.color = "gray60",
     #main = paste("Local Tree at position:", target_pos, "\n(Estimated Node States)")
     )

# --- 图例 ---
legend_demes <- unique(round(seq(1, MAXDEME, length.out = 5)))
legend("topright", title = "Deme", legend = legend_demes, 
       fill = pal[legend_demes], bty = "n", cex = 1.2, inset = 0.02)

#dev.off()












