# 安装缺失的包（如果需要）：
# install.packages(c("sf", "rnaturalearth", "rnaturalearthdata", "dplyr", "spdep", "ggplot2"))

library(sf)
library(rnaturalearth)
library(dplyr)
library(spdep)
library(ggplot2)
library(igraph)

############################ map of the americas ############################

# ==========================================
# 1. 获取地图并进行初步过滤
# ==========================================
# 获取中等精度的全球地图
world <- ne_countries(scale = "medium", returnclass = "sf")

# 筛选南北美洲，并明确剔除格陵兰岛
americas <- world %>%
  filter(continent %in% c("North America", "South America")) %>%
  filter(name != "Greenland")

# 截断东半球部分 (0 ~ 180E)
# 创建一个仅包含西半球的边界框 (xmin = -180, xmax = 0)
# suppressWarnings 是为了避免 sf 提示 attribute spatially constant 的常规警告
bbox_west <- st_bbox(c(xmin = -180, xmax = 0, ymin = -90, ymax = 90), crs = st_crs(4326))
americas_clipped <- st_crop(americas, bbox_west)

# ==========================================
# 2. 投影转换与轮廓合并
# ==========================================
# 使用以美洲为中心的 Lambert Azimuthal Equal Area (LAEA) 投影
# 确保在计算 300km 网格时，物理距离不受纬度变形影响
proj_laea <- "+proj=laea +lon_0=-80 +lat_0=0 +datum=WGS84 +units=m +no_defs"
americas_proj <- st_transform(americas_clipped, proj_laea)

# 将所有国家的多边形合并成一个完整的大陆轮廓
americas_union <- st_union(americas_proj)

# ==========================================
# 3. 生成 300km 边长的正六边形网格
# ==========================================
# 对于正六边形，若边长为 a，那么平行的两条边之间的距离（即 st_make_grid 的 cellsize）为 a * sqrt(3)
a_meters <- 300000 
cell_size <- a_meters * sqrt(3)

# 生成覆盖整个范围的六边形网格 (square = FALSE 表示六边形)
hex_grid <- st_make_grid(americas_union, cellsize = cell_size, square = FALSE)

# 转换为 sf 对象以便后续操作
hex_sf <- st_as_sf(data.frame(id = 1:length(hex_grid)), geometry = hex_grid)

# ==========================================
# 4. 海岸线截断 (Clipping)
# ==========================================
# 只保留与大陆相交的网格，并沿着海岸线精准截断
hex_clipped <- st_intersection(hex_sf, americas_union) %>%
  st_cast("MULTIPOLYGON") # 统一几何类型，避免出现意外的点或线段

# ==========================================
# 5. 确保 100% 连通并过滤孤立岛屿
# ==========================================
# 计算当前所有网格之间的邻居关系 (基于边界接触)
nb <- poly2nb(hex_clipped)

# 查找空间图中的独立连通分支
comp <- n.comp.nb(nb)

# 找到包含网格数量最多的那个连通分支（即美洲主大陆）
mainland_comp_id <- as.integer(names(which.max(table(comp$comp.id))))

# 过滤掉所有不连通的网格（离岛）
hex_mainland <- hex_clipped[comp$comp.id == mainland_comp_id, ]

# ==========================================
# 6. 坐标系转换与重新编号
# ==========================================
# 转换回 WGS84 经纬度坐标系。这一步将自然产生你要求的“网格随纬度拉伸”的二维视觉效果
hex_final <- st_transform(hex_mainland, 4326) %>%
  mutate(cell_id = 1:n()) %>%
  select(cell_id, geometry) # 整理列，保持清爽

# ==========================================
# 7. 可视化检查
# ==========================================
p <- ggplot() +
  geom_sf(data = hex_final, fill = "transparent", color = "black", size = 0.3) +
  theme_minimal() +
  labs(title = "Americas Connected Land Grid",
       subtitle = "300km Edge, Coastline Clipped, Mainland Only")
print(p)

# ==========================================
# 8. 保存 GeoPackage (GPKG)
# ==========================================
st_write(hex_final, "data/geo/landgrid_americas_300km.gpkg", append = FALSE)


# ==========================================
# 9. 提取并保存 Adjacency Matrix
# ==========================================
# 基于最终的网格重新计算完全连通的邻接关系
nb_final <- poly2nb(hex_final)

# 转换为矩阵格式 (style="B" 表示二进制的 0/1 矩阵，zero.policy=TRUE 防止无邻居报错虽然这里不会发生)
adj_matrix <- nb2mat(nb_final, style = "B", zero.policy = TRUE)

# 设置行列名为真实的 cell_id
rownames(adj_matrix) <- hex_final$cell_id
colnames(adj_matrix) <- hex_final$cell_id

# 保存为 CSV
write.csv(adj_matrix, "data/geo/landgrid-adjmat-americas.csv", row.names = TRUE)



############################ map of the americas and eurasia ############################
# 关闭 S2 球面几何引擎，避免在处理跨越日期线的多边形时出现拓扑报错
sf_use_s2(FALSE)

# ==========================================
# 1. 获取地图并进行“太平洋中心”重构
# ==========================================
# 获取全球地图
world <- ne_countries(scale = "medium", returnclass = "sf")

# 核心技巧：将地图向右平移，使得地图中心从 0° 移至 150°W (太平洋)
# 这样大西洋将成为地图的边缘，而白令海峡将处于地图中间位置
shift <- 150
world_pacific <- world %>%
  st_break_antimeridian(lon_0 = shift)

# 筛选大洲，剔除格陵兰岛等碎岛
target_continents <- c("North America", "South America", "Europe", "Asia")
islands_to_remove <- c("Greenland", "United Kingdom", "Ireland", "Japan", "Iceland", 
                       "Madagascar", "Taiwan", "Philippines", "Indonesia", "Sri Lanka")

land_mass <- world_pacific %>%
  filter(continent %in% target_continents) %>%
  filter(!name %in% islands_to_remove)

# ==========================================
# 2. 投影转换：使用等面积投影
# ==========================================
# 为了计算 300km 网格，我们需要一个平面投影。
# 既然要看白令海峡，使用 Lambert Azimuthal Equal Area (LAEA)，中心设在太平洋
proj_pacific <- "+proj=laea +lat_0=30 +lon_0=-150 +datum=WGS84 +units=m +no_defs"
land_proj <- st_transform(land_mass, proj_pacific)

# 合并轮廓并修复拓扑
land_union <- st_make_valid(st_union(land_proj))

# ==========================================
# 3. 生成 300km 六边形网格
# ==========================================
a_meters <- 300000 
cell_size <- a_meters * sqrt(3)

hex_grid <- st_make_grid(land_union, cellsize = cell_size, square = FALSE)
hex_sf <- st_as_sf(data.frame(id = 1:length(hex_grid)), geometry = hex_grid)

# ==========================================
# 4. 海岸线截断 (提取 Polygon 避免 Point/Line 报错)
# ==========================================
cat("正在沿海岸线截断网格（太平洋中心模式）...\n")
hex_clipped <- st_intersection(st_make_valid(hex_sf), land_union) %>%
  st_collection_extract(type = "POLYGON") %>%
  st_cast("MULTIPOLYGON")

# ==========================================
# 5. 过滤岛屿并分配 ID
# ==========================================
nb <- poly2nb(hex_clipped)
comp <- n.comp.nb(nb)
comp_sizes <- table(comp$comp.id)

# 阈值调低：确保通过白令陆桥连接后的超级大陆被保留
main_comp_ids <- as.integer(names(comp_sizes[comp_sizes > 20]))
hex_mainland <- hex_clipped[comp$comp.id %in% main_comp_ids, ] %>%
  mutate(cell_id = 1:n())

# ==========================================
# 6. 提取邻接矩阵 (关键点：此时白令海峡已连通)
# ==========================================
nb_final <- poly2nb(hex_mainland)
adj_matrix <- nb2mat(nb_final, style = "B", zero.policy = TRUE)
rownames(adj_matrix) <- hex_mainland$cell_id
colnames(adj_matrix) <- hex_mainland$cell_id

# 保存邻接矩阵
write.csv(adj_matrix, "data/geo/landgrid-adjmat-global.csv", row.names = TRUE)

# ==========================================
# 7. 可视化与保存
# ==========================================
# 为了可视化，我们转换到一个专门展示太平洋中心的 Robinson 投影
p <- ggplot() +
  geom_sf(data = hex_mainland, fill = "gray90", color = "darkred", linewidth = 0.1) +
  theme_minimal() +
  labs(title = "Pacific-Centric Land Grid (300km)",
       subtitle = "Connected via Beringia | Africa-Eurasia-Americas")

print(p)

# 保存 GPKG
st_write(hex_mainland, "data/geo/landgrid_global_300km.gpkg", append = FALSE)


##### from adjmat to distmat
# 1. 从邻接矩阵构建图网络对象
# mode = "undirected" 因为网格连通性是双向的
# diag = FALSE 忽略自己到自己的连接（对角线为0）
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)

# 2. 计算最短路径距离矩阵 (Cost Matrix)
# weights = NA 非常关键！它强制算法忽略任何潜在的权重，将每一次跨越格子严格计为 1 步 (整数步)
cost_matrix <- distances(g, weights = NA)

# 3. 整理输出格式，以匹配你下游的 MPR 代码
# 确保列名和行名与网格 ID 一致（或者根据你的需求清空它们）
# dimnames(cost_matrix) <- NULL # 如果你的后续代码不需要 dimnames，可以取消这行的注释

# 保存为你需要的 costmat 文件
write.csv(cost_matrix, "data/geo/landgrid-distmat-steps-global.csv", row.names = FALSE)


############################ map of africa and eurasia ############################
world <- ne_countries(scale = "medium", returnclass = "sf")

# 安全裁切框：避开白令海峡 180度经线
bbox_safe <- st_bbox(c(xmin = -20, xmax = 170, ymin = -40, ymax = 90), crs = st_crs(4326))
world_safe <- suppressWarnings(st_crop(world, bbox_safe))

# ==========================================
# 修改点 1：定义纯净大陆提取函数，从根源剔除所有离岛
# ==========================================
# 该函数会将拼合后的大陆打散成独立的多边形，并强制只保留面积最大的那一块。
# 这会完美剔除印尼群岛、日本、英国、马达加斯加以及北极的所有破碎岛屿，但保留马六甲半岛等陆上延伸。
get_mainland_polygon <- function(sf_obj) {
  sf_obj %>%
    st_union() %>%
    st_make_valid() %>% 
    st_cast("POLYGON") %>%
    st_as_sf() %>%
    mutate(area = as.numeric(st_area(x))) %>%
    arrange(desc(area)) %>%
    slice(1) %>%
    st_geometry()
}

# 精准分离出两块纯净的“主大陆”轮廓
africa_outline <- world_safe %>% filter(continent == "Africa") %>% get_mainland_polygon()
eurasia_outline <- world_safe %>% filter(continent %in% c("Europe", "Asia")) %>% get_mainland_polygon()

# ==========================================
# 2. 投影转换与生成通用底层网格 (Universal Grid)
# ==========================================
proj_laea <- "+proj=laea +lon_0=50 +lat_0=30 +datum=WGS84 +units=m +no_defs"

africa_laea <- st_transform(africa_outline, proj_laea)
eurasia_laea <- st_transform(eurasia_outline, proj_laea)

# 合并以获取范围
combined_laea <- st_union(africa_laea, eurasia_laea)

# 生成一个覆盖全体的“统一制式”网格
a_meters <- 200000 
cell_size <- a_meters * sqrt(3)
universal_grid <- st_make_grid(combined_laea, cellsize = cell_size, square = FALSE)
universal_sf <- st_as_sf(data.frame(id = 1:length(universal_grid)), geometry = universal_grid)

# ==========================================
# 3. 独立裁切 
# ==========================================
africa_hex <- suppressWarnings(st_intersection(universal_sf, africa_laea)) %>% 
  st_cast("MULTIPOLYGON")
eurasia_hex <- suppressWarnings(st_intersection(universal_sf, eurasia_laea)) %>% 
  st_cast("MULTIPOLYGON")

# ==========================================
# 4. 独立过滤（作为二次保险，清洗切边可能产生的小碎块）
# ==========================================
filter_mainland <- function(hex_sf) {
  nb <- poly2nb(hex_sf)
  comp <- n.comp.nb(nb)
  mainland_id <- as.integer(names(which.max(table(comp$comp.id))))
  return(hex_sf[comp$comp.id == mainland_id, ])
}

africa_mainland <- filter_mainland(africa_hex)
eurasia_mainland <- filter_mainland(eurasia_hex)

# ==========================================
# 5. 合并组装与最终清洗
# ==========================================
afroeurasia_combined <- bind_rows(
  africa_mainland %>% mutate(continent = "Africa"),
  eurasia_mainland %>% mutate(continent = "Eurasia")
)

# 转回 4326 并分配最终唯一的 cell_id
afroeurasia_final <- st_transform(afroeurasia_combined, 4326) %>%
  mutate(cell_id = 1:n()) %>%
  select(cell_id, continent) 

# ==========================================
# 6. 可视化检查
# ==========================================
p <- ggplot() +
  geom_sf(data = afroeurasia_final, aes(fill = continent), color = "black", linewidth = 0.3, alpha = 0.5) +
  scale_fill_manual(values = c("Africa" = "#E69F00", "Eurasia" = "#56B4E9")) +
  theme_minimal() +
  labs(title = "Afro-Eurasia Landgrid (Islands Removed)",
       subtitle = "Gibraltar & Bab-el-Mandeb connected in Adjacency Matrix")
print(p)

# ==========================================
# 7. 计算高精度的拓扑邻接矩阵 (Adjacency Matrix)
# ==========================================
# 修改点 2：放宽容差，使得海峡两岸的网格可以强行握手
# 将 snap 设定为 35000 (35km)。
# 曼德海峡(~26km) 和 直布罗陀海峡(~14km) 将会被判定为邻居。
# 而红海主体部分(>200km) 远大于这个容差，绝对不会发生错误连通。
afroeurasia_laea_final <- st_transform(afroeurasia_final, proj_laea)
nb_final <- poly2nb(afroeurasia_laea_final, snap = 35000)

adj_matrix <- nb2mat(nb_final, style = "B", zero.policy = TRUE)
rownames(adj_matrix) <- afroeurasia_final$cell_id
colnames(adj_matrix) <- afroeurasia_final$cell_id

g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
cost_matrix <- distances(g)

# 保存文件
st_write(afroeurasia_final, "data/geo/landgrid_afroeurasia_200km.gpkg", append = FALSE)
write.csv(adj_matrix, "data/geo/landgrid-adjmat-afroeurasia_200km.csv", row.names = TRUE)
write.csv(adj_matrix, "data/geo/landgrid-adjmat-afroeurasia_200km_norownames.csv", row.names = FALSE)
write.csv(cost_matrix, "data/geo/landgrid-distmat-afroeurasia_200km.csv", row.names = TRUE)





############################ map of southern asia ############################
world <- ne_countries(scale = "medium", returnclass = "sf")

# ==========================================
# 1. 提取整个亚洲主大陆 (到里海/乌拉尔山一线)
# ==========================================
# 定义提纯函数：打碎所有不相连的区块，按面积倒序排列，只保留最大的那一块纯净大陆
get_mainland_polygon <- function(sf_obj) {
  sf_obj %>%
    st_union() %>%
    st_make_valid() %>% 
    st_cast("POLYGON") %>%
    st_as_sf() %>%
    mutate(area = as.numeric(st_area(x))) %>%
    arrange(desc(area)) %>%
    slice(1) %>%
    st_geometry()
}

# 使用边界框精准控制西侧边界。
# xmin = 45 刚好切在里海和高加索山脉。
# 【注】：如果你后续的模拟需要中东作为连接非洲的“走廊”（如苏伊士地峡），请将 xmin 改为 30。
bbox_asia <- st_bbox(c(xmin = 55, xmax = 125, ymin = -15, ymax = 50), crs = st_crs(4326))

world_cropped <- suppressWarnings(st_crop(world, bbox_asia))

# 提取并提纯（这会把日本、印尼、马来群岛等所有没有陆地相连的岛屿全部物理剔除）
asia_outline <- world_cropped %>% 
  filter(continent %in% c("Asia", "Europe")) %>% 
  get_mainland_polygon()

# ==========================================
# 2. 投影转换与生成底层网格
# ==========================================
# 将投影中心点移动到亚洲腹地 (约新疆/中亚)，让大范围网格变形最小化
proj_laea <- "+proj=laea +lon_0=90 +lat_0=45 +datum=WGS84 +units=m +no_defs"

asia_laea <- st_transform(asia_outline, proj_laea)

a_meters <- 150000    # 150km
cell_size <- a_meters * sqrt(3)
universal_grid <- st_make_grid(asia_laea, cellsize = cell_size, square = FALSE)
universal_sf <- st_as_sf(data.frame(id = 1:length(universal_grid)), geometry = universal_grid)

# ==========================================
# 3. 独立裁切 (已彻底解决 POINT 报错)
# ==========================================
asia_inter <- suppressWarnings(st_intersection(universal_sf, asia_laea))

# 使用 st_collection_extract 专门抽取出所有的 Polygon（面），丢弃相切产生的点和线
asia_hex <- st_collection_extract(asia_inter, "POLYGON") %>% 
  st_cast("MULTIPOLYGON")

# ==========================================
# 4. 清洗组装与格式化
# ==========================================
asia_final <- st_transform(asia_hex, 4326) %>%
  mutate(
    cell_id = 1:n(),
    region = "Asia Mainland" 
  ) %>%
  select(cell_id, region) 

# ==========================================
# 5. 可视化检查
# ==========================================
p <- ggplot() +
  geom_sf(data = asia_final, aes(fill = region), color = "black", linewidth = 0.3, alpha = 0.5) +
  scale_fill_manual(values = c("Asia Mainland" = "#D55E00")) + 
  theme_minimal() +
  labs(title = "Asian Mainland Grid (150km)",
       subtitle = "Extended westward to the Caspian Sea, all islands removed")
print(p)

# ==========================================
# 6. 计算高精度的拓扑邻接矩阵和距离矩阵
# ==========================================
asia_laea_final <- st_transform(asia_final, proj_laea)

# 因为所有海岛都去掉了，网格之间全是硬连接，snap 保持 1000（1公里容差）防微小裂缝即可
nb_final <- poly2nb(asia_laea_final, snap = 1000)

adj_matrix <- nb2mat(nb_final, style = "B", zero.policy = TRUE)
rownames(adj_matrix) <- asia_final$cell_id
colnames(adj_matrix) <- asia_final$cell_id

g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
cost_matrix <- distances(g)

# ==========================================
# 7. 导出数据
# ==========================================
st_write(asia_final, "data/geo/landgrid_asia_150km.gpkg", append = FALSE)
write.csv(adj_matrix, "data/geo/landgrid-adjmat-asia_150km.csv", row.names = TRUE)
write.csv(adj_matrix, "data/geo/landgrid-adjmat-asia_150km_norownames.csv", row.names = FALSE)
write.csv(cost_matrix, "data/geo/landgrid-distmat-asia_150km.csv", row.names = TRUE)
