#target_crs <- "+proj=robin"    # Robinson projection
target_crs <- "+proj=laea +lat_0=90 +lon_0=0 +datum=WGS84 +units=m +no_defs"    # arctic projection

# target ethnic groups
target_names <- c(
  "Karitiana", "Maya", "Pima", "Surui", "Chane", "Mayan", "Mixe", 
  "Mixtec", "Piapoco", "Quechua", "Tlingit", "Zapotec",  
  "Aleut", "Chukchi", "Itelman", "Eskimo_Chaplin", "Eskimo_Naukan", 
  "Eskimo_Sireniki", "Even", "Yakut", "Ulchi", "Altaian", "Mansi", 
  "Tubalar", "Oroqen", "Hezhen",
  "Han", "CHB", "Mongolian", "Mongola", "Daur", "Xibo",
  "Afanasievo", "Denisovan", "Altai", "Chagyrskaya", "Vindija" #no comtemporary samples
)
# from America: Karitiana, Maya, Pima, Surui, Chane, Mayan, Mixe, Mixtec, Piapoco, Quechua, Zapotec
# from CentralAsiaSiberia: Tlingit, Aleut, Chukchi, Itelman, Eskimo_Chaplin, Eskimo_Naukan, Eskimo_Sireniki, Even, Yakut, Ulchi, Altaian, Mansi, Tubalar, Mongola, Altai
# from EAST_ASIA: Oroqen, Hezhen, Han, Mongolian, Daur, Xibo
# not identified (and no concurrent samples): CHB, Afanasievo, Denisovan, Chagyrskaya, Vindija


set.seed(1123581321L)

# full tree sequence after pre-processing
D <- process_genetic_data_asiaamericas(18L)
print("summary information of full tree sequence:")
print_summary_stats(D)
# inspect_tree_nodes(D, n_show = 15)
#p_full <- plot_spatial_distribution(D)

# example of a subsample tree sequence
D_single_sub <- data_subsample(D)
print("summary information of a subsample tree sequence:")
print_summary_stats(D_single_sub)
# inspect_tree_nodes(D_single_sub, n_show = 15)
#p_sub <- plot_spatial_distribution(D_single_sub)

treeseq_write(D$ts, "data/trees/empirical_tree_asia-americas.trees")
write.csv(D$data, file="data/genetics/sample_states_asia-americas.csv", row.names=FALSE)
#write.csv(D$sample.coords, file="data/genetics/sample_coords_asia-americas.csv", row.names=FALSE)

# generate 100 subsets
num_subsets = 100L
for (i in 1:num_subsets)
{
  D2 = data_subsample(D)
  treeseq_write(D2$ts, sprintf("data/trees/empirical_tree_asia-americas_subset_%d.trees", i))
  write.csv(D2$data, file=sprintf("data/genetics/subsets/sample_states_asia-americas_subset_%d.csv", i), row.names=FALSE)
  #write.csv(D2$sample.coords, file=sprintf("data/genetics/subsets/sample_coords_asia-americas_subset_%d.csv", i), row.names=FALSE)
  
  #if (i %% 10 == 0) cat(sprintf("finished: %d / %d\n", i, num_subsets))
}

