##### parameters #####
WORLD = my_world
TIMEDEPTH = 7000

data = read.csv(sprintf("data/genetics/sample_states_%s.csv", WORLD))
pops_to_sample = sort(unique(data[,2])-1L)
landgrid = st_read("data/geo/landgrid_afro-eurasia.gpkg", quiet=TRUE)

# extract estimated location of oldest ancestors for all sampled replications
y_friction <- get_y(world=WORLD, map="friction", timedepth=TIMEDEPTH)
y_naive <- get_y(world=WORLD, map="naive", timedepth=TIMEDEPTH)

coords <- st_coordinates(st_centroid(landgrid$geom))

# count frequency
count_friction <- tabulate(y_friction[,2], nrow(coords))
count_naive <- tabulate(y_naive[,2], nrow(coords))

# plot
output_file <- sprintf("output/figures/compare_ancestor_estimates_%s.png", WORLD)
pdf(output_file, width = 7, height = 6)

par(mar=c(1,1,2,1))
plot(landgrid$geom, border=8, lwd=0.5, main="ancestor estimates")
coords <- st_coordinates(st_centroid(landgrid$geom))
# compute estimation error based on cost matrix
x_range <- range(coords[,1])
offset <- (x_range[2] - x_range[1]) * 0.005   # for display
coords_left  <- coords
coords_right <- coords
coords_left[,1]  <- coords_left[,1]  - offset   # naive on the left
coords_right[,1] <- coords_right[,1] + offset   # friction on the right
scale_factor <- 0.2 

points(coords_left,
       cex = scale_factor *  sqrt(count_naive),
       pch = 21,
       bg = "#66c2a5",
       col = NA)

points(coords_right,
       cex = scale_factor * sqrt(count_friction),
       pch = 21,
       bg = "#fc8d62",
       col = NA)

origin_idx <- 58
plot(st_geometry(landgrid)[origin_idx],
     col = rgb(1,0,0,0.2),
     border = "red",
     lwd = 2,
     add = TRUE)

legend("topright",
       legend=c("naive", "friction", "true origin"),
       pch=c(21,21,NA),
       pt.bg=c("#66c2a5", "#fc8d62", NA),
       col=c(NA,NA,"red"),
       lty=c(NA,NA,1),
       lwd=c(NA,NA,2),
       pt.cex=c(1.5,1.5,NA),
       bty="n",
       cex=0.9)

dev.off()
message(sprintf("saved figure to: %s", output_file))

cat("naive:", count_naive[58], " out of ", sum(count_naive), " on target.\n")
cat("friction:", count_friction[58], " out of ", sum(count_friction), " on target.\n")

