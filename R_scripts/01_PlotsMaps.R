# Colour Polymorphism Manuscript
# 1. Plotting maps

rm(list=ls())

# load libraries
library(openxlsx)
library(dplyr)
library(tibble)
library(geosphere)
library(ggplot2)
library(tmap)
library(sf)
library(raster)
library(rnaturalearth)
library(RColorBrewer)
library(ggforce)
library(scatterpie)
library(ggiraph)
library(PieGlyph)


# set working directory
setwd("C:/Users/feiner/Dropbox/DataAnalysis")

############### Fig. 1B - P. muralis distribution map

world <- ne_countries(scale = "large", returnclass = "sf") #change to large for final plot!
sf_use_s2(FALSE)
europe_cropped <- st_crop(world, xmin = -10.00, xmax = 35.00, ymin = 35.00, ymax = 52.00)

#read in distribution map of P. muralis from IUCN
Pmur_dist <- st_read("ToSubmit/Pmuralis_distribution/data_0.shp")
#pdf("./Plots/DistributionMap.pdf", useDingbats = F)
ggplot(data = europe_cropped) + 
  geom_sf(data = europe_cropped, fill= "white", color = "grey", size = 0.1) +
  geom_sf(data = Pmur_dist, fill = "lightgrey", color = NA) +
  coord_sf(xlim = c(-10.00, 35.00), ylim = c(35.00, 52.00), expand = F) +
  annotate("rect", xmin = 7.57, xmax = 16.44, ymin = 39.78, ymax = 46.37,color = "black", fill=NA) +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.background = element_rect(fill = "#EBF9FA"), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
#dev.off()


############### Fig. 1C - Distribution of morph frequencies in Italy

# merging locations for plotting morph frequency since some locations 
# are too close to each other for plotting

# Read in data set
freq_combined <- read.csv("ToSubmit/DataS1_220_Populations.csv") 

# Calculate pairwise distances (in meters)
dist_matrix <- distm(freq_combined[, c("Longitude", "Latitude")], fun = distHaversine)

# Convert to kilometers
dist_matrix <- dist_matrix / 1000

# Perform hierarchical clustering
hc <- hclust(as.dist(dist_matrix), method = "complete")

# Cut tree to form clusters. Using a 5 km threshold produces good results.
clusters <- cutree(hc, h = 5)

# Add clusters to the data frame
freq_combined$Cluster <- clusters

#summarize data for each cluster
freq_combined_sum <- freq_combined %>%
  group_by(Cluster) %>%
  summarise(
    Combined_Abbpop = paste(abbpop, collapse = ", "),  # Combine all 'abbpop' values into a string
    cluster_Latitude = mean(Latitude, na.rm = TRUE),
    cluster_Longitude = mean(Longitude, na.rm = TRUE),
    cluster_white = sum(white, na.rm = TRUE),
    cluster_orange = sum(orange, na.rm = TRUE) + 0.5 * sum(orange_yellow, na.rm = TRUE),
    cluster_yellow = sum(yellow, na.rm = TRUE) + 0.5 * sum(orange_yellow, na.rm = TRUE),
    cluster_totalN = sum(totalN, na.rm = TRUE)
  )

# Plot the morph frequencies as pie charts on map
#pdf("./Plots/MorphFreqsPie.pdf", useDingbats = F)
ggplot() + 
  geom_sf(data = europe_cropped, fill= "white", color = "grey", size = 0.1) +
  geom_sf(data = Pmur_dist, fill = "lightgrey", color = NA) +
  coord_sf(xlim = c(7.57, 16.44), ylim = c(39.78, 46.37), expand = F) +
  xlab("Longitude") + ylab("Latitude") + labs(type="Morph") + theme_void() +
  geom_pie_glyph(aes(y = cluster_Latitude, x = cluster_Longitude),
                 data = freq_combined_sum, colour = 'black', alpha = 0.8, radius = 0.25, linewidth = 0.05,
                 slices = c('cluster_white','cluster_orange','cluster_yellow'))+
  scale_fill_manual(values = c('white','orange', 'yellow'), name = 'Morph')+
  theme(panel.background = element_rect(fill = "#EBF9FA"), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
#dev.off()


############### Fig. 2A - Greenness across Italy

# Average greenness exists for 48 locations (each with at least 5 males with greenness data and lineage assignment). 
# Read in data set
data <- read.csv("ToSubmit/DataS2_148_Populations_ForGreenness.csv") 
Italy_cropped <- st_crop(europe_cropped, xmin = 7.57, xmax = 16.44, ymin = 39.78, ymax = 46.37)
g <- st_transform(Italy_cropped, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

coordinates(data) <- ~Longitude+Latitude
proj4string(data) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84")
sp_mydata <- st_as_sf(data, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

#extent(data) #this gives the extreme coordinates of all populations, used in next step (plus margin)
grd <- expand.grid(x = seq(from = 7.67, to = 13.63, length.out = 10000),
                   y = seq(from = 41.28, to = 44.96,  length.out = 10000)) # change this for modifying the frame around data points
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(data) <- proj4string(grd)
Krig = gstat::idw(Green~1,data,newdata=grd, idp=4.0)

# Convert to raster object the interpolation and clip it to the Italian raster
r       <- raster(Krig)
r.m     <- mask(r,g)
cbp1 <- c("#996633","#996633","#666600","#669900","#33CC00","#339900","#336600")# Manual set of the palette.WHAT IS IMPORTANT HERE IS TO BE SURE THAT IF WE CHANGE ANY OF THE PALETTES HERE WE ALSO NEED TO CHANGE THEM IN THE BARPLOT FOR THE ADMIXTURE
 #cbp1 <- c("#4d3e1d","#675326","#816830","#66CC33","#339900","#336600")
crs(r.m) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
 
 # Plot
green <- tm_shape(r.m, raster.downsample = FALSE) + #plot interpolation in raster
   tm_raster(n=10, palette = cbp1, interpolate=TRUE, #plot using the palet and interpolating across map
             legend.show = TRUE, title = "") + # add 'style = "cont"' to get the continuous scale
   tm_shape(Italy_cropped) + # Add a shape for the land outline
   tm_borders(col = "grey", lwd = 0.1) + # Grey outline for the land
  tm_shape(sp_mydata) +
   tm_symbols(size = 0.95, shape = 21, col = "snow", fill = "snow") +
   tm_symbols(size=0.85, shape=21, col="Green", palette = cbp1, border.col = "snow") +
   tm_legend(position = c("left", "bottom"), frame = FALSE, legend.text.size = 3)+ #position of the legend
   tm_layout(bg.color = "#EBF9FA", inner.margins = c(0, 0, 0, 0))+  #lay out of the graph
   tm_scale_bar(breaks = c(0, 50, 100), # Specify 100-km scale bar
    position = c("center", "bottom"), # Place scale bar on the bottom right
    text.size = 0.8) +
   tm_layout(frame = T)

 tmap_save(green, "./Plots/GreennessMap.pdf")
