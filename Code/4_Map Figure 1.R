###########################################
#                                         #
#     Code for creating Figure 1 map      #     
#                                         #   
# #########################################

# Author: Andrew Du

# NB: The following code requires the GRAY_HR_SR.tif
# file in the Datasets folder, but the file is too 
# large. File can be obtained upon request. 


# read in packages
library(terra)
library(giscoR)
library(tidyverse)
library(tidyterra)
library(ggspatial)
library(ggnewscale)
library(patchwork)

# read in files for creating map
dem <- rast("C:/Users/adu/OneDrive - Colostate/Desktop/Manuscripts/Du & Alemseged 2019_Paranthropus afar/GitHub/TaxonAbsence/Datasets/GRAY_HR_SR.tif") # world DEM
lakes <- vect("C:/Users/adu/OneDrive - Colostate/Desktop/Manuscripts/Du & Alemseged 2019_Paranthropus afar/GitHub/TaxonAbsence/Datasets/ne_10m_lakes.shp") # shapefile of world's lakes
sites <- read.csv("C:/Users/adu/OneDrive - Colostate/Desktop/Manuscripts/Du & Alemseged 2019_Paranthropus afar/GitHub/TaxonAbsence/Datasets/NISP data.csv")

# get other necessary geographic data
coast <- gisco_get_coastallines( # world coastlines for masking
  resolution = "20",
  epsg = "4326", # WGS84
  year = "2016"
)
africa <- gisco_get_countries( # country borders for Africa
  region = "Africa",
  year = "2016",
  epsg = "4326" # WGS84
)

# define eastern African extent using decimal degrees and crop objects
extent <- ext(c(32, 53, -12, 17))
EA <- crop(dem, extent) # crop world DEM
EA.lakes <- crop(lakes, extent) # crop world lakes 
EA.countries <- crop(vect(africa), extent) # crop coastlines

# mask area of interest for DEM (get rid of Indian Ocean)
EA1 <- mask(EA, vect(coast))

# transform into a data frame for plotting
EA.df <- as.data.frame(EA1, xy = TRUE)
colnames(EA.df)[3] <- "elevation"

# wrangles sites dataset for plotting purposes
sites1 <- sites %>%
  group_by(Site_Formation) %>%
  summarize(
    mean_lat = mean(Latitude),
    mean_long = mean(Longitude),
    paran = sum(Paran_nisp > 0) / n()) # figuring out how proportion of a site's time slices occupied by Paranthropus  

sites1$paran[sites1$paran > 0 & sites1$paran < 1] <- 0.5 # code all "sometimes" occupation as 0.5


# make Figure 1 

## main figure
p1 <- ggplot() +
  # DEM
  geom_raster( 
    data = EA.df,
    aes(
      x = x, 
      y = y, 
      fill = elevation)
  ) +
  theme_bw() + 
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.position = c(0.685, 0.2),
    legend.text = element_text(size = 10.5),
    legend.text.align = 0) +
  guides(fill = "none") +
  # make DEM gray-scale
  scale_fill_distiller(
    palette = "Greys", 
    limits = c(20, 250), # manually set to make DEM lighter
    direction = -1) + 
  # country borders
  geom_spatvector( 
    data = EA.countries,
    fill = NA,
    color = "gray40",
    lwd = 0.45
  ) +
  # lake borders
  geom_spatvector( 
    data = EA.lakes, 
    fill = "lightskyblue", 
    color = "dodgerblue4", 
    lwd = 0.45) +
  # required b/c fill is already used by DEM
  new_scale_fill() + 
  # plot sites
  geom_point( 
    data = sites1,
    mapping = aes(
      x = mean_long,
      y = mean_lat,
      fill = factor(paran)),
    size = 5,
    shape = 21) + 
  # fill color for sites
  scale_fill_manual(
    values = c(
      "0" = "firebrick", 
      "0.5" = "gold",
      "1" = "palegreen3"),
  # legend text
    labels = c(
      expression(italic("Paranthropus") ~ "never present"),
      expression(italic("Paranthropus") ~ "sometimes present"),
      expression(italic("Paranthropus") ~ "always present")),
    name = "") +
  # scale bar
  annotation_scale( 
    location = "br",
    text_cex = 0.85) +
  # site labels
  annotate(geom = "text", 
           x = 40.5,
           y = 12.5,
           label = "Mille-Logya") +
  annotate(geom = "text", 
           x = 43.6,
           y = 11.5,
           label = "Ledi-Geraru") +
  annotate(geom = "text", 
           x = 38.9,
           y = 11.3,
           label = "Gona") +
  annotate(geom = "text", 
           x = 42.1,
           y = 10.7,
           label = "Hadar") +
  annotate(geom = "text", 
           x = 41.7,
           y = 9.7,
           label = "Bouri") +
  annotate(geom = "text", 
           x = 36.2,
           y = 9.3,
           label = "Melka Kunture") +
  annotate(geom = "text", 
           x = 39,
           y = 5.2,
           label = "Konso") +
  annotate(geom = "text", 
           x = 34,
           y = 5.7,
           label = "Shungura") +
  annotate(geom = "text", 
           x = 38.8,
           y = 4,
           label = "Koobi Fora") +
  annotate(geom = "text", 
           x = 33.7,
           y = 4.1,
           label = "Nachukui") +
  annotate(geom = "text", 
           x = 38.3,
           y = -0.5,
           label = "Chemeron") +
  annotate(geom = "text", 
           x = 36.6,
           y = 0.5,
           label = "Kanjera South") +
  annotate(geom = "text", 
           x = 34,
           y = -1,
           label = "Nyayanga") +
  annotate(geom = "text", 
           x = 37.3,
           y = -2.6,
           label = "Peninj") +
  annotate(geom = "text", 
           x = 33.6,
           y = -2.5,
           label = "Olduvai") +
  annotate(geom = "text", 
           x = 35,
           y = -4,
           label = "Laetoli") +
  annotate(geom = "text", 
           x = 33.9,
           y = -9,
           label = "Malema")
  
## Inset map of Africa
p2 <- ggplot(africa) + 
  geom_sf(fill = "gray30") + # color of Africa
  coord_sf(
    xlim = c(-15, 55), # define extent of Africa
    ylim = c(-33, 35)
  ) + 
  theme( # get rid of frills
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    panel.grid.major = element_blank(),
    panel.background = element_rect("white"),
    panel.border = element_rect(
      fill = NA, 
      linewidth = 0.8),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0)
  ) +
  annotate( # add bounding rectangle
    "rect",
    xmin = ext(EA1)[1],
    xmax = ext(EA1)[2],
    ymin = ext(EA1)[3],
    ymax = ext(EA1)[4],
    color = "black",
    linewidth = 0.8,
    fill = "white",
    alpha = 0.8) +
  geom_text( # add "Study Area" text
    x = mean(ext(EA1)[1:2]),
    y = mean(ext(EA1)[3:4]),
    label = "Study\nArea",
    size = 5,
    fontface = "italic"
  ) 

# Combine main plot and inset together
p1 + inset_element(p2, 
                   left = 0.555, 
                   bottom = 0, 
                   right = 0.955, 
                   top = 1,
                   align_to = "full")

# Save the figure
ggsave(filename = "C:/Users/adu/OneDrive - Colostate/Desktop/Manuscripts/Du & Alemseged 2019_Paranthropus afar/GitHub/TaxonAbsence/Figures/Fig 1.pdf")