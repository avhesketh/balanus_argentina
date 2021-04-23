## creating maps for Figure S1 in Supplement 1

require(ggplot2)
require(ggmap)
require(maps)
require(mapproj)
require(mapdata)
require(rgeos)
require(maptools)
require(sp)
require(raster)
require(rgdal)
require(dismo)

## ggmap

library(ggmap)
library(ggsn)

## need valid Google key here to extract stamen map data

lat_bc <- c(48.5,52) 
long_bc <- c(-126,-123)

bc_inset_base <- get_map(location=c(-125.27,48.88), zoom=11, source = "stamen",
                         maptype = "terrain-background")

map_bp_inset <- ggmap(bc_inset_base) +
  geom_point(aes(x=-125.16477,y=48.82099), size = 6, colour = "darkblue", alpha = 0.3) +
  geom_point(aes(x=-125.16477,y=48.82099), size = 6, pch = 21, colour = "black") +
  annotate(geom = "text", y=48.82099, x=-125.27, label="Bluestone Pt.",
           color = "darkblue", fontface = "bold", size = 8, check_overlap = FALSE) +
  geom_point(aes(y = 48.8248, x = -125.1358), size = 6, colour = "black", alpha = 0.4) +
  annotate(geom = "text", y= 48.8, x=-125.122, label="Bamfield",
           color = "black", fontface = "bold", size = 8, check_overlap = FALSE) +
  annotate(geom = "text", y = 48.9, x = -125.2, label = "Barkley Sound", angle = 49, fontface = "italic", 
           color = "black", size = 8, check_overlap = FALSE) +
  theme(axis.ticks.length=unit(-0.25, "cm")) +
  theme(axis.text.y = element_text(margin = margin(l = 20, r = -45), size = 20)) +
  theme(axis.text.x = element_text(vjust = 8, size = 20))+
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(panel.border = element_rect(fill = NA, colour= "black", size = 2))
map_bp_inset

base_bc_large <- get_map(location=c(-122.5,51), zoom=5,
                         source = "stamen", maptype = "terrain-background")

map_bp_large <- ggmap(base_bc_large) +
  geom_point(aes(x=-125.16477,y=48.82099), size = 1, colour = "darkblue", alpha = 0.4) +
  geom_point(aes(x=-125.16477,y=48.82099), size = 1, pch = 21, colour = "black") +
  geom_point(aes(x=-125.275133, y = 48.904970), size = 2, pch = 22, colour = "black") +
  annotate(geom = "text", y=47.2, x=-128.5, label="Bluestone Point",
           color = "darkblue", fontface = "bold", size = 3, check_overlap = FALSE) +
  annotate(geom = "segment", y= 47.6, yend = 48.82, x = -128.5, xend=-125.1, colour = "darkblue") +
  annotate(geom = "text", y=49, x=-131.4, label="Barkley Sound",
          color = "black", fontface = "italic", size = 3, check_overlap = FALSE) +
  annotate(geom = "segment", y= 48.90, yend = 48.90, x = -126.9, xend=-125.55, colour = "black") +
  annotate(geom = "text", y=44, x=-130.5, label="Pacific Ocean",
           color = "black", fontface = "italic", size = 3, check_overlap = FALSE) +
  annotate(geom = "text", y =55, x = -127, label = "Canada", fontface = "italic",
           color = "black", size = 3, check_overlap = FALSE) +
  annotate(geom = "text", y =45, x = -115, label = "USA", fontface = "italic", 
           color = "black", size = 3, check_overlap = FALSE) +
  xlab("Longitude (˚W)") +
  ylab("Latitude (˚N)") +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(panel.border = element_rect(fill = NA, colour= "black", size = 1))

north2(map_bp_large, y = 0.9, x = 0.3, symbol = 10)
dev.off()

tiff("figures/S1a.tiff", units = "in", width = 3, height = 3, res = 600)

map_bp_large

base_arg_large <- get_map(location = c(-63.5, -47), zoom = 5,
                          source = "stamen", maptype = "terrain-background")

map_arg_large <- ggmap(base_arg_large) +
  geom_point(aes(y = -42.61972222, x = -64.873889), colour = "darkred", size = 1, alpha = 0.4) +
  geom_point(aes(y = -42.61972222, x = -64.873889), colour = "black", pch = 21, size = 1) + 
  geom_point(aes(y = -42.73, x = -64.58), pch = 22, colour = "black", size = 2) +
  annotate(geom = "text", y = -42.6, x = -55.5, label = "Punta Ameghino",
           colour = "darkred", fontface = "bold", size = 3, check_overlap = FALSE) +
  annotate(geom = "segment", y = -42.6, yend = -42.6, xend = -64.7, x = -61.2, colour = "darkred") +
  annotate(geom = "text", y = -44, x = -58, label = "Nuevo Gulf", size = 3, fontface = "italic") +
  annotate(geom = "segment", y = -44, yend = -43, xend = -64.2, x = -61.7, colour = "black") +
  annotate(geom = "text", y=-40, x=-55, label="Atlantic Ocean",
           color = "black", fontface = "italic", size = 3, check_overlap = FALSE) +
  annotate(geom = "text", y= -45, x=-69.5, label="Argentina", angle = 55,
           color = "black", fontface = "italic", size = 3, check_overlap = FALSE) +
  xlab("Longitude (˚W)") +
  ylab("Latitude (˚S)") +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(panel.border = element_rect(fill = NA, colour= "black", size = 1))

north2(map_arg_large, y = 0.9, x = 0.3, symbol = 10)
dev.off()

tiff("./figures/S1b.tiff", units = "in", width = 3, height = 3, res = 600)


base_arg_inset <- get_map(location = c(-64.95, -42.7), zoom = 11,
                          source = "stamen", maptype = "terrain-background")

map_arg_inset <- ggmap(base_arg_inset) +
  geom_point(aes(y = -42.61972222, x = -64.873889), colour = "darkred", size = 6, alpha = 0.3) +
  geom_point(aes(y = -42.61972222, x = -64.873889), colour = "black", pch = 21, size = 6) + 
  annotate(geom = "text", y = -42.593, x = -64.9, label = "Punta Ameghino",
           colour = "darkred", fontface = "bold", size = 8, check_overlap = FALSE) +
  annotate(geom = "text", y = -42.79, x = -64.97, label = "Puerto Madryn",
           colour = "black", fontface = "bold", size = 8, check_overlap = FALSE) +
  annotate(geom = "text", y = -42.71, x = -64.85, label = "Nuevo Gulf", size = 8, fontface = "italic",
           colour = "black") +
  geom_point(aes(y = -42.766, x = -65.030), colour = "black", size = 6, alpha = 0.3) +
  theme(axis.ticks.length=unit(-0.25, "cm")) +
  theme(axis.text.y = element_text(margin = margin(l = 20, r = -50), size = 20)) +
  theme(axis.text.x = element_text(vjust = 8, size = 20))
map_arg_inset

