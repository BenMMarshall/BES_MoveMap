
library(dplyr)
library(lubridate)
library(move)
# library(adehabitatHR) #  not needed if not plotting OC
library(spatialEco)

library(ggplot2)
library(ggspatial)
library(magick)
library(ggtext)
library(ggforce)
library(scico)
library(ggpubr)

# having to use ragg to fix transparency rendering issues
# install.packages("ragg")
library(ragg)

ophaData <- read.csv("./Data/Ophiophagus hannah 2014_2018.csv")
nestData <- read.csv("./Data/ApproxNestLocations.csv")
names(nestData)[1] <- "year"

ophaData <- ophaData %>% 
  filter(id == "AF017") %>% 
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%dT%H:%M:%S",
                               tz = "Asia/Bangkok"),
         gps.accuracy = ifelse(is.na(gps.accuracy) | gps.accuracy == 0,
                               mean(gps.accuracy, na.rm = TRUE),
                               gps.accuracy)) %>% 
  arrange(datetime)

moveObj <- move(x = ophaData$x, y = ophaData$y, time = ophaData$datetime, data = ophaData,
                proj = CRS("+proj=utm +zone=47 +datum=WGS84 +units=m +no_defs"),
                animal = ophaData$id)

dbbmmObj <- brownian.bridge.dyn(moveObj,
                                location.error = ophaData$gps.accuracy,
                                margin = 5, window = 25,
                                ext = 1, dimSize = 1000)

## DROPPED AS NOT PLOTTING OC
# dbbmmSP <- as(dbbmmObj, "SpatialPixelsDataFrame")
# dbbmmSP <- new("estUD", dbbmmSP)
# dbbmmSP@vol = FALSE; dbbmmSP@h$meth = "dBBMM"
# dbbmmUD <- getvolumeUD(dbbmmSP, standardize = TRUE) # convert UD values to volume
# 
# plot(dbbmmUD)
# 
# polyCon <- getverticeshr(dbbmmUD, percent = 99)
# polyDf <- fortify(polyCon)
# 
# dbbmmUDdf <- as.data.frame.estUD(dbbmmUD)
# names(dbbmmUDdf) <- c("od", "x", "y")

ophaData$moVar <- getMotionVariance(dbbmmObj)

coreArea <- readOGR("./Data/SBRcore.shp")
klong <- readOGR("./Data/Klong_2.shp")
klong <- spTransform(klong, CRS("+proj=utm +zone=47 +datum=WGS84 +units=m +no_defs"))
plot(coreArea)
plot(klong, add = TRUE)
points(ophaData$x, ophaData$y)

ophaDataPlot <- point.in.poly(SpatialPointsDataFrame(ophaData[,5:6],
                                                     proj4string = CRS("+proj=utm +zone=47 +datum=WGS84 +units=m +no_defs"),
                                                     data = ophaData),
                              coreArea)
ophaDataPlot <- as.data.frame(ophaDataPlot)

# Time to plot ------------------------------------------------------------

nestImg <- image_read("OPHA017_nest.png" , depth = 16)
nestRaster <- as.raster(nestImg)
hideImg <- image_read("OPHA017_hide.png" , depth = 16)
hideRaster <- as.raster(hideImg)

ophaDataPlot <- ophaDataPlot %>% 
  mutate(yearday = yday(datetime),
         year = year(datetime),
         nesting = ifelse(yearday  > 140 & yearday < 220,
                          "Nesting", "Non-nesting")) %>% 
  dplyr::select(-ogc_fid, -objectid, -shape_area, -coords.x1, -coords.x2, -shape_leng)

yearcolours <- scico(7, palette = "roma", begin = 0.1, end = 0.9)[c(1,2,3)]

protectedcolours <- scico(7, palette = "roma", begin = 0.1, end = 0.9)[6:7]

scalelocy <- 1608200

(map <- ophaDataPlot %>% 
    ggplot() +
    geom_spatial_polygon(data = coreArea, aes(x = long, y = lat, group = group),
                         alpha = 0.5, fill = NA,
                         colour = protectedcolours[2], size = 2,
                         linetype = 1) +
    geom_spatial_polygon(data = coreArea, aes(x = long+20, y = lat, group = group),
                         alpha = 0.5, fill = NA, colour = "grey45", size = 1,
                         linetype = 1) +
    geom_spatial_polygon(data = klong, aes(x = long, y = lat, group = group),
                         alpha = 0.75, fill = protectedcolours[1]) +
    geom_path(aes(x = x, y = y, colour = as.factor(year)), size = 0.5, alpha = 0.75) +
    geom_point(aes(x = x, y = y, colour = as.factor(year)), size = 0.5, alpha = 0.75) +
    # cap and death points
    geom_point(data = ophaDataPlot[c(1, dim(ophaDataPlot)[1]),],
              aes(x = x, y = y), size = 2, alpha = 1,
              colour = "black") +
    geom_point(data = ophaDataPlot[c(dim(ophaDataPlot)[1]),],
              aes(x = x, y = y), size = 1.1, alpha = 1,
              colour = "white", pch = 4, stroke = 0.75) +
    geom_point(data = ophaDataPlot[1,],
              aes(x = x, y = y), size = 1.1, alpha = 1,
              colour = "white", pch = 6, stroke = 0.75) +
    geom_richtext(data = ophaDataPlot[c(1, dim(ophaDataPlot)[1]),],
                  aes(x = x +c(20,-22), y = y+c(18,-20), label = c("Capture", "Death"), hjust = c(0,1)),
                  fill = NA, label.color = NA,
                  label.padding = unit(rep(0, 4), "pt"),
                  size = 3, fontface = 3) +
    scale_colour_manual(values = yearcolours, na.value = "black") +
    ## protect area boundary labels
    annotate("text", x = 817710, y = 1607700, label = "Agriculture",
             lineheight = 0.85, fontface = 4, vjust = 1, hjust = 0,
             colour = "grey45") +
    annotate("segment",
             x = 817710, xend = 817710 +550,
             y = 1607720, yend = 1607720,
             size = 1,
             colour = "grey45", arrow = arrow(angle = 25, length = unit(3, "mm"),
                                              type = "closed")) +
    annotate("text", x = 817560, y = 1607700, label = "Protected\nForest",
             lineheight = 0.85, fontface = 4, vjust = 1, hjust = 1,
             colour = protectedcolours[2]) +
    annotate("segment",
             xend = 817560 -550, x = 817560,
             y = 1607720, yend = 1607720,
             size = 1,
             colour = protectedcolours[2], arrow = arrow(angle = 25, length = unit(3, "mm"),
                                              type = "closed")) +
    ## nesting location points
    geom_point(data = nestData, aes(x = Easting, y = Northing, colour = as.factor(year)),
               size = 4, fill = NA, pch = 21, stroke = 2) +
    geom_point(data = nestData, aes(x = Easting, y = Northing),
               colour = "black",
               size = 5, fill = NA, pch = 21) +
    geom_curve(data = nestData, aes(xend = Easting, x = 817000,
                                    yend = Northing, y = 1607000),
               arrow = arrow(angle = 25, length = unit(2, "mm"),
                             type = "closed"), size = 1,
               curvature = 0.3) +
    ## nesting locations label
    annotation_raster(nestRaster,
                      xmin = 816950-250, xmax = 816950+250,
                      ymin = 1607100-250, ymax = 1607100+250) +
    geom_circle(data = data.frame(x = c(816950), y = c(1607100),
                                  r = 250),
                aes(x0 = x, y0 = y, r = r), size = 1.5) +
    annotate("label", x = 816950, y = 1607100+250,
             label = "Forested\nnest locations",
             lineheight = 0.85, fontface = 2, vjust = 0.5, hjust = 0.5,
             label.padding = unit(2.5, "mm"),
             label.r = unit(0.5, "lines"),
             label.size = 1.2) + 
    ## usual HR arrows
    geom_curve(data = data.frame(x = c(818950, 819250, 819500),
                                 xend = rep(819450, 3),
                                 y = c(1607380, 1607550, 1607650),
                                 yend = rep(1606750, 3)),
               aes(xend = x, x = xend,
                   yend = y, y = yend),
               arrow = arrow(angle = 25, length = unit(2, "mm"),
                             type = "closed"), size = 1,
               curvature = 0.3) +
    ## usual HR label
    annotation_raster(hideRaster,
                      xmin = 819450-250, xmax = 819450+250,
                      ymin = 1606750-250, ymax = 1606750+250) +
    geom_circle(data = data.frame(x = c(819450), y = c(1606750),
                                  r = 250),
                aes(x0 = x, y0 = y, r = r), size = 1.5) +
    annotate("label", x = 819450, y = 1606750-250,
             label = "Usual home range along\nmain irrigation canal",
             lineheight = 0.85, fontface = 2, vjust = 0.5, hjust = 0.5,
             label.padding = unit(2.5, "mm"),
             label.r = unit(0.5, "lines"),
             label.size = 1.2) + 
    ## journey label
    geom_curve(data = data.frame(x = c(817850,817850),
                                 xend = c(817800,818000),
                                 y = c(1606320,1606345),
                                 yend = c(1606880,1606140)),
               aes(xend = xend, x = x,
                   yend = yend, y = y),
               arrow = arrow(angle = 25, length = unit(2, "mm"),
                             type = "closed"), size = 1,
               curvature = -0.3) +
    annotate("label", x = 817830, y = 1606325,
             label = "Annual journey\nto/from the forest",
             lineheight = 0.85, fontface = 2, vjust = 0, hjust = 0.5,
             label.padding = unit(2.5, "mm"),
             label.r = unit(0.5, "lines"),
             label.size = 1.2) + 
    ## scale bar and N
    annotate("segment", x = 818900, xend = 819900, y = scalelocy, yend = scalelocy,
             size = 1.5) +
    annotate("text", x = 819400, y = scalelocy-90, vjust = 0, label = "1 km",
             fontface = 2) +
    annotate("segment", x = 819850, xend = 819850, y = scalelocy-80, yend = scalelocy-50,
             size = 1, arrow = arrow(angle = 25, length = unit(5, "mm"),
                                     type = "closed")) +
    annotate("text", x = 819850, y = scalelocy-100, colour = "white",
             label = "N", fontface = 2, hjust = 0.5) +
    ## title
    annotate("richtext", x = 816470, y = 1608280,
             label = "**Nesting movements of a King Cobra**<br><i>(Ophiophagus hannah)</i>",
             hjust = 0, vjust = 1, lineheight = 0.85,
             fill = alpha("white", 0.9), label.color = NA,
             label.padding = unit(rep(8, 4), "pt"),
             size = 5) +
    ## theming and look extras
    coord_equal(xlim = range(ophaData$x)+c(-250,250),
                ylim = range(ophaData$y)+c(-250,+250),
                expand = 0, clip = "off") +
    theme_void() +
    theme(legend.position = "none",
          panel.background = element_rect(fill = NA, colour = NA),
          plot.background = element_rect(fill = NA, colour = NA)
          )
    )

dateBreakData <- data.frame(yday = 1:365,
                            month = month(1:365 - 1 + as.Date("1990-01-01"),
                                          label = TRUE, abbr = TRUE)) %>% 
  group_by(month) %>% 
  summarise(startDate = min(yday),
            length = n())

timeseriesdata <- ophaDataPlot %>% 
  mutate(moVarTrans = ifelse(year == 2016, moVar + 165,
                             ifelse(year == 2015, moVar + 340,
                                    moVar)),
         yearmod = ifelse(!is.na(name), NA, year))

capdeathdf <- rbind(
  timeseriesdata %>% 
    filter(!is.na(moVarTrans)) %>% 
    filter(year == 2015) %>% 
    filter(yearday == min(yearday)) %>% 
    dplyr::select(yearday, moVarTrans),
  timeseriesdata %>% 
    filter(!is.na(moVarTrans)) %>% 
    filter(year == 2017) %>% 
    filter(yearday == max(yearday)) %>% 
    dplyr::select(yearday, moVarTrans))
capdeathdf <- capdeathdf[c(1,4),]
capdeathdf$lab <- c("Capture", "Death")
capdeathdf$hjust <- c(1,0)

(timeplot <- timeseriesdata %>% 
    ggplot() +
    geom_hline(yintercept = c(0, 165, 340), linetype = 2, size = 0.25,
               alpha = 0.5, colour = "black") +
    ## year annotations
    geom_richtext(data = data.frame(year = c(paste0("<span style='color:", yearcolours[3], "'>**2017**</span>"),
                                             paste0("<span style='color:", yearcolours[2], "'>**2016**</span>"),
                                             paste0("<span style='color:", yearcolours[1], "'>**2015**</span>")),
                                y = c(135, 310, 495),
                                x = rep(360,3)),
              aes(x = x, y = y, label = year), fontface = 2,
              fill = NA, label.color = NA,
              label.padding = grid::unit(rep(0, 4), "pt"),
              vjust = 1, hjust = 1) +
    ## highlighting nesting periods
    geom_rect(data = timeseriesdata %>% 
                group_by(year) %>% 
                filter(!is.na(name), !yearday == 100) %>% 
                summarise(ymin = min(moVarTrans, na.rm = TRUE),
                          ymax = max(moVarTrans, na.rm = TRUE),
                          xmin = min(yearday),
                          xmax = max(yearday)),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = c(495, 310, 135), group = year),
              alpha = 0.15, fill = protectedcolours[2]) +
    geom_path(aes(x = yearday, y = moVarTrans, colour = as.factor(yearmod), group = year)) +
    scale_colour_manual(values = yearcolours, na.value = protectedcolours[2]) +
    ## capture and death annotations
    geom_point(data = capdeathdf,
               aes(x = yearday, y = moVarTrans),
               size = 2) +
    geom_point(data = capdeathdf,
               aes(x = yearday, y = moVarTrans, shape = c("Capture", "Death")),
               size = 1.1, stroke = 0.75, colour = "white") +
    scale_shape_manual(values = c(6, 4)) +
    geom_richtext(data = capdeathdf,
               aes(x = yearday +c(-2,3), y = moVarTrans +c(0,-2), label = lab, hjust = hjust),
               fill = NA, label.color = NA,
               label.padding = unit(rep(0, 4), "pt"),
               size = 3, fontface = 3) +
    ## motion variance arrows and text
    geom_segment(data = data.frame(x = rep(-10, 3),
                                   xend = rep(-10, 3),
                                   y = c(0,165,340),
                                   yend = c(135,310,480)),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 size = 1, arrow = arrow(angle = 25, length = unit(1.5, "mm"),
                                           type = "closed")) +
    annotate("richtext", x = -5, y = 485, label = "**\u03C3\u00B2m**<br><i>(higher = more movement)</i>",
             hjust = 0, vjust = 1, lineheight = 0.85,
             fill = NA, label.color = NA,
             label.padding = unit(rep(0, 4), "pt"),
             size = 3) +
    ## other scales and formatting extas
    scale_y_continuous(limits = c(-20, 525), expand = c(0,0)) +
    scale_x_continuous(breaks = dateBreakData %>% 
                         pull(startDate),
                       labels = dateBreakData %>% 
                         pull(month),
                       limits = c(-20, 365), expand = c(0,0)) +
    labs(x = "Month") +
    theme_void() +
    theme(
      legend.position = "none",
      axis.title.x = element_text(hjust = 1, vjust = 1, face = 2),
      axis.line.x = element_line(size = 1.2, lineend = "square"),
      axis.ticks.x = element_line(size = 1.2, lineend = "square"),
      axis.ticks.length.x = unit(5, "mm"),
      axis.text.x = element_text(hjust = -0.2, vjust = 5, face = 2),
      panel.background = element_rect(fill = NA, colour = NA),
      plot.background = element_rect(fill = alpha("white", 1), colour = NA),
      )
)

agg_png(filename = "Figure - combined.png", width = 160, height = 160,
        res = 300, units = "mm")
ggarrange(plotlist = list(map, timeplot),
          nrow = 2, ncol = 1,
          heights = c(1, 0.65))
dev.off()

ggsave("Figure - combined.pdf", width = 180, height = 160,
       units = "mm")



