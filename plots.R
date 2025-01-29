library(tidyverse)
library(zoo)
library(RColorBrewer)
library(sp)
library(ggrepel)
library(patchwork)
library(reticulate)
library(ggplot2)
library(ncdf4)
library(plotly)

##### Read in Data #####

tuned2000 <- read_csv(file = 'data/DBSCAN_output_esmm__NLDAS_tuned_2000.csv') %>%
  dplyr::select(!`...1`)%>%
  mutate(date = as_date(parse_date_time(x = paste(year, day), orders = "yj")))

list_of_files <- list.files(path = "output/DBSCAN/ESMM/NLDAS",
                            recursive =FALSE,
                            pattern = "\\.csv$",
                            full.names = TRUE)

svdi_clusters <- readr::read_csv(list_of_files, id = "file_name",col_types = cols(...1 = col_skip())) %>%
  dplyr::select(!file_name)


land_coords <- read_csv("data/land_coordinates.csv")
NLDAS_landcoordinates = read_csv("NLDASlandcoordinates.csv")

np <- import("numpy")
# data reading
autocorr_2003 <- np$load("../bilevelclustering/data/autocorrelation_2003.npy")
autocorr_r = py_to_r(autocorr_2003)

lonlats <- nc_open('data/NLDAS_lon_lat.nc')

lon <- ncvar_get(lonlats, "lon")
lat <- ncvar_get(lonlats, "lat")
lon_lats = expand_grid(lat,lon)

nc_data <- nc_open('data/NLDAS/NLDAS_SVDI_2003.nc') 
svdi.array <- ncvar_get(nc_data, "SVDI")
fillvalue <- ncatt_get(nc_data, "SVDI", "_FillValue")
svdi.array[svdi.array == fillvalue$value] <- NA


lagged_mean <- read_csv('../bilevelclustering/data/laggedtimeseries_mean.csv') %>%
  select(!...1) 

standard_mean <- read_csv('../bilevelclustering/data/standardtimeseries_mean.csv') %>%
  select(!...1) 

cluster_assignment = read_csv('../bilevelclustering/data/clusterassignment_laggedstandard_withlongs.csv') %>%
  select(!...1) %>%
  pivot_longer(cols = contains('_'),
               names_to = 'kmeanstype',
               names_prefix = 'ClusterID_',
               values_to = 'Cluster ID') %>%
  mutate( kmeanstype = case_when(
    kmeanstype == 'Standard' ~ "Standard",
    kmeanstype == 'Lagged' ~ "Space-time")) 

ls_hist <- read_csv('../bilevelclustering/histogram_laggedstandard.csv')%>%
  select(!...1) %>%
  rename(Standard = labels_standard,
         `Space-time` = labels_lagged)

ls_hist_fixed <- ls_hist %>%
  pivot_longer(cols = !SVDI,
               names_to = 'kmeanstype',
               values_to = 'Cluster ID') %>%
  mutate(`Cluster ID` = as.factor(`Cluster ID`))

lagged_fixed <-    lagged_mean %>%
  pivot_longer(cols = everything(),
               names_to = c("Cluster",'kmeanstype'),
               names_sep = '_',
               values_to = 'Mean SVDI') %>%
  mutate(Day = rep(seq(1,365), each = 5),
         kmeanstype = case_when(
           kmeanstype == 'standard' ~ "Standard",
           kmeanstype == 'lagged' ~ "Space-time"))


standard_fixed <- standard_mean %>%
  pivot_longer(cols = everything(),
               names_to = c("Cluster",'kmeanstype'),
               names_sep = '_',
               values_to = 'Mean SVDI')%>%
  mutate(Day = rep(seq(1,365), each = 5),
         kmeanstype = case_when(
           kmeanstype == 'standard' ~ "Standard",
           kmeanstype == 'lagged' ~ "Space-time"))


all_centroids = rbind(lagged_fixed, standard_fixed)


### Functions #####
pentad_id = function(day){
  i = 1
  counter = 1
  while (i <= 365){
    if(day >= i & day < i+5 ){
      return(counter)
    }
    else{
      i = i + 5 
      counter = counter + 1
    }
  }
  
}


###### Creating SVDI categories #####
svdi_clusters <- svdi_clusters %>%
  mutate(SVDI_severity = case_when(
    SVDI <= -3  ~ 'Extreme Wetness',
    SVDI > -3 & SVDI <= -2 ~ 'Severe Wetness',
    SVDI > -2 & SVDI <= -1 ~ 'Moderate Wetness',
    SVDI > -1 & SVDI <= -0.5 ~ 'Mild Wetness',
    SVDI > -0.5 & SVDI <= 0.5 ~ 'Neutral',
    SVDI > 0.5 & SVDI <= 1 ~ 'Mild Dryness',
    SVDI > 1 & SVDI <= 1.5 ~ 'Moderate Dryness',
    SVDI > 1.5 & SVDI <= 2 ~ 'Severe Dryness',
    SVDI > 2 ~ 'Extreme Dryness'))
svdi_clusters <- svdi_clusters %>%
  mutate(SVDI_severity = factor(SVDI_severity, 
                                levels =  rev(c("Extreme Wetness",
                                                'Severe Wetness',
                                                'Moderate Wetness',
                                                'Mild Wetness',
                                                'Neutral',
                                                'Mild Dryness',
                                                'Moderate Dryness',
                                                'Severe Dryness',
                                                'Extreme Dryness'))),
         drought_severity = factor(drought_severity,
                                   levels = c("mild",'moderate','severe','extreme')),
         cluster_id = as_factor(cluster_id),
         #day starts in 0 from python
         day = day + 1 )

svdi_clusters <- svdi_clusters %>% 
  rowwise() %>%
  mutate(pentad_num = pentad_id(day))

day_to_season <- tibble(day = seq(0,365)) %>%
  mutate(season = case_when(
    day < 79 ~ 'Winter',
    day >= 79 & day < 171 ~ 'Spring',
    day >= 171 & day < 264 ~ 'Summer',
    day >= 264 & day < 355 ~ 'Fall',
    day >= 355 ~ 'Winter'))

svdi_clusters <- svdi_clusters %>%
  left_join(day_to_season)

svdi_clusters<- svdi_clusters %>%
  mutate(date = parse_date_time(x = paste(year, day), orders = "yj"))

#### Plotting Begins #####
timeseries = all_centroids %>%
  ggplot(aes(x = Day, y = `Mean SVDI`, color = Cluster)) +
  geom_line(linewidth = 1 ) +
  facet_grid(rows = vars(kmeanstype)) +
  scale_color_viridis(discrete=TRUE,
                      guide = guide_legend(reverse = TRUE)) +
  theme_minimal() 

histmean_labels = ls_hist_fixed %>%
  group_by(kmeanstype, `Cluster ID`) %>%
  summarise(mean = mean(SVDI)) %>%
  mutate(SVDI = rep(2.5,5),
         y_pos = seq(2500,6000,length.out = length(levels(ls_hist_fixed$`Cluster ID`))),
         mean = signif(mean, 2),
         mean_label = paste0('Cluster ',`Cluster ID`,' mean: ', mean)) %>%
  arrange(`Cluster ID`)

histograms <- ls_hist_fixed %>% 
  ggplot() +
  geom_histogram(aes(x = SVDI,fill = `Cluster ID`),bins =200) +
  geom_text(data = histmean_labels,  aes(label = mean_label,x = SVDI, y = y_pos)) +
  scale_fill_viridis(discrete = TRUE) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(-3,4)) +
  facet_grid(cols = vars(kmeanstype)) +
  labs( y = 'Count',
        fill = 'Cluster') +
  theme_minimal()

histogramdifferences_fig03 = histograms / timeseries  + plot_annotation(tag_levels = 'a', tag_prefix = '(',tag_suffix = ')')
ggsave(file = '~/fig03_histogramdifferences.png' , plot = histogramdifferences_fig03 ,dpi = 300, 
       width = 25, height = 12.5, units = 'cm')


### change axis to coordiantes by everyt10
map = cluster_assignment %>%
  filter(Day == 153, !is.na(`Cluster ID`)) %>%
  ggplot(aes(x = Long, y = Lat, fill = as_factor(`Cluster ID`))) +
  geom_raster() +
  scale_fill_viridis(discrete = TRUE) +
  scale_x_continuous(name = 'Longitude',
                     breaks = seq(-125,-65,5),
                     labels = ~ paste0(.x*-1, "\u00B0", " W")) +
  scale_y_continuous(name = 'Latitude', 
                     breaks = seq(20,60,5),
                     expand = c(0,0),
                     labels = ~ paste0(.x, "\u00B0", ' N')) + 
  coord_fixed() + 
  labs(fill = 'Cluster ID') + 
  #annotate('rect', xmin = 14.5, xmax = 15.5, ymin = 9.5, ymax = 10.5, color = 'red', fill = NA) + 
  theme_minimal() + 
  theme(legend.position = 'bottom',
        title = element_text(size = 10),
        axis.text.x = element_text(angle = 90))

assignment_timeline <- cluster_assignment %>%
  filter(X == 15, Y == 10) %>%
  ggplot(aes(x = Day, y = `Cluster ID`)) +
  scale_x_continuous(breaks = seq(0,365,50)) +
  scale_color_viridis(discrete = TRUE) +
  geom_line(alpha = .5) +
  geom_point(aes(color = as.factor(`Cluster ID`))) +
  facet_wrap(~kmeanstype,ncol = 1) +
  theme_minimal()+
  labs(color = 'Cluster ID') + 
  theme(panel.grid.minor = element_line(color = NULL),
        legend.position = 'bottom')

assignmenttimeline_fig02 = map   + assignment_timeline  + plot_annotation(tag_levels = 'a', tag_prefix = '(',tag_suffix = ')')  
ggsave(file = '~/fig02_assignmenttimeline.png' , plot = assignmenttimeline_fig02 ,dpi = 300, 
       width = 25, height = 12.5, units = 'cm')





summer_autocorr_2003 = t(autocorr_r[9:12,])
summer_autocorr_2003 <- as_tibble(summer_autocorr_2003) %>%
  rename('Autocorrelation at Lag 1' = V1,
         'Autocorrelation at Lag 7' = V2,
         'Autocorrelation at Lag 14' = V3,
         'Autocorrelation at Lag 21' = V4) %>%
  mutate(lats = lon_lats$lat,
         long = lon_lats$lon) %>%
  pivot_longer(cols = contains('Lag'), names_to = 'Lag Time', values_to='Autocorrelation') 


summer_autocorr_2003$LagTimeF= factor(summer_autocorr_2003$`Lag Time`,
                                      levels=c('Autocorrelation at Lag 1',
                                               'Autocorrelation at Lag 7',
                                               'Autocorrelation at Lag 14',
                                               'Autocorrelation at Lag 21'))
autocorrelation_fig05 = summer_autocorr_2003 %>%
  ggplot() +
  geom_tile(aes(x = long, y = lats, fill = Autocorrelation)) +
  theme_minimal() +
  facet_wrap(~LagTimeF) +
  scale_x_continuous(name = 'Longitude',
                     labels = ~paste0(.x,"\u00B0", " W"),
                     expand = c(0,0)) +
  scale_y_continuous(name = 'Latitude',
                     breaks = seq(20,60,5),
                     labels = ~paste0(.x,"\u00B0"," N"),
                     transform = 'reverse',
                     expand = c(0,0)) +
  scale_fill_continuous(type = 'viridis',
                        na.value = NA) +
  
  coord_fixed()

autocorrelation_fig05
ggsave(file = 'fig05_autocorrelation.png' , plot = autocorrelation_fig05,dpi = 300, 
       width = 25, height = 12.5, units = 'cm')

## Drought Progression ####
propin4 = cluster_assignment %>%
  filter(Day >= 213 & Day <= 245, kmeanstype == 'Space-time') %>%
  group_by(Long,Lat, `Cluster ID`) %>%
  tally() %>%
  mutate(prop_days_in_cluster = n/(245-213)) %>%
  filter(`Cluster ID`  == 4, .preserve = TRUE)

base_map = cluster_assignment %>% filter( Day == 1, kmeanstype == 'Space-time', !is.na(`Cluster ID`))%>%
  ggplot(aes(x = Long, y = Lat)) +
  geom_raster(fill = 'black')


drought_Area = cluster_assignment %>%
  filter(X <= 31, X>= 17, Y<= 21, Y >=7) 

cluster4_assignment = base_map + 
  geom_raster(data = propin4, aes(x = Long, y = Lat, fill = prop_days_in_cluster)) +
  scale_fill_viridis(
    name = 'Proportion of time \n assigned to driest cluster',
    option = "magma") +
  scale_x_continuous(name = 'Longitude',
                     breaks = seq(-125,-65,5),
                     labels = ~ paste0(.x*-1, "\u00B0", " W")) +
  scale_y_continuous(name = 'Latitude', 
                     breaks = seq(20,60,5),
                     expand = c(0,0),
                     labels = ~ paste0(.x, "\u00B0", ' N')) + 
  coord_fixed() + 
  annotate('rect', xmin = min(drought_Area$Long), xmax = max(drought_Area$Long), 
           ymin = min(drought_Area$Lat), ymax = max(drought_Area$Lat), color = 'red', fill = NA, size = 1)  +
  theme_minimal() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 90))   


progression_fig06 = cluster_assignment %>%
  filter(Day %in% c(185, 200, 215,230,245, 260), !is.na(`Cluster ID`)) %>%
  ggplot(aes(x = Long, y = Lat, fill = as_factor(`Cluster ID`))) +
  geom_raster() +
  scale_fill_viridis(discrete = TRUE) +
  scale_x_continuous(name = 'Longitude',
                     breaks = seq(-125,-65,25),
                     labels = ~ paste0(.x*-1, "\u00B0", " W")) +
  scale_y_continuous(name = 'Latitude', 
                     breaks = seq(20,60,10),
                     expand = c(0,0),
                     labels = ~ paste0(.x, "\u00B0", ' N')) + 
  coord_fixed() + 
  labs(fill = 'Cluster ID',
       x = NULL,
       y = NULL) + 
  annotate('rect', xmin = min(drought_Area$Long), xmax = max(drought_Area$Long), 
           ymin = min(drought_Area$Lat), ymax = max(drought_Area$Lat), color = 'red', fill = NA, size = 1)  +
  facet_wrap(~Day, 
             nrow = 2,
             labeller = labeller(Day = c("185" = "4 July 2003",
                                         "200" = "19 July 2003",
                                         "215" = "3 August 2003",
                                         "230" = "18 August 2003",
                                         "245" = "2 September 2003",
                                         "260" = "17 September 2003 "))) + 
  theme_minimal() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 0))   +
  guides(colour = guide_legend(nrow = 1))

ggsave(file = '~/fig06_progression.png' , plot = progression_fig06,dpi = 300, 
       width = 25, height = 12.5, units = 'cm')




# Save the print(nc) dump to a text file



###### GP PLOTS #####
august_svdi <- svdi.array[ , ,213:245]
august_svdi_avg <- apply(august_svdi,c(1,2), mean)


r_august <- raster(t(august_svdi_avg), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
                   crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
r_august <- flip(r_august, direction ='y')

august_avg = gplot(r_august) +
  geom_tile(aes(fill=value)) +
  coord_fixed() +
  scale_fill_viridis(option = 'turbo',
                     name = 'SVDI',
                     limits = c(-1.5,2.5),
                     breaks = seq(-1.5,2.5,1)) +
  scale_x_continuous(name = 'Longitude',
                     breaks = seq(-125,-65,5),
                     labels = ~ paste0(.x, "\u00B0")
  ) +
  scale_y_continuous(name = 'Latitude', 
                     breaks = seq(20,60,5),
                     expand = c(0,0),
                     labels = ~ paste0(.x, "\u00B0")) + 
  geom_rect(xmin =-102, xmax = -87,
            ymin = 35, ymax = 50, fill = NA, color ='black') +
  theme(
    panel.background = element_rect(fill = 'gray50'),
    legend.position = 'bottom',
    panel.grid = element_line(colour = 'gray50'),
    axis.ticks = element_line(colour = 'gray50'))




august_avg+cluster4_assignment + plot_annotation(tag_levels = 'A')  





##### Publication Plots #######
#Create a custom color scale
myFill <- brewer.pal(9,"RdBu")
myFill[5] = "#FFFFFF"
names(myFill) <- levels(svdi_clusters$SVDI_severity)
fillScale <- scale_fill_manual(name = "SVDI Severity",values = myFill)



##### Create base map ####
base_map = 
  NLDAS_landcoordinates %>%
  ggplot(aes(x = long, y = lat)) +
  coord_equal() +
  geom_tile(fill = 'gray70') +
  scale_x_continuous(name = 'Longitude',
                     breaks = seq(-125,-65,5),
                     labels = ~ paste0(.x*-1, "\u00B0", " W")) +
  scale_y_continuous(name = 'Latitude', 
                     breaks = seq(20,60,5),
                     expand = c(0,0),
                     labels = ~ paste0(.x, "\u00B0", ' N')) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

### Drought Centroid Tracking #####

centroids_5 <- svdi_clusters %>%
  filter(year == 2003, drought_severity == 'moderate', cluster_id == 5) %>%
  group_by(pentad_num,cluster_id) %>%
  reframe( lat = mean(lat),
           long = mean(long)) 


centroids_6 <- svdi_clusters%>%
  filter(year == 2003, drought_severity == 'moderate', cluster_id == 6) %>%
  group_by(pentad_num,cluster_id) %>%
  reframe( lat = mean(lat),
           long = mean(long))

small_date_5 = centroids_5 %>% filter(row_number() %% 3 == 1)
small_date_6 = centroids_6 %>% filter(row_number() %% 5 == 1)
combined_small_dates = rbind(small_date_5,small_date_6)
centroid_6_small =  centroids_6 %>% filter(pentad_num %in% seq(24,32))
combined_centroids = rbind(centroids_5,centroid_6_small)
endpoints_5 = rbind(centroids_5[2:9,], centroids_5[9,]) %>% rename(lat_end = lat, long_end = long)
endpoints_6 = rbind(centroid_6_small[2:7,],centroid_6_small[7,]) %>% rename(lat_end = lat, long_end = long)
combined_endpoints = rbind(endpoints_5,endpoints_6)

centroidtracking_fig09 <- base_map + 
  geom_segment(data = combined_centroids, 
               aes(x = long,
                   y = lat,
                   xend = combined_endpoints$long_end,
                   yend = combined_endpoints$lat_end,
                   group = cluster_id), 
               color = 'white', 
               arrow = arrow(length=unit(.2,"cm"),
                             type = 'closed')) +
  geom_point(data = combined_centroids, aes(x = long, y = lat, group = cluster_id, color = cluster_id)) +
  geom_text_repel(data=combined_small_dates %>% filter(pentad_num %in% seq(24,32)), aes(label =pentad_num),min.segment.length = .05, color = 'black')+
  scale_color_manual(values = viridis::viridis(n = 2),
                     labels = c("Moderate Cluster 5", 'Moderate Cluster 6')) +
  theme_minimal() + 
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90))

ggsave("fig09_centroidtracking.png", centroidtracking_fig09, dpi = 300,
       width = 25, height = 12.5, units = 'cm')

#### Average Drought Duration ####
svdi_clusters_w_area %>%
  group_by(year,drought_severity,cluster_id) %>%
  summarise(min_day = min(day),
            max_day = max(day)) %>%
  mutate(drought_duration = max_day - min_day) %>%
  ungroup() %>%
  group_by(drought_severity,year) %>%
  summarise(avg_drought_duration = mean(drought_duration)) %>% 
  ggplot(aes(x = year, y = avg_drought_duration,color = drought_severity))  +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0,400), breaks = seq(0,365,100), expand = c(0,0)) +
  scale_x_continuous(limits = c(1979,2022), breaks = seq(1980,2020,5), expand = c(0,0)) +
  scale_color_manual(values = viridis::viridis(n = 4),
                     labels = c("Mild Drought", 'Moderate Drought','Severe Drought', 'Extreme Drought')) +
  facet_wrap(~drought_severity,nrow = 4) +  
  labs(y = 'Average Drought Duration (Days)',
       color = 'Drought Severity',
       x = 'Year') +
  theme_minimal() +
  theme(axis.ticks = element_blank(),
        legend.title = element_blank(),
        panel.grid.minor.x= element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.position = 'bottom',
        strip.background = element_blank(),
        strip.text.x = element_blank()) 



#### Mean Intensity #####
mean_intensity = svdi_clusters_w_area %>%
  group_by(year,drought_severity,cluster_id) %>%
  summarise(intensity = mean(SVDI)) %>%
  ungroup() %>%
  group_by(drought_severity,year) %>%
  summarise(avg_drought_intensity = mean(intensity)) %>% 
  ggplot(aes(x = year, y = avg_drought_intensity, color = drought_severity)) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(0,3), breaks = seq(0,3,.5), expand = c(0,0)) +
  scale_x_continuous(limits = c(1979,2022), breaks = seq(1980,2020,5), expand = c(0,0)) + 
  scale_color_manual(values = viridis::viridis(n = 4),
                     labels = c("Mild Drought", 'Moderate Drought','Severe Drought', 'Extreme Drought')) +
  labs(y = 'Mean SVDI',
       x = 'Year') +
  theme_minimal() +
  theme(axis.ticks = element_blank(),
        legend.title = element_blank(),
        panel.grid.minor.x= element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.position = 'bottom') 


# Example of stacked Histogram Plot ##### 
svdi_clusters_w_area  %>%
  filter(year == 2003, drought_severity == 'moderate') %>%
  group_by(year,pentad_num,SVDI_severity, cluster_id) %>%
  summarise(n = n())   %>%
  ggplot(aes(x = pentad_num, y = n, fill = SVDI_severity, color = cluster_id)) + 
  geom_col( position = 'stack') + 
  geom_vline(xintercept = 43, lty = 2) +
  geom_vline(xintercept = 49, lty = 2)+
  annotate(geom = 'text', x = 15, y = 1000, label = 'Black lines delineate the timing \n of 2003 flash drought') +
  scale_y_continuous(expand = c(0,0),  limits = c(0,1300), name = 'Pixel Count') +
  scale_x_continuous(expand = c(0,0), name = 'Pentad', breaks = seq(0,75,5)) +
  scale_color_manual(name = 'Cluster ID', values = viridis::viridis(3)) +
  #fillScale +
  scale_fill_manual(name = ' SVDI Severity', values = viridis::viridis(7, option = 'magma', direction = 1)) + 
  theme_minimal() +
  #scale_color_manual(values = cluster_colors) + 
  theme(legend.position = 'bottom',
        axis.ticks.y= element_blank())


#### Drought Spatial Extent Snapshots #####

drought_2003 = svdi_clusters_w_area %>%
  filter(year == 2003,
         drought_severity == 'moderate') %>%
  mutate(week = lubridate::week(date),
         start_of_week = format(lubridate::ymd(lubridate::floor_date(date, unit="week")),"%d %B %Y")) 


#format(lubridate::ymd(unique(drought_2003$start_of_week)), "%d %B %Y")
snapshots_fig11 <- base_map + geom_tile(data = drought_2003,
                                        aes(x = long, y = lat, fill = cluster_id)) +
  #theme(legend.position = "none") +
  facet_wrap(vars(start_of_week)) +
  scale_x_continuous(name = ' Longitude',
                     labels = ~paste0(.x,'\u00B0'," W")) +
  scale_y_continuous(name = ' Latitude',
                     labels = ~paste0(.x,'\u00B0'," N")) +
  scale_fill_viridis(discrete = TRUE,
                     labels=c("Moderate Cluster 0 ",
                              "Moderate Cluster 5 ",
                              "Moderate Cluster 6")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggsave(file = 'fig11_snapshotsv1.png',plot = snapshots_fig11, dpi = 300,
       width = 27, height = 25, units = 'cm')

drought_2000 = svdi_clusters_w_area %>%
  filter(year == 2000,
         cluster_id == 11, drought_severity == 'severe') %>%
  mutate(week = lubridate::week(date),
         start_of_week = format(lubridate::ymd(lubridate::floor_date(date, unit="week")),"%d %B %Y")) 


original_map <- base_map + geom_tile(data = drought_2000,
                                     aes(x = long, y = lat, fill = cluster_id)) +
  scale_fill_manual(
                     name = NULL,
                     values = viridis::viridis(1,direction = -1) + 
                     labels=c("Severe Cluster 11 ")) +
  facet_wrap(vars(start_of_week)) +
  theme_minimal() + 
  theme(legend.position = 'bottom')




tuned_drought_2000 = tuned2000 %>%
  filter(year == 2000,
         cluster_id == 17, drought_severity == 'extreme') %>%
  mutate(week = lubridate::week(date),
         start_of_week = format(lubridate::ymd(lubridate::floor_date(date, unit="week")),"%d %B %Y")) 


tuned_map = base_map + geom_tile(data = tuned_drought_2000,
                                 aes(x = long, y = lat, fill =as.factor( cluster_id))) +
  scale_fill_manual(name = NULL,
                    values = viridis:viridis(1),
                    labels = 'Extreme Cluster 17') +
  facet_wrap(vars(start_of_week)) +
  theme(legend.position = 'bottom')


ggsave(file = 'fig13a_tunedmap.png',plot = tuned_map, dpi = 300,
       width = 27, height = 25, units = 'cm')

ggsave(file = 'fig13b_originalmap.png',plot = original_map, dpi = 300,
       width = 27, height = 25, units = 'cm')




test = read_csv('../bilevelclustering/data/DBSCAN_output_2003_full.csv') %>%
  mutate(cluster_id = as.factor(cluster_id))

moderate_drought_cluster <- plot_ly(data = test %>% filter(drought_severity =='moderate'), 
                                    type = 'scatter3d', 
                                    mode = "markers",
                                    x =  ~long , 
                                    z = ~lat, 
                                    y = ~day,
                                    #colors =  c('#7E2954','#337538','#DCCD7D'),
                                    colorscale='Viridis',
                                    color = ~cluster_id,
                                    opacity = .85,
                                    size = 15,
                                    scene ='scene1') %>%
  layout(legend=list(title=list(text='Cluster ID')),
         scene = list(xaxis=list(title = 'Longitude (\u00B0 W)',
                                 autorange = "reversed"),
                      yaxis=list(title = 'Day of Year'),
                      zaxis=list(title = 'Latitude (\u00B0 N)')))
moderate_drought_cluster


mild_drought_cluster <- plot_ly(data = test %>% filter(drought_severity =='mild'), 
                                type = 'scatter3d', 
                                mode = "markers",
                                x =  ~long , 
                                z = ~lat, 
                                y = ~day,
                                # colors =  c('#9F4A96'),
                                colorscale='Viridis',
                                color = ~cluster_id,
                                opacity = .15,
                                size = 20, 
                                scene ='scene2') %>%
  layout(legend=list(title=list(text='Cluster ID')),
         scene2 = list(xaxis=list(title = 'Longitude (\u00B0 W)',
                                  autorange = "reversed"),
                       yaxis=list(title = 'Day'),
                       zaxis=list(title = 'Latitude (\u00B0 N)')))
mild_drought_cluster

plotly::subplot(moderate_drought_cluster, mild_drought_cluster) %>%
  layout(legend=list(title=list(text='Cluster ID')),
         scene = list(xaxis=list(title = 'Longitude (\u00B0 W)',
                                 autorange = "reversed"),
                      yaxis=list(title = 'Day of Year'),
                      zaxis=list(title = 'Latitude (\u00B0 N)')),
         scene2 = list(xaxis=list(title = 'Longitude (\u00B0 W)',
                                  autorange = "reversed"),
                       yaxis=list(title = 'Day of Year'),
                       zaxis=list(title = 'Latitude (\u00B0 N)')))