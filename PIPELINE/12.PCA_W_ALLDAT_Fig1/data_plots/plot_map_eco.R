### Map and ecology figure.
library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)
library(reshape2)
library(rnaturalearth)
library(rnaturalearthdata)
library(nasapower)

### Plot map
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

parents = 
data.frame(
  pop= c("VT","SK"),
  long= c(-72.471669,-62.809850),
  lat= c(44.362721,17.373190)
    )

ggplot(data = world) +
  geom_sf(fill= "antiquewhite", alpha = 0.8) +
  coord_sf(xlim = c(-60, -83), ylim = c(13.00, 47.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue")) +
  geom_point(data = parents,
             aes(x=long, y=lat, fill = pop),size = 5, shape = 21) +
  scale_fill_manual(values = c("firebrick", "steelblue")) 

#### Weather
VT_weather <-
get_power(
  community = "ag",
  lonlat = c(-72.47167, 44.36272),
  pars = c("RH2M", "T2M", "PRECTOTCORR"),
  dates = c("2015-01-01", 
            paste(2025, "-01-01", sep = "")), # gets data from earliest occurence - summer 2024
  temporal_api = "hourly",
  time_standard = "UTC"
)

SK_weather <-
  get_power(
    community = "ag",
    lonlat = c(-62.80985, 17.37319),
    pars = c("RH2M", "T2M", "PRECTOTCORR"),
    dates = c("2015-01-01", 
              paste(2025, "-01-01", sep = "")), # gets data from earliest occurence - summer 2024
    temporal_api = "hourly",
    time_standard = "UTC"
  )

VT_weather %<>% mutate(pop = "VT")
SK_weather %<>% mutate(pop = "SK")

rbind(VT_weather, SK_weather) -> sites.weather
save(sites.weather, file = "/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/12.PCA_W_ALLDAT_Fig1/data_plots/sites.weather.nasa.Rdata")

load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/12.PCA_W_ALLDAT_Fig1/data_plots/sites.weather.nasa.Rdata")
sites.weather

sites.weather %>%
  melt(id = c("LON", "LAT", "YEAR", "MO", "DY", "HR", "pop")) %>%
  group_by(YEAR, MO, pop, variable) %>%
  summarise(m.val = mean(value)) %>%
  ggplot(aes(
    x=factor((MO), levels = 1:12),
    y=m.val,
    col=pop
  )) + geom_boxplot() +
  theme_bw() +
  facet_wrap(~variable, scales = "free_y")

