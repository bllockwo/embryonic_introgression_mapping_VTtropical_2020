### Map and ecology figure.
library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)
library(reshape2)
library(rnaturalearth)
library(rnaturalearthdata)
library(nasapower)
library(gmodels)

### Plot map
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

parents = 
data.frame(
  pop= c("VT","SK", "VA"),
  long= c(-72.471669,-62.809850,-78.4897),
  lat= c(44.362721,17.373190, 37.979)
    )

ggplot(data = world) +
  geom_sf(fill= "antiquewhite", alpha = 0.8) +
  coord_sf(xlim = c(-60, -83), ylim = c(13.00, 47.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue")) +
  geom_point(data = parents,
             aes(x=long, y=lat, fill = pop),size = 5, shape = 21) +
  scale_fill_manual(values = c("firebrick", "springgreen", 
                               "steelblue")) -> map.plot

ggsave(map.plot, 
       file = "/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/20.WEATHER_PLOT/map.plot.pdf",
       w= 4, h=6)

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
	
VA_weather <-
  get_power(
    community = "ag",
    lonlat = c(-78.4897, 37.979),
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
VA_weather %<>% mutate(pop = "VA")

rbind(VT_weather, SK_weather) -> sites.weather
rbind(sites.weather, VA_weather) -> sites.weather
save(sites.weather, file = "/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/12.PCA_W_ALLDAT_Fig1/data_plots/sites.weather.nasa.Rdata")

load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/12.PCA_W_ALLDAT_Fig1/data_plots/sites.weather.nasa.Rdata")
sites.weather

sites.weather %>%
  melt(id = c("LON", "LAT", "YEAR", "MO", "DY", "HR", "pop")) %>%
  group_by(MO, pop, variable) %>%
  summarise(m.val = ci(value)[1],
            m.lci = ci(value)[2],
            m.uci = ci(value)[3],
            ) %>%
  ggplot(aes(
    x=as.numeric((MO), levels = 1:12),
    y=m.val,
    ymin=m.lci,
    ymax=m.uci,
    col=pop
  )) + #geom_boxplot() +
  geom_vline(xintercept = 6) +
  geom_vline(xintercept = 11) +
  geom_line() +
  geom_errorbar() +
  geom_point() +
  theme_classic() +
  scale_x_continuous(name="Month", limits=c(1, 12),
                     breaks=seq(1,12, by =2 )) +
  facet_wrap(~variable, 
             scales = "free_y") -> weather.trace

ggsave(weather.trace, 
       file = "/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/20.WEATHER_PLOT/weather.trace.pdf",
       w= 9, h=3)


sites.weather %>%
  melt(id = c("LON", "LAT", "YEAR", "MO", "DY", "HR", "pop")) %>%
  group_by(YEAR, MO, pop, variable) %>%
  summarise(m.val = ci(value)[1]) %>%
  filter(MO > 6 & MO < 11) %>%
  ggplot(aes(
    x=m.val,
    fill=pop
  )) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  facet_wrap(~variable, scales = "free_x") -> weather.dens

ggsave(weather.dens, 
       file = "/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/20.WEATHER_PLOT/weather.dens.pdf",
       w= 9, h=3)

sites.weather %>%
  melt(id = c("LON", "LAT", "YEAR", "MO", "DY", "HR", "pop")) %>%
  group_by(YEAR, MO, pop, variable) %>%
  filter(MO > 6 & MO < 11) %>%
  group_by(pop, variable) %>%
  summarise(m.val = ci(value)[1],
            sd.val = sd(value))
