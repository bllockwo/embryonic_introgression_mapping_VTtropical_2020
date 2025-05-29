## Precipitation
library(nasapower)
library(tidyverse)
library(magrittr)

#### Weather
VT_weather <- get_power(
  community = "ag",
    lonlat = c(-72.47167, 44.36272),
    pars = c("RH2M", "T2M", "PRECTOTCORR"),
    dates = c("2015-01-01", 
              paste(2025, "-01-01", sep = "")),
    temporal_api = "hourly"
)

VA_weather <-
  get_power(
    community = "ag",
    lonlat = c(-78.4897, 37.979),
    pars = c("RH2M", "T2M", "PRECTOTCORR"),
    dates = c("2015-01-01", 
              paste(2025, "-01-01", sep = "")),
    temporal_api = "hourly"
  )

SK_weather <-
  get_power(
    community = "ag",
    lonlat = c(-62.80985, 17.37319),
    pars = c("RH2M", "T2M", "PRECTOTCORR"),
    dates = c("2015-01-01", 
              paste(2025, "-01-01", sep = "")),
    temporal_api = "hourly"
  )

VT_weather <- get_power(
  community = "ag",
  lonlat = c(-72.47167, 44.36272),
  pars = c("RH2M", "T2M", "PRECTOTCORR"),
  dates = c("2015-01-01", 
            paste(2025, "-01-01", sep = "")),
  temporal_api = "hourly"
)

VT_weather %<>% mutate(pop = "VT")
SK_weather %<>% mutate(pop = "SK")

cor.test(
  filter(SK_weather, MO >= 1 & MO <= 12)$PRECTOTCORR,
  filter(SK_weather, MO >= 1 & MO <= 12)$T2M )

cor.test(
  filter(VT_weather, MO >= 1 & MO <= 12)$PRECTOTCORR,
  filter(VT_weather, MO >= 1 & MO <= 12)$T2M )


rbind(VT_weather, SK_weather) %>%
  filter(MO >= 7 & MO <= 12) %>%
  group_by(MO, pop) %>%
  summarize(MaxPrecip = sum(PRECTOTCORR > 10)) %>%
  ggplot(aes(
    x=MO, y=MaxPrecip, color = pop
  )) + geom_point() -> 
  precipt


ggsave(precipt, file = "precipt.pdf")

