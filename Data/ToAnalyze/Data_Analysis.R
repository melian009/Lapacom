rm(list = ls())

# 📥 Importación y limpieza de datos
library(readxl)
library(janitor)
library(tidyverse)

# 📊 Estadística descriptiva y análisis
library(moments)
library(DescTools)
library(descr)
library(effsize)
library(vcd)

# 📈 Modelos estadísticos
library(car)
library(multcomp)
library(lmtest)
library(conover.test)
library(pROC)
library(moderndive)

# 🖼️ Visualización avanzada
library(ggthemes)
library(hrbrthemes)
library(viridis)
library(ggridges)
library(gridExtra)
library(flextable)

# 🗺️ Datos espaciales
library(sf)
library(ggspatial)
library(ggmap)
library(geodata)  # Nuevo: reemplaza a raster::getData

# 🧰 Auxiliares
library(stringr)

# ───────────────────────────────────────────────
# DATOS

Mad_RAW <- read_excel("Madeira/BD_LIMPETS_MAD_1996-2018.xlsx",
                      sheet = "Data", range = "A1:T69707")
Mad_RAW_Clean <- Mad_RAW %>% clean_names()

Mad_selected <- Mad_RAW_Clean %>%
  dplyr::select(species, year, month, total_length_mm, total_length_class_mm,
                weight_g, mature_imature, sampling_site, lat, long,
                protective_regime, proximity_human_settlements,
                accessibility, age_lt, age_months, age_class) %>%
  mutate(across(where(is.character), as.factor))

# Filtrar por especies objetivo
Mad_2sp <- Mad_selected %>%
  filter(species %in% c("Patella ordinaria", "Patella aspera"))

Mad_2sp_site <- Mad_2sp %>%
  filter(!is.na(sampling_site) & trimws(sampling_site) != "")


# ───────────────────────────────────────────────
# COORDENADAS Y MAPA

# 1. Función para convertir DMS a decimal
convert_dms_to_decimal <- function(dms) {
  parts <- str_match(dms, "(\\d+)°(\\d+)'(\\d+\\.?\\d*)\"?([NSEW])")
  deg <- as.numeric(parts[, 2])
  min <- as.numeric(parts[, 3])
  sec <- as.numeric(parts[, 4])
  dir <- parts[, 5]
  decimal <- deg + min / 60 + sec / 3600
  ifelse(dir %in% c("S", "W"), -decimal, decimal)
}

# 2. Crear data frame con coordenadas únicas y convertir a decimal
coords <- Mad_2sp_site %>%
  dplyr::select(sampling_site, lat, long) %>%
  distinct() %>%
  drop_na(lat, long) %>%
  mutate(
    lat_dd = convert_dms_to_decimal(lat),
    long_dd = convert_dms_to_decimal(long)
  )

# 3. Convertir a objeto sf
coords_sf <- st_as_sf(coords, coords = c("long_dd", "lat_dd"), crs = 4326)

# 4. Obtener mapa detallado de Portugal con GADM desde geodata
portugal_gadm <- geodata::gadm(country = "PRT", level = 1, path = tempdir())
portugal_sf <- st_as_sf(portugal_gadm)

# Filtrar región de Madeira
madeira_sf <- portugal_sf %>% filter(NAME_1 == "Madeira")

# 5. Graficar con mapa detallado
ggplot() +
  geom_sf(data = madeira_sf, fill = "gray90", color = "gray60") +
  geom_sf(data = coords_sf, aes(color = sampling_site), size = 3) +
  coord_sf(xlim = c(-17.6, -16.2), ylim = c(32.4, 33.2), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  theme_minimal() +
  labs(title = "Sampling Sites in Madeira Archipelago",
       x = "Longitude", y = "Latitude", color = "Sampling Site") +
  theme(legend.position = "bottom")
