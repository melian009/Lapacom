rm(list = ls())

#Importación y limpieza de datos
library(readxl)
library(janitor)
library(tidyverse)
library(tidyr)

#Estadística descriptiva y análisis
library(moments)
library(DescTools)
library(descr)
library(effsize)
library(vcd)
library(tibble)

# Modelos estadísticos
library(car)
library(multcomp)
library(lmtest)
library(conover.test)
library(pROC)
library(moderndive)

#️ Visualización avanzada
library(ggthemes)
library(hrbrthemes)
library(viridis)
library(ggridges)
library(gridExtra)
library(flextable)
library(cowplot)

#️ Datos espaciales
library(sf)
library(ggspatial)
library(ggmap)
library(geodata)
library(ggpubr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(grid) 

# Auxiliares
library(stringr)

# ───────────────────────────────────────────────
# DATOS

Mad_RAW <- read_excel("E:/Investigación/Repositorios/Lapacom/Data/ToAnalyze/Madeira/BD_LIMPETS_MAD_1996-2018.xlsx",
                      sheet = "Data", range = "A1:T69707")
Mad_RAW_Clean <- Mad_RAW %>% clean_names()

Mad_selected <- Mad_RAW_Clean %>%
  dplyr::select(species, year, month, total_length_mm, total_length_class_mm,
                weight_g, mature_imature, sampling_site, lat, long,
                protective_regime, proximity_human_settlements,
                accessibility, age_lt, age_months, age_class)
# Filtrar por especies objetivo
Mad_2sp <- Mad_selected %>%
  filter(species %in% c("Patella ordinaria", "Patella aspera"))

Mad_2sp_maps <- Mad_2sp %>%
  filter(!is.na(sampling_site) & 
           trimws(sampling_site) != "" & 
           sampling_site %in% c("Porto Moniz", "Paúl do Mar", "Funchal", "Desertas", 
                                "Caniçal", "Santa Cruz", "Ribeira Brava", "São Vicente")) %>% 
  mutate(across(where(is.character), as.factor))


# ───────────────────────────────────────────────
# COORDENADAS Y MAPAs
# ───────────────────────────────────────────────
# MAPAS
# 1. Mapa regional
world <- ne_countries(scale = 10, returnclass = "sf")
portugal_gadm <- geodata::gadm("PRT", level = 1, path = tempdir()) %>% st_as_sf()
portugal_mainland <- portugal_gadm %>% filter(!NAME_1 %in% c("Madeira", "Açores"))
madeira_sf <- portugal_gadm %>% filter(NAME_1 == "Madeira")
azores_sf <- portugal_gadm %>% filter(NAME_1 == "Açores")

etiquetas <- tibble::tibble(
  lugar = c("Portugal", "Selvagens", "Azores", "Canary Islands", "Morocco", "Gibraltar", "Spain"),
  lon = c(-10.5, -16.1, -26.0, -16.5, -8.5, -5.4, -6.5),
  lat = c(39.5, 31.0, 38.6, 29.2, 31.0, 36.1, 43.0)
)

map_region <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "darkgrey") +
  geom_sf(data = portugal_mainland, fill = "darkgrey") +
  geom_sf(data = madeira_sf, fill = "darkgrey") +
  geom_sf(data = azores_sf, fill = "darkgrey") +
  geom_text(data = etiquetas, aes(x = lon, y = lat, label = lugar),
            size = 3.2, fontface = "italic") +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotate("rect", xmin = -17.6, xmax = -16.2, ymin = 32.4, ymax = 33.2,
           fill = "grey", alpha = 0.3, color = "grey20", size = 0.7) +
  annotate("text", x = -20, y = 36, label = "Atlantic Ocean", angle = 45, size = 7,
           fontface = "italic", color = "gray30") +
  annotate("text", x = -17.0, y = 33.5, label = "Madeira", size = 3.2, fontface = "bold",
           color = "black") +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  coord_sf(xlim = c(-30, -5), ylim = c(27, 44), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "darkgrey", size = 0.2),
    axis.title = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

# 2. Coordenadas en decimal
convert_dms_to_decimal <- function(dms) {
  parts <- str_match(dms, "(\\d+)°(\\d+)'(\\d+\\.?\\d*)\"?([NSEW])")
  deg <- as.numeric(parts[, 2])
  min <- as.numeric(parts[, 3])
  sec <- as.numeric(parts[, 4])
  dir <- parts[, 5]
  decimal <- deg + min / 60 + sec / 3600
  ifelse(dir %in% c("S", "W"), -decimal, decimal)
}

coords <- Mad_2sp %>%
  dplyr::select(sampling_site, lat, long) %>%
  distinct() %>%
  drop_na() %>%
  mutate(
    lat_dd = convert_dms_to_decimal(lat),
    long_dd = convert_dms_to_decimal(long)
  )

coords_sf <- st_as_sf(coords, coords = c("long_dd", "lat_dd"), crs = 4326)

# 3. Mapa Madeira
map_madeira <- ggplot() +
  geom_sf(data = madeira_sf, fill = "gray90", color = "darkgrey") +
  geom_sf(data = coords_sf, aes(color = sampling_site), size = 3) +
  coord_sf(xlim = c(-17.6, -16.2), ylim = c(32.4, 33.2), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  labs(x = "Longitude", y = "Latitude", color = "Sampling Site:") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.box = "horizontal",
    plot.margin = margin(5, 5, 5, 5)
  )

# 4. Extraer leyenda como grob
legend_b <- ggplotGrob(map_madeira) %>%
  gtable::gtable_filter("guide-box")

# 5. Quitar leyenda del mapa B
map_madeira_noleg <- map_madeira + theme(legend.position = "none")

# 6. Combinar mapas A y B
panel_maps <- plot_grid(
  map_region, map_madeira_noleg,
  labels = c("A)", "B)"),
  label_size = 14,
  ncol = 2,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1, 1.2)
)

# 7. Añadir leyenda debajo, centrada
final_plot <- plot_grid(
  panel_maps,
  legend_b,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

# Mostrar resultado final
print(final_plot)

