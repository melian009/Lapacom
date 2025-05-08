rm(list = ls())

#Importación y limpieza de datos
library(readxl)
library(janitor)
library(tidyverse)

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

#️ Datos espaciales
library(sf)
library(ggspatial)
library(ggmap)
library(geodata)  # Nuevo: reemplaza a raster::getData
library(ggpubr)

# Auxiliares
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
  filter(species %in% c("Patella ordinaria", "Patella aspera")) %>% 
  filter(sampling_site %in% c("",))

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
       x = "Longitude", y = "Latitude", shape = "Sampling Site") +
  theme(legend.position = "bottom")

# ───────────────────────────────────────────────
# Group Analysis MPA (Control) vs Full access (Exploited)
summary(Mad_2sp_site)

# Función para análisis por especie
analisis_por_regimen <- function(data, especie) {
  cat("\n Análisis para:", especie, "\n")
  
  df <- data %>% filter(species == especie) %>%
    filter(!is.na(total_length_mm) & !is.na(protective_regime))
  
  n_groups <- nlevels(df$protective_regime)
  
  # Boxplot
  p <- ggplot(df, aes(x = protective_regime, y = total_length_mm, fill = protective_regime)) +
    geom_boxplot() +
    labs(title = paste("Total Length by Protective Regime:", especie),
         x = "Protective Regime", y = "Total Length (mm)") +
    theme_minimal()
  print(p)
  
  # Levene Test
  levene <- car::leveneTest(total_length_mm ~ protective_regime, data = df)
  print(levene)
  
  # Inicializar
  resumen <- tibble(
    Especie = especie,
    Test = NA,
    P_value = NA,
    Post_Hoc = NA,
    N_Grupos = n_groups,
    Efecto = NA,
    Magnitud = NA
  )
  modelo <- NULL
  post_hoc <- NULL
  
  if (levene$`Pr(>F)`[1] > 0.05) {
    # ANOVA
    cat("\n Varianzas homogéneas → ANOVA\n")
    modelo <- aov(total_length_mm ~ protective_regime, data = df)
    res_anova <- summary(modelo)
    print(res_anova)
    
    resumen$Test <- "ANOVA"
    resumen$P_value <- res_anova[[1]]$`Pr(>F)`[1]
    
    if (resumen$P_value < 0.05) {
      if (n_groups > 2) {
        cat("Post Hoc (Tukey):\n")
        post_hoc <- TukeyHSD(modelo)
        print(post_hoc)
        resumen$Post_Hoc <- "Tukey"
      } else {
        # Tamaño del efecto con Cohen's d
        d <- effsize::cohen.d(total_length_mm ~ protective_regime, data = df)
        resumen$Post_Hoc <- "No necesario"
        resumen$Efecto <- round(d$estimate, 3)
        resumen$Magnitud <- d$magnitude
        print(d)
      }
    }
  } else {
    # Kruskal-Wallis
    cat("\n Varianzas NO homogéneas → Kruskal-Wallis\n")
    modelo <- kruskal.test(total_length_mm ~ protective_regime, data = df)
    print(modelo)
    
    resumen$Test <- "Kruskal-Wallis"
    resumen$P_value <- modelo$p.value
    
    if (modelo$p.value < 0.05) {
      if (n_groups > 2) {
        cat("Post Hoc (Conover-Iman):\n")
        post_hoc <- conover.test::conover.test(df$total_length_mm, df$protective_regime, method = "holm")
        print(post_hoc)
        resumen$Post_Hoc <- "Conover"
      } else {
        # Tamaño del efecto con Cliff's delta
        d <- effsize::cliff.delta(total_length_mm ~ protective_regime, data = df)
        resumen$Post_Hoc <- "No necesario"
        resumen$Efecto <- round(d$estimate, 3)
        resumen$Magnitud <- d$magnitude
        print(d)
      }
    }
  }
  
  return(list(
    resumen = resumen,
    modelo = modelo,
    post_hoc = post_hoc
  ))
}


res_asp <- analisis_por_regimen(Mad_2sp_site, "Patella aspera")
res_ord <- analisis_por_regimen(Mad_2sp_site, "Patella ordinaria")

# Tabla combinada
tabla_resumen <- bind_rows(res_asp$resumen, res_ord$resumen)
print(tabla_resumen)

# Usar modelos posteriormente:
summary(res_asp$modelo)  # 



