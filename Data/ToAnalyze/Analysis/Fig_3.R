# Data import and manipulation
library(readxl)
library(janitor)
library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)

#Descriptive statistics 
library(moments)
library(DescTools)
library(descr)
library(effsize)
library(vcd)
library(tibble)


# Data wrangling & plotting
library(ggplot2)
library(scales)
library(cowplot)
library(flextable)
library(ggthemes)
library(hrbrthemes)
library(viridis)
library(ggridges)
library(gridExtra)
library(flextable)
library(cowplot)

# Mapping plots
library(sf)
library(ggspatial)
library(ggmap)
library(geodata)
library(ggpubr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(grid) 

# Statistical modeling and tests
library(car)
library(broom)
library(emmeans)
library(conover.test)
library(MASS)
library(multcomp)
library(lmtest)
library(pROC)
library(moderndive)
library(stringr)

#Data

Mad_RAW <- read_excel("E:/Investigación/Repositorios/Lapacom/Data/ToAnalyze/Madeira/BD_LIMPETS_MAD_1996-2018.xlsx", sheet = "Data", range = "A1:T69707")
Mad_RAW_Clean <- Mad_RAW %>% clean_names()

Mad_selected <- Mad_RAW_Clean %>%
  dplyr::select(species, year, month, total_length_mm, total_length_class_mm,
                weight_g, mature_imature, sampling_site, lat, long,
                protective_regime, proximity_human_settlements,
                accessibility, age_lt, age_months, age_class)

Mad_2sp <- Mad_selected %>%
  filter(species %in% c("Patella ordinaria", "Patella aspera")) %>%
  filter(!is.na(sampling_site) & trimws(sampling_site) != "") %>%
  mutate(across(where(is.character), as.factor),
         regulation_period = if_else(year < 2007, "NR", "R") %>% as.factor())

#For mapping
Mad_2sp_maps <- Mad_2sp %>%
  filter(!is.na(sampling_site) & 
           trimws(sampling_site) != "" & 
           sampling_site %in% c("Porto Moniz", "Paúl do Mar", "Funchal", "Desertas", 
                                "Caniçal", "Santa Cruz", "Ribeira Brava", "São Vicente")) %>% 
  mutate(across(where(is.character), as.factor))

#For analisis
Mad_2sp_analisis <- Mad_2sp %>%
  filter(!is.na(sampling_site) & 
           trimws(sampling_site) != "" & 
           sampling_site %in% c("Porto Moniz", "Paúl do Mar", "Funchal", "Desertas", 
                                "Caniçal", "Santa Cruz", "Ribeira Brava", "São Vicente"))%>% 
  mutate(across(where(is.character), as.factor))

summary(Mad_2sp_analisis)


# Monthly proportion of mature individuals by protection regime
madurez_prop_mensual_pr <- Mad_2sp_analisis %>%
  filter(!is.na(mature_imature)) %>%
  group_by(year, month, species, sampling_site, accessibility, protective_regime) %>%
  summarise(
    total = n(),
    maduros = sum(mature_imature == "Mature"),
    proporcion_madura = maduros / total,
    .groups = "drop"
  ) %>%
  mutate(
    grupo = interaction(species, protective_regime, sep = " - "),
    fecha = lubridate::make_date(year, month, 1)  # For a continuous time axis
  )

# Visual style for protection regime plot
colores <- c(
  "Patella aspera - Full acess" = "#D39C00",
  "Patella ordinaria - Full acess" = "#0072B2",
  "Patella aspera - MPA" = "#D39C00",
  "Patella ordinaria - MPA" = "#0072B2"
)

lineas <- c(
  "Patella aspera - Full acess" = "solid",
  "Patella ordinaria - Full acess" = "solid",
  "Patella aspera - MPA" = "dotted",
  "Patella ordinaria - MPA" = "dotted"
)

# A plot: mature proportion by protection regime
MATURE_PROP_protreg <- ggplot(
  data = madurez_prop_mensual_pr, 
  aes(x = fecha,
      y = proporcion_madura,
      color = grupo,
      linetype = grupo)) +
  
  geom_point(aes(color = grupo)) +
  
  geom_line(linewidth = 0.1) +
  
  geom_vline(
    xintercept = as.Date("2007-01-01"),
    linetype = "dashed",
    color = "black",
    linewidth = 0.7) +
  
  
  
  scale_color_manual(
    values = colores) +
  
  scale_linetype_manual(
    values = lineas) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)) +
  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_date(
    date_breaks = "2 year",
    date_labels = "%Y"
  ) +
  
  labs(
    x = "Date",
    y = "Proportion of mature adults (%)") +
  
  theme_minimal() +
  
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# Monthly proportion of mature individuals by accessibility
madurez_prop_mensual_ac <- Mad_2sp_analisis %>%
  filter(!is.na(mature_imature),!is.na(accessibility)) %>%
  group_by(year, month, species, sampling_site, accessibility, protective_regime) %>%
  summarise(
    total = n(),
    maduros = sum(mature_imature == "Mature"),
    proporcion_madura = maduros / total,
    .groups = "drop"
  ) %>%
  mutate(
    grupo = interaction(species, accessibility, sep = " - "),
    fecha = lubridate::make_date(year, month, 1)  # For a continuous time axis
  )

# Visual style for accessibility plot
colores <- c(
  "Patella aspera - North" = "#D39C00",
  "Patella ordinaria - North" = "#0072B2",
  "Patella aspera - South" = "#D39C00",
  "Patella ordinaria - South" = "#0072B2"
)

lineas <- c(
  "Patella aspera - South" = "solid",
  "Patella ordinaria - South" = "solid",
  "Patella aspera - North" = "dotted",
  "Patella ordinaria - North" = "dotted"
)

# B plot: mature proportion by accessibility
MATURE_PROP_ACC <- ggplot(
  madurez_prop_mensual_ac,
  aes(x = fecha,
      y = proporcion_madura,
      color = grupo, 
      linetype = grupo)) +
  
  geom_point(aes(color = grupo)) +
  
  geom_line(
    linewidth = 0.1) +
  
  geom_vline(
    xintercept = as.Date("2007-01-01"),
    linetype = "dashed",
    color = "black",
    linewidth = 0.7) +
  
  scale_color_manual(
    values = colores) +
  
  scale_linetype_manual(
    values = lineas) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)) +
  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_date(
    date_breaks = "2 year",
    date_labels = "%Y"
  ) +
  
  labs(
    x = "Date",
    y = "Proportion of mature adults (%)") +
  
  theme_minimal() +
  
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# Plot both panels (A and B)
plot_grid(
  MATURE_PROP_protreg, MATURE_PROP_ACC,
  labels = c("A", "B"),
  label_size = 12,
  label_x = 0.05,   # Horizontal alignment (closer to edge)
  label_y = 1.001,    # Vertical alignment
  hjust = 0,        # Left justified
  vjust = 1,        # Top justified
  ncol = 1,
  rel_heights = c(1, 1)
)


