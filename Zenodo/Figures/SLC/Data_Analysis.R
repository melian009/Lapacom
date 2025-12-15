rm(list = ls())
# Data import and manipulation
library(readxl)
library(janitor)
library(dplyr)
library(tidyr)
library(stringr)

# Data wrangling & plotting
library(ggplot2)
library(scales)
library(cowplot)
library(flextable)

# Statistical modeling and tests
library(car)
library(broom)
library(emmeans)
library(conover.test)
library(MASS)


# Data import ans manipulation

Mad_RAW <- read_excel("Madeira/BD_LIMPETS_MAD_1996-2018.xlsx",
                      sheet = "Data", range = "A1:T69707")
Mad_RAW_Clean <- Mad_RAW %>% clean_names()

Mad_selected <- Mad_RAW_Clean %>%
  dplyr::select(species, year, month, total_length_mm, total_length_class_mm,
                weight_g, mature_imature, sampling_site, lat, long,
                protective_regime, proximity_human_settlements,
                accessibility, age_lt, age_months, age_class)

# Filter by studied species
Mad_2sp <- Mad_selected %>%
  filter(species %in% c("Patella ordinaria", "Patella aspera"))


Mad_2sp_analisis <- Mad_2sp %>%
  filter(!is.na(sampling_site) & 
           trimws(sampling_site) != "" & 
           sampling_site %in% c("Porto Moniz", "Paúl do Mar", "Funchal", "Desertas", 
                                "Caniçal", "Santa Cruz", "Ribeira Brava", "São Vicente"))%>% 
  mutate(across(where(is.character), as.factor))

# Grouping of samples before and after the 2007 regulations

Mad_2sp_analisis <- Mad_2sp_analisis %>%
  mutate(regulation_period = if_else(year < 2007, "NR", "R") %>% as.factor())



# ───────────────────────────────────────────────
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

# Final plot: mature proportion by protection regime
MATURE_PROP_protreg <- ggplot(madurez_prop_mensual_pr, aes(x = fecha, y = proporcion_madura,
                                                           color = grupo, linetype = grupo)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = as.Date("2007-01-01"), linetype = "dashed", color = "black", linewidth = 0.7) +
  scale_color_manual(values = colores) +
  scale_linetype_manual(values = lineas) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Date", y = "Mature Proportion (%)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# Monthly proportion of mature individuals by accessibility
madurez_prop_mensual_ac <- Mad_2sp_analisis %>%
  filter(!is.na(mature_imature)) %>%
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

# Final plot: mature proportion by accessibility
MATURE_PROP_ACC <- ggplot(madurez_prop_mensual_ac, aes(x = fecha, y = proporcion_madura,
                                                       color = grupo, linetype = grupo)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = as.Date("2007-01-01"), linetype = "dashed", color = "black", linewidth = 0.7) +
  scale_color_manual(values = colores) +
  scale_linetype_manual(values = lineas) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Date", y = "Mature proportion (%)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# Plot both panels (A and B)
plot_grid(
  MATURE_PROP_protreg, MATURE_PROP_ACC,
  labels = c("A)", "B)"),
  label_size = 16,
  label_x = 0.02,   # Horizontal alignment (closer to edge)
  label_y = 0.2,    # Vertical alignment
  hjust = 0,        # Left justified
  vjust = 1,        # Top justified
  ncol = 1,
  rel_heights = c(1, 1)
)



# Diagnostic plot function
plot_glm_diagnostics <- function(model, species) {
  par(mfrow = c(2, 2))
  plot(model, main = paste("GLM Diagnostic -", species))
  par(mfrow = c(1, 1))
}



# Patella aspera analysis

mad_aspera <- Mad_2sp_analisis %>%
  filter(species == "Patella aspera" & !is.na(mature_imature)) %>%
  mutate(mature_binary = ifelse(mature_imature == "Mature", 1, 0))

mad_aspera_modelo <- mad_aspera %>%
  group_by(year, protective_regime, regulation_period) %>%
  summarise(
    total = n(),
    maduros = sum(mature_binary),
    .groups = "drop"
  )

modelo_aspera <- glm(cbind(maduros, total - maduros) ~ protective_regime * regulation_period,
                     family = binomial, data = mad_aspera_modelo)
modelo_opt_aspera <- stepAIC(modelo_aspera, direction = "both", trace = FALSE)

# Store model summary table
resultado_glm_aspera <- tidy(modelo_opt_aspera) %>%
  mutate(Species = "Patella aspera",
         p.value = formatC(p.value, format = "e", digits = 2))

# Store emmeans contrasts table
resultado_emm_aspera <- emmeans(modelo_opt_aspera, pairwise ~ protective_regime * regulation_period, type = "response")$contrasts %>%
  as.data.frame() %>%
  mutate(Species = "Patella aspera",
         p.value = formatC(p.value, format = "e", digits = 2))

# Diagnostics
plot_glm_diagnostics(modelo_opt_aspera, "Patella aspera")



# Patella ordinaria analysis

mad_ordinaria <- Mad_2sp_analisis %>%
  filter(species == "Patella ordinaria" & !is.na(mature_imature)) %>%
  mutate(mature_binary = ifelse(mature_imature == "Mature", 1, 0))

mad_ordinaria_modelo <- mad_ordinaria %>%
  group_by(year, protective_regime, regulation_period) %>%
  summarise(
    total = n(),
    maduros = sum(mature_binary),
    .groups = "drop"
  )

modelo_ordinaria <- glm(cbind(maduros, total - maduros) ~ protective_regime * regulation_period,
                        family = binomial, data = mad_ordinaria_modelo)
modelo_opt_ordinaria <- stepAIC(modelo_ordinaria, direction = "both", trace = FALSE)

# Store model summary table
resultado_glm_ordinaria <- tidy(modelo_opt_ordinaria) %>%
  mutate(Species = "Patella ordinaria",
         p.value = formatC(p.value, format = "e", digits = 2))

# Store emmeans contrasts table
resultado_emm_ordinaria <- emmeans(modelo_opt_ordinaria, pairwise ~ protective_regime * regulation_period, type = "response")$contrasts %>%
  as.data.frame() %>%
  mutate(Species = "Patella ordinaria",
         p.value = formatC(p.value, format = "e", digits = 2))

# Diagnostics
plot_glm_diagnostics(modelo_opt_ordinaria, "Patella ordinaria")



# Present tables with flextable

# GLM results
flextable(resultado_glm_aspera %>% mutate(across(where(is.numeric), round, 3))) %>%
  set_caption("GLM Summary - Patella aspera") %>% autofit()

flextable(resultado_glm_ordinaria %>% mutate(across(where(is.numeric), round, 3))) %>%
  set_caption("GLM Summary - Patella ordinaria") %>% autofit()

# EMMEANS contrasts
flextable(resultado_emm_aspera %>% mutate(across(where(is.numeric), round, 3))) %>%
  set_caption("Marginal Effects (EMMEANS) - Patella aspera") %>% autofit()

flextable(resultado_emm_ordinaria %>% mutate(across(where(is.numeric), round, 3))) %>%
  set_caption("Marginal Effects (EMMEANS) - Patella ordinaria") %>% autofit()
