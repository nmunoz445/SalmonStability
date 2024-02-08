#### Salmon stability - Muñoz et al ####
### 05: tables and figures ###

# SET UP -------------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(dabestr)
library(magrittr)
library(lemon)

# FIG 1: ABUNDANCE, REG. STABILITY, SPATIAL STB. --------------------------
# 1c - abundance in North Pacific
Abundance_all <- data.long %>% # all salmon. by year
  filter(type == "all") %>% 
  group_by(year) %>% 
  summarise(Tot.Abundance = sum(abundance))

Abundance_sp <- data.long %>% # all salmon, by year and species
  filter(type == "all") %>% 
  group_by(year, species) %>% 
  summarise(Tot.Abundance = sum(abundance))

fig1c <- Abundance_all %>% 
  mutate(species = "Combined") %>% # Dummy variable, to merge with species data
  full_join(Abundance_sp) %>% 
  ggplot(aes(year, Tot.Abundance, color = species)) +
  geom_line(aes(size = species)) +
  scale_y_continuous(name = "Abundance (millions of fish)",
                     breaks = seq(0, 800000000, 200000000),
                     labels = c("0"="0", "200000000"="200", "400000000"="400", "600000000"="600", "800000000"="800")) +
  scale_x_continuous(name = element_blank(),
                     breaks = seq(1960, 2010, 20)) +
  labs(title = "North Pacific",
       tag = "(c)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = c(.32, .78),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(3.4, "mm"),
        plot.title = element_text(size = 9, face = "bold"),
        plot.tag = element_text(face = "bold", size = 9),
        plot.margin = unit(c(1,1,1,0.1), "mm")) +
  scale_color_manual(breaks = c("Combined", "pink", "chum", "sockeye"),
                     values = c("#39568CFF", "#CC79A7", "#440154FF", "#e31a1c")) +
  scale_size_manual(values = c(1.3, 0.8, 0.8, 0.8),
                    breaks = c("Combined", "pink", "chum", "sockeye")) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# 1f - abundance in northern BC
Abundance_allbc <- data_bc.long %>% # all salmon. by year
  filter(type == "run") %>% 
  group_by(year) %>% 
  summarise(Tot.Abundance = sum(abundance))

Abundance_spbc <- data_bc.long %>% # all salmon, by year and species
  filter(type == "run") %>% 
  group_by(year, species) %>% 
  summarise(Tot.Abundance = sum(abundance))

fig1f <- Abundance_allbc %>% 
  mutate(species = "Combined") %>% # Dummy variable, to merge with species data
  full_join(Abundance_spbc) %>% 
  ggplot(aes(year, Tot.Abundance, color = species)) +
  geom_line(aes(size = species)) +
  scale_y_continuous(name = "Abundance (millions of fish)",
                     breaks = seq(0, 50000000, 10000000),
                     labels = c("0"="0", "10000000"="10", "20000000"="20", "30000000"="30", 
                                "40000000"="40", "50000000"="50")) +
  coord_cartesian(xlim = c(1952, 2015)) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2010, 20)) +
  labs(title = "Northern BC",
       tag = "(f)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=6), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 9, face = "bold"),
        plot.tag = element_text(face = "bold", size = 9),
        plot.margin = unit(c(1,1,1,0.1), "mm")) +
  scale_color_manual(breaks = c("Combined", "pink", "chum", "sockeye"),
                     values = c("#39568CFF", "#CC79A7", "#440154FF", "#e31a1c")) +
  scale_size_manual(values = c(0.8, 1.3, 0.8, 0.8))

# 1d - regional stability in the NP
fig1d <- ggplot(reg.stb.np, (aes(year, gamma))) +
  geom_line(color = "#5ec962", size = 1.3) +
  ylab(expression(paste("Regional stability (", gamma["stb"],")"))) +
  scale_y_continuous(breaks = seq(2, 12, 2)) +
  coord_cartesian(xlim = c(1956, 2010),
                  ylim = c(1.8, 13)) +
  scale_x_continuous(name = element_blank(),
                     breaks = seq(1960, 2010, 10)) +
  labs(tag = "(d)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 9),
        plot.margin = unit(c(1,1,1,1), "mm"))

# 1g - regional stability in BC
fig1g <- ggplot(reg.stb.bc, (aes(year, gamma))) +
  geom_line(color = "#5ec962", size = 1.3) +
  ylab(expression(paste("Regional stability (", gamma["stb"],")"))) +
  scale_y_continuous(breaks = seq(2, 12, 2)) +
  coord_cartesian(xlim = c(1956, 2010),
                  ylim = c(1.8, 13)) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2010, 10)) +
  labs(tag = "(g)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=6), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 9),
        plot.margin = unit(c(1,1,1,1), "mm"))

# 1e - spatial stabilization in the NP
fig1e <- ggplot(reg.stb.np, (aes(year, spatial_stb))) +
  geom_line(color = "#D55E00", size = 1.3) +
  ylab(expression(paste("Spatial stabilization (", tau,")"))) +
  scale_y_continuous(breaks = seq(1, 5, 1)) +
  coord_cartesian(xlim = c(1956, 2010),
                  ylim = c(1, 5)) +
  scale_x_continuous(name = element_blank(),
                     breaks = seq(1960, 2010, 10)) +
  labs(tag = "(e)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 9),
        plot.margin = unit(c(1,0.1,1,1), "mm"))

# 1h - spatial stabilization in BC
fig1h <- ggplot(reg.stb.bc, (aes(year, spatial_stb))) +
  geom_line(color = "#D55E00", size = 1.3) +
  ylab(expression(paste("Spatial stabilization (", tau,")"))) +
  scale_y_continuous(breaks = seq(1, 5, 1)) +
  coord_cartesian(xlim = c(1956, 2010),
                  ylim = c(1, 5)) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2010, 10)) +
  labs(tag = "(h)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=6), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 9),
        plot.margin = unit(c(1,0.1,1,1), "mm"))

# Fig 2: compile panels
fig1 <- (fig1c | fig1d | fig1e) / (fig1f | fig1g | fig1h)
fig1
ggsave("fig1.tiff", plot = fig1, path = dir.out, dpi = 600, width = 180, height = 100, units = "mm")
ggsave("fig1_lowq.tiff", plot = fig1, path = dir.out, dpi = 96, width = 180, height = 100, units = "mm")


# FIG S1: 8 AND 12-YEAR WINDOWS -----------------------------------------
# regional stability from both spatial scales, using components derived from total abundance from 8 and 12 year windows
stb.np8 <- stb.all8$gamma_w
spt.stb.np8 <- stb.all8$spatial_stb
year.np8 <- 1955:2011    # represents the middle year of each 10-year rolling window
window.8 <- "8 years"
reg.stb.np8 <- data.frame(stb.np8, spt.stb.np8, year.np8, window.8)
reg.stb.np8 <- rename(reg.stb.np8, 
                      gamma = stb.np8, spatial_stb = spt.stb.np8, year = year.np8, windows = window.8)

stb.bc8 <- stb.run8$gamma_w
spt.stb.bc8 <- stb.run8$spatial_stb
year.bc8 <- 1963:2008
reg.stb.bc8 <- data.frame(stb.bc8, spt.stb.bc8, year.bc8, window.8)
reg.stb.bc8 <- rename(reg.stb.bc8, 
                      gamma = stb.bc8, spatial_stb = spt.stb.bc8, year = year.bc8, windows = window.8)

stb.np12 <- stb.all12$gamma_w
spt.stb.np12 <- stb.all12$spatial_stb
year.np12 <- 1957:2009
window.12 <- "12 years"
reg.stb.np12 <- data.frame(stb.np12, spt.stb.np12, year.np12, window.12)
reg.stb.np12 <- rename(reg.stb.np12, 
                       gamma = stb.np12, spatial_stb = spt.stb.np12, year = year.np12, windows = window.12)

stb.bc12 <- stb.run12$gamma_w
spt.stb.bc12 <- stb.run12$spatial_stb
year.bc12 <- 1965:2006
reg.stb.bc12 <- data.frame(stb.bc12, spt.stb.bc12, year.bc12, window.12)
reg.stb.bc12 <- rename(reg.stb.bc12, 
                       gamma = stb.bc12, spatial_stb = spt.stb.bc12, year = year.bc12, windows = window.12)

# format 10-year window dataframe
window.10 <- "10 years"
reg.stb.np10 <- data.frame(stb.np, spt.stb.np, year.np, window.10)
reg.stb.np10 <- rename(reg.stb.np10, 
                       gamma = stb.np, spatial_stb = spt.stb.np, year = year.np, windows = window.10)

reg.stb.bc10 <- data.frame(stb.bc, spt.stb.bc, year.bc, window.10)
reg.stb.bc10 <- rename(reg.stb.bc10, 
                       gamma = stb.bc, spatial_stb = spt.stb.bc, year = year.bc, windows = window.10)

# combine into dataframes
windows.np <- reg.stb.np8 %>% 
  full_join(reg.stb.np10) %>% 
  full_join(reg.stb.np12)

windows.bc <- reg.stb.bc8 %>% 
  full_join(reg.stb.bc10) %>% 
  full_join(reg.stb.bc12)

# S1a regional stability in the NP among windows
figs1a <- ggplot(windows.np, (aes(year, gamma, color = windows))) +
  geom_line(size = 1.3, aes(linetype = windows)) +
  ylab(expression(paste("Regional stability (", gamma["stb"],")"))) +
  scale_y_continuous(breaks = seq(2, 12, 2)) +
  coord_cartesian(xlim = c(1956, 2010),
                  ylim = c(1.8, 13.5)) +
  scale_x_continuous(name = element_blank(),
                     breaks = seq(1960, 2010, 10)) +
  labs(title = "North Pacific",
       tag = "(a)") +
  scale_linetype_manual(values = c("solid", "dashed", "twodash")) +
  scale_color_manual(values= c("#5ec962", "grey65", "grey40")) +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 10),
        axis.title.y = element_text(margin = margin(r=1.5), size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = c(0.45, 0.8),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(2.5, "mm"),
        plot.title = element_text(size = 9, face = "bold"),
        plot.tag = element_text(face = "bold", size = 9),
        plot.margin = unit(c(1,1,1,1), "mm"))

# S1b spatial stabilization in the NP among windows
figs1b <- ggplot(windows.np, (aes(year, spatial_stb, color = windows))) +
  geom_line(size = 1.3, aes(linetype = windows)) +
  ylab(expression(paste("Spatial stabilization (", tau,")"))) +
  scale_y_continuous(breaks = seq(1, 5, 1)) +
  coord_cartesian(xlim = c(1956, 2010),
                  ylim = c(1, 5)) +
  scale_x_continuous(name = element_blank(),
                     breaks = seq(1960, 2010, 10)) +
  labs(tag = "(b)") +
  scale_linetype_manual(values = c("solid", "dashed", "twodash")) +
  scale_color_manual(values= c("#D55E00", "grey60", "grey40")) +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 10),
        axis.title.y = element_text(margin = margin(r=1.5), size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.45, 0.8),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 10),
        plot.tag = element_text(face = "bold", size = 9),
        plot.margin = unit(c(1,0.1,1,1), "mm"))

# S1c regional stability in BC among windows
figs1c <- ggplot(windows.bc, (aes(year, gamma, color = windows))) +
  geom_line(size = 1.3, aes(linetype = windows)) +
  ylab(expression(paste("Regional stability (", gamma["stb"],")"))) +
  scale_y_continuous(breaks = seq(2, 12, 2)) +
  coord_cartesian(xlim = c(1956, 2010),
                  ylim = c(1.8, 13.5)) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2010, 10)) +
  labs(title = "Northern BC",
       tag = "(c)") +
  scale_linetype_manual(values = c("solid", "dashed", "twodash")) +
  scale_color_manual(values= c("#5ec962", "grey65", "grey40")) +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 10),
        axis.title.y = element_text(margin = margin(r=1.5), size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 9, face = "bold"),
        plot.tag = element_text(face = "bold", size = 9),
        plot.margin = unit(c(1,1,1,1), "mm"))

# S1d spatial stabilization in BC among windows
figs1d <- ggplot(windows.bc, (aes(year, spatial_stb, color = windows))) +
  geom_line(size = 1.3, aes(linetype = windows)) +
  ylab(expression(paste("Spatial stabilization (", tau,")"))) +
  scale_y_continuous(breaks = seq(1, 5, 1)) +
  coord_cartesian(xlim = c(1956, 2010),
                  ylim = c(1, 5)) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2010, 10)) +
  labs(tag = "(d)") +
  scale_linetype_manual(values = c("solid", "dashed", "twodash")) +
  scale_color_manual(values= c("#D55E00", "grey60", "grey40")) +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 10),
        axis.title.y = element_text(margin = margin(r=1.5), size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 9),
        plot.margin = unit(c(1,0.1,1,1), "mm"))

figs1 <- (figs1a | figs1b) / (figs1c | figs1d)
figs1
ggsave("figs1.tiff", plot = figs1, path = dir.out, dpi = 600, width = 134, height = 100, units = "mm")
ggsave("figs1_lowq.tiff", plot = figs1, path = dir.out, dpi = 96, width = 134, height = 100, units = "mm")

# FIG S2: EFFECT SIZE PLOTS -----------------------------------------------
# For paired comparisons of reg. stability and spat. stabilization across spatial scales
figs2a <- dabest_plot(paired.gamma.diff,
               swarm_label = "Regional stability (γstb)",
               raw_marker_size = 0.7,
               swarm_y_text = 12,
               contrast_y_text = 12)
figs2b <- dabest_plot(paired.spat.diff,
               swarm_label = "Spatial stabilization (τ)",
               raw_marker_size = 0.7,
               swarm_y_text = 12,
               contrast_y_text = 12)

# Compile figures
figs2 <- figs2a / figs2b
figs2
ggsave("figs2.tiff", plot = figs2, path = dir.out, dpi = 600, width = 141, height = 140, units = "mm")
ggsave("figs2_lowq.tiff", plot = figs2, path = dir.out, dpi = 96, width = 141, height = 140, units = "mm")

# FIG 2: INFLUENCE OF LOCAL STABILITY AND SPATIAL SYNCHRONY --------------------------------------------
# at the NP scale
df_np <- cbind.data.frame(stb.np, alpha.np, phi.np)

fig2a <- ggplot(df_np, (aes(alpha.np, stb.np))) +
  geom_point(color = "#56B4E9", size = 3) +
  ylab(expression(paste("Regional stability (", gamma["stb"],")"))) +
  xlab(expression(paste("Local stability (", alpha["stb"],")"))) +
  scale_y_continuous(breaks = seq(4, 12, 2), limits=c(3.8,13.1)) +
  scale_x_continuous(breaks=seq(1.8,2.8,0.2), limits=c(1.8,2.85)) +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=4), size = 11),
        axis.title.y = element_text(margin = margin(r=4), size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid = element_blank(),
        plot.margin = margin(t = 12, r = 1, b = 2, l = 8, unit = "mm"))

fig2b <- ggplot(df_np, (aes(phi.np, stb.np))) +
  geom_point(color = "#E69F00", size = 3) +
  ylab(element_blank()) +
  xlab(expression(paste("Spatial synchrony (", phi,")"))) +
  scale_y_continuous(breaks = seq(4, 12, 2), limits=c(3.8,13.1)) +
  scale_x_continuous(breaks=seq(0.1,0.3,0.1), limits = c(0.03,0.39)) +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=4), size = 11),
        axis.title.y = element_text(margin = margin(r=4), size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid = element_blank(),
        plot.margin = margin(t = 12, r = 0.5, b = 2, l = 2, unit = "mm"))

# at the BC scale
df_bc <- cbind.data.frame(stb.bc, alpha.nbc, phi.nbc)

fig2c <- ggplot(df_bc, (aes(alpha.nbc, stb.bc))) +
  geom_point(color = "#56B4E9", size = 3) +
  ylab(expression(paste("Regional stability (", gamma["stb"],")"))) +
  xlab(expression(paste("Local stability (", alpha["stb"],")"))) +
  scale_y_continuous(breaks = seq(2, 5, 1), limits = c(1.7,5.2)) +
  scale_x_continuous(breaks=seq(1,2,0.2), limits = c(1.2,2.05)) +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=4), size = 11),
        axis.title.y = element_text(margin = margin(r=4), size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid = element_blank(),
        plot.margin = margin(t = 10, r = 1, b = 0.3, l = 8, unit = "mm"))

fig2d <- ggplot(df_bc, (aes(phi.nbc, stb.bc))) +
  geom_point(color = "#E69F00", size = 3) +
  ylab(element_blank()) +
  xlab(expression(paste("Spatial synchrony (", phi,")"))) +
  scale_y_continuous(breaks = seq(2, 5, 1), limits = c(1.7,5.2)) +
  scale_x_continuous(breaks=seq(0.2,0.6,0.2), limits = c(0.15,0.7)) +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=4), size = 11),
        axis.title.y = element_text(margin = margin(r=4), size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid = element_blank(),
        plot.margin = margin(t = 10, r = 0.5, b = 0.3, l = 2, unit = "mm"))
#compile plots
fig2 <- (fig2a | fig2b) / (fig2c | fig2d)
fig2
ggsave("fig2.tiff", plot = fig2, path = dir.out, dpi = 600, width = 141, height = 140, units = "mm")
ggsave("fig2_lowq.tiff", plot = fig2, path = dir.out, dpi = 96, width = 141, height = 140, units = "mm")

# note: bar displaying variance in reg. stab. by component added externally

# FIG 3: CONTRIBUTION OF SPECIES TO REGIONAL STABILITY ## -----------------------------------
# 3a: regional stab in NP
# create dataframe with all species contributions
con.s.gamma.all <- as_tibble(con.s.gamma.all) %>% 
  mutate(species = "sockeye", year = 1956:2010) 

con.c.gamma.all <- as_tibble(con.c.gamma.all) %>% 
  mutate(species = "chum", year = 1956:2010)

con.p.gamma.all <- as_tibble(con.p.gamma.all) %>% 
  mutate(species = "pink", year = 1956:2010)

fig3a <- con.p.gamma.all %>% 
  full_join(con.c.gamma.all) %>% 
  full_join(con.s.gamma.all) %>%
  rename(gamma = value) %>% 
  ggplot(aes(year, gamma, color = species)) +
  geom_line(size = 1.3) +
  scale_y_continuous(name = "Contrib. regional stability (%)",
                     breaks = seq(-200, 50, 50)) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2000, 20)) +
  coord_cartesian(xlim = c(1956, 2010),
                  ylim = c(-200, 50)) +
  labs(title = "North Pacific",
       tag = "(a)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=4), size = 10),
        axis.title.y = element_text(margin = margin(r=4), size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = c(.25, .2),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.key.size = unit(2, "mm"),
        plot.title = element_text(size = 9, face = "bold"),
        plot.tag = element_text(size = 10, face = "bold"),
        plot.margin = unit(c(1,1,1,1), "mm")) +
  scale_color_manual(values = c("#CC79A7", "#440154FF", "#e31a1c"), 
                     breaks = c("pink", "chum", "sockeye")) +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  geom_hline(yintercept = 0, size = 1.3, color = "black", linetype = "dashed")
fig3a

# 3b: regional stab in BC
# create dataframe with all species contributions
con.bc.s.gamma <- as_tibble(con.bc.s.gamma) %>% 
  mutate(species = "sockeye", year = 1964:2007) 

con.bc.c.gamma <- as_tibble(con.bc.c.gamma) %>% 
  mutate(species = "chum", year = 1964:2007)

con.bc.p.gamma <- as_tibble(con.bc.p.gamma) %>% 
  mutate(species = "pink", year = 1964:2007)

fig3b <- con.bc.p.gamma %>% 
  full_join(con.bc.c.gamma) %>% 
  full_join(con.bc.s.gamma) %>%
  rename(gamma = value) %>% 
  ggplot(aes(year, gamma, color = species)) +
  geom_line(size = 1.3) +
  scale_y_continuous(name = "Dummy",
                     breaks = seq(-200, 50, 50)) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2000, 20)) +
  coord_cartesian(xlim = c(1956, 2010),
                  ylim = c(-200, 50)) +
  labs(title = "Northern BC",
       tag = "(b)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=4), size = 10),
        axis.title.y = element_text(margin = margin(t=1), size = 10, colour = "white"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 9, face = "bold"),
        plot.tag = element_text(size = 10, face = "bold"),
        plot.margin = unit(c(1,1,0.1,0.1), "mm")) +
  scale_color_manual(values = c("#CC79A7", "#440154FF", "#e31a1c"), 
                     breaks = c("pink", "chum", "sockeye")) +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  geom_hline(yintercept = 0, size = 1.3, color = "black", linetype = "dashed")
fig3b

fig3 <- fig3a + fig3b
fig3
ggsave("fig3.tiff", plot = fig3, dpi = 600, width = 141, height = 58, units = "mm")
ggsave("fig3_lowq.tiff", plot = fig5b, dpi = 96, width = 141, height = 58, units = "mm")


# FIG S4: DRIVERS OF REG STABILITY AMONG SPECIES IN NP ### --------------
## chum - alpha
chum_a_lab <- paste0("~R^{2} == 0.09") # r2 value

con.c.alpha.gamma <- as_tibble(con.c.alpha.all) %>% 
  mutate(gamma = gamma_all)

figs4a <- ggplot(con.c.alpha.gamma, aes(x = value, y = gamma)) +
  geom_point(color = "#440154FF", size = 2) +
  stat_smooth(method = 'lm', color = "#440154FF", size = 1.3) +
  geom_text(aes(x = 11.5, y = 8, label = chum_a_lab), 
                stat = "unique", size = 3.6, color = "black",
                parse = TRUE) +
  geom_text(aes(x = 14.25, y = 11.5, label = "chum"), 
                stat = "unique", size = 3.6, color = "black") +
  ylab(expression(paste("Regional stability (", gamma["stb"],")"))) +
  scale_y_continuous(breaks = seq(4, 12, 2),
                     limits = c(3.8, 12.8)) +
  scale_x_continuous(name = "Contrib. local stability (%)",
                     breaks = seq(10, 18, 2)) +
  labs(title = "North Pacific",
       tag = "(a)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.title = element_text(size = 9, face = "bold"),
        plot.tag = element_text(size = 9, face = "bold"))
figs4a

# chum - spatial stabilization
con.c.spat.gamma <- as_tibble(con.c.spat.all) %>% 
  mutate(gamma = gamma_all)

figs4d <- ggplot(con.c.spat.gamma, aes(x = value, y = gamma)) +
  geom_point(color = "#440154FF", size = 2) +
  ylab(expression(paste("Regional stability (", gamma["stb"],")"))) +
  scale_y_continuous(breaks = seq(4, 12, 2),
                     limits = c(3.8, 12.8)) +
  scale_x_continuous(name = "Contrib. spatial stabilization (%)",
                     breaks = seq(-10, 15, 5)) +
  labs(tag = "(d)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.tag = element_text(size = 9, face = "bold"))
figs4d

# sockeye - alpha
sock_a_lab <- paste0("~R^{2} == 0.08") # r2 value

con.s.alpha.gamma <- as_tibble(con.s.alpha.all) %>% 
  mutate(gamma = gamma_all)

figs4b <- ggplot(con.s.alpha.gamma, aes(x = value, y = gamma)) +
  geom_point(color = "#e31a1c", size = 2) +
  stat_smooth(method = 'lm', color = "#e31a1c", size = 1.3) +
  geom_text(aes(x = 10.7, y = 8, label = sock_a_lab), 
                stat = "unique", size = 3.6, color = "black",
                parse = TRUE) +
  geom_text(aes(x = 6.25, y = 11.5, label = "sockeye"), 
                stat = "unique", size = 3.6, color = "black") +
  scale_y_continuous(name = element_blank(), 
                     breaks = seq(4, 12, 2),
                     limits = c(3.8, 12.8)) +
  scale_x_continuous(name = "Contrib. local stability (%)") +
  labs(tag = "(b)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.tag = element_text(size = 9, face = "bold"))
figs4b

# sockeye - spatial stabilization
con.s.spat.gamma <- as_tibble(con.s.spat.all) %>% 
  mutate(gamma = gamma_all)

figs4e <- ggplot(con.s.spat.gamma, aes(x = value, y = gamma)) +
  geom_point(color = "#e31a1c", size = 2) +
  scale_y_continuous(name = element_blank(),
                     breaks = seq(4, 12, 2),
                     limits = c(3.8, 12.8)) +
  scale_x_continuous(name = "Contrib. spatial stabilization (%)",
                     breaks = seq(-15, 10, 5)) +
  labs(tag = "(e)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.tag = element_text(size = 9, face = "bold"))
figs4e

# pink - alpha
pink_a_lab <- paste0("~R^{2} == 0.08") # r2 value

con.p.alpha.gamma <- as_tibble(con.p.alpha.all) %>% 
  mutate(gamma = gamma_all)

figs4c <- ggplot(con.p.alpha.gamma, aes(x = value, y = gamma)) +
  geom_point(color = "#CC79A7", size = 2) +
  stat_smooth(method = 'lm', color = "#CC79A7", size = 1.3) +
  geom_text(aes(x = -65, y = 8, label = pink_a_lab), 
                stat = "unique", size = 3.6, color = "black",
                parse = TRUE) +
  geom_text(aes(x = -35, y = 11.5, label = "pink"), 
                stat = "unique", size = 3.6, color = "black") +
  scale_y_continuous(name = element_blank(), 
                     breaks = seq(4, 12, 2),
                     limits = c(3.8, 12.8)) +
  labs(tag = "(c)") +
  scale_x_continuous(name = "Contrib. local stability (%)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.tag = element_text(size = 9, face = "bold"))
figs4c

# pink - spatial stabilization
pink_s_lab <- paste0("~R^{2} == 0.39") # r2 value

con.p.spat.gamma <- as_tibble(con.p.spat.all) %>% 
  mutate(gamma = gamma_all)

figs4f <- ggplot(con.p.spat.gamma, aes(x = value, y = gamma)) +
  geom_point(color = "#CC79A7", size = 2) +
  stat_smooth(method = 'lm', color = "#CC79A7", size = 1.3) +
  geom_text(aes(x = -53, y = 8, label = pink_s_lab), 
                stat = "unique", size = 3.6, color = "black",
                parse = TRUE) +
  scale_y_continuous(name = element_blank(), 
                     breaks = seq(4, 12, 2),
                     limits = c(3.8, 12.8)) +
  scale_x_continuous(name = "Contrib. spatial stabilization (%)",
                     breaks = seq(-75, 50, 25)) +
  labs(tag = "(f)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.title = element_text(size = 9, face = "bold"),
        plot.tag = element_text(size = 9, face = "bold"))
figs4f

# compile plots
figs4 <- (figs4a | figs4b | figs4c) / (figs4d | figs4e | figs4f)
figs4

ggsave("figs4.tiff", plot = figs4, path = dir.out, dpi = 600, width = 180, height = 100, units = "mm")
ggsave("figs4_lowq.tiff", plot = figs4, path = dir.out, dpi = 96, width = 180, height = 100, units = "mm")

# FIG S5: DRIVERS OF REG STABILITY AMONG SPECIES IN BC ### --------------
## chum - alpha
con.bc.c.alpha.gamma <- as_tibble(con.bc.c.alpha) %>% 
  mutate(gamma = gamma_run)

figs5a <- ggplot(con.bc.c.alpha.gamma, aes(x = value, y = gamma)) +
  geom_point(color = "#440154FF", size = 2) +
  geom_text(aes(x = 8.2, y = 4.7, label = "chum"), 
                stat = "unique", size = 3.6, color = "black") +
  ylab(expression(paste("Regional stability (", gamma["stb"],")"))) +
  scale_y_continuous(breaks = seq(2, 5, 1),
                     limits = c(1.7, 5)) +
  scale_x_continuous(name = "Contrib. local stability (%)",
                     breaks = seq(4, 14, 2)) +
  labs(title = "Northern BC",
       tag = "(a)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.title = element_text(size = 9, face = "bold"),
        plot.tag = element_text(size = 9, face = "bold"))
figs5a

# chum - spatial stabilization
chum_s_labbc <- paste0("~R^{2} == 0.11") # r2 value

con.bc.c.stab.gamma <- as_tibble(con.bc.c.stab) %>% 
  mutate(gamma = gamma_run)

figs5d <- ggplot(con.bc.c.stab.gamma, aes(x = value, y = gamma)) +
  geom_point(color = "#440154FF", size = 2) +
  stat_smooth(method = 'lm', color = "#440154FF", size = 1.3) +
  geom_text(aes(x = -4, y = 4, label = chum_s_labbc), 
                stat = "unique", size = 3.6, color = "black",
                parse = TRUE) +
  ylab(expression(paste("Regional stability (", gamma["stb"],")"))) +
  scale_y_continuous(breaks = seq(2, 5, 1),
                     limits = c(1.7, 5)) +
  scale_x_continuous(name = "Contrib. spatial stabilization (%)") +
  labs(tag = "(d)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.tag = element_text(size = 9, face = "bold"))
figs5d

# sockeye - alpha
con.bc.s.alpha.gamma <- as_tibble(con.bc.s.alpha) %>% 
  mutate(gamma = gamma_run)

figs5b <- ggplot(con.bc.s.alpha.gamma, aes(x = value, y = gamma)) +
  geom_point(color = "#e31a1c", size = 2) +
  geom_text(aes(x = 12.25, y = 4.7, label = "sockeye"), 
                stat = "unique", size = 3.6, color = "black") +
  scale_y_continuous(name = element_blank(), 
                     breaks = seq(2, 5, 1),
                     limits = c(1.7, 5)) +
  scale_x_continuous(name = "Contrib. local stability (%)") +
  labs(tag = "(b)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.tag = element_text(size = 9, face = "bold"))
figs5b

# sockeye - spatial stabilization
sock_s_labbc <- paste0("~R^{2} == 0.19") # r2 value

con.bc.s.stab.gamma <- as_tibble(con.bc.s.stab) %>% 
  mutate(gamma = gamma_run)

figs5e <- ggplot(con.bc.s.stab.gamma, aes(x = value, y = gamma)) +
  geom_point(color = "#e31a1c", size = 2) +
  stat_smooth(method = 'lm', color = "#e31a1c", size = 1.3) +
  geom_text(aes(x = -3.5, y = 4, label = sock_s_labbc), 
                stat = "unique", size = 3.6, color = "black",
            parse = TRUE) +
  scale_y_continuous(name = element_blank(),
                     breaks = seq(2, 5, 1),
                     limits = c(1.7, 5)) +
  scale_x_continuous(name = "Contrib. spatial stabilization (%)",
                     breaks = seq(-5, 15, 5)) +
  labs(tag = "(e)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.tag = element_text(size = 9, face = "bold"))
figs5e

# pink - alpha
pink_a_labbc <- paste0("~R^{2} == 0.10") # r2 value

con.bc.p.alpha.gamma <- as_tibble(con.bc.p.alpha) %>% 
  mutate(gamma = gamma_run)

figs5c <- ggplot(con.bc.p.alpha.gamma, aes(x = value, y = gamma)) +
  geom_point(color = "#CC79A7", size = 2) +
  stat_smooth(method = 'lm', color = "#CC79A7", size = 1.3) +
  geom_text(aes(x = -62, y = 4, label = pink_a_labbc), 
                stat = "unique", size = 3.6, color = "black",
                parse = TRUE) +
  geom_text(aes(x = -36, y = 4.7, label = "pink"), 
                stat = "unique", size = 3.6, color = "black") +
  scale_y_continuous(name = element_blank(), 
                     breaks = seq(2, 5, 1),
                     limits = c(1.7, 5)) +
  scale_x_continuous(name = "Contrib. local stability (%)") +
  labs(tag = "(c)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.title = element_text(size = 9, face = "bold"),
        plot.tag = element_text(size = 9, face = "bold"))
figs5c

# pink - spatial stabilization
pink_s_labbc <- paste0("~R^{2} == 0.37") # r2 value

con.bc.p.stab.gamma <- as_tibble(con.bc.p.stab) %>% 
  mutate(gamma = gamma_run)

figs5f <- ggplot(con.bc.p.stab.gamma, aes(x = value, y = gamma)) +
  geom_point(color = "#CC79A7", size = 2) +
  stat_smooth(method = 'lm', color = "#CC79A7", size = 1.3) +
  geom_text(aes(x = -44, y = 4, label = pink_s_labbc), 
                stat = "unique", size = 3.6, color = "black",
                parse = TRUE) +
  scale_y_continuous(name = element_blank(), 
                     breaks = seq(2, 5, 1),
                     limits = c(1.7, 5)) +
  scale_x_continuous(name = "Contrib. spatial stabilization (%)") +
  labs(tag = "(f)") +
  theme_bw(base_line_size = 0.5) +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 9),
        axis.title.y = element_text(margin = margin(r=1.5), size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.tag = element_text(size = 9, face = "bold"))
figs5f

# compile plots
figs5 <- (figs5a | figs5b | figs5c) / (figs5d | figs5e | figs5f)
figs5

ggsave("figs5.tiff", plot = figs5, path = dir.out, dpi = 600, width = 180, height = 100, units = "mm")
ggsave("figs5_lowq.tiff", plot = figs5, path = dir.out, dpi = 96, width = 180, height = 100, units = "mm")

# FIG S6: CONTRIBUTION OF LOCAL COMMUNITIES TO REG STABILITY IN NP --------
con.1.gamma.all <- con.1.gamma %>% 
  as.tibble() %>% 
  mutate(Year = 1956:2010,
         Region = "Washington & Oregon") %>% 
  rename(contrib = value)

con.2.gamma.all <- con.2.gamma %>% 
  as.tibble() %>% 
  mutate(Year = 1956:2010,
         Region = "Southern BC") %>% 
  rename(contrib = value)

con.3.gamma.all <- con.3.gamma %>% 
  as.tibble() %>% 
  mutate(Year = 1956:2010,
         Region = "Northern BC") %>% 
  rename(contrib = value)

con.4.gamma.all <- con.4.gamma %>% 
  as.tibble() %>% 
  mutate(Year = 1956:2010,
         Region = "Cook Inlet Alaska") %>% 
  rename(contrib = value)

con.5.gamma.all <- con.5.gamma %>% 
  as.tibble() %>% 
  mutate(Year = 1956:2010,
         Region = "Kodiak Island Alaska") %>% 
  rename(contrib = value)

con.6.gamma.all <- con.6.gamma %>% 
  as.tibble() %>% 
  mutate(Year = 1956:2010,
         Region = "PW Sound Alaska") %>% 
  rename(contrib = value)

con.7.gamma.all <- con.7.gamma %>% 
  as.tibble() %>% 
  mutate(Year = 1956:2010,
         Region = "South-Central Alaska") %>% 
  rename(contrib = value)

con.8.gamma.all <- con.8.gamma %>% 
  as.tibble() %>% 
  mutate(Year = 1956:2010,
         Region = "Southeast Alaska") %>% 
  rename(contrib = value)

con.9.gamma.all <- con.9.gamma %>% 
  as.tibble() %>% 
  mutate(Year = 1956:2010,
         Region = "Western Alaska") %>% 
  rename(contrib = value)

con.10.gamma.all <- con.10.gamma %>% 
  as.tibble() %>% 
  mutate(Year = 1956:2010,
         Region = "Japan") %>% 
  rename(contrib = value)

con.11.gamma.all <- con.11.gamma %>% 
  as.tibble() %>% 
  mutate(Year = 1956:2010,
         Region = "Eastern Kamchatka") %>% 
  rename(contrib = value)

con.12.gamma.all <- con.12.gamma %>% 
  as.tibble() %>% 
  mutate(Year = 1956:2010,
         Region = "Russia mainland") %>% 
  rename(contrib = value)

con.13.gamma.all <- con.13.gamma %>% 
  as.tibble() %>% 
  mutate(Year = 1956:2010,
         Region = "Western Kamchatka") %>% 
  rename(contrib = value)

con.lc.gamma.all <- con.1.gamma.all %>% 
  full_join(con.2.gamma.all) %>% 
  full_join(con.3.gamma.all) %>%
  full_join(con.4.gamma.all) %>% 
  full_join(con.5.gamma.all) %>% 
  full_join(con.6.gamma.all) %>%
  full_join(con.7.gamma.all) %>%
  full_join(con.8.gamma.all) %>% 
  full_join(con.9.gamma.all) %>% 
  full_join(con.10.gamma.all) %>%
  full_join(con.11.gamma.all) %>% 
  full_join(con.12.gamma.all) %>%
  full_join(con.13.gamma.all)

con.lc.gamma.all$Region <- factor(con.lc.gamma.all$Region, 
                                  levels = c("Washington & Oregon", "Southern BC", "Northern BC", "Southeast Alaska",
                                             "PW Sound Alaska", "Cook Inlet Alaska", "Kodiak Island Alaska", 
                                             "South-Central Alaska", "Western Alaska", "Eastern Kamchatka", "Western Kamchatka",
                                             "Russia mainland", "Japan"))

figs6 <- ggplot(con.lc.gamma.all, aes(Year, contrib)) +
  geom_line(size = 1.3, color = "#5ec962") +
  scale_y_continuous(name = "Contrib. regional stability (%)",
                     breaks = seq(-150, 50, 50)) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2000, 20)) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 8.3), 
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t=6), size = 10),
        axis.title.y = element_text(margin = margin(r=4), size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(size = 0.4),
        panel.grid.minor.x = element_line(size = 0.4),
        legend.position = "none") +
  geom_hline(yintercept = 0, size = 1.3, color = "black", linetype = "dashed") +
  facet_wrap(~Region, nrow = 3)
figs6

ggsave("figs6.tiff", plot = figs6, path = dir.out, dpi = 600, width = 180, height = 158, units = "mm")
ggsave("figs6_lowq.tiff", plot = figs5, path = dir.out, dpi = 96, width = 180, height = 158, units = "mm")

# FIG S7: CONTRIBUTION OF LOCAL COMMUNITIES TO REG STABILITY IN BC --------
con.bc.1.g.all <- as.tibble(con.bc.1.g) %>% 
  mutate(area = "Area 1", year = 1964:2007) %>% 
  rename(Contrib = value)

con.bc.2.g.all <- as.tibble(con.bc.2.g) %>% 
  mutate(area = "Area 2", year = 1964:2007) %>% 
  rename(Contrib = value)

con.bc.3.g.all <- as.tibble(con.bc.3.g) %>% 
  mutate(area = "Area 3", year = 1964:2007) %>% 
  rename(Contrib = value)

con.bc.4.g.all <- as.tibble(con.bc.4.g) %>% 
  mutate(area = "Area 4", year = 1964:2007) %>% 
  rename(Contrib = value)

con.bc.5.g.all <- as.tibble(con.bc.5.g) %>% 
  mutate(area = "Area 5", year = 1964:2007) %>% 
  rename(Contrib = value)

con.bc.6.g.all <- as.tibble(con.bc.6.g) %>% 
  mutate(area = "Area 6", year = 1964:2007) %>% 
  rename(Contrib = value)

con.bc.7.g.all <- as.tibble(con.bc.7.g) %>% 
  mutate(area = "Area 7", year = 1964:2007) %>% 
  rename(Contrib = value)

con.bc.8.g.all <- as.tibble(con.bc.8.g) %>% 
  mutate(area = "Area 8", year = 1964:2007) %>% 
  rename(Contrib = value)

con.bc.9.g.all <- as.tibble(con.bc.9.g) %>% 
  mutate(area = "Area 9", year = 1964:2007) %>% 
  rename(Contrib = value)

con.bc.10.g.all <- as.tibble(con.bc.10.g) %>% 
  mutate(area = "Area 10", year = 1964:2007) %>% 
  rename(Contrib = value)

con.st.g.all <- con.bc.1.g.all %>% 
  full_join(con.bc.2.g.all) %>% 
  full_join(con.bc.3.g.all) %>%
  full_join(con.bc.4.g.all) %>% 
  full_join(con.bc.5.g.all) %>% 
  full_join(con.bc.6.g.all) %>%
  full_join(con.bc.7.g.all) %>%
  full_join(con.bc.8.g.all) %>% 
  full_join(con.bc.9.g.all) %>% 
  full_join(con.bc.10.g.all)

con.st.g.all$area <- factor(con.st.g.all$area, 
                                levels = c("Area 1", "Area 2", "Area 3", "Area 4", "Area 5", 
                                           "Area 6", "Area 7", "Area 8", "Area 9", "Area 10"))

figs7 <- ggplot(con.st.g.all, aes(year, Contrib)) +
  geom_line(size = 1.3, color = "#5ec962") +
  scale_y_continuous(name = "Contrib. regional stability (%)",
                     breaks = seq(-40, 40, 20)) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2000, 20),
                     limits = c(1960, 2012)) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 8.5), 
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t=6), size = 10),
        axis.title.y = element_text(margin = margin(r=4), size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(size = 0.4),
        panel.grid.minor.x = element_line(size = 0.4),
        legend.position = "none") +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  geom_hline(yintercept = 0, size = 1.3, color = "black", linetype = "dashed") +
  facet_wrap(~area, nrow = 2)
figs7

ggsave("figs7.tiff", plot = figs7, path = dir.out, dpi = 600, width = 180, height = 82, units = "mm")
ggsave("figs7_lowq.tiff", plot = figs7, path = dir.out, dpi = 96, width = 180, height = 82, units = "mm")

# FIG S8 AND S9: ABUNDANCE BY SPECIES AND LOCAL COMMUNITY -----------------
# Fig S8 in the North Pacific
Abundance_sp_lc <- data.long %>% # all salmon, by year, species and local community
  filter(type == "all") %>% 
  group_by(year, species, region) %>% 
  summarise(Tot.Abundance = sum(abundance))

Abundance_sp_lc$region <- factor(Abundance_sp_lc$region, 
                                 levels = c("wa.or", "bc.s", "bc.n", "ak.se", "ak.pws", "ak.cook",
                                            "ak.kod", "ak.sc", "ak.west", "russ.e.kam", "russ.w.kam",
                                            "russ.main", "japan"),
                                 labels = c("Washington & Oregon", "Southern BC", "Northern BC", "Southeast Alaska",
                                            "PW Sound Alaska", "Cook Inlet Alaska", "Kodiak Island Alaska", 
                                            "South-central Alaska", "Western Alaska", "Eastern Kamchatka", "Western Kamchatka",
                                            "Russia Mainland", "Japan"))

figs8 <- ggplot(Abundance_sp_lc, aes(year, Tot.Abundance, color = species)) +
  geom_line(size = 0.75) +
  scale_y_continuous(name = "Annual abundance (millions of fish)",
                     labels = c("0"="0", "50000000"="50", "100000000"="100", "150000000"="150","
                                200000000"="200", "250000000"="250")) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2000, 20)) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 8), 
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t=6), size = 10),
        axis.title.y = element_text(margin = margin(r=4), size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(size = 0.4),
        panel.grid.minor.x = element_line(size = 0.4),
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values = c("#CC79A7", "#440154FF", "#e31a1c"), 
                     breaks = c("pink", "chum", "sockeye")) +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  facet_wrap(~region, nrow = 3)

figs8save <- reposition_legend(figs8, 'center', panel = 'panel-2-3')

ggsave("figs8.tiff", plot = figs8save, path = dir.out, dpi = 600, width = 180, height = 158, units = "mm")
ggsave("figs8_lowq.tiff", plot = figs7save, path = dir.out, dpi = 96, width = 180, height = 158, units = "mm")

# Fig S9 in BC
Abundance_sp_lcbc <- data_bc.long %>% # all salmon, by year, species and local community
  filter(type == "run") %>% 
  group_by(year, species, stat.area) %>% 
  summarise(Tot.Abundance = sum(abundance))

Abundance_sp_lcbc$stat.area <- factor(Abundance_sp_lcbc$stat.area,
                                      labels = c("Area 1", "Area 2", "Area 3", "Area 4", "Area 5",
                                                 "Area 6", "Area 7", "Area 8", "Area 9", "Area 10"))

figs9 <- ggplot(Abundance_sp_lcbc, aes(year, Tot.Abundance, color = species)) +
  geom_line(size = 0.75) +
  scale_y_continuous(name = "Annual abundance (millions of fish)",
                     labels = c("0"="0", "10000000"="10", "20000000"="20", "30000000"="30")) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2000, 20)) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 8.5), 
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t=6), size = 10),
        axis.title.y = element_text(margin = margin(r=4), size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(size = 0.4),
        panel.grid.minor.x = element_line(size = 0.4),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = c(0.1, 0.9),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(2.5, "mm")) +
  scale_color_manual(values = c("#CC79A7", "#440154FF", "#e31a1c"), 
                     breaks = c("pink", "chum", "sockeye")) +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  facet_wrap(~stat.area, nrow = 2)
figs9
ggsave("figs9.tiff", plot = figs9, path = dir.out, dpi = 600, width = 180, height = 106, units = "mm")
ggsave("figs9_lowq.tiff", plot = figs9, path = dir.out, dpi = 96, width = 180, height = 106, units = "mm")

# FIG 4: EFFECT OF HATCHERIES ON REGIONAL AND LOCAL STABILITY --------------------------
## Fig 4a: annual abundance of hatchery salmon
fig4a <- Abundance_or %>% 
  filter(type == "hatch") %>% 
  ggplot(aes(year, Tot.Abundance)) + 
  geom_line(size = 1.3, color = "#39568CFF") +
  scale_y_continuous(name = "Abundance (millions of fish)",
                     breaks = seq(0, 200000000, 50000000),
                     labels = c("0"="0", "50000000"="50", "100000000"="100", "150000000"="150", "200000000"="200")) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2000, 20)) +
  theme_bw() +
  theme(axis.title.x = element_text(margin = margin(t=6), size = 10),
        axis.title.y = element_text(margin = margin(r=4), size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(size = 0.4),
        panel.grid.minor.x = element_line(size = 0.4),
        legend.position = "none")
fig4a

ggsave("fig4a.tiff", plot = fig4a, path = dir.out, dpi = 600, width = 82, height = 59, units = "mm")
ggsave("figs4a_lowq.tiff", plot = figs4a, path = dir.out, dpi = 96, width = 82, height = 59, units = "mm")

## Fig 4b: contribution of hatchery salmon to regional stability
fig4b <- as_tibble(con.hatch.gamma) %>% 
  mutate(year = year.np) %>% 
  rename(contrib = value) %>% 
  ggplot(aes(year, contrib)) +
  geom_line(size = 1.3, color = "#5ec962") +
  scale_y_continuous(name = "Contrib. regional stability (%)",
                     breaks = seq(-0, 25, 5)) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2000, 20),
                     limits = c(1956, 2010)) +
  theme_bw() +
  theme(axis.title.x = element_text(margin = margin(t=6), size = 10),
        axis.title.y = element_text(margin = margin(r=4), size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(size = 0.4),
        panel.grid.minor.x = element_line(size = 0.4),
        legend.position = "none") +
  geom_hline(yintercept = 0, size = 1.3, color = "black", linetype = "dashed")
fig4b

ggsave("fig4b.tiff", plot = fig4b, path = dir.out, dpi = 600, width = 82, height = 59, units = "mm")
ggsave("fig4b_lowq.tiff", plot = fig4b, path = dir.out, dpi = 96, width = 82, height = 59, units = "mm")

## Fig 4c: effect of hatchery salmon on local stability of each local community
# convert matrix of alpha_i into a long-form data frame
alpha_i_natdf <- as_tibble(stb.nat$alpha_i) %>% 
  set_colnames(year.np) %>% 
  mutate(Region = regions) %>% 
  pivot_longer(cols = '1956':'2010',
               names_to = "Year") %>% 
  rename(alpha_i = value) %>% 
  mutate(Origin = "nat") 

alpha_i_df <- as_tibble(stb.all$alpha_i) %>% 
  set_colnames(year.np) %>% 
  mutate(Region = regions) %>% 
  pivot_longer(cols = '1956':'2010',
               names_to = "Year") %>% 
  rename(alpha_i = value) %>% 
  mutate(Origin = "all") %>% 
  full_join(alpha_i_natdf)

# new data frame for the plot
alpha_i_plot <- alpha_i_df
alpha_i_plot$Year <- as.integer(alpha_i_plot$Year)
alpha_i_plot$Region <- factor(alpha_i_plot$Region, 
                              levels = c("wa.or", "bc.s", "bc.n", "ak.se", "ak.pws", "ak.cook",
                                         "ak.kod", "ak.sc", "ak.west", "russ.e.kam", "russ.w.kam",
                                         "russ.main", "japan"),
                              labels = c("Washington & Oregon", "Southern BC", "Northern BC", "Southeast Alaska",
                                         "PW Sound Alaska", "Cook Inlet Alaska", "Kodiak Island Alaska", 
                                         "South-central Alaska", "Western Alaska", "Eastern Kamchatka", "Western Kamchatka",
                                         "Russia Mainland", "Japan")) # change the order and labels for the plot

fig4c <- ggplot(alpha_i_plot, aes(Year, alpha_i, color = Origin)) +
  geom_line(size = 1.3) +
  ylab(expression(paste("Local community stability (", alpha["stb,i"],")"))) +
  scale_y_continuous(breaks = seq(2, 8, 2)) +
  scale_x_continuous(breaks = seq(1960, 2000, 20)) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 8.5), 
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t=6), size = 10),
        axis.title.y = element_text(margin = margin(r=4), size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(size = 0.4),
        panel.grid.minor.x = element_line(size = 0.4),
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values = c("black", "#009E73"),
                     breaks = c("all", "nat"),
                     labels = c("Natural + Hatchery", "Natural")) +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  facet_wrap(~Region, nrow = 3)

fig4csave <- reposition_legend(fig4c, 'center', panel = c('panel-2-3', 'panel-5-3'))
fig4csave

ggsave("fig4c.tiff", plot = fig4csave, path = dir.out, dpi = 600, width = 180, height = 158, units = "mm")
ggsave("fig4c_lowq.tiff", plot = fig4csave, path = dir.out, dpi = 96, width = 180, height = 158, units = "mm")

# FIG S10 : HATCHERY- AND NATURAL-ORIGIN ABUNDANCES ---------------------------
Abundance_or_lc <- data.long %>% # by origin (hatchery vs natural), year and local community
  filter(type %in% c("nat","hatch")) %>% 
  group_by(year, type, region) %>% 
  summarise(Tot.Abundance = sum(abundance)) %>% 
  mutate(log.abundance = log(Tot.Abundance +1))

# new data frame for the plot
Abundance_or_lc.plot <- Abundance_or_lc
Abundance_or_lc.plot$region <- factor(Abundance_or_lc.plot$region, 
                                      levels = c("wa.or", "bc.s", "bc.n", "ak.se", "ak.pws", "ak.cook",
                                                 "ak.kod", "ak.sc", "ak.west", "russ.e.kam", "russ.w.kam",
                                                 "russ.main", "japan"),
                                      labels = c("Washington & Oregon", "Southern BC", "Northern BC", "Southeast Alaska",
                                                 "PW Sound Alaska", "Cook Inlet Alaska", "Kodiak Island Alaska", 
                                                 "South-central Alaska", "Western Alaska", "Eastern Kamchatka", "Western Kamchatka",
                                                 "Russia Mainland", "Japan")) # change the order and labels for the plot

figs10 <- ggplot(Abundance_or_lc.plot, aes(year, Tot.Abundance, fill = type)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_y_continuous(name = "Annual abundance (millions of fish)",
                     labels = c("0"="0", "100000000"="100", "200000000"="200", "300000000"="300")) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2000, 20)) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 8.3), 
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t=6), size = 10),
        axis.title.y = element_text(margin = margin(r=4), size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(size = 0.4),
        panel.grid.minor.x = element_line(size = 0.4),
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#999999", "#E69F00"),
                    breaks = c("nat", "hatch"),
                    labels = c("Natural", "Hatchery")) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  facet_wrap(~region, nrow = 3)

figs10save <- reposition_legend(figs10, 'center', panel = 'panel-2-3')

ggsave("figs10.tiff", plot = figs10save, path = dir.out, dpi = 600, width = 180, height = 158, units = "mm")
ggsave("figs10_lowq.tiff", plot = figs10save, path = dir.out, dpi = 96, width = 180, height = 158, units = "mm")

# FIG 5: STABILITY OF ESCAPEMENT AND CATCH IN BC ------------------------------------------
# 5a: abundances
Abundance_esc <- data_bc.long %>% 
  filter(type == "escapement") %>% 
  group_by(year) %>% 
  summarise(Tot.Abundance = sum(abundance)) %>% 
  mutate(type = "escapement")

Abundance_catch <- data_bc.long %>% 
  filter(type == "harvest") %>% 
  group_by(year) %>% 
  summarise(Tot.Abundance = sum(abundance)) %>% 
  mutate(type = "catch") %>% 
  full_join(Abundance_esc)


fig5a <- ggplot(Abundance_catch, aes(year, Tot.Abundance, color = type)) +
  geom_line(size = 1) +
  scale_y_continuous(name = "Abundance (millions of fish)",
                     labels = c("0"="0", "10000000"="10", "20000000"="20", "30000000"="30")) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2010, 20)) +
  labs(tag = "(a)") +
  theme_bw() +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 10),
        axis.title.y = element_text(margin = margin(r=4), size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 9),
        plot.margin = unit(c(1,1,1,1), "mm")) +
  scale_color_manual(values = c("#0072B2", "#D55E00"),
                     breaks = c("escapement","catch"))

# 5b: regional stability
fig5b <- ggplot(reg.stb.catch.esc, (aes(year, gamma, color = type))) +
  geom_line(size = 1.3) +
  ylab(expression(paste("Regional stability (", gamma["stb"],")"))) +
  scale_x_continuous(name = "Year",
                     breaks = seq(1960, 2010, 10)) +
  labs(tag = "(b)") +
  theme_bw() +
  theme(axis.title.x = element_text(margin = margin(t=2), size = 10),
        axis.title.y = element_text(margin = margin(r=3), size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size = 1.1, fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = c(.68, .91),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.key.size = unit(3.4, "mm"),
        legend.background = element_rect(fill = "transparent"),
        plot.tag = element_text(face = "bold", size = 9),
        plot.margin = unit(c(1,1,1,1), "mm")) +
  scale_color_manual(values = c("#0072B2", "#D55E00"),
                     breaks = c("escapement","catch")) + 
  guides(color = guide_legend(override.aes = list(size = 3)))

#compile plots
fig5 <- fig5a | fig5b
fig5
ggsave("fig5.tiff", plot = fig5, path = dir.out, dpi = 600, width = 110, height = 59, units = "mm")
ggsave("fig5_lowq.tiff", plot = fig5, path = dir.out, dpi = 96, width = 110, height = 59, units = "mm")

# FIG S11: CATCH VS ESCAPEMENT EFFECT SIZE ----------------------------------
figs11 <- dabest_plot(catch.diff,
               swarm_label = "Regional stability (γstb)",
               raw_marker_size = 0.7,
               swarm_y_text = 12,
               contrast_y_text = 12)
figs11

ggsave("figs11.tiff", plot = figs11, path = dir.out, dpi = 600, width = 141, height = 70, units = "mm")
ggsave("figs11_lowq.tiff", plot = figs12, path = dir.out, dpi = 96, width = 141, height = 70, units = "mm")