#### Salmon stability - Muñoz et al ####
### 03: comparing stability components across spatial scales ###

## SET UP -------------------------------------------------------------------
#library(tidyverse)
library(dabestr)
library(vegan)


## REGIONAL STABILITY AND SPATIAL STABILIZATION ----------------------------
# regional stability from both spatial scales, using components derived from total abundance
stb.np <- stb.all$gamma_w
spt.stb.np <- stb.all$spatial_stb
year.np <- 1956:2010    # represents the middle year of each 10-year rolling window
scale.np <- "North Pacific"        # add spatial scale label
reg.stb.np <- data.frame(stb.np, spt.stb.np, year.np, scale.np)

stb.bc <- stb.run$gamma_w
spt.stb.bc <- stb.run$spatial_stb
year.bc <- 1964:2007
scale.bc <- "Northern BC"
reg.stb.bc <- data.frame(stb.bc, spt.stb.bc, year.bc, scale.bc)

# rename variables
reg.stb.np <- rename(reg.stb.np, 
                     gamma = stb.np, spatial_stb = spt.stb.np, year = year.np, scale = scale.np)

reg.stb.bc <- rename(reg.stb.bc, 
                     gamma = stb.bc, spatial_stb = spt.stb.bc, year = year.bc, scale = scale.bc)

# combine data, only for the overlapping years
reg.stb <- reg.stb.np %>% 
  filter(year %in% c(1964:2007)) %>% 
  full_join(reg.stb.bc)

# mean differences and effect size using DABEST
# regional stability
paired.gamma <- load(reg.stb, x = scale, y = gamma,
         idx = c("Northern BC", "North Pacific"),
         paired = "sequential", id_col = year)

paired.gamma.diff <- mean_diff(paired.gamma)
paired.gamma.diff

# spatial stabilization
paired.spat <- load(reg.stb, x = scale, y = spatial_stb,
                     idx = c("Northern BC", "North Pacific"),
                     paired = "sequential", id_col = year)

paired.spat.diff <- mean_diff(paired.spat)
paired.spat.diff

## VARIANCE PARTITIONING ---------------------------------------------------
# North Pacific scale
alpha.np <- stb.all$alpha_w
phi.np <- stb.all$phi

var_np <- varpart(stb.np, ~alpha.np, ~phi.np)
var_np

# northern BC scale
alpha.nbc <- stb.run$alpha_w
phi.nbc <- stb.run$phi

var_nbc <- varpart(stb.bc, ~alpha.nbc, ~phi.nbc)
var_nbc

## ESCAPEMENT VS. CATCH ----------------------------------------------------
esc.gamma <- stb.esc$gamma_w
esc.type <- "escapement"
reg.stb.esc <- data.frame(year.bc, esc.gamma, esc.type)

catch.gamma <- stb.catch$gamma_w
catch.type <- "catch"
reg.stb.catch <- data.frame(year.bc, catch.gamma, catch.type)

# rename variables
reg.stb.esc <- rename(reg.stb.esc,
                      year = year.bc, gamma = esc.gamma, type = esc.type)

reg.stb.catch <- rename(reg.stb.catch,
                        year = year.bc, gamma = catch.gamma, type = catch.type)

# combine data
reg.stb.catch.esc <- reg.stb.esc %>% 
  full_join(reg.stb.catch)

# mean differences and effect size using DABEST
paired.gamma.catch <- reg.stb.catch.esc %>%
  load(x = type, y = gamma,
       idx = c("catch", "escapement"),
       paired = "sequential", id_col = year)

catch.diff <- mean_diff(paired.gamma.catch)
catch.diff