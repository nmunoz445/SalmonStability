#### Salmon stability - Muñoz et al ####
### 04: contribution calculations ###

# SET UP -------------------------------------------------------------------
library(tidyverse)
library(psych)

## note: for the contributions of species and local communities below, each species/community was removed manually.
## These calculations should be automated to reduce the extensive amount of code.

# ### Species contributions in North Pacific Ocean ### -------------------------------------------
# Sockeye contribution ----------------------------------------------------
data.soc <- data[!(data$species== "sockeye"),] # remove sockeye from dataset

# Reshape dataset for use in function
data.long.soc <- data.soc %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))
# Function to calculate stability metrics for each moving window across time series
stb.fcn.s <- function(data.soc) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.soc[data.soc$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1

    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}
## USE FUNCTION
# Stability values for total run
stb.all.soc <- stb.fcn.s(data.long.soc[data.long.soc$type == "all",])
stb.all.soc

# Contributions
con.s.gamma.all <- ((stb.all$gamma_w - stb.all.soc$gamma_w) / stb.all$gamma_w) * 100
con.s.gamma.all

con.s.alpha.all <- ((stb.all$alpha_w - stb.all.soc$alpha_w) / stb.all$alpha_w) * 100
con.s.alpha.all

con.s.spat.all <- ((stb.all$spatial_stb - stb.all.soc$spatial_stb) / stb.all$spatial_stb) * 100
con.s.spat.all

# Pink contribution -------------------------------------------------------
data.p <- data[!(data$species== "pink"),]

# Reshape dataset for use in function
data.long.p <- data.p %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))

stb.fcn.p <- function(data.p) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.p[data.p$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}

## USE FUNCTION
# Stability values for total run
stb.all.p <- stb.fcn.p(data.long.p[data.long.p$type == "all",])
stb.all.p

# Contributions
con.p.gamma.all <- ((stb.all$gamma_w - stb.all.p$gamma_w) / stb.all$gamma_w) * 100
con.p.gamma.all

con.p.alpha.all <- ((stb.all$alpha_w - stb.all.p$alpha_w) / stb.all$alpha_w) * 100
con.p.alpha.all

con.p.spat.all <- ((stb.all$spatial_stb - stb.all.p$spatial_stb) / stb.all$spatial_stb) * 100
con.p.spat.all

# Chum contribution -------------------------------------------------------
data.c <- data[!(data$species== "chum"),]

# Reshape dataset for use in function
data.long.c <- data.c %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))

stb.fcn.c <- function(data.c) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.c[data.c$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}
## USE FUNCTION
# Stability values for total run
stb.all.c <- stb.fcn.c(data.long.c[data.long.c$type == "all",])
stb.all.c

# Contributions
con.c.gamma.all <- ((stb.all$gamma_w - stb.all.c$gamma_w) / stb.all$gamma_w) * 100
con.c.gamma.all

con.c.alpha.all <- ((stb.all$alpha_w - stb.all.c$alpha_w) / stb.all$alpha_w) * 100
con.c.alpha.all

con.c.spat.all <- ((stb.all$spatial_stb - stb.all.c$spatial_stb) / stb.all$spatial_stb) * 100
con.c.spat.all


# ### Species contributions in northern BC ### ----------------------------
# Sockeye contribution ----------------------------------------------------
data_bc.soc <- data_bc[!(data_bc$species== "sockeye"),]

# Reshape data_bcset for use in function
data_bc.soc.long <- data_bc.soc %>% gather(type, abundance, c(escapement, total.harvest, total.run)) %>% 
  mutate(type = str_replace_all(type, c("escapement" = "escapement", "total.harvest" = "harvest", 
                                        "total.run" = "run")))
data_bc.soc.long$cdn.harvest <- NULL
data_bc.soc.long$noncdn.harvest <- NULL

## STABILITY FUNCTION
# Function to calculate stability metrics for each moving window across time series
stb.fcn.bc.s <- function(data_bc.soc) {
  for(y in 0:(length(years.bc) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window - 1))
    temp <- data_bc.soc[data_bc.soc$year %in% temp.years.bc,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m.bc, nrow = length(temp.years.bc))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m.bc) {   # for each community
      
      for(t in 1:length(temp.years.bc)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$stat.area == stat.areas[i] & temp$year == temp.years.bc[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M.bc[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m.bc
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m.bc
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M.bc[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M.bc[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE.bc[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m.bc^2) * sum(w_i_j) + (1/m.bc) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi.bc[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi.bc)
    
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi.bc
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w.bc[y+1] <- u_M.bc[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])
    
    }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, phi = phi.bc, spatial_stb = spatial_stb.bc, 
              u_M = u_M.bc, PE = PE.bc))
}
## USE FUNCTION
# Stability values for total run
stb.run.soc <- stb.fcn.bc.s(data_bc.soc.long[data_bc.soc.long$type == "run",])
stb.run.soc

# Contributions
con.bc.s.gamma <- ((stb.run$gamma_w - stb.run.soc$gamma_w) / stb.run$gamma_w) * 100
con.bc.s.gamma

con.bc.s.alpha <- ((stb.run$alpha_w - stb.run.soc$alpha_w) / stb.run$alpha_w) * 100
con.bc.s.alpha

con.bc.s.stab <- ((stb.run$spatial_stb - stb.run.soc$spatial_stb) / stb.run$spatial_stb) * 100
con.bc.s.stab

# Pink contribution ----------------------------------------------------
data_bc.p <- data_bc[!(data_bc$species== "pink"),]

# Reshape data_bcset for use in function
data_bc.p.long <- data_bc.p %>% gather(type, abundance, c(escapement, total.harvest, total.run)) %>% 
  mutate(type = str_replace_all(type, c("escapement" = "escapement", "total.harvest" = "harvest", 
                                        "total.run" = "run")))
data_bc.p.long$cdn.harvest <- NULL
data_bc.p.long$noncdn.harvest <- NULL

## STABILITY FUNCTION
# Function to calculate stability metrics for each moving window across time series
stb.fcn.bc.p <- function(data_bc.p) {
  for(y in 0:(length(years.bc) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window - 1))
    temp <- data_bc.p[data_bc.p$year %in% temp.years.bc,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m.bc, nrow = length(temp.years.bc))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m.bc) {   # for each community
      
      for(t in 1:length(temp.years.bc)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$stat.area == stat.areas[i] & temp$year == temp.years.bc[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M.bc[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m.bc
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m.bc
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M.bc[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M.bc[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE.bc[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m.bc^2) * sum(w_i_j) + (1/m.bc) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi.bc[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi.bc)
    
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi.bc
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w.bc[y+1] <- u_M.bc[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])
    
  }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, phi = phi.bc, spatial_stb = spatial_stb.bc, 
              u_M = u_M.bc, PE = PE.bc))
}
## USE FUNCTION
# Stability values for total run
stb.run.p <- stb.fcn.bc.p(data_bc.p.long[data_bc.p.long$type == "run",])
stb.run.p

# Contributions
con.bc.p.gamma <- ((stb.run$gamma_w - stb.run.p$gamma_w) / stb.run$gamma_w) * 100
con.bc.p.gamma

con.bc.p.alpha <- ((stb.run$alpha_w - stb.run.p$alpha_w) / stb.run$alpha_w) * 100
con.bc.p.alpha

con.bc.p.stab <- ((stb.run$spatial_stb - stb.run.p$spatial_stb) / stb.run$spatial_stb) * 100
con.bc.p.stab

# Chum contribution ----------------------------------------------------
data_bc.c <- data_bc[!(data_bc$species== "chum"),]

# Reshape data_bcset for use in function
data_bc.c.long <- data_bc.c %>% gather(type, abundance, c(escapement, total.harvest, total.run)) %>% 
  mutate(type = str_replace_all(type, c("escapement" = "escapement", "total.harvest" = "harvest", 
                                        "total.run" = "run")))
data_bc.c.long$cdn.harvest <- NULL
data_bc.c.long$noncdn.harvest <- NULL

## STABILITY FUNCTION 
# Function to calculate stability metrics for each moving window across time series
stb.fcn.bc.c <- function(data_bc.c) {
  for(y in 0:(length(years.bc) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window - 1))
    temp <- data_bc.c[data_bc.c$year %in% temp.years.bc,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m.bc, nrow = length(temp.years.bc))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m.bc) {   # for each community
      
      for(t in 1:length(temp.years.bc)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$stat.area == stat.areas[i] & temp$year == temp.years.bc[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M.bc[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m.bc
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m.bc
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M.bc[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M.bc[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE.bc[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m.bc^2) * sum(w_i_j) + (1/m.bc) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi.bc[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi.bc)
    
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi.bc
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w.bc[y+1] <- u_M.bc[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])
    
  }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, phi = phi.bc, spatial_stb = spatial_stb.bc, 
              u_M = u_M.bc, PE = PE.bc))
}
## USE FUNCTION
# Stability values for total run
stb.run.c <- stb.fcn.bc.c(data_bc.c.long[data_bc.c.long$type == "run",])
stb.run.c

# Contributions
con.bc.c.gamma <- ((stb.run$gamma_w - stb.run.c$gamma_w) / stb.run$gamma_w) * 100
con.bc.c.gamma

con.bc.c.alpha <- ((stb.run$alpha_w - stb.run.c$alpha_w) / stb.run$alpha_w) * 100
con.bc.c.alpha

con.bc.c.stab <- ((stb.run$spatial_stb - stb.run.c$spatial_stb) / stb.run$spatial_stb) * 100
con.bc.c.stab

# ### Partitioning variance in regional stability among species contributions ----------------------
## NP
gamma_all <- stb.all$gamma_w

# gamma vs. contributions to alpha stability and spatial stabilization
pink_a <- lm(gamma_all ~ con.p.alpha.all)
summary(pink_a) #r2 = 0.08
pink_s <- lm(gamma_all ~ con.p.spat.all)
summary(pink_s) #r2 = 0.39

chum_a <- lm(gamma_all ~ con.c.alpha.all)
summary(chum_a) #r2 = 0.09
chum_s <- lm(gamma_all ~ con.c.spat.all)
summary(chum_s)

sock_a <- lm(gamma_all ~ con.s.alpha.all)
summary(sock_a) # r2 = 0.08
sock_s <- lm(gamma_all ~ con.s.spat.all)
summary(sock_s)

## BC
gamma_run <- stb.run$gamma_w

# gamma vs. contributions to alpha stability and spatial stabilization
pinkbc_a <- lm(gamma_run ~ con.bc.p.alpha)
summary(pinkbc_a) #r2 = 0.10
pinkbc_s <- lm(gamma_run ~ con.bc.p.stab)
summary(pinkbc_s) #r2 = 0.37

chumbc_a <- lm(gamma_run ~ con.bc.c.alpha)
summary(chumbc_a)
chumbc_s <- lm(gamma_run ~ con.bc.c.stab)
summary(chumbc_s) #r2 = 0.11

sockbc_a <- lm(gamma_run ~ con.bc.s.alpha)
summary(sockbc_a)
sockbc_s <- lm(gamma_run ~ con.bc.s.stab)
summary(sockbc_s) #r2 = 0.19


# ### Contributions of local communities to reg stability in NP ### -------
# Washington/Oregon contribution -------------------------------------------------------
data.1 <- data[!(data$region=="wa.or"),]

# Store for later
regions.1 <- unique(data.1$region)

# Number of communities in the metacommunity
m_con <- length(unique(data.1$region))

# Reshape dataset for use in function
data.long.1 <- data.1 %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))

# Function to calculate stability metrics for each moving window across time series
stb.fcn.1 <- function(data.1) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.1[data.1$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m_con, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m_con) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions.1[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m_con
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1

    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
    # Confirm equal
    # stb.all$gamma_w / stb.all$alpha_w
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}

## USE FUNCTION
stb.all.1 <- stb.fcn.1(data.long.1[data.long.1$type == "all",])
stb.all.1

# Contributions
con.1.gamma <- ((stb.all$gamma_w - stb.all.1$gamma_w) / stb.all$gamma_w) * 100
con.1.gamma

# Southern BC contribution -------------------------------------------------------
data.2 <- data[!(data$region=="bc.s"),]

# Store for later
regions.2 <- unique(data.2$region)

# Reshape dataset for use in function
data.long.2 <- data.2 %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))

# Function to calculate stability metrics for each moving window across time series
stb.fcn.2 <- function(data.2) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.2[data.2$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m_con, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m_con) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions.2[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m_con
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
    # Confirm equal
    # stb.all$gamma_w / stb.all$alpha_w
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}

## USE FUNCTION
stb.all.2 <- stb.fcn.2(data.long.2[data.long.2$type == "all",])
stb.all.2

# Contributions
con.2.gamma <- ((stb.all$gamma_w - stb.all.2$gamma_w) / stb.all$gamma_w) * 100
con.2.gamma

# Northern BC contribution -------------------------------------------------------
data.3 <- data[!(data$region=="bc.n"),]

# Store for later
regions.3 <- unique(data.3$region)

# Reshape dataset for use in function
data.long.3 <- data.3 %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))

# Function to calculate stability metrics for each moving window across time series
stb.fcn.3 <- function(data.3) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.3[data.3$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m_con, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m_con) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions.3[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m_con
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
    # Confirm equal
    # stb.all$gamma_w / stb.all$alpha_w
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}

## USE FUNCTION
stb.all.3 <- stb.fcn.3(data.long.3[data.long.3$type == "all",])
stb.all.3

# Contributions
con.3.gamma <- ((stb.all$gamma_w - stb.all.3$gamma_w) / stb.all$gamma_w) * 100
con.3.gamma

# Alaska Cook Inlet contribution -------------------------------------------------------
data.4 <- data[!(data$region=="ak.cook"),]

# Store for later
regions.4 <- unique(data.4$region)

# Reshape dataset for use in function
data.long.4 <- data.4 %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))

# Function to calculate stability metrics for each moving window across time series
stb.fcn.4 <- function(data.4) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.4[data.4$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m_con, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m_con) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions.4[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m_con
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
    # Confirm equal
    # stb.all$gamma_w / stb.all$alpha_w
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}

## USE FUNCTION
stb.all.4 <- stb.fcn.4(data.long.4[data.long.4$type == "all",])
stb.all.4

# Contributions
con.4.gamma <- ((stb.all$gamma_w - stb.all.4$gamma_w) / stb.all$gamma_w) * 100
con.4.gamma

# Alaska Kodiak Island contribution -------------------------------------------------------
data.5 <- data[!(data$region=="ak.kod"),]

# Store for later
regions.5 <- unique(data.5$region)

# Reshape dataset for use in function
data.long.5 <- data.5 %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))

# Function to calculate stability metrics for each moving window across time series
stb.fcn.5 <- function(data.5) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.5[data.5$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m_con, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m_con) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions.5[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m_con
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
    # Confirm equal
    # stb.all$gamma_w / stb.all$alpha_w
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}

## USE FUNCTION
stb.all.5 <- stb.fcn.5(data.long.5[data.long.5$type == "all",])
stb.all.5

# Contributions
con.5.gamma <- ((stb.all$gamma_w - stb.all.5$gamma_w) / stb.all$gamma_w) * 100
con.5.gamma

# Alaska Prince William Sound contribution -------------------------------------------------------
data.6 <- data[!(data$region=="ak.pws"),]

# Store for later
regions.6 <- unique(data.6$region)

# Reshape dataset for use in function
data.long.6 <- data.6 %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))

# Function to calculate stability metrics for each moving window across time series
stb.fcn.6 <- function(data.6) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.6[data.6$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m_con, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m_con) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions.6[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m_con
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
    # Confirm equal
    # stb.all$gamma_w / stb.all$alpha_w
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}

## USE FUNCTION
stb.all.6 <- stb.fcn.6(data.long.6[data.long.6$type == "all",])
stb.all.6

# Contributions
con.6.gamma <- ((stb.all$gamma_w - stb.all.6$gamma_w) / stb.all$gamma_w) * 100
con.6.gamma

# Alaska South Central contribution -------------------------------------------------------
data.7 <- data[!(data$region=="ak.sc"),]

# Store for later
regions.7 <- unique(data.7$region)

# Reshape dataset for use in function
data.long.7 <- data.7 %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))

# Function to calculate stability metrics for each moving window across time series
stb.fcn.7 <- function(data.7) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.7[data.7$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m_con, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m_con) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions.7[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m_con
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
    # Confirm equal
    # stb.all$gamma_w / stb.all$alpha_w
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}

## USE FUNCTION
stb.all.7 <- stb.fcn.7(data.long.7[data.long.7$type == "all",])
stb.all.7

# Contributions
con.7.gamma <- ((stb.all$gamma_w - stb.all.7$gamma_w) / stb.all$gamma_w) * 100
con.7.gamma

# Alaska Southeast contribution -------------------------------------------------------
data.8 <- data[!(data$region=="ak.se"),]

# Store for later
regions.8 <- unique(data.8$region)

# Reshape dataset for use in function
data.long.8 <- data.8 %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))

# Function to calculate stability metrics for each moving window across time series
stb.fcn.8 <- function(data.8) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.8[data.8$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m_con, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m_con) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions.8[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m_con
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
    # Confirm equal
    # stb.all$gamma_w / stb.all$alpha_w
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}

## USE FUNCTION
stb.all.8 <- stb.fcn.8(data.long.8[data.long.8$type == "all",])
stb.all.8

# Contributions
con.8.gamma <- ((stb.all$gamma_w - stb.all.8$gamma_w) / stb.all$gamma_w) * 100
con.8.gamma

# Alaska west contribution -------------------------------------------------------
data.9 <- data[!(data$region=="ak.west"),]

# Store for later
regions.9 <- unique(data.9$region)

# Reshape dataset for use in function
data.long.9 <- data.9 %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))

# Function to calculate stability metrics for each moving window across time series
stb.fcn.9 <- function(data.9) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.9[data.9$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m_con, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m_con) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions.9[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m_con
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
    # Confirm equal
    # stb.all$gamma_w / stb.all$alpha_w
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}

## USE FUNCTION
stb.all.9 <- stb.fcn.9(data.long.9[data.long.9$type == "all",])
stb.all.9

# Contributions
con.9.gamma <- ((stb.all$gamma_w - stb.all.9$gamma_w) / stb.all$gamma_w) * 100
con.9.gamma

# Japan contribution -------------------------------------------------------
data.10 <- data[!(data$region=="japan"),]

# Store for later
regions.10 <- unique(data.10$region)

# Reshape dataset for use in function
data.long.10 <- data.10 %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))

# Function to calculate stability metrics for each moving window across time series
stb.fcn.10 <- function(data.10) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.10[data.10$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m_con, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m_con) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions.10[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m_con
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
    # Confirm equal
    # stb.all$gamma_w / stb.all$alpha_w
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}

## USE FUNCTION
stb.all.10 <- stb.fcn.10(data.long.10[data.long.10$type == "all",])
stb.all.10

# Contributions
con.10.gamma <- ((stb.all$gamma_w - stb.all.10$gamma_w) / stb.all$gamma_w) * 100
con.10.gamma

# Russia eastern Kamchatka contribution -------------------------------------------------------
data.11 <- data[!(data$region=="russ.e.kam"),]

# Store for later
regions.11 <- unique(data.11$region)

# Reshape dataset for use in function
data.long.11 <- data.11 %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))

# Function to calculate stability metrics for each moving window across time series
stb.fcn.11 <- function(data.11) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.11[data.11$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m_con, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m_con) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions.11[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m_con
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
    # Confirm equal
    # stb.all$gamma_w / stb.all$alpha_w
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}

## USE FUNCTION
stb.all.11 <- stb.fcn.11(data.long.11[data.long.11$type == "all",])
stb.all.11

# Contributions
con.11.gamma <- ((stb.all$gamma_w - stb.all.11$gamma_w) / stb.all$gamma_w) * 100
con.11.gamma

# Russia mainland and islands contribution -------------------------------------------------------
data.12 <- data[!(data$region=="russ.main"),]

# Store for later
regions.12 <- unique(data.12$region)

# Reshape dataset for use in function
data.long.12 <- data.12 %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))

# Function to calculate stability metrics for each moving window across time series
stb.fcn.12 <- function(data.12) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.12[data.12$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m_con, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m_con) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions.12[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m_con
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
    # Confirm equal
    # stb.all$gamma_w / stb.all$alpha_w
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}

## USE FUNCTION
stb.all.12 <- stb.fcn.12(data.long.12[data.long.12$type == "all",])
stb.all.12

# Contributions
con.12.gamma <- ((stb.all$gamma_w - stb.all.12$gamma_w) / stb.all$gamma_w) * 100
con.12.gamma

# Russia western kamchatka contribution -------------------------------------------------------
data.13 <- data[!(data$region=="russ.w.kam"),]

# Store for later
regions.13 <- unique(data.13$region)

# Reshape dataset for use in function
data.long.13 <- data.13 %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))

# Function to calculate stability metrics for each moving window across time series
stb.fcn.13 <- function(data.13) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data.13[data.13$year %in% temp.years,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m_con, nrow = length(temp.years))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m_con) {   # for each community
      
      for(t in 1:length(temp.years)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$region == regions.13[i] & temp$year == temp.years[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m_con
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m^2) * sum(w_i_j) + (1/m) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w[y+1] <- u_M[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w[y+1] <- sum((u_i/u_M[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
    # Confirm equal
    # stb.all$gamma_w / stb.all$alpha_w
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, phi = phi, spatial_stb = spatial_stb, u_M = u_M, PE = PE))
}

## USE FUNCTION
stb.all.13 <- stb.fcn.13(data.long.13[data.long.13$type == "all",])
stb.all.13

# Contributions
con.13.gamma <- ((stb.all$gamma_w - stb.all.13$gamma_w) / stb.all$gamma_w) * 100
con.13.gamma


# ### Contributions of local communities to reg stability in BC ### -------
# SA 1 contribution -------------------------------------------------------
data_bc.1 <- data_bc[!(data_bc$stat.area==1),]

# Store for later
stat.areas1 <- unique(data_bc.1$stat.area)

# Number of communities in the metacommunity
m.bc_con <- length(unique(data_bc.1$stat.area))

# Reshape dataset for use in function
data_bc.long.1 <- data_bc.1 %>% gather(type, abundance, c(escapement, total.harvest, total.run)) %>% 
  mutate(type = str_replace_all(type, c("escapement" = "escapement", "total.harvest" = "harvest", 
                                        "total.run" = "run")))
data_bc.long.1$cdn.harvest <- NULL
data_bc.long.1$noncdn.harvest <- NULL

# Function to calculate stability metrics for each moving window across time series
stb.fcn.bc.1 <- function(data_bc.1) {
  for(y in 0:(length(years.bc) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window - 1))
    temp <- data_bc.1[data_bc.1$year %in% temp.years.bc,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m.bc_con, nrow = length(temp.years.bc))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m.bc_con) {   # for each community
      
      for(t in 1:length(temp.years.bc)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$stat.area == stat.areas1[i] & temp$year == temp.years.bc[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M.bc[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m.bc_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m.bc_con
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M.bc[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M.bc[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE.bc[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m.bc_con^2) * sum(w_i_j) + (1/m.bc_con) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2

    # Spatial synchrony among communities
    phi.bc[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi.bc)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi.bc
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w.bc[y+1] <- u_M.bc[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])
    
  }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, phi = phi.bc, spatial_stb = spatial_stb.bc, u_M = u_M.bc, PE = PE.bc))
}

## USE FUNCTION
# Stability values for total run
stb.run.1 <- stb.fcn.bc.1(data_bc.long.1[data_bc.long.1$type == "run",])
stb.run.1

# Contributions
con.bc.1.g <- ((stb.run$gamma_w - stb.run.1$gamma_w) / stb.run$gamma_w) * 100
con.bc.1.g

# SA 2 contribution -------------------------------------------------------
data_bc.2 <- data_bc[!(data_bc$stat.area==2),]

# Store for later
stat.areas2 <- unique(data_bc.2$stat.area)

# Reshape dataset for use in function
data_bc.long.2 <- data_bc.2 %>% gather(type, abundance, c(escapement, total.harvest, total.run)) %>% 
  mutate(type = str_replace_all(type, c("escapement" = "escapement", "total.harvest" = "harvest", 
                                        "total.run" = "run")))
data_bc.long.2$cdn.harvest <- NULL
data_bc.long.2$noncdn.harvest <- NULL

# Function to calculate stability metrics for each moving window across time series
stb.fcn.bc.2 <- function(data_bc.2) {
  for(y in 0:(length(years.bc) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window - 1))
    temp <- data_bc.2[data_bc.2$year %in% temp.years.bc,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m.bc_con, nrow = length(temp.years.bc))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m.bc_con) {   # for each community
      
      for(t in 1:length(temp.years.bc)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$stat.area == stat.areas2[i] & temp$year == temp.years.bc[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M.bc[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m.bc_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m.bc_con
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M.bc[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M.bc[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE.bc[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m.bc_con^2) * sum(w_i_j) + (1/m.bc_con) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi.bc[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi.bc)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi.bc
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w.bc[y+1] <- u_M.bc[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])
    
  }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, phi = phi.bc, spatial_stb = spatial_stb.bc, u_M = u_M.bc, PE = PE.bc))
}

## USE FUNCTION
# Stability values for total run
stb.run.2 <- stb.fcn.bc.2(data_bc.long.2[data_bc.long.2$type == "run",])
stb.run.2

# Contributions
con.bc.2.g <- ((stb.run$gamma_w - stb.run.2$gamma_w) / stb.run$gamma_w) * 100
con.bc.2.g

# SA 3 contribution -------------------------------------------------------
data_bc.3 <- data_bc[!(data_bc$stat.area==3),]

# Store for later
stat.areas3 <- unique(data_bc.3$stat.area)

# Reshape dataset for use in function
data_bc.long.3 <- data_bc.3 %>% gather(type, abundance, c(escapement, total.harvest, total.run)) %>% 
  mutate(type = str_replace_all(type, c("escapement" = "escapement", "total.harvest" = "harvest", 
                                        "total.run" = "run")))
data_bc.long.3$cdn.harvest <- NULL
data_bc.long.3$noncdn.harvest <- NULL

# Function to calculate stability metrics for each moving window across time series
stb.fcn.bc.3 <- function(data_bc.3) {
  for(y in 0:(length(years.bc) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window - 1))
    temp <- data_bc.3[data_bc.3$year %in% temp.years.bc,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m.bc_con, nrow = length(temp.years.bc))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m.bc_con) {   # for each community
      
      for(t in 1:length(temp.years.bc)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$stat.area == stat.areas3[i] & temp$year == temp.years.bc[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M.bc[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m.bc_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m.bc_con
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M.bc[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M.bc[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE.bc[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m.bc_con^2) * sum(w_i_j) + (1/m.bc_con) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi.bc[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi.bc)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi.bc
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w.bc[y+1] <- u_M.bc[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])
    
  }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, phi = phi.bc, spatial_stb = spatial_stb.bc, u_M = u_M.bc, PE = PE.bc))
}

## USE FUNCTION
# Stability values for total run
stb.run.3 <- stb.fcn.bc.3(data_bc.long.3[data_bc.long.3$type == "run",])
stb.run.3

# Contributions
con.bc.3.g <- ((stb.run$gamma_w - stb.run.3$gamma_w) / stb.run$gamma_w) * 100
con.bc.3.g

# SA 4 contribution -------------------------------------------------------
data_bc.4 <- data_bc[!(data_bc$stat.area==4),]

# Store for later
stat.areas4 <- unique(data_bc.4$stat.area)

# Reshape dataset for use in function
data_bc.long.4 <- data_bc.4 %>% gather(type, abundance, c(escapement, total.harvest, total.run)) %>% 
  mutate(type = str_replace_all(type, c("escapement" = "escapement", "total.harvest" = "harvest", 
                                        "total.run" = "run")))
data_bc.long.4$cdn.harvest <- NULL
data_bc.long.4$noncdn.harvest <- NULL

# Function to calculate stability metrics for each moving window across time series
stb.fcn.bc.4 <- function(data_bc.4) {
  for(y in 0:(length(years.bc) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window - 1))
    temp <- data_bc.4[data_bc.4$year %in% temp.years.bc,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m.bc_con, nrow = length(temp.years.bc))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m.bc_con) {   # for each community
      
      for(t in 1:length(temp.years.bc)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$stat.area == stat.areas4[i] & temp$year == temp.years.bc[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M.bc[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m.bc_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m.bc_con
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M.bc[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M.bc[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE.bc[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m.bc_con^2) * sum(w_i_j) + (1/m.bc_con) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi.bc[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi.bc)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi.bc
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w.bc[y+1] <- u_M.bc[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])
    
  }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, phi = phi.bc, spatial_stb = spatial_stb.bc, u_M = u_M.bc, PE = PE.bc))
}

## USE FUNCTION
# Stability values for total run
stb.run.4 <- stb.fcn.bc.4(data_bc.long.4[data_bc.long.4$type == "run",])
stb.run.4

# Contributions
con.bc.4.g <- ((stb.run$gamma_w - stb.run.4$gamma_w) / stb.run$gamma_w) * 100
con.bc.4.g

# SA 5 contribution -------------------------------------------------------
data_bc.5 <- data_bc[!(data_bc$stat.area==5),]

# Store for later
stat.areas5 <- unique(data_bc.5$stat.area)

# Reshape dataset for use in function
data_bc.long.5 <- data_bc.5 %>% gather(type, abundance, c(escapement, total.harvest, total.run)) %>% 
  mutate(type = str_replace_all(type, c("escapement" = "escapement", "total.harvest" = "harvest", 
                                        "total.run" = "run")))
data_bc.long.5$cdn.harvest <- NULL
data_bc.long.5$noncdn.harvest <- NULL

# Function to calculate stability metrics for each moving window across time series
stb.fcn.bc.5 <- function(data_bc.5) {
  for(y in 0:(length(years.bc) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window - 1))
    temp <- data_bc.5[data_bc.5$year %in% temp.years.bc,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m.bc_con, nrow = length(temp.years.bc))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m.bc_con) {   # for each community
      
      for(t in 1:length(temp.years.bc)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$stat.area == stat.areas5[i] & temp$year == temp.years.bc[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M.bc[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m.bc_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m.bc_con
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M.bc[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M.bc[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE.bc[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m.bc_con^2) * sum(w_i_j) + (1/m.bc_con) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi.bc[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi.bc)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi.bc
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w.bc[y+1] <- u_M.bc[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])
    
  }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, phi = phi.bc, spatial_stb = spatial_stb.bc, u_M = u_M.bc, PE = PE.bc))
}

## USE FUNCTION
# Stability values for total run
stb.run.5 <- stb.fcn.bc.5(data_bc.long.5[data_bc.long.5$type == "run",])
stb.run.5

# Contributions
con.bc.5.g <- ((stb.run$gamma_w - stb.run.5$gamma_w) / stb.run$gamma_w) * 100
con.bc.5.g

# SA 6 contribution -------------------------------------------------------
data_bc.6 <- data_bc[!(data_bc$stat.area==6),]

# Store for later
stat.areas6 <- unique(data_bc.6$stat.area)

# Reshape dataset for use in function
data_bc.long.6 <- data_bc.6 %>% gather(type, abundance, c(escapement, total.harvest, total.run)) %>% 
  mutate(type = str_replace_all(type, c("escapement" = "escapement", "total.harvest" = "harvest", 
                                        "total.run" = "run")))
data_bc.long.6$cdn.harvest <- NULL
data_bc.long.6$noncdn.harvest <- NULL

# Function to calculate stability metrics for each moving window across time series
stb.fcn.bc.6 <- function(data_bc.6) {
  for(y in 0:(length(years.bc) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window - 1))
    temp <- data_bc.6[data_bc.6$year %in% temp.years.bc,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m.bc_con, nrow = length(temp.years.bc))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m.bc_con) {   # for each community
      
      for(t in 1:length(temp.years.bc)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$stat.area == stat.areas6[i] & temp$year == temp.years.bc[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M.bc[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m.bc_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m.bc_con
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M.bc[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M.bc[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE.bc[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m.bc_con^2) * sum(w_i_j) + (1/m.bc_con) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi.bc[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi.bc)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi.bc
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w.bc[y+1] <- u_M.bc[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])
    
  }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, phi = phi.bc, spatial_stb = spatial_stb.bc, u_M = u_M.bc, PE = PE.bc))
}

## USE FUNCTION
# Stability values for total run
stb.run.6 <- stb.fcn.bc.6(data_bc.long.6[data_bc.long.6$type == "run",])
stb.run.6

# Contributions
con.bc.6.g <- ((stb.run$gamma_w - stb.run.6$gamma_w) / stb.run$gamma_w) * 100
con.bc.6.g

# SA 7 contribution -------------------------------------------------------
data_bc.7 <- data_bc[!(data_bc$stat.area==7),]

# Store for later
stat.areas7 <- unique(data_bc.7$stat.area)

# Reshape dataset for use in function
data_bc.long.7 <- data_bc.7 %>% gather(type, abundance, c(escapement, total.harvest, total.run)) %>% 
  mutate(type = str_replace_all(type, c("escapement" = "escapement", "total.harvest" = "harvest", 
                                        "total.run" = "run")))
data_bc.long.7$cdn.harvest <- NULL
data_bc.long.7$noncdn.harvest <- NULL

# Function to calculate stability metrics for each moving window across time series
stb.fcn.bc.7 <- function(data_bc.7) {
  for(y in 0:(length(years.bc) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window - 1))
    temp <- data_bc.7[data_bc.7$year %in% temp.years.bc,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m.bc_con, nrow = length(temp.years.bc))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m.bc_con) {   # for each community
      
      for(t in 1:length(temp.years.bc)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$stat.area == stat.areas7[i] & temp$year == temp.years.bc[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M.bc[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m.bc_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m.bc_con
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M.bc[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M.bc[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE.bc[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m.bc_con^2) * sum(w_i_j) + (1/m.bc_con) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi.bc[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi.bc)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi.bc
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w.bc[y+1] <- u_M.bc[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])
    
  }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, phi = phi.bc, spatial_stb = spatial_stb.bc, u_M = u_M.bc, PE = PE.bc))
}

## USE FUNCTION
# Stability values for total run
stb.run.7 <- stb.fcn.bc.7(data_bc.long.7[data_bc.long.7$type == "run",])
stb.run.7

# Contributions
con.bc.7.g <- ((stb.run$gamma_w - stb.run.7$gamma_w) / stb.run$gamma_w) * 100
con.bc.7.g

# SA 8 contribution -------------------------------------------------------
data_bc.8 <- data_bc[!(data_bc$stat.area==8),]

# Store for later
stat.areas8 <- unique(data_bc.8$stat.area)

# Reshape dataset for use in function
data_bc.long.8 <- data_bc.8 %>% gather(type, abundance, c(escapement, total.harvest, total.run)) %>% 
  mutate(type = str_replace_all(type, c("escapement" = "escapement", "total.harvest" = "harvest", 
                                        "total.run" = "run")))
data_bc.long.8$cdn.harvest <- NULL
data_bc.long.8$noncdn.harvest <- NULL

# Function to calculate stability metrics for each moving window across time series
stb.fcn.bc.8 <- function(data_bc.8) {
  for(y in 0:(length(years.bc) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window - 1))
    temp <- data_bc.8[data_bc.8$year %in% temp.years.bc,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m.bc_con, nrow = length(temp.years.bc))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m.bc_con) {   # for each community
      
      for(t in 1:length(temp.years.bc)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$stat.area == stat.areas8[i] & temp$year == temp.years.bc[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M.bc[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m.bc_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m.bc_con
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M.bc[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M.bc[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE.bc[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m.bc_con^2) * sum(w_i_j) + (1/m.bc_con) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi.bc[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi.bc)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi.bc
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w.bc[y+1] <- u_M.bc[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])
    
  }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, phi = phi.bc, spatial_stb = spatial_stb.bc, u_M = u_M.bc, PE = PE.bc))
}

## USE FUNCTION
# Stability values for total run
stb.run.8 <- stb.fcn.bc.8(data_bc.long.8[data_bc.long.8$type == "run",])
stb.run.8

# Contributions
con.bc.8.g <- ((stb.run$gamma_w - stb.run.8$gamma_w) / stb.run$gamma_w) * 100
con.bc.8.g

# SA 9 contribution -------------------------------------------------------
data_bc.9 <- data_bc[!(data_bc$stat.area==9),]

# Store for later
stat.areas9 <- unique(data_bc.9$stat.area)

# Reshape dataset for use in function
data_bc.long.9 <- data_bc.9 %>% gather(type, abundance, c(escapement, total.harvest, total.run)) %>% 
  mutate(type = str_replace_all(type, c("escapement" = "escapement", "total.harvest" = "harvest", 
                                        "total.run" = "run")))
data_bc.long.9$cdn.harvest <- NULL
data_bc.long.9$noncdn.harvest <- NULL

# Function to calculate stability metrics for each moving window across time series
stb.fcn.bc.9 <- function(data_bc.9) {
  for(y in 0:(length(years.bc) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window - 1))
    temp <- data_bc.9[data_bc.9$year %in% temp.years.bc,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m.bc_con, nrow = length(temp.years.bc))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m.bc_con) {   # for each community
      
      for(t in 1:length(temp.years.bc)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$stat.area == stat.areas9[i] & temp$year == temp.years.bc[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M.bc[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m.bc_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m.bc_con
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M.bc[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M.bc[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE.bc[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m.bc_con^2) * sum(w_i_j) + (1/m.bc_con) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi.bc[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi.bc)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi.bc
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w.bc[y+1] <- u_M.bc[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])
    
  }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, phi = phi.bc, spatial_stb = spatial_stb.bc, u_M = u_M.bc, PE = PE.bc))
}

## USE FUNCTION
# Stability values for total run
stb.run.9 <- stb.fcn.bc.9(data_bc.long.9[data_bc.long.9$type == "run",])
stb.run.9

# Contributions
con.bc.9.g <- ((stb.run$gamma_w - stb.run.9$gamma_w) / stb.run$gamma_w) * 100
con.bc.9.g

# SA 10 contribution -------------------------------------------------------
data_bc.10 <- data_bc[!(data_bc$stat.area==10),]

# Store for later
stat.areas10 <- unique(data_bc.10$stat.area)

# Reshape dataset for use in function
data_bc.long.10 <- data_bc.10 %>% gather(type, abundance, c(escapement, total.harvest, total.run)) %>% 
  mutate(type = str_replace_all(type, c("escapement" = "escapement", "total.harvest" = "harvest", 
                                        "total.run" = "run")))
data_bc.long.10$cdn.harvest <- NULL
data_bc.long.10$noncdn.harvest <- NULL

# Function to calculate stability metrics for each moving window across time series
stb.fcn.bc.10 <- function(data_bc.10) {
  for(y in 0:(length(years.bc) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window - 1))
    temp <- data_bc.10[data_bc.10$year %in% temp.years.bc,]
    
    # Community biomass in community i at time t (sum across species)
    N_i_t <- matrix(ncol = m.bc_con, nrow = length(temp.years.bc))
    
    # Temporal mean of community biomass in community i (temporal mean = mean across the time period)
    u_i <- NULL
    
    w_i_i <- NULL
    
    for(i in 1:m.bc_con) {   # for each community
      
      for(t in 1:length(temp.years.bc)) {   # for each year
        N_i_t[t,i] <- sum(temp$abundance[temp$stat.area == stat.areas10[i] & temp$year == temp.years.bc[t]])
      }
      
      u_i[i] <- mean(N_i_t[,i])
      w_i_i[i] <- var(N_i_t[,i])
    }
    
    # Temporal covariance of community biomass between communities i and j
    w_i_j <- cov(N_i_t)
    
    # Temporal mean of metacommunity biomass
    u_M.bc[y+1] <- sum(u_i)
    
    # Regional average of temporal means of community biomasses
    u_mean <- sum(u_i) / m.bc_con
    
    # Regional average of temporal standard deviations of community biomasses
    sd_mean <- sum(sqrt(w_i_i)) / m.bc_con
    
    
    ## Covariation, spatial variation, and synchrony --------
    # still sticking with notation in W&L 2014 for most part, except spatial variability terms (CV^2 in W&L) now called SV
    
    # Coefficient of temporal variation of community biomass in community i
    CV_i <- sqrt(w_i_i) / u_i
    
    # Weighted average of coefficients of temporal variation across communities
    CV_L <- sum(sqrt(w_i_i)) / u_M.bc[y+1]
    
    # Coefficient of temporal variation of metacommunity biomass
    CV_M <- sqrt(sum(w_i_j)) / u_M.bc[y+1]
    
    # Portfolio effect
    # See Anderson et al. 2013 MEE
    PE.bc[y+1] <- mean(CV_i) / CV_M
    
    # Spatial variability: squared coefficient of spatial variation
    SV <- (sd_mean^2 - (1/m.bc_con^2) * sum(w_i_j) + (1/m.bc_con) * sum((u_i - u_mean)^2 + (sqrt(w_i_i) - sd_mean)^2)) / u_mean^2
    
    # Spatial synchrony among communities
    phi.bc[y+1] <- sum(w_i_j) / sum(sqrt(w_i_i))^2
    
    # Spatial variability that is related to spatial asynchrony (but also affected by spatial unevenness)
    SV_asyn <- CV_L^2 * (1 - phi.bc)
    
    ## Alpha, beta, and gamma variability --------
    # Alpha variability: temporal variability at community level
    alpha <- CV_L^2
    
    # Multiplicative beta variability: spatial asynchrony or the reciprocal of spatial synchrony
    beta_1 <- 1 / phi.bc
    
    # Additive beta variability: asynchrony-related spatial variability
    beta_2 <- SV_asyn
    
    # Gamma variability: temporal variability at the metacommunity scale
    gamma <- CV_M^2
    
    # Checking numbers (gamma + beta_2 = alpha; phi = CV_M^2 / CV_L^2)
    gamma + beta_2
    CV_M^2/CV_L^2
    
    ## Gamma and alpha stability --------
    # Gamma stability: stability of the metacommunity
    gamma_w.bc[y+1] <- u_M.bc[y+1] / sqrt(sum(w_i_j))
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])
    
  }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, phi = phi.bc, spatial_stb = spatial_stb.bc, u_M = u_M.bc, PE = PE.bc))
}

## USE FUNCTION
# Stability values for total run
stb.run.10 <- stb.fcn.bc.10(data_bc.long.10[data_bc.long.10$type == "run",])
stb.run.10

# Contributions
con.bc.10.g <- ((stb.run$gamma_w - stb.run.10$gamma_w) / stb.run$gamma_w) * 100
con.bc.10.g

# ### CONTRIBUTION OF HATCHERY SALMON ### -----------------------------------------
con.hatch.gamma <- ((stb.all$gamma_w - stb.nat$gamma_w) / stb.all$gamma_w) * 100

con.hatch.alpha.loc <- ((stb.all$alpha_i - stb.nat$alpha_i) / stb.all$alpha_i) * 100

# contributions of hatchery salmon to stability in each local community
con.hatch.alpha.loc <- as_tibble(con.hatch.alpha.loc) %>% 
  set_colnames(year.np) %>% 
  mutate(Region = regions) %>% 
  pivot_longer(cols = '1956':'2010',
               names_to = "Year") %>% 
  rename(contrib.alpha_i = value)

