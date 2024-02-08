#### Salmon stability - Muñoz et al ####
### 01: stability components for the North Pacific Ocean metacommunity

## START ---------
rm(list=ls())
graphics.off()


## INPUTS ---------
dir.gen <- "C:/Users/nmuno/Documents/Manuscripts/Salmon stability"
dir.out <- "C:/Users/nmuno/Documents/Manuscripts/Salmon stability/plots"         # plot directory
datafile <- "data/abundance_returns_GEN.csv"    # datafile
window <- 10                              # number of years in rolling window


## SET UP --------
# Note: we're using the notation from Wang and Loreau 2014 ...
# ... but what they call "patches" we call "communities" as in Wilcox et al. 2017
# ... so the "metacommunity" = ocean basin, and "local community" = country/region

library(codyn)
library(tidyverse)

# Read in data
setwd(dir.gen)
data <- read.csv(datafile, header=T, stringsAsFactors=F)

# change abundance from millions of individuals to individuals
data$abundance.all <- data$abundance.all*1000000
data$abundance.nat <- data$abundance.nat*1000000
data$abundance.hatch <- data$abundance.hatch*1000000

# Remove Korea from dataset
data <- data[data$region != "korea",]

# Store for later
regions <- unique(data$region)
years <- unique(data$year)
species <- unique(data$species)

# Number of local communities in the metacommunity
m <- length(unique(data$region))

# Empty objects for function
alpha_i <- matrix(ncol = length(years) - window + 1, nrow = length(regions))
alpha_w <- NULL
gamma_w <- NULL
phi <- NULL
spatial_stb <- NULL
phi_pop_s <- matrix(ncol = length(years) - window + 1, nrow = length(species))
phi_sp_M <- NULL
stb_sp_M <- NULL
phi_pop_M <- NULL
u_M  <- NULL
PE <- NULL

# Reshape dataset for use in function
data.long <- data %>% gather(type, abundance, c(abundance.all, abundance.hatch, abundance.nat)) %>% 
  mutate(type = str_replace_all(type, c("abundance.all" = "all", "abundance.nat" = "nat", 
                                        "abundance.hatch" = "hatch")))


## STABILITY FUNCTION --------
# Function to calculate stability metrics for each moving window across time series
stb.fcn <- function(data) {
  for(y in 0:(length(years) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window - 1))
    temp <- data[data$year %in% temp.years,]
    
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
    
    # Alpha stability for each local community
    for(i in 1:m) {   # for each local community
      alpha_i[i,y+1] <- (sqrt(w_i_i[i]) / u_i[i])^-1
    }
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
    ## Species synchrony and stability --------
    # Species synchrony within each community - different species responding differently through time.
    # Using Loreau and Mazancourt's (2008) metric (also W&L and Wilcox), calculated using the codyn package (Hallatt et al. 2016 MEE)
    phi_sp_i <- synchrony(df = temp, time.var = "year", species.var = "species", abundance.var = "abundance", replicate.var = "region",
                          metric = "Loreau")[,2]
    
    # Species synchrony at the metacommunity level - weighted average across communities
    phi_sp_M[y+1] <- sum(phi_sp_i * u_i / sum(u_i))
    
    
    ## Species stability --------
    # Biomass for species s in community i at time t
    # Order is chum > pink > sockeye
    N_i_t_sp <- list()
    N_i_t_sp[[1]] <- N_i_t_sp[[2]] <- N_i_t_sp[[3]] <- matrix(ncol = m, nrow = length(temp.years))
    
    # Temporal mean and variance of biomass for species s in patch i
    u_i_sp <- matrix(ncol = m, nrow = length(species))
    w_i_i_sp <- matrix(ncol = m, nrow = length(species))
    
    # Stability for species s in patch i
    stb_sp_s_i <- matrix(ncol = m, nrow=length(species))
    
    for(s in 1:length(species)) {
      for(i in 1:m) {   # for each patch
        
        for(t in 1:length(temp.years)) {   # for each year 
          
          N_i_t_sp[[s]][t,i] <- sum(data$abundance[data$region == regions[i] & data$year == temp.years[t] & data$species == species[s]])
        }
        
        u_i_sp[s,i] <- mean(N_i_t_sp[[s]][,i])
        w_i_i_sp[s,i] <- var(N_i_t_sp[[s]][,i])
        
        # From equation 5 in Wilcox et al.: weighted average of species stability across communities
        stb_sp_s_i[s,i] <- sum((u_i_sp[s,i]/u_i[i]) * (sqrt(w_i_i_sp[s,i])/u_i_sp[s,i]))
      }
    }
    
    # Convert NA for Japan's sockeye to 0
    stb_sp_s_i[is.na(stb_sp_s_i)] <- 0
    
    # Species stability for each community
    stb_sp_i <- colSums(stb_sp_s_i)^-1
    stb_sp_i
    
    # Species stability at metacommunity level - weighted average across communities
    stb_sp_M[y+1] <- sum(stb_sp_i * u_i / sum(u_i))
    
    
    ## Population synchrony --------
    # Population synchrony for each species: the synchrony of "populations" within a species
    # i.e., The degree that a species' abundance through time within one region aligns with its abundance in different regions
    phi_pop_s[,y+1] <- synchrony(df = temp, time.var = "year", 
                                 species.var = "region", 
                                 abundance.var = "abundance", 
                                 replicate.var = "species", 
                                 metric = "Loreau")[,2]
    
    phi_pop_M[y+1] <- sum(phi_pop_s[,y+1] * rowSums(u_i_sp) / sum(u_i))
    
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, alpha_i = alpha_i, phi = phi, spatial_stb = spatial_stb, 
              phi_sp_M = phi_sp_M, stb_sp_M = stb_sp_M, phi_pop_s = phi_pop_s, phi_pop_M = phi_pop_M, u_M = u_M, PE = PE))
}


## USE FUNCTION --------
# Stability values for all fish (hatchery and natural)
stb.all <- stb.fcn(data.long[data.long$type == "all",])
stb.all

# Stability values for only natural fish
stb.nat <- stb.fcn(data.long[data.long$type == "nat",])
stb.nat

# Stability values for only hatchery fish
# Much of this can't be calculated because so many zeros
stb.hatch <- stb.fcn(data.long[data.long$type == "hatch",])
stb.hatch

### 8 and 12 year windows ### -----------------------------------------------
## 8 year window
window8 <- 8

alpha_i8 <- matrix(ncol = length(years) - window8 + 1, nrow = length(regions))
phi_pop_s8 <- matrix(ncol = length(years) - window8 + 1, nrow = length(species))

stb.fcn8 <- function(data) {
  for(y in 0:(length(years) - window8)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window8 - 1))
    temp <- data[data$year %in% temp.years,]
    
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
    
    # Alpha stability for each local community
    for(i in 1:m) {   # for each local community
      alpha_i8[i,y+1] <- (sqrt(w_i_i[i]) / u_i[i])^-1
    }
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])
    
    ## Species synchrony and stability --------
    # Species synchrony within each community - different species responding differently through time.
    # Using Loreau and Mazancourt's (2008) metric (also W&L and Wilcox), calculated using the codyn package (Hallatt et al. 2016 MEE)
    phi_sp_i <- synchrony(df = temp, time.var = "year", species.var = "species", abundance.var = "abundance", replicate.var = "region",
                          metric = "Loreau")[,2]
    
    # Species synchrony at the metacommunity level - weighted average across communities
    phi_sp_M[y+1] <- sum(phi_sp_i * u_i / sum(u_i))
    
    
    ## Species stability --------
    # Biomass for species s in community i at time t
    # Order is chum > pink > sockeye
    N_i_t_sp <- list()
    N_i_t_sp[[1]] <- N_i_t_sp[[2]] <- N_i_t_sp[[3]] <- matrix(ncol = m, nrow = length(temp.years))
    
    # Temporal mean and variance of biomass for species s in patch i
    u_i_sp <- matrix(ncol = m, nrow = length(species))
    w_i_i_sp <- matrix(ncol = m, nrow = length(species))
    
    # Stability for species s in patch i
    stb_sp_s_i <- matrix(ncol = m, nrow=length(species))
    
    for(s in 1:length(species)) {
      for(i in 1:m) {   # for each patch
        
        for(t in 1:length(temp.years)) {   # for each year 
          
          N_i_t_sp[[s]][t,i] <- sum(data$abundance[data$region == regions[i] & data$year == temp.years[t] & data$species == species[s]])
        }
        
        u_i_sp[s,i] <- mean(N_i_t_sp[[s]][,i])
        w_i_i_sp[s,i] <- var(N_i_t_sp[[s]][,i])
        
        # From equation 5 in Wilcox et al.: weighted average of species stability across communities
        stb_sp_s_i[s,i] <- sum((u_i_sp[s,i]/u_i[i]) * (sqrt(w_i_i_sp[s,i])/u_i_sp[s,i]))
      }
    }
    
    # Convert NA for Japan's sockeye to 0
    stb_sp_s_i[is.na(stb_sp_s_i)] <- 0
    
    # Species stability for each community
    stb_sp_i <- colSums(stb_sp_s_i)^-1
    stb_sp_i
    
    # Species stability at metacommunity level - weighted average across communities
    stb_sp_M[y+1] <- sum(stb_sp_i * u_i / sum(u_i))
    
    
    ## Population synchrony --------
    # Population synchrony for each species: the synchrony of "populations" within a species
    # i.e., The degree that a species' abundance through time within one region aligns with its abundance in different regions
    phi_pop_s8[,y+1] <- synchrony(df = temp, time.var = "year", 
                                  species.var = "region", 
                                  abundance.var = "abundance", 
                                  replicate.var = "species", 
                                  metric = "Loreau")[,2]
    
    phi_pop_M[y+1] <- sum(phi_pop_s8[,y+1] * rowSums(u_i_sp) / sum(u_i))
    
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, alpha_i = alpha_i8, phi = phi, spatial_stb = spatial_stb, 
              phi_sp_M = phi_sp_M, stb_sp_M = stb_sp_M, phi_pop_s = phi_pop_s8, phi_pop_M = phi_pop_M, u_M = u_M, PE = PE))
}


# USE FUNCTION
# Stability values for all fish (hatchery and natural)
stb.all8 <- stb.fcn8(data.long[data.long$type == "all",])
stb.all8

## 12 year window
window12 <- 12

alpha_i12 <- matrix(ncol = length(years) - window12 + 1, nrow = length(regions))
phi_pop_s12 <- matrix(ncol = length(years) - window12 + 1, nrow = length(species))

stb.fcn12 <- function(data) {
  for(y in 0:(length(years) - window12)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years <- c((min(years) + y):(min(years) + y + window12 - 1))
    temp <- data[data$year %in% temp.years,]
    
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
    
    # Alpha stability for each local community
    for(i in 1:m) {   # for each local community
      alpha_i12[i,y+1] <- (sqrt(w_i_i[i]) / u_i[i])^-1
    }
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb[y+1] <- sqrt(1/phi[y+1])

    ## Species synchrony and stability --------
    # Species synchrony within each community - different species responding differently through time.
    # Using Loreau and Mazancourt's (2008) metric (also W&L and Wilcox), calculated using the codyn package (Hallatt et al. 2016 MEE)
    phi_sp_i <- synchrony(df = temp, time.var = "year", species.var = "species", abundance.var = "abundance", replicate.var = "region",
                          metric = "Loreau")[,2]
    
    # Species synchrony at the metacommunity level - weighted average across communities
    phi_sp_M[y+1] <- sum(phi_sp_i * u_i / sum(u_i))
    
    
    ## Species stability --------
    # Biomass for species s in community i at time t
    # Order is chum > pink > sockeye
    N_i_t_sp <- list()
    N_i_t_sp[[1]] <- N_i_t_sp[[2]] <- N_i_t_sp[[3]] <- matrix(ncol = m, nrow = length(temp.years))
    
    # Temporal mean and variance of biomass for species s in patch i
    u_i_sp <- matrix(ncol = m, nrow = length(species))
    w_i_i_sp <- matrix(ncol = m, nrow = length(species))
    
    # Stability for species s in patch i
    stb_sp_s_i <- matrix(ncol = m, nrow=length(species))
    
    for(s in 1:length(species)) {
      for(i in 1:m) {   # for each patch
        
        for(t in 1:length(temp.years)) {   # for each year 
          
          N_i_t_sp[[s]][t,i] <- sum(data$abundance[data$region == regions[i] & data$year == temp.years[t] & data$species == species[s]])
        }
        
        u_i_sp[s,i] <- mean(N_i_t_sp[[s]][,i])
        w_i_i_sp[s,i] <- var(N_i_t_sp[[s]][,i])
        
        # From equation 5 in Wilcox et al.: weighted average of species stability across communities
        stb_sp_s_i[s,i] <- sum((u_i_sp[s,i]/u_i[i]) * (sqrt(w_i_i_sp[s,i])/u_i_sp[s,i]))
      }
    }
    
    # Convert NA for Japan's sockeye to 0
    stb_sp_s_i[is.na(stb_sp_s_i)] <- 0
    
    # Species stability for each community
    stb_sp_i <- colSums(stb_sp_s_i)^-1
    stb_sp_i
    
    # Species stability at metacommunity level - weighted average across communities
    stb_sp_M[y+1] <- sum(stb_sp_i * u_i / sum(u_i))
    
    
    ## Population synchrony --------
    # Population synchrony for each species: the synchrony of "populations" within a species
    # i.e., The degree that a species' abundance through time within one region aligns with its abundance in different regions
    phi_pop_s12[,y+1] <- synchrony(df = temp, time.var = "year", 
                                   species.var = "region", 
                                   abundance.var = "abundance", 
                                   replicate.var = "species", 
                                   metric = "Loreau")[,2]
    
    phi_pop_M[y+1] <- sum(phi_pop_s12[,y+1] * rowSums(u_i_sp) / sum(u_i))
    
    
  }
  
  return(list(gamma_w = gamma_w, alpha_w = alpha_w, alpha_i = alpha_i12, phi = phi, spatial_stb = spatial_stb, 
              phi_sp_M = phi_sp_M, stb_sp_M = stb_sp_M, phi_pop_s = phi_pop_s12, phi_pop_M = phi_pop_M, u_M = u_M, PE = PE))
}


## USE FUNCTION
# Stability values for all fish (hatchery and natural)
stb.all12 <- stb.fcn12(data.long[data.long$type == "all",])
stb.all12

### SUMMARY STATISTICS ------------------------------------------------------
# annual abundance


Abundance_or <- data.long %>% # by origin (hatchery vs natural) and year
  filter(type %in% c("nat","hatch")) %>% 
  group_by(year, type) %>% 
  summarise(Tot.Abundance = sum(abundance))

Abundance_or_sp <- data.long %>% # by origin (hatchery vs natural), year and species
  filter(type %in% c("nat","hatch")) %>% 
  group_by(year, type, species) %>% 
  summarise(Tot.Abundance = sum(abundance))


# Summary statistics: stability components
mean(stb.all$gamma_w) # metacommunity stability
sd(stb.all$gamma_w)

mean(stb.all$spatial_stb) # spatial stabilization
sd(stb.all$spatial_stb)
max(stb.all$spatial_stb)
min(stb.all$spatial_stb)
