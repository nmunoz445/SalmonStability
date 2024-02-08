#### Salmon stability - Muñoz et al ####
### 02: stability components for the Northern British Columbia metacommunity

## INPUTS ---------
datafile_bc <- "data/bc_salmondata.csv"    

## SET UP --------
# Note: we're using the notation from Wang and Loreau 2014 ...
# ... but what they call "patches" we call "communities" as in Wilcox et al. 2017
# ... so the "metacommunity" = north/central BC coast, and "local community" = stat area

#library(codyn)
#library(tidyverse)

# Read in data
data_bc <- read.csv(datafile_bc, header=T, stringsAsFactors=F, fileEncoding="UTF-8-BOM")

# Removing incomplete stat areas and years
data_bc <- data_bc[!(data_bc$year== 1954),]
data_bc <- data_bc[!(data_bc$year== 1955),]
data_bc <- data_bc[!(data_bc$year== 1956),]
data_bc <- data_bc[!(data_bc$year== 1957),]
data_bc <- data_bc[!(data_bc$year== 1958),]
data_bc <- data_bc[!(data_bc$year== 1959),]
data_bc <- data_bc[!(data_bc$year== 2013),]
data_bc <- data_bc[!(data_bc$year== 2014),]
data_bc <- data_bc[!(data_bc$year== 2015),]
data_bc <- data_bc[!(data_bc$year== 2016),]
data_bc <- data_bc[!(data_bc$year== 2017),]
data_bc <- data_bc[!(data_bc$species== "sockeye" & data_bc$stat.area ==3),] # removing sockeye in stat area 3
data_bc <- data_bc[!(data_bc$species== "pink" & data_bc$stat.area ==10),] # removing pink in stat area 10

# Store for later
stat.areas <- unique(data_bc$stat.area)
years.bc <- unique(data_bc$year)
species.bc <- unique(data_bc$species)

# Number of communities in the metacommunity
m.bc <- length(unique(data_bc$stat.area))

# Empty objects for function
alpha_i.bc <- matrix(ncol = length(years.bc) - window + 1, nrow = length(stat.areas))
alpha_w.bc <- NULL
gamma_w.bc <- NULL
phi.bc <- NULL
spatial_stb.bc <- NULL
phi_pop_s.bc <- matrix(ncol = length(years.bc) - window + 1, nrow = length(species.bc))
phi_sp_M.bc <- NULL
stb_sp_M.bc <- NULL
phi_pop_M.bc <- NULL
u_M.bc  <- NULL
PE.bc <- NULL

# Reshape data_bcset for use in function
data_bc.long <- data_bc %>% gather(type, abundance, c(escapement, total.harvest, total.run)) %>% 
  mutate(type = str_replace_all(type, c("escapement" = "escapement", "total.harvest" = "harvest", 
                                        "total.run" = "run")))
data_bc.long$cdn.harvest <- NULL
data_bc.long$noncdn.harvest <- NULL


## STABILITY FUNCTION --------
# Function to calculate stability metrics for each moving window across time series
stb.fcn.bc <- function(data_bc) {
  for(y in 0:(length(years.bc) - window)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window - 1))
    temp <- data_bc[data_bc$year %in% temp.years.bc,]
    
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
    
    # Alpha stability for each local community
    for(i in 1:m.bc) {   # for each local community
      alpha_i.bc[i,y+1] <- (sqrt(w_i_i[i]) / u_i[i])^-1
    }
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])

    ## Species synchrony and stability --------
    # Species synchrony within each community - different species responding differently through time.
    # Using Loreau and Mazancourt's (2008) metric (also W&L and Wilcox), calculated using the codyn package (Hallatt et al. 2016 MEE)
    phi_sp_i <- synchrony(df = temp, time.var = "year", species.var = "species", abundance.var = "abundance", replicate.var = "stat.area",
                          metric = "Loreau")[,2]
    
    # Species synchrony at the metacommunity level - weighted average across communities
    phi_sp_M.bc[y+1] <- sum(phi_sp_i * u_i / sum(u_i))
    
    
    ## Species stability --------
    # Biomass for species s in community i at time t
    # Order is chum > pink > sockeye
    N_i_t_sp <- list()
    N_i_t_sp[[1]] <- N_i_t_sp[[2]] <- N_i_t_sp[[3]] <- matrix(ncol = m.bc, nrow = length(temp.years.bc))
    
    # Temporal mean and variance of biomass for species s in patch i
    u_i_sp <- matrix(ncol = m.bc, nrow = length(species.bc))
    w_i_i_sp <- matrix(ncol = m.bc, nrow = length(species.bc))
    
    # Stability for species s in patch i
    stb_sp_s_i <- matrix(ncol = m.bc, nrow=length(species.bc))
    
    for(s in 1:length(species.bc)) {
      for(i in 1:m.bc) {   # for each patch
        
        for(t in 1:length(temp.years.bc)) {   # for each year 
          
          N_i_t_sp[[s]][t,i] <- sum(data_bc$abundance[data_bc$stat.area == stat.areas[i] & data_bc$year == temp.years.bc[t] & data_bc$species == species.bc[s]])
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
    stb_sp_M.bc[y+1] <- sum(stb_sp_i * u_i / sum(u_i))
    
    
    ## Population synchrony --------
    # Population synchrony for each species: the synchrony of "populations" within a species
    # i.e., The degree that a species' abundance through time within one region aligns with its abundance in different regions
    phi_pop_s.bc[,y+1] <- synchrony(df = temp, time.var = "year", 
                                    species.var = "stat.area", 
                                    abundance.var = "abundance", 
                                    replicate.var = "species", 
                                    metric = "Loreau")[,2]
    
    phi_pop_M.bc[y+1] <- sum(phi_pop_s.bc[,y+1] * rowSums(u_i_sp) / sum(u_i))
    
    
  }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, alpha_i = alpha_i.bc, phi = phi.bc, spatial_stb = spatial_stb.bc, 
              phi_sp_M = phi_sp_M.bc, stb_sp_M = stb_sp_M.bc, phi_pop_s = phi_pop_s.bc, phi_pop_M = phi_pop_M.bc, u_M = u_M.bc, PE = PE.bc))
}


## USE FUNCTION --------
# Stability values for total run
stb.run <- stb.fcn.bc(data_bc.long[data_bc.long$type == "run",])
stb.run

# Stability values for escapement
stb.esc <- stb.fcn.bc(data_bc.long[data_bc.long$type == "escapement",])
stb.esc

# Stability values for catch
# zero values (in later years) make some of these values Nan
stb.catch <- stb.fcn.bc(data_bc.long[data_bc.long$type == "harvest",])
stb.catch


### 8 AND 12 YEAR ROLLING WINDOWS -------------------------------------------
# 8 year window
alpha_i.bc8 <- matrix(ncol = length(years.bc) - window8 + 1, nrow = length(stat.areas))
phi_pop_s.bc8 <- matrix(ncol = length(years.bc) - window8 + 1, nrow = length(species.bc))

stb.fcn.bc8 <- function(data_bc) {
  for(y in 0:(length(years.bc) - window8)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window8 - 1))
    temp <- data_bc[data_bc$year %in% temp.years.bc,]
    
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
    
    # Alpha stability for each local community
    for(i in 1:m.bc) {   # for each local community
      alpha_i.bc8[i,y+1] <- (sqrt(w_i_i[i]) / u_i[i])^-1
    }
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])
    
    ## Species synchrony and stability --------
    # Species synchrony within each community - different species responding differently through time.
    # Using Loreau and Mazancourt's (2008) metric (also W&L and Wilcox), calculated using the codyn package (Hallatt et al. 2016 MEE)
    phi_sp_i <- synchrony(df = temp, time.var = "year", species.var = "species", abundance.var = "abundance", replicate.var = "stat.area",
                          metric = "Loreau")[,2]
    
    # Species synchrony at the metacommunity level - weighted average across communities
    phi_sp_M.bc[y+1] <- sum(phi_sp_i * u_i / sum(u_i))
    
    
    ## Species stability --------
    # Biomass for species s in community i at time t
    # Order is chum > pink > sockeye
    N_i_t_sp <- list()
    N_i_t_sp[[1]] <- N_i_t_sp[[2]] <- N_i_t_sp[[3]] <- matrix(ncol = m.bc, nrow = length(temp.years.bc))
    
    # Temporal mean and variance of biomass for species s in patch i
    u_i_sp <- matrix(ncol = m.bc, nrow = length(species.bc))
    w_i_i_sp <- matrix(ncol = m.bc, nrow = length(species.bc))
    
    # Stability for species s in patch i
    stb_sp_s_i <- matrix(ncol = m.bc, nrow=length(species.bc))
    
    for(s in 1:length(species.bc)) {
      for(i in 1:m.bc) {   # for each patch
        
        for(t in 1:length(temp.years.bc)) {   # for each year 
          
          N_i_t_sp[[s]][t,i] <- sum(data_bc$abundance[data_bc$stat.area == stat.areas[i] & data_bc$year == temp.years.bc[t] & data_bc$species == species.bc[s]])
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
    stb_sp_M.bc[y+1] <- sum(stb_sp_i * u_i / sum(u_i))
    
    
    ## Population synchrony --------
    # Population synchrony for each species: the synchrony of "populations" within a species
    # i.e., The degree that a species' abundance through time within one region aligns with its abundance in different regions
    phi_pop_s.bc8[,y+1] <- synchrony(df = temp, time.var = "year", 
                                     species.var = "stat.area", 
                                     abundance.var = "abundance", 
                                     replicate.var = "species", 
                                     metric = "Loreau")[,2]
    
    phi_pop_M.bc[y+1] <- sum(phi_pop_s.bc8[,y+1] * rowSums(u_i_sp) / sum(u_i))
    
    
  }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, alpha_i = alpha_i.bc8, phi = phi.bc, spatial_stb = spatial_stb.bc, 
              phi_sp_M = phi_sp_M.bc, stb_sp_M = stb_sp_M.bc, phi_pop_s = phi_pop_s.bc8, phi_pop_M = phi_pop_M.bc, u_M = u_M.bc, PE = PE.bc))
}


## USE FUNCTION 
# Stability values for total run
stb.run8 <- stb.fcn.bc8(data_bc.long[data_bc.long$type == "run",])
stb.run8

# 12 year window
alpha_i.bc12 <- matrix(ncol = length(years.bc) - window12 + 1, nrow = length(stat.areas))
phi_pop_s.bc12 <- matrix(ncol = length(years.bc) - window12 + 1, nrow = length(species.bc))

stb.fcn.bc12 <- function(data_bc) {
  for(y in 0:(length(years.bc) - window12)) {   # for each window
    ## Initial parameters --------
    # Subset data to just include years in this window
    temp.years.bc <- c((min(years.bc) + y):(min(years.bc) + y + window12 - 1))
    temp <- data_bc[data_bc$year %in% temp.years.bc,]
    
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
    
    # Alpha stability for each local community
    for(i in 1:m.bc) {   # for each local community
      alpha_i.bc12[i,y+1] <- (sqrt(w_i_i[i]) / u_i[i])^-1
    }
    
    # Alpha stability: averaged stability of the local communities
    # If any element in u_i == 0, need to remove from this otherwise
    alpha_w.bc[y+1] <- sum((u_i/u_M.bc[y+1]) * ((sqrt(w_i_i))/u_i))^-1
    
    ## Spatial synchrony and stabilization --------
    # Spatial synchrony (phi) - similarity of temporal fluctuations of different communities - same as W&L
    # phi[y+1]
    
    # Spatial stabilization - the amount of stability enhanced when moving from the community to metacommunity level
    spatial_stb.bc[y+1] <- sqrt(1/phi.bc[y+1])
    
    ## Species synchrony and stability --------
    # Species synchrony within each community - different species responding differently through time.
    # Using Loreau and Mazancourt's (2008) metric (also W&L and Wilcox), calculated using the codyn package (Hallatt et al. 2016 MEE)
    phi_sp_i <- synchrony(df = temp, time.var = "year", species.var = "species", abundance.var = "abundance", replicate.var = "stat.area",
                          metric = "Loreau")[,2]
    
    # Species synchrony at the metacommunity level - weighted average across communities
    phi_sp_M.bc[y+1] <- sum(phi_sp_i * u_i / sum(u_i))
    
    
    ## Species stability --------
    # Biomass for species s in community i at time t
    # Order is chum > pink > sockeye
    N_i_t_sp <- list()
    N_i_t_sp[[1]] <- N_i_t_sp[[2]] <- N_i_t_sp[[3]] <- matrix(ncol = m.bc, nrow = length(temp.years.bc))
    
    # Temporal mean and variance of biomass for species s in patch i
    u_i_sp <- matrix(ncol = m.bc, nrow = length(species.bc))
    w_i_i_sp <- matrix(ncol = m.bc, nrow = length(species.bc))
    
    # Stability for species s in patch i
    stb_sp_s_i <- matrix(ncol = m.bc, nrow=length(species.bc))
    
    for(s in 1:length(species.bc)) {
      for(i in 1:m.bc) {   # for each patch
        
        for(t in 1:length(temp.years.bc)) {   # for each year 
          
          N_i_t_sp[[s]][t,i] <- sum(data_bc$abundance[data_bc$stat.area == stat.areas[i] & data_bc$year == temp.years.bc[t] & data_bc$species == species.bc[s]])
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
    stb_sp_M.bc[y+1] <- sum(stb_sp_i * u_i / sum(u_i))
    
    
    ## Population synchrony --------
    # Population synchrony for each species: the synchrony of "populations" within a species
    # i.e., The degree that a species' abundance through time within one region aligns with its abundance in different regions
    phi_pop_s.bc12[,y+1] <- synchrony(df = temp, time.var = "year", 
                                      species.var = "stat.area", 
                                      abundance.var = "abundance", 
                                      replicate.var = "species", 
                                      metric = "Loreau")[,2]
    
    phi_pop_M.bc[y+1] <- sum(phi_pop_s.bc12[,y+1] * rowSums(u_i_sp) / sum(u_i))
    
    
  }
  
  return(list(gamma_w = gamma_w.bc, alpha_w = alpha_w.bc, alpha_i = alpha_i.bc12, phi = phi.bc, spatial_stb = spatial_stb.bc, 
              phi_sp_M = phi_sp_M.bc, stb_sp_M = stb_sp_M.bc, phi_pop_s = phi_pop_s.bc12, phi_pop_M = phi_pop_M.bc, u_M = u_M.bc, PE = PE.bc))
}


## USE FUNCTION
# Stability values for total run
stb.run12 <- stb.fcn.bc12(data_bc.long[data_bc.long$type == "run",])
stb.run12

## SUMMARY STATISTICS -----------------------------------
# Annual abundance
Abundance_allbc <- data_bc.long %>% # all salmon. by year
  filter(type == "run") %>% 
  group_by(year) %>% 
  summarise(Tot.Abundance = sum(abundance))

Abundance_spbc <- data_bc.long %>% # all salmon, by year and species
  filter(type == "run") %>% 
  group_by(year, species) %>% 
  summarise(Tot.Abundance = sum(abundance))

# Stability components: summary statistics
mean(stb.run$gamma_w) # metacommunity stability
sd(stb.run$gamma_w)

mean(stb.run$spatial_stb) # spatial stabilization
sd(stb.run$spatial_stb)
max(stb.run$spatial_stb)
min(stb.run$spatial_stb)