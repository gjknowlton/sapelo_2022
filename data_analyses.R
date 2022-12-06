#### Data Analysis for Sapelo Course #### 
# Author: Garrett J. Knowlton - gjknowlton@wisc.edu

####  load packages ####

#install.packages("pacman")
library(pacman)
p_load(tidyverse, lme4, car, vegan, ggpubr)

####  read & tidy data, calculate diversity indices ####

covariates <- read.csv("G:\\My Drive\\Classes\\Sapelo\\Data\\covariates_all.csv")
comm.data <-  read.csv("G:\\My Drive\\Classes\\Sapelo\\diversity_calcs_data.csv")
total_cover <- read.csv("G:\\My Drive\\Classes\\Sapelo\\total_pct_cover.csv")

# calc shannon diversity for all 100 quads
shannon_dune <- diversity(comm.data, index = "shannon", groups = comm.data$Site_ID, equalize.groups = F)

model.df <- covariates %>% 
  mutate(shannon_diversity = shannon_dune,
         pct_cover = total_cover$total_cover_pct)

model.df.north <- model.df %>% 
  filter(n_s == "N")

model.df.south <- model.df %>% 
  filter(n_s == "S") 

 ####  question 1 ####

## Which abiotic characteristics drive plant community diversity and cover? 
## Does this vary b/w North and South ends of Nannygoat Beach? 

for(i in c(8:11)){model.df[,i] <- (model.df[,i] - mean(model.df[,i]))/sd(model.df[,i])} 

#model.df.north %>% 


#model.df.south %>% 
  

dunecommrand <- lmer(shannon_diversity ~ dist.to.strand..m. + specific_conductance 
                     + soil.moisture.pct + pct.OM + (1 | Site), data = model.df)

anova(dunecommrand)

library(lme4)
library(nlme)
mod1 <- lme(shannon_diversity ~ dist.to.strand..m. + specific_conductance 
            + soil.moisture.pct + pct.OM, random= ~1 | Site, method='ML', data = model.df)

summary(dunecommrand)

###

confint <- as.data.frame(confint(dunecommrand)) 
confint <- confint[-c(1:3),]

estimate_mod <- as.data.frame(coef(summary(dunecommrand)))
estimate_mod <- estimate_mod[-1,]

coef <- c("dist.to.strand..m.","specific_conductance", "soil.moisture.pct", "pct.OM")

confint <- cbind(confint, coef)
colnames(confint) <- c("lowerci", "upperci", "coef")
confint <- confint %>% 
  mutate(estimate = estimate_mod$Estimate,
         class = c(1:4))

#Standardized coefficients for model
#Add TotTreeCones when running PICO
#Remove Elev_m when running PIEN

ggplot(confint, aes(estimate, reorder(coef, -class))) + 
  geom_errorbar(aes(xmin = lowerci, xmax = upperci), width = 0, size = 1) +
  geom_point(aes(color = coef), size = 7) + xlim(-0.02, 0.1) + labs(x="Coefficient estimate", y=element_blank()) +
  #scale_y_discrete(labels = c('Seed source Distance\nby Northeasterliness', 'Slope (degrees)', #'Number of PICO\nwith Cones',
   #                           'Seed Source Distance (m)',
    #                          'Seed Source\nNortheasterliness', 'TMI', 'Seed Source Slope')) +#, 'Elevation (m)')) +
  theme(panel.border = element_rect(color = "black", fill="NA"), legend.position = "none", 
        axis.text=element_text(size = 23), axis.title = element_text(size = 28),
        panel.background = element_blank(), panel.grid.major = element_line(color = "grey")) +
  geom_vline(xintercept = 0, color = "grey", size = 1.5, linetype="dashed")

###

dunecommrand_north <- lmer(shannon_diversity ~ dist.to.strand..m. + specific_conductance 
                           + soil.moisture.pct + pct.OM +(1 | Site), data = model.df.north)

summary(dunecommrand_north)

dunecommrand_south <- lmer(shannon_diversity ~ dist.to.strand..m. + specific_conductance 
                           + soil.moisture.pct + pct.OM +(1 | Site), data = model.df.south)

summary(dunecommrand_south)



####  question 2 ####

## Do plant community diversity and cover or dune/soil abiotic characteristics vary b/w N & S sites?

# salinity

var.test(specific_conductance ~ n_s, covariates, 
         alternative = "two.sided")

t.test(model.df.north$specific_conductance, model.df.south$specific_conductance, var.equal = F)

t.test(log(model.df.north$specific_conductance), log(model.df.south$specific_conductance), var.equal = F)

boxplot(covariates$specific_conductance ~ covariates$n_s, outline =F)

# soil moisture - NOT CORRECT, NEED ALTERNATE TEST FOR NON-CONTINUOUS, NON-NORMAL DATA

var.test(soil.moisture.pct ~ n_s, covariates, 
         alternative = "two.sided")

t.test(model.df.north$soil.moisture.pct, model.df.south$soil.moisture.pct, var.equal = F)

t.test(log(model.df.north$specific_conductance), log(model.df.south$specific_conductance), var.equal = F)

hist(model.df.$soil.moisture.pct)

# organic material - NOT CORRECT, NEED ALTERNATE TEST FOR NON-CONTINUOUS, NON-NORMAL DATA

var.test(soil.moisture.pct ~ n_s, covariates, 
         alternative = "two.sided")

t.test(model.df.north$pct.OM, model.df.south$pct.OM, var.equal = F)

t.test(log(model.df.north$specific_conductance), log(model.df.south$specific_conductance), var.equal = F)



####  question 3 ####

## How is dune soil salinity related to distance to strand line? Does this differ b/w N & S sites? 

# all sites

ggplot(model.df, aes(x = dist.to.dune..m., y = specific_conductance)) + 
  geom_point() +
  stat_smooth(method = "lm") +
  stat_cor(label.x = 60, label.y = 1700) +
  stat_regline_equation(label.x = 60, label.y = 1625)

summary(lm(model.df$specific_conductance ~ model.df$dist.to.dune..m.))

# north site

ggplot(model.df.north, aes(x = dist.to.dune..m., y = specific_conductance)) + 
  geom_point() +
  stat_smooth(method = "lm") +
  stat_cor(label.x = 60, label.y = 400) +
  stat_regline_equation(label.x = 60, label.y = 380)

summary(lm(model.df.north$specific_conductance ~ model.df.north$dist.to.dune..m.))

# south site

ggplot(model.df.south, aes(x = dist.to.dune..m., y = specific_conductance)) + 
  geom_point() +
  stat_smooth(method = "lm") +
  stat_cor(label.x = 60, label.y = 1600) +
  stat_regline_equation(label.x = 60, label.y = 1545)

summary(lm(model.df.south$specific_conductance ~ model.df.south$dist.to.dune..m.))




#### test analyses ####

# compare diversity b/w north and south sites

ns_groups <- read.csv("G:\\My Drive\\Classes\\Sapelo\\site_df.csv")

shannon_ns_comp <- diversity(comm.data, index = "shannon", groups = comm.data$Site, equalize.groups = F)

ns_groups <- ns_groups %>% 
  mutate(shannon_index = shannon_ns_comp)

ns_groups_n <- ns_groups %>% 
  filter(n_s == "N")

ns_groups_s <- ns_groups %>% 
  filter(n_s == "S")

t.test(ns_groups_n$shannon_index, ns_groups_s$shannon_index, var.equal = F)

# cca analysis

duneCCA <- cca(comm.data ~ dist.to.strand..m. + specific_conductance + soil.moisture.pct + pct.OM, data = model.df)
duneCCA
summary(duneCCA)

duneCCAplot <- plot(duneCCA)

# is CCA significant?
anova(duneCCA)

# are axes significant?
anova(duneCCA, by = "axis")

# are vectors significant?
anova(duneCCA, by = "term", permutations = 999)
