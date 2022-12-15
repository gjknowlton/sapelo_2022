#### Data Analysis for Sapelo Course #### 
# Author: Garrett J. Knowlton - gjknowlton@wisc.edu

####  load packages ####

#install.packages("compensateR")
library(pacman)
p_load(tidyverse, lme4, car, vegan, ggpubr, compensateR)

library(parameters)
library(see)
library(sjPlot)
library(ggplot2)
library(report)


####  read & tidy data, calculate diversity indices ####

covariates <- read.csv("data/covariates_all.csv")
comm.data <-  read.csv("data/diversity_calcs_data.csv")
total_cover <- read.csv("data/total_pct_cover.csv")

# Convert specific conductance to salinity (ppt)

salinity <- function(salinity){
  sc <- c(salinity)
  
  R <- sc/53087
  
  k1 <- 0.0120
  k2 <- -0.2174
  k3 <- 25.3283
  k4 <- 13.7714
  k5 <- -6.4788
  k6 <- 2.5842
  
  salinity <- k1 + (k2 * (R^(1/2))) + (k3 * R) + (k4 * (R^(3/2))) + (k5 *( R^2)) + (k6 * (R^(5/2)))
  
  return(salinity)
}

covariates <- covariates %>% 
  mutate(salinity = (salinity(covariates$salinity)))

comm.data.north <- comm.data %>% 
  filter(Site_ID > 600)

comm.data.south <- comm.data %>% 
  filter(Site_ID < 600)

# calc shannon diversity for all 100 quads
shannon_dune <- diversity(comm.data, index = "shannon", groups = comm.data$Site_ID, equalize.groups = F)

shannon_dune.s <- diversity(comm.data.north, index = "shannon", groups = comm.data.north$Site_ID, equalize.groups = F)
shannon_dune.n <- diversity(comm.data.south, index = "shannon", groups = comm.data.south$Site_ID, equalize.groups = F)


model.df <- covariates %>% 
  mutate(shannon_diversity = shannon_dune,
         pct_cover = total_cover$total_cover_pct)

for(i in c(8:11)){model.df[,i] <- (model.df[,i] - mean(model.df[,i]))/sd(model.df[,i])} 


model.df.north <- model.df %>% 
  filter(n_s == "N")

model.df.south <- model.df %>% 
  filter(n_s == "S") 

 ####  question 1 ####

## Which abiotic characteristics drive plant community diversity and cover? 
## Does this vary b/w North and South ends of Nannygoat Beach? 


#model.df.north %>% 


#model.df.south %>% 
  

dunecommrand <- lmer(shannon_diversity ~ salinity + dist.to.strand..m. 
                     + soil.moisture.pct + pct.OM + (1 | Site), data = model.df)

report.dunemod <- report(dunecommrand)

#anova(dunecommrand)
#
#library(lme4)
#library(nlme)
#mod1 <- lme(shannon_diversity ~ dist.to.strand..m. + specific_conductance 
#            + soil.moisture.pct + pct.OM, random= ~1 | Site, method='ML', data = model.df)
#
#summary(dunecommrand)

###

confint <- as.data.frame(confint(dunecommrand)) 
confint <- confint[-c(1:3),]

estimate_mod <- as.data.frame(coef(summary(dunecommrand)))
estimate_mod <- estimate_mod[-1,]

coef <- c("dist.to.strand..m.","salinity", "soil.moisture.pct", "pct.OM")

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

dunecommrand_north <- lmer(shannon_diversity ~ dist.to.strand..m. + salinity 
                           + soil.moisture.pct + pct.OM +(1 | Site), data = model.df.north)

summary(dunecommrand_north)

dunecommrand_south <- lmer(shannon_diversity ~ dist.to.strand..m. + salinity 
                           + soil.moisture.pct + pct.OM +(1 | Site), data = model.df.south)

summary(dunecommrand_south)



####  question 2 ####

## Do plant community diversity and cover or dune/soil abiotic characteristics vary b/w N & S sites?

# salinity

var.test(salinity ~ n_s, covariates, 
         alternative = "two.sided")

ttest_salt <- t.test(model.df.north$salinity, model.df.south$salinity, var.equal = F)
report(ttest_salt)
t.test(log(model.df.north$salinity), log(model.df.south$salinity), var.equal = F)

boxplot(covariates$salinity ~ covariates$n_s, outline =F)
boxplot(model.df$salinity ~ model.df$n_s, outline =F)


ttestplot_salt<-ggplot(covariates, aes(x=n_s, y=salinity, fill=n_s)) +
  geom_boxplot() +
  labs(x="Site Location",
       y ="Practical Salinity Unit")
ttestplot_salt

# soil moisture - NOT CORRECT, NEED ALTERNATE TEST FOR NON-CONTINUOUS, NON-NORMAL DATA

var.test(soil.moisture.pct ~ n_s, covariates, 
         alternative = "two.sided")

ttest_sm <- t.test(model.df.north$soil.moisture.pct, model.df.south$soil.moisture.pct, var.equal = F)
report(ttest_sm)
t.test(log(model.df.north$salinity), log(model.df.south$salinity), var.equal = F)

hist(model.df.$soil.moisture.pct)

boxplot(model.df$soil.moisture.pct ~ model.df$n_s, outline =T)

ttestplot_SM<-ggplot(covariates, aes(x=n_s, y=soil.moisture.pct, fill=n_s)) +
  geom_boxplot() +
  labs(x="Site Location",
       y ="Soil Moisture %")
ttestplot_SM

# organic material - NOT CORRECT, NEED ALTERNATE TEST FOR NON-CONTINUOUS, NON-NORMAL DATA

var.test(soil.moisture.pct ~ n_s, covariates, 
         alternative = "two.sided")

ttest_om <- t.test(model.df.north$pct.OM, model.df.south$pct.OM, var.equal = F)
report(ttest_om)

boxplot(model.df$pct.OM ~ model.df$n_s, outline =F)

ttestplot_OM<-ggplot(covariates, aes(x=n_s, y=pct.OM, fill=n_s)) +
  geom_boxplot() +
  labs(x="Site Location",
       y ="Organic Matter %")
ttestplot_OM


####  question 3 ####

## How is dune soil salinity related to distance to strand line? Does this differ b/w N & S sites? 

# all sites

ggplot(model.df, aes(x = dist.to.dune..m., y = salinity)) + 
  geom_point() +
  stat_smooth(method = "lm") +
  stat_cor(label.x = 60, label.y = 1) +
  stat_regline_equation(label.x = 60, label.y = 1.1) +
  labs(x = "Distance to Dune (m)",
       y = "Practical Salinity Unit")

lm_all <- lm(model.df$salinity ~ model.df$dist.to.dune..m.)
report(lm_all)

summary(lm(model.df$salinity ~ model.df$dist.to.dune..m.))

# north site

ggplot(model.df.north, aes(x = dist.to.dune..m., y = salinity)) + 
  geom_point() +
  stat_smooth(method = "lm") +
  stat_cor(label.x = 60, label.y = 0.3) +
  stat_regline_equation(label.x = 60, label.y = 0.32) +
  labs(x = "Distance to Dune (m)",
       y = "Practical Salinity Unit")

lm_n <- lm(model.df.north$salinity ~ model.df.north$dist.to.dune..m.)
report(lm_n)
summary(lm(model.df.north$salinity ~ model.df.north$dist.to.dune..m.))

# south site

ggplot(model.df.south, aes(x = dist.to.dune..m., y = salinity)) + 
  geom_point() +
  stat_smooth(method = "lm") +
  stat_cor(label.x = 60, label.y = 0.7) +
  stat_regline_equation(label.x = 60, label.y = 0.75) +
  labs(x = "Distance to Dune (m)",
       y = "Practical Salinity Unit")

  

lm_s <- lm(model.df.south$salinity ~ model.df.south$dist.to.dune..m.)
report(lm_s)
summary(lm(model.df.south$salinity ~ model.df.south$dist.to.dune..m.))




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

report(t.test(ns_groups_n$shannon_index, ns_groups_s$shannon_index, var.equal = F))

t.test(ns_groups_n$shannon_index, ns_groups_s$shannon_index, var.equal = T)

ttestplot_div<-ggplot(ns_groups, aes(x=n_s, y=shannon_index, fill=n_s)) +
  geom_boxplot() +
  labs(x = "Site Location",
       y = "Shannon Diversity Index")
ttestplot_div

boxplot(ns_groups$shannon_index ~ ns_groups$n_s)
boxplot(covariates$salinity ~ covariates$n_s, outline =T)




# cca analysis

duneCCA <- cca(comm.data ~ dist.to.strand..m. + salinity + soil.moisture.pct + pct.OM, data = model.df)
duneCCA
summary(duneCCA)

duneCCAplot <- plot(duneCCA)
duneCCAplot

x <- plot(duneCCA)
geom_text_repel(data = sp, mapping = aes(x = CCA1, y = CCA2, label = Label))

sp <- fortify(duneCCA, display = "sp")


plot(duneCCA, display=c("sp","lc","cn"), scaling=2, xlab="CCA1 ()", ylab="CCA2 ()")



# is CCA significant?
anova(duneCCA)

# are axes significant?
anova(duneCCA, by = "axis")

# are vectors significant?
anova(duneCCA, by = "term", permutations = 999)

#cca north



duneCCA.north <- cca(comm.data.north ~ dist.to.strand..m. + salinity + soil.moisture.pct + pct.OM, data = model.df.north)
duneCCA.north
summary(duneCCA.north)

duneCCAplot.north <- plot(duneCCA.north)

# is CCA significant?
anova(duneCCA.north)

# are axes significant?
anova(duneCCA.north, by = "axis")

# are vectors significant?
anova(duneCCA.north, by = "term", permutations = 999)

#cca south


duneCCA.south<- cca(comm.data.south ~ dist.to.strand..m. + salinity + soil.moisture.pct + pct.OM, data = model.df.south)
duneCCA.south
summary(duneCCA.north)

duneCCAplot.south <- plot(duneCCA.south)

# is CCA significant?
anova(duneCCA.south)

# are axes significant?
anova(duneCCA.south, by = "axis")

# are vectors significant?
anova(duneCCA.south, by = "term", permutations = 999)


### testing lmer for N & S

dunecommrand_north <- lmer(shannon_diversity ~ dist.to.strand..m. + salinity 
                           + soil.moisture.pct + pct.OM +(1 | Site), data = model.df.north)
plot(parameters(dunecommrand_north))

summary(dunecommrand_north)

results.n <- psycho::analyze(dunecommrand_north, CI = 95)

dunecommrand_south <- lmer(shannon_diversity ~ salinity + dist.to.strand..m. 
                           + soil.moisture.pct + pct.OM +(1 | Site), data = model.df.south)
plot(parameters(dunecommrand_south))

plot(parameters(dunecommrand))



report.s <- report(dunecommrand_south)

report(dunecommrand)

confint <- as.data.frame(confint(dunecommrand_south)) 
confint <- confint[-c(1:3),]

estimate_mod <- as.data.frame(coef(summary(dunecommrand_south)))
estimate_mod <- estimate_mod[-1,]

coef <- c("dist.to.strand..m.","salinity", "soil.moisture.pct", "pct.OM")

confint <- cbind(confint, coef)
colnames(confint) <- c("lowerci", "upperci", "coef")
confint <- confint %>% 
  mutate(estimate = estimate_mod$Estimate,
         class = c(1:4))

plot(parameters(dunecommrand))
report(dunecommrand) %>% 
  summary()
summary(as.data.frame(report_all))


# plant community plots 
plant_comms_vis <- read.csv("data/plant_communities.csv")


plant_comm_plot <- ggplot(plant_comms_vis, aes(fill=Species, y=Occurrence, x=Site)) + 
  geom_bar(position="stack", stat="identity")
plant_comm_plot


# Stacked
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="stack", stat="identity")


plot(parameters(dunecommrand_north))
plot(parameters(dunecommrand_south))

report(dunecommrand) %>% 
  summary()

sjPlot:: tab_model(dunecommrand)


sjPlot::tab_model(dunecommrand, 
                  show.re.var= TRUE, 
                  pred.labels =c("(Intercept)", "Salinity", "Dist. to Stand (m)", "Soil Moisture (%)", "Org. Material (%)"),
                  dv.labels= "Effects of Abiotic Factors on Shannon Diversity")

report(dunecommrand_north) %>% 
  summary()

sjPlot::tab_model(dunecommrand_north)


report(dunecommrand_south) %>% 
  summary()
sjPlot::tab_model(dunecommrand_south)


