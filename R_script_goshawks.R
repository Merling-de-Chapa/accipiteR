## Script written by Manuela Merling, IZW Berlin

library(readxl)
library(spaMM)
library(DHARMa)
library(lme4) 
library(car)
library(purrr)
library(binom)
library(ggplot2)
library(maps)
library(mapdata)
library(dplyr)
library(doSNOW)
spaMM.options(nb_cores = 3, separation_max = 1)
boot.repl <- 1000
nb_cores <- 3


############################# behavioural flexibility #####################################################

table_per_nest <- read_excel("./source_data/goshawk_data_nest.xlsx")

goshawk_nest1 <- subset(table_per_nest, Year == "2015" | Year == "2016")
goshawk_nest <- goshawk_nest1[, c("Year", "Location", "Age_youngest", "No_nestlings",
                                  "Laying_begin_day", "Territory", "Habitat", 
                                  "reaction_female", "Rainfall")]
goshawk_nest$Year <- as.factor(goshawk_nest$Year)
goshawk_nest$Location <- as.factor(goshawk_nest$Location)
goshawk_nest$Territory <- as.factor(goshawk_nest$Territory)
goshawk_nest$Habitat <- as.factor(goshawk_nest$Habitat)
goshawk_nest$reaction_female <- factor(goshawk_nest$reaction_female)
goshawk_nest <- as.data.frame(goshawk_nest)
goshawk_nest <- droplevels(na.omit(goshawk_nest)) 

goshawk_nest$reaction <- goshawk_nest$reaction_female != "no reaction"
goshawk_nest$Age <- ifelse(goshawk_nest$Age_youngest <= 16,"young","old")
goshawk_nest$Age <- as.factor(goshawk_nest$Age)

str(goshawk_nest)

# subset Stadt/Land
urban <- subset(goshawk_nest, Habitat == "urban")
rural <- subset(goshawk_nest, Habitat == "rural")

test_react_int <- fitme(reaction ~ Habitat*Age + No_nestlings + Laying_begin_day + 
                        Year + Rainfall + (1|Location/Territory),
                        family = binomial(link = "logit"), data = goshawk_nest, 
                        method = "PQL/L")

test_react_0 <- fitme(reaction~1+(1|Location/Territory),
                      family = binomial(link = "logit"), data = goshawk_nest, 
                      method = "PQL/L")
anova(test_react_int, test_react_0)

## effects
test_no_Hab_age <- fitme(reaction ~ No_nestlings +  Laying_begin_day + Year + 
                         Rainfall + (1|Location/Territory), 
                         family = binomial(link = "logit"), data = goshawk_nest, 
                         method = "PQL/L")
anova(test_react_int, test_no_Hab_age)

test_noHabitat <- fitme(reaction ~ Age + No_nestlings +  Laying_begin_day + Year + 
                        Rainfall + (1|Location/Territory),
                        family = binomial(link = "logit"), data = goshawk_nest, 
                        method = "PQL/L")
anova(test_react_int, test_noHabitat)

test_nonest <- fitme(reaction ~ Habitat*Age + Laying_begin_day + Year +  Rainfall 
                     + (1|Location/Territory), family = binomial(link = "logit"), 
                     data = goshawk_nest, method = "PQL/L")
anova(test_react_int, test_nonest)

test_noYear <- fitme(reaction ~ Habitat*Age + No_nestlings + Laying_begin_day + 
                     Rainfall + (1|Location/Territory), 
                     family = binomial(link = "logit"), data = goshawk_nest, 
                     method = "PQL/L")
anova(test_react_int, test_noYear)

test_noDay <- fitme(reaction ~ Habitat*Age + No_nestlings +  Year + Rainfall 
                    + (1|Location/Territory), family = binomial(link = "logit"), 
                    data = goshawk_nest, method = "PQL/L")
anova(test_react_int, test_noDay)

test_noRain <- fitme(reaction ~ Habitat*Age + No_nestlings + Laying_begin_day + 
                     Year + (1|Location/Territory), family = binomial(link = "logit"),
                     data = goshawk_nest, method = "PQL/L")
anova(test_react_int, test_noRain)

test_noYoung <- fitme(reaction ~ Habitat + Rainfall + No_nestlings + 
                      Laying_begin_day + Year + (1|Location/Territory),
                      family = binomial(link = "logit"), data = goshawk_nest, 
                      method = "PQL/L")
anova(test_react_int, test_noYoung)

#### Habitat
ODD <- exp(test_react_int$fixef["Habitaturban"])
ODD 

exp(test_react_int$fixef)

## ODDS 
rural_young <- data.frame(Habitat = "rural", Age = "young", No_nestlings = 0, Laying_begin_day = 0,
                          Year = "2015", Rainfall = 0, Location = "new", Territory = "new") 
urban_young <- data.frame(Habitat = "urban", Age = "young", No_nestlings = 0, Laying_begin_day = 0,
                          Year = "2015", Rainfall = 0, Location = "new", Territory = "new") 
rural_old   <- data.frame(Habitat = "rural", Age = "old",   No_nestlings = 0, Laying_begin_day = 0,
                          Year = "2015", Rainfall = 0, Location = "new", Territory = "new") 
urban_old   <- data.frame(Habitat = "urban", Age = "old",   No_nestlings = 0, Laying_begin_day = 0,
                          Year = "2015", Rainfall = 0, Location = "new", Territory = "new") 
prob_rural_young <- predict(test_react_int, newdata = rural_young, re.form = NA)[[1]]
prob_urban_young <- predict(test_react_int, newdata = urban_young, re.form = NA)[[1]]
prob_rural_old   <- predict(test_react_int, newdata = rural_old, re.form = NA)[[1]]
prob_urban_old   <- predict(test_react_int, newdata = urban_old, re.form = NA)[[1]]

OR <- function(p1, p2) (p1/(1 - p1))/(p2/(1 - p2))

odd_ratio_urban_young_vs_rural_young <- round(OR(prob_urban_young, prob_rural_young), 2)
odd_ratio_urban_old_vs_rural_old <- round(OR(prob_urban_old, prob_rural_old), 2)
odd_ratio_urban_young_vs_urban_old <- round(OR(prob_urban_young, prob_urban_old), 2)
odd_ratio_rural_young_vs_rural_old <- round(OR(prob_rural_young, prob_rural_old), 2)

#### bootstrap
anova(test_react_int, test_noHabitat, boot.repl = boot.repl, nb_cores = nb_cores)
anova(test_react_int, test_nonest, boot.repl = boot.repl, nb_cores = nb_cores)
anova(test_react_int, test_noYear, boot.repl = boot.repl, nb_cores = nb_cores)
anova(test_react_int, test_noDay, boot.repl = boot.repl, nb_cores = nb_cores)
anova(test_react_int, test_noRain, boot.repl = boot.repl, nb_cores = nb_cores)
anova(test_react_int, test_noYoung, boot.repl = boot.repl, nb_cores = nb_cores)

##### test without interaction ###############

test_no_int <- fitme(reaction ~ Habitat + Age + No_nestlings + Laying_begin_day + 
                     Year + Rainfall + (1|Location/Territory),
                     family = binomial(link = "logit"), data = goshawk_nest, 
                     method = "PQL/L")
round(exp(fixef(test_no_int)["Habitaturban"]), 2)

test_no_Habitat <- fitme(reaction ~ Age + No_nestlings + Laying_begin_day + Year + 
                         Rainfall + (1|Location/Territory),
                         family = binomial(link = "logit"), data = goshawk_nest, 
                         method = "PQL/L")
test_no_Age <- fitme(reaction ~ Habitat + No_nestlings + Laying_begin_day + Year + 
                     Rainfall + (1|Location/Territory),
                     family = binomial(link = "logit"), data = goshawk_nest, 
                     method = "PQL/L")


anova(test_react_int, test_no_int, boot.repl = boot.repl, nb_cores =nb_cores)
anova(test_no_int, test_no_Age, boot.repl = boot.repl, nb_cores = nb_cores)
anova(test_no_int, test_no_Habitat, boot.repl = boot.repl, nb_cores = nb_cores)

#### figures #####################
#### barplot raw data

str(goshawk_nest$reaction_female)

total_rural <- nrow(goshawk_nest[goshawk_nest$Habitat  == "rural", ])
total_urban <- nrow(goshawk_nest[goshawk_nest$Habitat  == "urban", ])

plot_data <- goshawk_nest %>% 
  group_by(reaction_female, Habitat) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(total = ifelse(Habitat == "rural", total_rural, total_urban), 
         prop = n/total, 
         CI_lwr = map2_dbl(n, total, ~ as.numeric(binom.confint( x = .x, n = .y)[5, 5])), ## method exact
         CI_upr = map2_dbl(n, total, ~ as.numeric(binom.confint( x = .x, n = .y)[5, 6]))) ## method exact


plot_data$reaction_female <- relevel(plot_data$reaction_female, ref = "no reaction")
plot_data$Habitat <- relevel(plot_data$Habitat, ref = "urban")
ggplot(plot_data, aes(x = reaction_female, y = prop, fill = Habitat, ymin = CI_lwr, ymax = CI_upr)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.8) +
  geom_errorbar(width = 0, position = position_dodge(preserve = "total", width = 0.9)) +
  scale_fill_grey() +
  scale_y_continuous("Proportion of females reacting") +
  scale_x_discrete("Reaction") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank()) 

ggsave(filename = "reaction_female_plot.pdf", width = 8, height =  7)


# predicite barplot Habitat/Age

pred_urban_young <- predict(test_no_int,
                            newdata = expand.grid(Laying_begin_day = median(goshawk_nest$Laying_begin_day, na.rm = TRUE),
                                                  Age = "young",
                                                  Rainfall = median(goshawk_nest$Rainfall, na.rm = TRUE),
                                                  No_nestlings = 3,
                                                  Year = "2015",
                                                  Habitat = "urban"),
                            re.form = NA, intervals = "predVar") # predVar = SD of Prediction^2
pred_1 <- data.frame(pred = pred_urban_young, attr(pred_urban_young, "intervals"), 
                     Habitat = "urban", Age = "young")
pred_1

pred_urban_old <- predict(test_no_int, newdata = expand.grid(Laying_begin_day = median(goshawk_nest$Laying_begin_day, na.rm = TRUE),
                                                             Age = "old",
                                                             Rainfall = median(goshawk_nest$Rainfall, na.rm = TRUE),
                                                             No_nestlings = 3,
                                                             Year = "2015",
                                                             Habitat = "urban"),
                          re.form = NA, intervals = "predVar") # predVar = SD of Prediction^2
pred_2 <- data.frame(pred = pred_urban_old, attr(pred_urban_old, "intervals"), 
                     Habitat = "urban", Age = "old")

pred_rural_young <- predict(test_no_int, newdata = expand.grid(Laying_begin_day = median(goshawk_nest$Laying_begin_day, na.rm = TRUE),
                                                               Age = "young",
                                                               Rainfall = median(goshawk_nest$Rainfall, na.rm = TRUE),
                                                               No_nestlings = 3,
                                                               Year = "2015",
                                                               Habitat = "rural"),
                            re.form = NA, intervals = "predVar") # predVar = SD of Prediction^2
pred_3 <- data.frame(pred = pred_rural_young, attr(pred_rural_young, "intervals"), 
                     Habitat = "rural", Age = "young")

pred_rural_old <- predict(test_no_int, newdata = expand.grid(Laying_begin_day = median(goshawk_nest$Laying_begin_day, na.rm = TRUE),
                                                             Age = "old",
                                                             Rainfall = median(goshawk_nest$Rainfall, na.rm = TRUE),
                                                             No_nestlings = 3,
                                                             Year = "2015",
                                                             Habitat = "rural"),
                          re.form = NA, intervals = "predVar") # predVar = SD of Prediction^2
pred_4 <- data.frame(pred = pred_rural_old, attr(pred_rural_old, "intervals"), 
                     Habitat = "rural", Age = "old")

pred_Habitat_age <- rbind(pred_1, pred_2, pred_3, pred_4)

ggplot(pred_Habitat_age, aes(x = Age, y = pred, fill = Habitat, ymin = predVar_0.025, ymax = predVar_0.975)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.8) +
  geom_errorbar(width = 0, position = position_dodge(preserve = "total", width = 0.9)) +
  scale_fill_grey() +
  scale_y_continuous("Proportion of females reacting") +
  scale_x_discrete("Age categorical") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank()) 

ggsave(filename = "./figures/reaction_Habitat.pdf", width = 8, height =  7)

# predicitve plot laying

predictor_lay <- seq(min(goshawk_nest$Laying_begin_day, na.rm = TRUE), max(goshawk_nest$Laying_begin_day, na.rm = TRUE), length = 30)
pred_lay_obj <- predict(test_no_int,
                        newdata = expand.grid(Laying_begin_day = predictor_lay,
                                              Age = "young",
                                              Rainfall = median(goshawk_nest$Rainfall, na.rm = TRUE),
                                              No_nestlings = 3,
                                              Year = "2015",
                                              Habitat = "urban"),
                        re.form = NA, intervals = "predVar") # predVar = SD of Prediction^2
pred_lay <- data.frame(pred = pred_lay_obj, attr(pred_lay_obj, "intervals"), lay = predictor_lay)

predictor_lay2 <- seq(min(goshawk_nest$Laying_begin_day, na.rm = TRUE),
                      max(goshawk_nest$Laying_begin_day, na.rm = TRUE),
                      length = 30)
pred_lay_obj2 <- predict(test_no_int,
                         newdata = expand.grid(Laying_begin_day = predictor_lay2,
                                               Age = "young",
                                               Rainfall = median(goshawk_nest$Rainfall, na.rm = TRUE), 
                                               No_nestlings = 3,
                                               Year = "2015",
                                               Habitat = "rural"),
                         re.form = NA, intervals = "predVar") # predVar = SD of Prediction^2
pred_lay2 <- data.frame(pred = pred_lay_obj2, attr(pred_lay_obj2, "intervals"), lay2 = predictor_lay2)


pdf("./figures/plot_prediction1.pdf", width = 10, height = 7)

plot(pred_lay$pred ~ pred_lay$lay, type = "l", lwd = 2, ylab = "Proportion of females reacting",
     las = 1, ylim = c(0, 1), xlab = "Laying date", axes = FALSE)
points(pred_lay$predVar_0.025~ pred_lay$lay, type = "l", lwd = 2, lty = 2)
points(pred_lay$predVar_0.975~ pred_lay$lay, type = "l", lwd = 2, lty = 2)
rug(x = jitter(urban$Laying_begin_day[!goshawk_nest$reaction]), side = 1)
rug(x = jitter(urban$Laying_begin_day[goshawk_nest$reaction]), side = 3)


lines(pred_lay2$pred ~ pred_lay2$lay2, type = "l", lwd = 2, ylab = "Proportion of females reacting",
      las = 1, ylim = c(0, 1), xlab = "Laying date",  col = "grey")
points(pred_lay2$predVar_0.025~ pred_lay2$lay2, type = "l", lwd = 2, lty = 2, col = "grey")
points(pred_lay2$predVar_0.975~ pred_lay2$lay2, type = "l", lwd = 2, lty = 2, col = "grey")
box()
rug(x = jitter(rural$Laying_begin_day[!goshawk_nest$reaction]), side = 1, col = "grey")
rug(x = jitter(rural$Laying_begin_day[goshawk_nest$reaction]), side = 3, col = "grey")
labelist <- as.vector(c("March 02","March 12", "March 22", "April 01", "April 11", "April 21", "May 01"))
axis(1, at = seq(60,120, by = 10), labels = labelist)
axis(2, las = 2)
dev.off()

###################### Diet composition ####################################################################
##### figure

Rupfungen_gesamt <- read_excel("./source_data/Goshawk_data_diet_plot.xlsx")
n_land <- nrow(Rupfungen_gesamt[Rupfungen_gesamt$Habitat == "Land", ])
n_stadt <- nrow(Rupfungen_gesamt[Rupfungen_gesamt$Habitat == "Stadt", ])

cut_off <- Rupfungen_gesamt %>% 
  group_by(species, Habitat) %>%
  summarise(n = n()) %>% 
  mutate(n_total = ifelse(Habitat == "Land", n_land, n_stadt), 
         prop = n / n_total) %>% 
  ungroup()

### table with all species even < 1%
cut_off <- cut_off  %>% 
  mutate(CI_lwr = map2_dbl(n, n_total, ~ as.numeric(binom.confint( x = .x, n = .y)[5, 5])), ## method exact
         CI_upr = map2_dbl(n, n_total, ~ as.numeric(binom.confint( x = .x, n = .y)[5, 6]))) ## method exact

## filter to keep only species > 1%
cut_off2 <- cut_off %>% 
  filter(prop > 0.01)

cut_off2$species <- relevel(as.factor(cut_off2$species), ref = "Feral pigeon")
cut_off2$species  <- factor(cut_off2$species,
                            levels = levels(cut_off2$species)[c(1, 22, 12, 20, 9, 3, 10, 14, 6, 4, 8, 5, 18, 15, 13, 21, 7, 16, 19, 17, 2, 11)])
cut_off2$Habitat <- factor(cut_off2$Habitat, levels = c("Stadt", "Land"))
cut_off2$Habitat <- relevel(as.factor(cut_off2$Habitat), ref = "Stadt")

ggplot(cut_off2, aes(x = species, y = prop, fill = Habitat, ymin = CI_lwr, ymax = CI_upr)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.8) +
  geom_errorbar(position = position_dodge(preserve = "total", width = 0.9), width = 0) +
  scale_fill_grey(labels = c("urban", "rural")) +
  scale_y_continuous("Proportion of prey species") +
  scale_x_discrete("Species") +
  geom_vline(xintercept = 4.5, linetype = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust =  1), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) 

ggsave(filename = "./figures/species_plot.pdf", width = 10, height =  7)

########## pigeons ###################
Rupfungen_per_nest_aktuell <- read_excel("./source_data/goshawk_data_diet.xlsx")

pigeon <- Rupfungen_per_nest_aktuell[, c("Location", "nest_side", "Habitat",
                                         "Species_richness_nest", "Pigeon", "Rest")]
pigeon$Location <- as.factor(pigeon$Location)
pigeon$nest_side <- as.factor(pigeon$nest_side)
pigeon$Habitat <- as.factor(pigeon$Habitat)
pigeon <- as.data.frame(pigeon)
str(pigeon)

GLMM_pigeon_SpaMM_loc <- fitme(cbind(Pigeon, Rest)~Habitat + (1|Location/nest_side),
                               family = binomial(link = "logit"),
                               data = pigeon,
                               method = "PQL/L")

GLMM_0_loc <- fitme(cbind(Pigeon, Rest) ~ 1 + (1|Location/nest_side),
                    family = binomial(link = "logit"),
                    data = pigeon,
                    method = "PQL/L")
anova(GLMM_0_loc, GLMM_pigeon_SpaMM_loc) 

# prediction
ODD <- exp(GLMM_pigeon_SpaMM_loc$fixef)
ODD #3.654

#### bootstrap ###
anova(GLMM_0_loc, GLMM_pigeon_SpaMM_loc, boot.repl = boot.repl, nb_cores = nb_cores) 

############ diversity ######################

Rupfungen_per_nest_aktuell <- read_excel("./source_data/goshawk_data_diet.xlsx")

diversity <- Rupfungen_per_nest_aktuell[, c("Location", "nest_side", "Habitat",
                                            "diversity_nest")]
diversity$Location <- as.factor(diversity$Location)
diversity$nest_side <- as.factor(diversity$nest_side)
diversity$Habitat <- as.factor(diversity$Habitat)
diversity$diversity_nest <- as.numeric(diversity$diversity_nest)
diversity <- na.omit(diversity)
diversity <- as.data.frame(diversity)
str(diversity)


lmer_diversity <- lmer(diversity_nest ~ Habitat + (1|Location), 
                       data = diversity, REML = FALSE)

spaMM_diversity_coor <- fitme(diversity_nest ~ Habitat +(1|Location), data = diversity)

bc_lmer <- powerTransform(lmer_diversity, family = "bcPower")

lmer_diversity2 <- update(lmer_diversity, bcPower(diversity_nest, bc_lmer$lambda) ~ .)
spaMM_diversity3 <- update(spaMM_diversity_coor, bcPower(diversity_nest, bc_lmer$lambda) ~ .)

logLik(lmer_diversity2)
logLik(spaMM_diversity3)

lmer_diversity_0 <- lmer(diversity_nest ~ 1 + (1|Location), data = diversity, 
                         REML = FALSE)
spaMM_diversity_0 <- fitme(diversity_nest ~ 1 + (1|Location), data = diversity)

lmer_diversity_upda <- update(lmer_diversity_0, bcPower(diversity_nest, bc_lmer$lambda) ~ .)
spaMM_diversity_upda_0 <- update(spaMM_diversity_0, bcPower(diversity_nest, bc_lmer$lambda) ~ .)

anova(spaMM_diversity3, spaMM_diversity_upda_0)

### effect ########
pred_trans_rural <- (spaMM_diversity3$fixef["(Intercept)"] * bc_lmer$lambda + 1)^(1/bc_lmer$lambda)
pred_trans_urban <- ((spaMM_diversity3$fixef["(Intercept)"] + spaMM_diversity3$fixef["HabitatUrban"]) * bc_lmer$lambda + 1)^(1/bc_lmer$lambda)
pred_trans_urban - pred_trans_rural

### bootstrap
anova(spaMM_diversity3, spaMM_diversity_upda_0, boot.repl = boot.repl, nb_cores = nb_cores)


################# species richness ############################################################################################

Rupfungen_per_nest_aktuell <- read_excel("./source_data/goshawk_data_diet.xlsx")

species <- Rupfungen_per_nest_aktuell[, c("Location","nest_side", "Habitat", 
                                          "Species_richness_nest")]

species$Location <- as.factor(species$Location)
species$nest_side <- as.factor(species$nest_side)
species$Habitat <- as.factor(species$Habitat)
species$Species_richness_nest <- as.numeric(species$Species_richness_nest)
species <- na.omit(species)
species <- as.data.frame(species)
str(species)

glm_species_spaMM <- fitme(Species_richness_nest ~ Habitat+(1|Location),
                           family = Tnegbin(),
                           data = species)

glm_0_spaMM <- fitme(Species_richness_nest ~ 1 + (1|Location), 
                     family = Tnegbin, data = species)
anova(glm_species_spaMM, glm_0_spaMM)


# Habitat Effekt
predict(glm_species_spaMM, newdata = data.frame(Habitat = c("Urban", "Rural")), 
        re.form = NA)

diff(predict(glm_species_spaMM, newdata = data.frame(Habitat = c("Urban", "Rural")),
             re.form = NA))

#### bootstrap 
anova(glm_species_spaMM, glm_0_spaMM, boot.repl = boot.repl, nb_cores = nb_cores)


########### Breeding performance #########################################################################################
### Laying date ###
table_per_nest <- read_excel("./source_data/goshawk_data_nest.xlsx")

goshawk_nest2 <- table_per_nest[, c("Year", "Location", "No_nestlings", "Laying_begin_day",
                                    "Territory",  "Habitat", "Temp_breeding_begin")]
goshawk_nest2$Year <- as.factor(goshawk_nest2$Year)
goshawk_nest2$Location <- as.factor(goshawk_nest2$Location)
goshawk_nest2$Territory <- as.factor(goshawk_nest2$Territory)
goshawk_nest2$Habitat <- as.factor(goshawk_nest2$Habitat)
goshawk_nest2 <- as.data.frame(goshawk_nest2)
goshawk_nest2 <- droplevels(na.omit(goshawk_nest2))
str(goshawk_nest2)

lmm_laying2 <- fitme(Laying_begin_day ~ Habitat + Temp_breeding_begin  + 
                     (1|Location/Territory), family = gaussian, data = goshawk_nest2)

# test against zero model
lmm_0 <- fitme(Laying_begin_day ~ 1 + (1|Location/Territory),
               family = gaussian,
               data = goshawk_nest2)
anova(lmm_laying2, lmm_0) 

# odds
lmm_laying2$fixef["Habitaturban"] 

### test  effects
lmm_laying_no_Habitat <- fitme(Laying_begin_day ~ Temp_breeding_begin + 
                              (1|Location/Territory), family = gaussian,
                              data = goshawk_nest2)
anova(lmm_laying2, lmm_laying_no_Habitat)

lmm_laying_no_Temp <- fitme(Laying_begin_day ~ Habitat + (1|Location/Territory),
                            family = gaussian, data = goshawk_nest2)
anova(lmm_laying2, lmm_laying_no_Temp)

### bootstrap
anova(lmm_laying2, lmm_laying_no_Habitat, boot.repl = boot.repl, nb_cores = nb_cores)
anova(lmm_laying2, lmm_laying_no_Temp, boot.repl = boot.repl, nb_cores = nb_cores)


########## reproductive output (number of nestlings) ############
table_per_nest <- read_excel("./source_data/goshawk_data_nest.xlsx")

goshawk_nest3 <- table_per_nest[, c("Year", "Location", "No_nestlings", "Laying_begin_day",
                                    "Territory", "Habitat", "Temp_breeding_begin")]
goshawk_nest3$Year <- as.factor(goshawk_nest3$Year)
goshawk_nest3$Location <- as.factor(goshawk_nest3$Location)
goshawk_nest3$Territory <- as.factor(goshawk_nest3$Territory)
goshawk_nest3$Habitat <- as.factor(goshawk_nest3$Habitat)
goshawk_nest3 <- as.data.frame(goshawk_nest3)
goshawk_nest3 <- droplevels(na.omit(goshawk_nest3))

goshawk_nest3$No_nestlings_binary <- goshawk_nest3$No_nestlings > 2
table(goshawk_nest3$No_nestlings_binary)

glmm_nestlings <- fitme(No_nestlings_binary ~ Habitat + Laying_begin_day + 
                        Temp_breeding_begin + (1|Location/Territory),
                        family = binomial(link = "logit"), data = goshawk_nest3, 
                        method = "PQL/L")

## 0-model
glmm_0_spaMM <- fitme(No_nestlings_binary ~ 1 + (1|Location/Territory),
                      family = binomial(link = "logit"),
                      data = goshawk_nest3,
                      method = "PQL/L")
anova(glmm_nestlings, glmm_0_spaMM)

# effects
glmm_nestlings_no_Habitat <- fitme(No_nestlings_binary ~ Laying_begin_day + 
                                   Temp_breeding_begin + (1|Location/Territory),
                                   family = binomial(link = "logit"), data = goshawk_nest3, 
                                   method = "PQL/L")
anova(glmm_nestlings, glmm_nestlings_no_Habitat)

glmm_nestlings_no_day <- fitme(No_nestlings_binary ~ Habitat + Temp_breeding_begin + 
                              (1|Location/Territory), family = binomial(link = "logit"), 
                              data = goshawk_nest3, method = "PQL/L")
anova(glmm_nestlings, glmm_nestlings_no_day)

glmm_nestlings_no_Temp <- fitme(No_nestlings_binary ~ Habitat + Laying_begin_day + 
                                (1|Location/Territory), family = binomial(link = "logit"), 
                                data = goshawk_nest3, method = "PQL/L")
anova(glmm_nestlings, glmm_nestlings_no_Temp)


# bootstraps
anova(glmm_nestlings, glmm_0_spaMM, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_nestlings, glmm_nestlings_no_Habitat, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_nestlings, glmm_nestlings_no_day, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_nestlings, glmm_nestlings_no_Temp, boot.repl = boot.repl, nb_cores = nb_cores)

##### predictive figures ####
# subset urban/rural
urban1 <- subset(goshawk_nest3, Habitat == "urban")
rural1 <- subset(goshawk_nest3, Habitat == "rural")

predictor_day <- seq(min(goshawk_nest3$Laying_begin_day, na.rm = TRUE),
                     max(goshawk_nest3$Laying_begin_day, na.rm = TRUE),
                     length = 30)
pred_day_obj <- predict(glmm_nestlings, newdata = expand.grid(Laying_begin_day = predictor_day,
                                                              Habitat = "urban",
                                                              Temp_breeding_begin = median(goshawk_nest3$Temp_breeding_begin, na.rm = TRUE)),
                        re.form = NA, intervals = "predVar") # predVar = SD of Prediction^2
pred_day <- data.frame(pred = pred_day_obj, attr(pred_day_obj, "intervals"), day = predictor_day)

predictor_day2 <- seq(min(goshawk_nest3$Laying_begin_day, na.rm = TRUE), max(goshawk_nest3$Laying_begin_day, na.rm = TRUE), length = 30)
pred_day_obj2 <- predict(glmm_nestlings, newdata = expand.grid(Laying_begin_day = predictor_day2,
                                                               Habitat = "rural",
                                                               Temp_breeding_begin = median(goshawk_nest3$Temp_breeding_begin, na.rm = TRUE)),
                         re.form = NA, intervals = "predVar") # predVar = SD of Prediction^2
pred_day2 <- data.frame(pred2 = pred_day_obj2, attr(pred_day_obj2, "intervals"), day2 = predictor_day2)

pdf("./figures/prediction_plot2.pdf", width = 10, height = 7)

plot(pred_day$pred ~ pred_day$day, type = "l", lwd = 2, ylab = "Proportion of large clutches",
     las = 1, ylim = c(0, 1), xlab = "Laying date", axes = FALSE)
points(pred_day$predVar_0.025~ pred_day$day, type = "l", lwd = 2, lty = 2)
points(pred_day$predVar_0.975~ pred_day$day, type = "l", lwd = 2, lty = 2)
rug(x = jitter(urban1$Laying_begin_day[!urban1$No_nestlings_binary]), side = 1)
rug(x = jitter(urban1$Laying_begin_day[urban1$No_nestlings_binary]), side = 3)

lines(pred_day2$pred2 ~ pred_day2$day2, type = "l", lwd = 2, ylab = "Proportion of large clutches",
      las = 1, ylim = c(0, 1), xlab = "Laying date", col = "grey")
points(pred_day2$predVar_0.025~ pred_day2$day2, type = "l", lwd = 2, lty = 2, col = "grey")
points(pred_day2$predVar_0.975~ pred_day2$day2, type = "l", lwd = 2, lty = 2, col = "grey")
box()
rug(x = jitter(rural1$Laying_begin_day[!rural1$No_nestlings_binary]), side = 1, col = "grey")
rug(x = jitter(rural1$Laying_begin_day[rural1$No_nestlings_binary]), side = 3, col = "grey")
labelist <- as.vector(c("March 02","March 12", "March 22", "April 01", "April 11", "April 21", "May 01"))
axis(1, at = seq(60,120, by = 10), labels = labelist)
axis(2, las = 2)

dev.off()

###### without laying start ##############################
glmm_nestlings2 <- fitme(No_nestlings_binary ~ Habitat + Temp_breeding_begin + 
                        (1|Location/Territory), family = binomial(link = "logit"),
                         data = goshawk_nest3, method = "PQL/L")

# 0-model
glmm_0_spaMM2 <- fitme(No_nestlings_binary ~ 1 + (1|Location/Territory),
                       family = binomial(link = "logit"),
                       data = goshawk_nest3, method = "PQL/L")
anova(glmm_nestlings2, glmm_0_spaMM2)

# Odd
exp(glmm_nestlings2$fixef["Habitaturban"]) 

#### testing effects########
glmm_nestlings_no_Habitat2 <- fitme(No_nestlings_binary ~ Temp_breeding_begin + 
                                    (1|Location/Territory), family = binomial(link = "logit"),
                                    data = goshawk_nest3, method = "PQL/L")
anova(glmm_nestlings2, glmm_nestlings_no_Habitat2)

glmm_nestlings_no_Temp2 <- fitme(No_nestlings_binary ~ Habitat + (1|Location/Territory),
                                 family = binomial(link = "logit"),
                                 data = goshawk_nest3,
                                 method = "PQL/L")
anova(glmm_nestlings2, glmm_nestlings_no_Temp2)

## bootstraps
anova(glmm_nestlings2, glmm_nestlings_no_Habitat2, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_nestlings2, glmm_nestlings_no_Temp2, boot.repl = boot.repl, nb_cores = nb_cores)


###################### Health status of goshawk nestlings ####################################

#### prevalence ####
Tabelle_Statistik_urban_ecology_aktuell <- read_excel("./source_data/Goshawk_data_nestlings.xlsx")


goshawk_urban <- Tabelle_Statistik_urban_ecology_aktuell[, c("Year", "Location","No_nestlings", "Sex", 
                                                             "Age", "Prevalence", "Territory", "Habitat", "laying_day",
                                                             "Temp_age", "Clinical_binary")]
goshawk_urban$Year <- as.factor(goshawk_urban$Year)
goshawk_urban$Location <- as.factor(goshawk_urban$Location)
goshawk_urban$Territory <- as.factor(goshawk_urban$Territory)
goshawk_urban$Sex <- as.factor(goshawk_urban$Sex)
goshawk_urban$Habitat <- as.factor(goshawk_urban$Habitat)
goshawk_urban$Clinical_binary <- as.numeric(goshawk_urban$Clinical_binary)
goshawk_urban <- as.data.frame(goshawk_urban)
goshawk_urban <- droplevels(na.omit(goshawk_urban))
str(goshawk_urban)

#### figure prevalence ###
### plot rawdata
plot_data <- goshawk_urban %>%
  group_by(Year, Habitat) %>%
  summarise(n = n(), 
            n_inf = sum(Prevalence == 1, na.rm = T), 
            n_not_inf = sum(Prevalence == 0, na.rm = T), 
            prop_inf = n_inf / n) %>% 
  ungroup() %>% 
  mutate(CI_lwr = map2_dbl(n_inf, n, ~ as.numeric(binom.confint( x = .x, n = .y)[5, 5])), ## method exact
         CI_upr = map2_dbl(n_inf, n, ~ as.numeric(binom.confint( x = .x, n = .y)[5, 6]))) ## method exact

plot_data$Habitat <- relevel(as.factor(plot_data$Habitat), ref = "urban")

ggplot(plot_data, aes(y = prop_inf, x = Year, ymin = CI_lwr, ymax = CI_upr, fill = Habitat)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.8) +
  geom_errorbar(width = 0, position =  position_dodge(width =  0.9)) +
  scale_fill_grey() +
  scale_y_continuous("Proportion of infected nestlings") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) 
ggsave(filename = "./figures/year_plot.pdf", width = 8, height = 8)

##
glmm_tricho_logit <- fitme(Prevalence ~ Habitat + Age + No_nestlings +  Sex + Temp_age + 
                           Year + laying_day +(1|Location/Territory),
                           family = binomial(link = "logit"),
                           data = goshawk_urban,
                           method = "PQL/L")

# 0-model
glmm_tricho_logit_0 <- fitme(Prevalence ~ 1 + (1|Location/Territory),
                             family = binomial(link = "logit"), 
                             data = goshawk_urban, method = "PQL/L") 

anova(glmm_tricho_logit, glmm_tricho_logit_0) 

# odds
exp(glmm_tricho_logit$fixef["Habitaturban"])

#test effects
glmm_tricho_logit_noHabitat <- fitme(Prevalence ~  Age + No_nestlings + Sex + Temp_age + 
                                     Year + laying_day + (1|Location/Territory),
                                     family = binomial(link = "logit"), 
                                     data = goshawk_urban, method = "PQL/L")
anova(glmm_tricho_logit, glmm_tricho_logit_noHabitat)

glmm_tricho_logit_nonest <- fitme(Prevalence ~ Habitat + Age + Sex + Temp_age + Year + 
                                  laying_day + (1|Location/Territory),
                                  family = binomial(link = "logit"), 
                                  data = goshawk_urban, method = "PQL/L")
anova(glmm_tricho_logit, glmm_tricho_logit_nonest)

glmm_tricho_logit_noage <- fitme(Prevalence ~ Habitat + No_nestlings + Sex + Temp_age + 
                                 Year + laying_day + (1|Location/Territory),
                                 family = binomial(link = "logit"), data = goshawk_urban, 
                                 method = "PQL/L")
anova(glmm_tricho_logit, glmm_tricho_logit_noage)

glmm_tricho_logit_noTemp <- fitme(Prevalence ~ Habitat + Age + No_nestlings + Sex + 
                                  Year + laying_day + (1|Location/Territory),
                                  family = binomial(link = "logit"), data = goshawk_urban, 
                                  method = "PQL/L")
anova(glmm_tricho_logit, glmm_tricho_logit_noTemp)

glmm_tricho_logit_noSex <- fitme(Prevalence ~ Habitat + Age + No_nestlings + Year + 
                                 laying_day + Temp_age + (1|Location/Territory),
                                 family = binomial(link = "logit"), data = goshawk_urban, 
                                 method = "PQL/L")
anova(glmm_tricho_logit, glmm_tricho_logit_noSex)

glmm_tricho_logit_noYear <- fitme(Prevalence ~ Habitat + Age + No_nestlings + Sex + 
                                  Temp_age + laying_day + (1|Location/Territory),
                                  family = binomial(link = "logit"), data = goshawk_urban, 
                                  method = "PQL/L")
anova(glmm_tricho_logit, glmm_tricho_logit_noYear)

glmm_tricho_logit_noDay <- fitme(Prevalence ~ Habitat + Age + No_nestlings + Sex + 
                                 Temp_age + Year + (1|Location/Territory),
                                 family = binomial(link = "logit"), data = goshawk_urban,
                                 method = "PQL/L")
anova(glmm_tricho_logit, glmm_tricho_logit_noDay)

### bootstraps
anova(glmm_tricho_logit, glmm_tricho_logit_noHabitat, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_tricho_logit, glmm_tricho_logit_nonest, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_tricho_logit, glmm_tricho_logit_noage, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_tricho_logit, glmm_tricho_logit_noTemp, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_tricho_logit, glmm_tricho_logit_noSex, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_tricho_logit, glmm_tricho_logit_noYear, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_tricho_logit, glmm_tricho_logit_noDay, boot.repl = boot.repl, nb_cores = nb_cores)


### clinical signs ##############

glmm_clinical <- fitme(Clinical_binary ~ Habitat + Age + Sex + Year + laying_day + 
                       No_nestlings + Temp_age + (1|Location/Territory),
                       family = binomial(link = "logit"), 
                       data = goshawk_urban, method = "PQL/L")

# 0-model
glmm_0 <- fitme(Clinical_binary ~ 1 + (1|Location/Territory),
                family = binomial(link = "logit"),
                data = goshawk_urban, method = "PQL/L")
anova(glmm_0, glmm_clinical)

# odds
exp(glmm_clinical$fixef["Habitaturban"])

#test effects
glmm_clinical_noHabitat <- fitme(Clinical_binary ~  Age + No_nestlings + Sex + Year + 
                                 laying_day + Temp_age + (1|Location/Territory),
                                 family = binomial(link = "logit"), data = goshawk_urban, 
                                 method = "PQL/L")
anova(glmm_clinical, glmm_clinical_noHabitat)

glmm_clinical_nonest <- fitme(Clinical_binary ~ Habitat + Age + Sex + Year + laying_day + 
                              Temp_age + (1|Location/Territory),
                              family = binomial(link = "logit"), data = goshawk_urban, 
                              method = "PQL/L")
anova(glmm_clinical, glmm_clinical_nonest)

glmm_clinical_noage <- fitme(Clinical_binary ~ Habitat + No_nestlings + Sex + Year + 
                             laying_day + Temp_age + (1|Location/Territory),
                             family = binomial(link = "logit"), data = goshawk_urban, 
                             method = "PQL/L")
anova(glmm_clinical, glmm_clinical_noage)

glmm_clinical_noSex <- fitme(Clinical_binary ~ Habitat + Age + No_nestlings + Year + 
                             laying_day + Temp_age + (1|Location/Territory),
                             family = binomial(link = "logit"), data = goshawk_urban, 
                             method = "PQL/L")
anova(glmm_clinical, glmm_clinical_noSex)

glmm_clinical_noTemp <- fitme(Clinical_binary ~ Habitat + Age + No_nestlings + Sex + 
                             Year + laying_day + (1|Location/Territory),
                             family = binomial(link = "logit"), data = goshawk_urban, 
                             method = "PQL/L")
anova(glmm_clinical, glmm_clinical_noTemp)

glmm_clinical_noYear <- fitme(Clinical_binary ~ Habitat + Age + No_nestlings + Sex + 
                              laying_day + Temp_age + (1|Location/Territory),
                              family = binomial(link = "logit"), data = goshawk_urban, 
                              method = "PQL/L")
anova(glmm_clinical, glmm_clinical_noYear)

glmm_clinical_noDay <- fitme(Clinical_binary ~ Habitat + Age + No_nestlings + Sex + Year + 
                             Temp_age + (1|Location/Territory),
                             family = binomial(link = "logit"), data = goshawk_urban, 
                             method = "PQL/L")
anova(glmm_clinical, glmm_clinical_noDay)

## bootstraps
anova(glmm_clinical, glmm_clinical_noHabitat, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_clinical, glmm_clinical_nonest, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_clinical, glmm_clinical_noage, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_clinical, glmm_clinical_noSex, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_clinical, glmm_clinical_noTemp, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_clinical, glmm_clinical_noYear, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_clinical, glmm_clinical_noDay, boot.repl = boot.repl, nb_cores = nb_cores)

#### predicitve plots clinical signs ######
# subset urban/rural
urban2 <- subset(goshawk_urban, Habitat == "urban")
rural2 <- subset(goshawk_urban, Habitat == "rural")

# Laying day
predictor_lay3 <- seq(min(goshawk_urban$laying_day, na.rm = TRUE),
                      max(goshawk_urban$laying_day, na.rm = TRUE),
                      length = 30)
pred_lay_obj3 <- predict(glmm_clinical, newdata = expand.grid(laying_day = predictor_lay3,
                                                              No_nestlings = 3,
                                                              Year = "2015",
                                                              Sex = "female",
                                                              Age = median(goshawk_urban$Age, na.rm = TRUE),
                                                              Habitat = "urban",
                                                              Temp_age = median(goshawk_urban$Temp_age, na.rm = TRUE)),
                         re.form = NA, intervals = "predVar") # predVar = SD of Prediction^2
pred_lay3 <- data.frame(pred = pred_lay_obj3, attr(pred_lay_obj3, "intervals"), lay = predictor_lay3)

predictor_lay4 <- seq(min(goshawk_urban$laying_day, na.rm = TRUE),
                      max(goshawk_urban$laying_day, na.rm = TRUE),
                      length = 30)
pred_lay_obj4 <- predict(glmm_clinical, newdata = expand.grid(laying_day = predictor_lay4,
                                                              No_nestlings = 3,
                                                              Year = "2015",
                                                              Sex = "female",
                                                              Age = median(goshawk_urban$Age, na.rm = TRUE),
                                                              Habitat = "rural",
                                                              Temp_age = median(goshawk_urban$Temp_age, na.rm = TRUE)),
                         re.form = NA, intervals = "predVar") # predVar = SD of Prediction^2
pred_lay4 <- data.frame(pred4 = pred_lay_obj4, attr(pred_lay_obj4, "intervals"), lay4 = predictor_lay4)


pdf("./figures/prediction_plot3.pdf", width = 10, height = 7) 

plot(pred_lay3$pred ~ pred_lay3$lay, type = "l", lwd = 2,
     ylab = "Proportion of nestlings presenting clinical signs", las = 1, ylim = c(0, 1),
     xlab = "Laying date", axes = FALSE)
points(pred_lay3$predVar_0.025~ pred_lay3$lay, type = "l", lwd = 2, lty = 2)
points(pred_lay3$predVar_0.975~ pred_lay3$lay, type = "l", lwd = 2, lty = 2)
rug(x = jitter(urban2$laying_day[urban2$Clinical_binary == "0"]), side = 1)
rug(x = jitter(urban2$laying_day[urban2$Clinical_binary == "1"]), side = 3)

lines(pred_lay4$pred4 ~ pred_lay4$lay4, type = "l", lwd = 2,
      ylab = "Proportion of nestlings presenting clinical signs",
      las = 1, ylim = c(0, 1), xlab = "Laying date", col = "grey")
points(pred_lay4$predVar_0.025~ pred_lay4$lay4, type = "l", lwd = 2, lty = 2, col = "grey")
points(pred_lay4$predVar_0.975~ pred_lay4$lay4,type = "l", lwd = 2, lty = 2, col = "grey")
box()
rug(x = jitter(rural2$laying_day[rural2$Clinical_binary == "0"]), side = 1, col = "grey")
rug(x = jitter(rural2$laying_day[rural2$Clinical_binary == "1"]), side = 3, col = "grey")
labelist <- as.vector(c("March 02","March 12", "March 22", "April 01", "April 11", "April 21", "May 01"))
axis(1, at = seq(60,120, by = 10), labels = labelist)
axis(2, las = 2)

dev.off()

#### Age
predictor_age <- seq(min(goshawk_urban$Age, na.rm = TRUE),
                     max(goshawk_urban$Age, na.rm = TRUE),
                     length = 30)
pred_age_obj <- predict(glmm_clinical, newdata = expand.grid(Age = predictor_age,
                                                             No_nestlings = 3,
                                                             Year = "2015",
                                                             Sex = "female",
                                                             laying_day = median(goshawk_urban$laying_day, na.rm = TRUE),
                                                             Temp_age = median(goshawk_urban$Temp_age, na.rm = TRUE),
                                                             Habitat = "urban"),
                        re.form = NA, intervals = "predVar") # predVar = SD of Prediction^2
pred_age <- data.frame(pred = pred_age_obj, attr(pred_age_obj, "intervals"), age = predictor_age)

predictor_age2 <- seq(min(goshawk_urban$Age, na.rm = TRUE),
                      max(goshawk_urban$Age, na.rm = TRUE),
                      length = 30)
pred_age_obj2 <- predict(glmm_clinical,
                         newdata = expand.grid(Age = predictor_age2,
                                               No_nestlings = 3,
                                               Year = "2015",
                                               Sex = "female",
                                               laying_day = median(goshawk_urban$laying_day, na.rm = TRUE),
                                               Temp_age = median(goshawk_urban$Temp_age, na.rm = TRUE),
                                               Habitat = "rural"),
                         re.form = NA, intervals = "predVar") # predVar = SD of Prediction^2
pred_age2 <- data.frame(pred2 = pred_age_obj2, attr(pred_age_obj2, "intervals"), age2 = predictor_age2)


pdf("./figures/prediction_plot4.pdf", width = 10, height = 7)

plot(pred_age$pred ~ pred_age$age, type = "l", lwd = 2,
     ylab = "Proportion of nestlings presenting clinical signs", las = 1, ylim = c(0, 1), xlab = "Age of nestlings")
points(pred_age$predVar_0.025~ pred_age$age, type = "l", lwd = 2, lty = 2)
points(pred_age$predVar_0.975~ pred_age$age, type = "l", lwd = 2, lty = 2)
rug(x = jitter(urban2$Age[urban2$Clinical_binary == "0"]), side = 1)
rug(x = jitter(urban2$Age[urban2$Clinical_binary == "1"]), side = 3)

lines(pred_age2$pred2 ~ pred_age2$age2, type = "l", lwd = 2,
      ylab = "Proportion of nestlings presenting clinical signs", las = 1, ylim = c(0, 1), xlab = "Age of nestlings", col = "grey")
points(pred_age2$predVar_0.025~ pred_age2$age2, type = "l", lwd = 2, lty = 2, col = "grey")
points(pred_age2$predVar_0.975~ pred_age2$age2, type = "l", lwd = 2, lty = 2, col = "grey")
rug(x = jitter(rural2$Age[rural2$Clinical_binary == "0"]), side = 1, col = "grey")
rug(x = jitter(rural2$Age[rural2$Clinical_binary == "1"]), side = 3, col = "grey")

dev.off()


###################### Causes of mortality ############################################################
###### tricho######
Todesursachen_Statistik <- read_excel("./source_data/goshawk_data_death.xlsx")

death <- Todesursachen_Statistik[, c("age", "location", "Trichomonosiasis",
                                     "Trauma (window)", "sex")]
death$location <- as.factor(death$location)
death$Trauma_window <- as.numeric(death$`Trauma (window)`)
death$age <- as.factor(death$age)
death$sex <- as.factor(death$sex)
death <- as.data.frame(death)
death <- droplevels(na.omit(death)) 
str(death)

glm_death <- fitme(Trichomonosiasis ~ location + age + sex,
                   family = binomial(link = "logit"),
                   data = death,
                   method = "PQL/L")


# test againt zero model
GLM_0 <- fitme(Trichomonosiasis ~ 1,
               family = binomial(link = "logit"),
               data = death,
               method = "PQL/L")
anova(GLM_0, glm_death)

# Model output
glm_death

# effects

glm_death_no_loc <- fitme(Trichomonosiasis ~ age + sex,
                          family = binomial(link = "logit"),
                          data = death, method = "PQL/L")
anova(glm_death_no_loc, glm_death)

glm_death_no_age <- fitme(Trichomonosiasis ~ location + sex,
                          family = binomial(link = "logit"),
                          data = death,
                          method = "PQL/L")
anova(glm_death_no_age, glm_death)

glm_death_no_sex <- fitme(Trichomonosiasis ~ location + age,
                          family = binomial(link = "logit"),
                          data = death,
                          method = "PQL/L")
anova(glm_death_no_sex, glm_death)

# odds ratio
exp(glm_death$fixef["locationUrban"]) 

## bootstraps
anova(glm_death_no_loc, glm_death, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glm_death_no_age, glm_death, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glm_death_no_sex, glm_death, boot.repl = boot.repl, nb_cores = nb_cores)

####### window####

glm_death2 <- fitme(Trauma_window ~ location + sex,
                    family = binomial(link = "logit"),
                    data = death,
                    method = "PQL/L")


# test againt zero model
GLM_0_2 <- fitme(Trauma_window ~ 1,
             family = binomial(link = "logit"),
             data = death,
             method = "PQL/L")
anova(GLM_0_2, glm_death2)

# Model output
glm_death2

# effects
glm_death_no_loc2 <- fitme(Trauma_window ~  sex,
                           family = binomial(link = "logit"),
                           data = death,
                           method = "PQL/L")
anova(glm_death_no_loc2, glm_death2)

glm_death_no_sex2 <- fitme(Trauma_window ~ location,
                           family = binomial(link = "logit"),
                           data = death,
                           method = "PQL/L")
anova(glm_death_no_sex2, glm_death2)

# odds ratio
exp(glm_death2$fixef["locationUrban"]) 
exp(glm_death2$fixef["sexm"]) 

## bootstraps
anova(glm_death_no_loc2, glm_death2, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glm_death_no_sex2, glm_death2, boot.repl = boot.repl, nb_cores = nb_cores)


##### figure ###

Todesursachen_Statistik <- read_excel("./source_data/goshawk_data_death.xlsx")
Todesursachen_Statistik <- Todesursachen_Statistik %>% 
  mutate(Cause_of_death = ifelse(Cause_of_death == "Unknwon", "Unknown", 
                                 ifelse(Cause_of_death == "Trichmononosis", "Trichomonosis", 
                                        ifelse(Cause_of_death == "Trauma (unknwon)" | 
                                                 Cause_of_death == "Trauma (unkown)",
                                               "Trauma (unknown)", as.character(Cause_of_death)))))  %>%
  mutate(Cause_of_death = as.factor(Cause_of_death))


levels(Todesursachen_Statistik$Cause_of_death)
Todesursachen_Statistik$Cause_of_death <- factor(Todesursachen_Statistik$Cause_of_death, 
                                                 levels = levels(Todesursachen_Statistik$Cause_of_death)[c(13, 11, 12, 14, 2,7,1, 9,3, 10,8,5,4, 6, 15)])


# variable proportion~causes_of_death
Todesursachen_Statistik$location <- relevel(as.factor(Todesursachen_Statistik$location), ref = "Urban")
total_urban <- nrow(Todesursachen_Statistik[Todesursachen_Statistik$location == "Urban", ])
total_rural <- nrow(Todesursachen_Statistik[Todesursachen_Statistik$location == "Rural", ])

data_plot2 <- Todesursachen_Statistik %>% 
  group_by(Cause_of_death, location) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(total = ifelse(location == "Urban", total_urban, total_rural), 
         prop = n / total) %>%  
  mutate(CI_lwr = map2_dbl(n, total, ~ as.numeric(binom.confint( x = .x, n = .y)[5, 5])), ## method exact
         CI_upr = map2_dbl(n, total, ~ as.numeric(binom.confint( x = .x, n = .y)[5, 6]))) ## method exact


ggplot(data_plot2, aes(y = prop, x = Cause_of_death, fill = location, ymin = CI_lwr, ymax = CI_upr))+ 
  geom_bar(stat = "identity", position =  position_dodge2(width = 0.9, preserve = "single"), alpha = 0.8) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_grey("Habitat", 
                  labels = c("urban", "rural")) +
  scale_y_continuous("Proportion") +
  scale_x_discrete("Cause of mortality") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust =  1), 
        panel.grid.minor.x = element_blank()) 
ggsave(filename = "./figures/cause_of_death_plot.pdf", width = 10, height = 8)




###################### Figures supporting information #####################################################
### Figure SI1

rural <- data.frame(long = c(6.12954, 9.563227, 13.699830, 8.5324708), 
                    lat = c(51.786726, 54.5239312, 52.857445, 52.0302285), 
                    names = c('Kleve', 'Schleswig', 'Barnim', 'Bielefeld'), stringsAsFactors = FALSE)
urban <- data.frame(long = c(13.413215, 6.958281, 9.993682), 
                    lat = c(52.521918, 50.941278, 53.551085), 
                    names = c('Berlin', 'Cologne', 'Hamburg'), stringsAsFactors = FALSE)

germany <- map_data("worldHires", "Germany")
Germ < -ggplot() + 
  geom_polygon(data = germany, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
  coord_fixed(1.5) +
  geom_point(data = rural, aes(x = long, y = lat), size = 3, shape = 17) +
  geom_point(data = urban, aes(x = long, y = lat), size = 3,shape = 15) + 
  geom_text(data = rural, aes(x = long + 0.1, y = lat + 0.45, label = names, fontface = 2),
            color = 'black', size = 3) +
  geom_text(data = urban, aes(x = long + 0.1, y = lat - 0.33, label = names, fontface = 2),
            color = 'black', size = 3) +
  labs(x = "Longitude") + labs(y = "Latitude") +
  theme_light() +
  theme(panel.grid.minor = element_blank())
Germ
ggsave(filename = "./figures/map_plot.pdf", width = 7, height = 10)


### Figures SI2

table_per_nest <- read_excel("./source_data/goshawk_data_nest.xlsx")

goshawk_imper <- subset(table_per_nest, Species == "Goshawk")
goshawk_imper <- as.data.frame(goshawk_imper)


ggplot(goshawk_imper, aes(x = Habitat, y = imperviousness2500)) +
  geom_boxplot() +
  theme_bw() +
  ylab("\n Imperviousness (%) \n")
ggsave(filename = "./figures/box_plot1.pdf", width = 8, height = 8)

