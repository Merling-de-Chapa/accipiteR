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
library(vegan)

set.seed(123L)
spaMM.options(nb_cores = 3, separation_max = 1)
boot.repl <- 1000 ## put 0 for not parametric bootstrap!
nb_cores <- 3

## map of study area
rural <- data.frame(long = c(6.12954, 9.563227, 13.699830, 8.5324708), 
                    lat = c(51.786726, 54.5239312, 52.857445, 52.0302285), 
                    names = c('Kleve', 'Schleswig', 'Barnim', 'Bielefeld'), stringsAsFactors = FALSE)
urban <- data.frame(long = c(13.413215, 6.958281, 9.993682), 
                    lat = c(52.521918, 50.941278, 53.551085), 
                    names = c('Berlin', 'Cologne', 'Hamburg'), stringsAsFactors = FALSE)

germany <- map_data("worldHires", "Germany")
Germ <- ggplot() + 
  geom_polygon(data = germany, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
  coord_fixed(1.5) +
  geom_point(data = rural, aes(x = long, y = lat), size = 3, shape = 17) +
  geom_point(data = urban, aes(x = long, y = lat), size = 3,shape = 15) + 
  geom_text(data = rural, aes(x = long + 0.1, y = lat + 0.45, label = names, fontface = 2), color = 'black', size = 3) +
  geom_text(data = urban, aes(x = long + 0.1, y = lat - 0.33, label = names, fontface = 2), color = 'black', size = 3) +
  labs(x = "Longitude") + labs(y = "Latitude") +
  theme_light() +
  theme(panel.grid.minor = element_blank())
Germ
ggsave(filename = "./figures/figure1.pdf", width = 7, height = 10)


############################# behavioural flexibility #####################################################

## behaviour model
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

# subset Urban/Rural
urban <- subset(goshawk_nest, Habitat == "urban")
rural <- subset(goshawk_nest, Habitat == "rural")

# median Age
median(urban$Age_youngest)
# 21
median(rural$Age_youngest)
# 20

nrow(urban)
# 77
nrow(rural)
# 74

wilcox.test(Age_youngest ~ Habitat, data = table_per_nest)
# Wilcoxon rank sum test with continuity correction
# 
# data:  Age_youngest by Habitat
# W = 4315.5, p-value = 0.2681
# alternative hypothesis: true location shift is not equal to 0

#GLMM
test_react_int <- fitme(reaction ~ Habitat*Age + No_nestlings + Laying_begin_day + 
                        Year + Rainfall + (1|Location/Territory),
                        family = binomial(link = "logit"), data = goshawk_nest, 
                        method = "PQL/L") ## FULL MODEL

test_no_Hab_age <- fitme(reaction ~ No_nestlings +  Laying_begin_day + Year + 
                         Rainfall + (1|Location/Territory), 
                         family = binomial(link = "logit"), data = goshawk_nest, 
                         method = "PQL/L")

test_noHabitat <- fitme(reaction ~ Age + No_nestlings +  Laying_begin_day + Year + 
                        Rainfall + (1|Location/Territory),
                        family = binomial(link = "logit"), data = goshawk_nest, 
                        method = "PQL/L")

test_nonest <- fitme(reaction ~ Habitat*Age + Laying_begin_day + Year +  Rainfall 
                     + (1|Location/Territory), family = binomial(link = "logit"), 
                     data = goshawk_nest, method = "PQL/L")

test_noYear <- fitme(reaction ~ Habitat*Age + No_nestlings + Laying_begin_day + 
                     Rainfall + (1|Location/Territory), 
                     family = binomial(link = "logit"), data = goshawk_nest, 
                     method = "PQL/L")

test_noDay <- fitme(reaction ~ Habitat*Age + No_nestlings +  Year + Rainfall 
                    + (1|Location/Territory), family = binomial(link = "logit"), 
                    data = goshawk_nest, method = "PQL/L")

test_noRain <- fitme(reaction ~ Habitat*Age + No_nestlings + Laying_begin_day + 
                     Year + (1|Location/Territory), family = binomial(link = "logit"),
                     data = goshawk_nest, method = "PQL/L")

test_noYoung <- fitme(reaction ~ Habitat + Rainfall + No_nestlings + 
                      Laying_begin_day + Year + (1|Location/Territory),
                      family = binomial(link = "logit"), data = goshawk_nest, 
                      method = "PQL/L")

#### Habitat

## Odds_ratioS 
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
odd_ratio_urban_young_vs_rural_young # [1] 10.81
odd_ratio_urban_old_vs_rural_old <- round(OR(prob_urban_old, prob_rural_old), 2)
odd_ratio_urban_old_vs_rural_old # [1] 23.46
odd_ratio_urban_young_vs_urban_old <- round(OR(prob_urban_young, prob_urban_old), 2)
odd_ratio_urban_young_vs_urban_old # [1] 2.87
odd_ratio_rural_young_vs_rural_old <- round(OR(prob_rural_young, prob_rural_old), 2)
odd_ratio_rural_young_vs_rural_old # [1] 6.24

#### bootstrap to get p-values of each effect
anova(test_react_int, test_noHabitat, boot.repl = boot.repl, nb_cores = nb_cores)
# bootstrap took 1907.8 s.
# chi2_LR df     p_value
# p_v 12.26368  2 0.002172577
# ======== Bootstrap: ========
#   Raw simulated p-value: 0.016
# Bartlett-corrected LR test:
#   chi2_LR df    p_value
# p_v 9.145138  2 0.01033138

anova(test_react_int, test_nonest, boot.repl = boot.repl, nb_cores = nb_cores)
# bootstrap took 3288.4 s.
# chi2_LR df   p_value
# p_v 1.036544  1 0.3086266
# ======== Bootstrap: ========
#   Raw simulated p-value: 0.337
# Bartlett-corrected LR test:
#   chi2_LR df   p_value
# p_v 0.8949581  1 0.3441372

anova(test_react_int, test_noYear, boot.repl = boot.repl, nb_cores = nb_cores)
# bootstrap took 2801.6 s.
# chi2_LR df    p_value
# p_v 8.014984  1 0.00463919
# ======== Bootstrap: ========
#   Raw simulated p-value: 0.011
# Bartlett-corrected LR test:
#   chi2_LR df     p_value
# p_v 7.454125  1 0.006329123

anova(test_react_int, test_noDay, boot.repl = boot.repl, nb_cores = nb_cores)
# bootstrap took 2297.4 s.
# chi2_LR df    p_value
# p_v 4.940328  1 0.02623706
# ======== Bootstrap: ========
#   Raw simulated p-value: 0.03
# Bartlett-corrected LR test:
#   chi2_LR df   p_value
# p_v 4.569614  1 0.0325439

anova(test_react_int, test_noRain, boot.repl = boot.repl, nb_cores = nb_cores)
# bootstrap took 1967.2 s.
# chi2_LR df   p_value
# p_v 0.6995435  1 0.4029371
# ======== Bootstrap: ========
#   Raw simulated p-value: 0.427
# Bartlett-corrected LR test:
#   chi2_LR df   p_value
# p_v 0.6404846  1 0.4235354

anova(test_react_int, test_noYoung, boot.repl = boot.repl, nb_cores = nb_cores)
# bootstrap took 2387 s.
# chi2_LR df    p_value
# p_v 5.979197  2 0.05030763
# ======== Bootstrap: ========
#   Raw simulated p-value: 0.0829
# Bartlett-corrected LR test:
#   chi2_LR df    p_value
# p_v 4.970473  2 0.08330584

## Odds ratio for Year effect
exp(fixef(test_react_int)["Year2016"])


##### test without interaction ###############

test_no_int <- fitme(reaction ~ Habitat + Age + No_nestlings + Laying_begin_day + 
                     Year + Rainfall + (1|Location/Territory),
                     family = binomial(link = "logit"), data = goshawk_nest, 
                     method = "PQL/L")
round(exp(fixef(test_no_int)["Habitaturban"]), 2)
# Habitaturban (Odds-ratio)
# 21.68
CI_behaviour <- confint(test_no_int, "Habitaturban")
exp(CI_behaviour$interval)

test_no_Habitat <- fitme(reaction ~ Age + No_nestlings + Laying_begin_day + Year + 
                         Rainfall + (1|Location/Territory),
                         family = binomial(link = "logit"), data = goshawk_nest, 
                         method = "PQL/L")

test_no_Age <- fitme(reaction ~ Habitat + No_nestlings + Laying_begin_day + Year + 
                     Rainfall + (1|Location/Territory),
                     family = binomial(link = "logit"), data = goshawk_nest, 
                     method = "PQL/L")

anova(test_react_int, test_no_int, boot.repl = boot.repl, nb_cores = nb_cores)
# bootstrap took 2789.7 s.
# chi2_LR df   p_value
# p_v 0.2841157  1 0.5940163
# ======== Bootstrap: ========
#   Raw simulated p-value: 0.736
# Bartlett-corrected LR test:
#   chi2_LR df   p_value
# p_v 0.2414185  1 0.6231834

anova(test_no_int, test_no_Age, boot.repl = boot.repl, nb_cores = nb_cores)
# bootstrap took 2018.7 s.
# chi2_LR df    p_value
# p_v 5.695081  1 0.01701252
# ======== Bootstrap: ========
#   Raw simulated p-value: 0.015
# Bartlett-corrected LR test:
#   chi2_LR df    p_value
# p_v 5.587429  1 0.01808983

anova(test_no_int, test_no_Habitat, boot.repl = boot.repl, nb_cores = nb_cores)
# bootstrap took 1943.1 s.
# chi2_LR df      p_value
# p_v 11.97957  1 0.0005378708
# ======== Bootstrap: ========
#   Raw simulated p-value: 0.005
# Bartlett-corrected LR test:
#   chi2_LR df     p_value
# p_v 7.658944  1 0.005649156

#### figures #####################
#### barplot raw data

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

ggsave(filename = "./figures/figure2.pdf", width = 8, height =  7)


# barplot predictions Habitat/Age

pred_urban_young <- predict(test_react_int,
                            newdata = expand.grid(Laying_begin_day = median(goshawk_nest$Laying_begin_day, na.rm = TRUE),
                                                  Age = "young",
                                                  Rainfall = median(goshawk_nest$Rainfall, na.rm = TRUE),
                                                  No_nestlings = 3,
                                                  Year = "2015",
                                                  Habitat = "urban"),
                            re.form = NA, intervals = "predVar") 
pred_1 <- data.frame(pred = pred_urban_young, attr(pred_urban_young, "intervals"), 
                     Habitat = "urban", Age = "young")

pred_urban_old <- predict(test_react_int, newdata = expand.grid(Laying_begin_day = median(goshawk_nest$Laying_begin_day, na.rm = TRUE),
                                                                Age = "old",
                                                                Rainfall = median(goshawk_nest$Rainfall, na.rm = TRUE),
                                                                No_nestlings = 3,
                                                                Year = "2015",
                                                                Habitat = "urban"),
                          re.form = NA, intervals = "predVar") 
pred_2 <- data.frame(pred = pred_urban_old, attr(pred_urban_old, "intervals"), 
                     Habitat = "urban", Age = "old")

pred_rural_young <- predict(test_react_int, newdata = expand.grid(Laying_begin_day = median(goshawk_nest$Laying_begin_day, na.rm = TRUE),
                                                                  Age = "young",
                                                                  Rainfall = median(goshawk_nest$Rainfall, na.rm = TRUE),
                                                                  No_nestlings = 3,
                                                                  Year = "2015",
                                                                  Habitat = "rural"),
                            re.form = NA, intervals = "predVar") 
pred_3 <- data.frame(pred = pred_rural_young, attr(pred_rural_young, "intervals"), 
                     Habitat = "rural", Age = "young")

pred_rural_old <- predict(test_react_int, newdata = expand.grid(Laying_begin_day = median(goshawk_nest$Laying_begin_day, na.rm = TRUE),
                                                                Age = "old",
                                                                Rainfall = median(goshawk_nest$Rainfall, na.rm = TRUE),
                                                                No_nestlings = 3,
                                                                Year = "2015",
                                                                Habitat = "rural"),
                          re.form = NA, intervals = "predVar") 
pred_4 <- data.frame(pred = pred_rural_old, attr(pred_rural_old, "intervals"), 
                     Habitat = "rural", Age = "old")

pred_Habitat_age <- rbind(pred_1, pred_2, pred_3, pred_4)
pred_Habitat_age$Habitat <- factor(pred_Habitat_age$Habitat, levels = c("urban", "rural"))
pred_Habitat_age$Age <- factor(pred_Habitat_age$Age, levels = c("young", "old"))

ggplot(pred_Habitat_age, aes(x = Age, y = pred, fill = Habitat, ymin = predVar_0.025, ymax = predVar_0.975)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.8) +
  geom_errorbar(width = 0, position = position_dodge(preserve = "total", width = 0.9)) +
  scale_fill_grey() +
  scale_y_continuous("Proportion of females reacting") +
  scale_x_discrete("Age categorical") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank()) 

ggsave(filename = "./figures/figure3.pdf", width = 8, height =  7)

# plot predictions laying

predictor_lay <- seq(min(goshawk_nest$Laying_begin_day, na.rm = TRUE), max(goshawk_nest$Laying_begin_day, na.rm = TRUE), length = 30)
pred_lay_obj <- predict(test_no_int,
                        newdata = expand.grid(Laying_begin_day = predictor_lay,
                                              Age = "young",
                                              Rainfall = median(goshawk_nest$Rainfall, na.rm = TRUE),
                                              No_nestlings = 3,
                                              Year = "2015",
                                              Habitat = "urban"),
                        re.form = NA, intervals = "predVar") 
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
                         re.form = NA, intervals = "predVar") 
pred_lay2 <- data.frame(pred = pred_lay_obj2, attr(pred_lay_obj2, "intervals"), lay2 = predictor_lay2)


pdf("./figures/figure4.pdf", width = 10, height = 7)
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

Rupfungen_gesamt <- read_excel("./source_data/goshawk_data_diet_plot.xlsx")
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

ggsave(filename = "./figures/figure5.pdf", width = 10, height =  7)

########## pigeons ###################
Rupfungen_per_nest_aktuell <- read_excel("./source_data/goshawk_data_diet.xlsx")

pigeon <- Rupfungen_per_nest_aktuell[, c("Location", "nest_side", "Habitat",
                                         "Species_richness_nest", "Pigeon", "Rest")]
pigeon$Location <- as.factor(pigeon$Location)
pigeon$nest_side <- as.factor(pigeon$nest_side)
pigeon$Habitat <- as.factor(pigeon$Habitat)
pigeon <- as.data.frame(pigeon)

GLMM_pigeon_spaMM_loc <- fitme(cbind(Pigeon, Rest)~Habitat + (1|Location/nest_side),
                               family = binomial(link = "logit"),
                               data = pigeon,
                               method = "PQL/L")

GLMM_0_loc <- fitme(cbind(Pigeon, Rest)~ 1 + (1|Location/nest_side),
                    family = binomial(link = "logit"),
                    data = pigeon,
                    method = "PQL/L")

# prediction
Odds_ratio <- exp(GLMM_pigeon_spaMM_loc$fixef["HabitatUrban"])
Odds_ratio
# HabitatUrban 
# 3.637754
CI_diet <- confint(GLMM_pigeon_spaMM_loc, "HabitatUrban")
exp(CI_diet$interval)
# lower HabitatUrban upper HabitatUrban 
# 2.436699           5.281086 

#### bootstrap ###
set.seed(1)
anova(GLMM_0_loc, GLMM_pigeon_spaMM_loc, boot.repl = boot.repl, nb_cores = nb_cores) 


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

## Here we use lme4 as an interim step to estimate the BoxCox transformation to 
## be applied in spaMM
lmer_diversity <- lmer(diversity_nest ~ Habitat + (1|Location), 
                       data = diversity, REML = FALSE)

spaMM_diversity_coor <- fitme(diversity_nest ~ Habitat +(1|Location), data = diversity)

bc_lmer <- powerTransform(lmer_diversity, family = "bcPower")

lmer_diversity2 <- update(lmer_diversity, bcPower(diversity_nest, bc_lmer$lambda) ~ .)
spaMM_diversity3 <- update(spaMM_diversity_coor, bcPower(diversity_nest, bc_lmer$lambda) ~ .)

logLik(lmer_diversity2)
logLik(spaMM_diversity3) ## same loglik across packages! GOOD

lmer_diversity_0 <- lmer(diversity_nest ~ 1 + (1|Location), data = diversity, 
                         REML = FALSE)
spaMM_diversity_0 <- fitme(diversity_nest ~ 1 + (1|Location), data = diversity)

lmer_diversity_upda <- update(lmer_diversity_0, bcPower(diversity_nest, bc_lmer$lambda) ~ .)
spaMM_diversity_upda_0 <- update(spaMM_diversity_0, bcPower(diversity_nest, bc_lmer$lambda) ~ .)

### effect ########
pred_trans_rural <- (spaMM_diversity3$fixef["(Intercept)"] * bc_lmer$lambda + 1)^(1/bc_lmer$lambda)
pred_trans_urban <- ((spaMM_diversity3$fixef["(Intercept)"] + spaMM_diversity3$fixef["HabitatUrban"]) * bc_lmer$lambda + 1)^(1/bc_lmer$lambda)
pred_trans_urban - pred_trans_rural
# (Intercept) 
# -0.1099321 

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

# Habitat Effect
predict(glm_species_spaMM, newdata = data.frame(Habitat = c("Urban", "Rural")), 
        re.form = NA)
# Point predictions:
# 4.017246 5.457141 

diff(predict(glm_species_spaMM, newdata = data.frame(Habitat = c("Urban", "Rural")),
             re.form = NA))
# Point predictions:
# 1.439895


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

lmm_laying2 <- fitme(Laying_begin_day ~ Habitat + Temp_breeding_begin  + 
                     (1|Location/Territory), family = gaussian, data = goshawk_nest2)

# odds
lmm_laying2$fixef["Habitaturban"] 
# Habitaturban 
# -12.49214
CI_laying <- confint(lmm_laying2, "Habitaturban")
CI_laying
#lower Habitaturban upper Habitaturban 
#-17.365852          -7.121119

### test  effects
lmm_laying_no_Habitat <- fitme(Laying_begin_day ~ Temp_breeding_begin + 
                              (1|Location/Territory), family = gaussian,
                              data = goshawk_nest2)

lmm_laying_no_Temp <- fitme(Laying_begin_day ~ Habitat + (1|Location/Territory),
                            family = gaussian, data = goshawk_nest2)

### bootstrap
anova(lmm_laying2, lmm_laying_no_Habitat, boot.repl = boot.repl, nb_cores = nb_cores)
anova(lmm_laying2, lmm_laying_no_Temp, boot.repl = boot.repl, nb_cores = nb_cores)


########## reproductive output (number of nestlings) ############
table_per_nest <- read_excel("./source_data/goshawk_data_nest.xlsx")

## wilcox Number of nestlings
wilcox.test(No_nestlings ~ Habitat, data = table_per_nest)
# Wilcoxon rank sum test with continuity correction
# 
# data:  No_nestlings by Habitat
# W = 3517.5, p-value = 0.0006975
# alternative hypothesis: true location shift is not equal to 0

## Figure Number of nestlings per nest

n_land <- nrow(table_per_nest[table_per_nest$Habitat == "rural", ])
n_stadt <- nrow(table_per_nest[table_per_nest$Habitat == "urban", ])

cut_off <- table_per_nest %>% 
  group_by(No_nestlings, Habitat) %>%
  summarise(n = n()) %>% 
  mutate(n_total = ifelse(Habitat == "urban", n_land, n_stadt), 
         prop = n / n_total) %>% 
  ungroup()

cut_off <- cut_off  %>% 
  mutate(CI_lwr = map2_dbl(n, n_total, ~ as.numeric(binom.confint( x = .x, n = .y)[5, 5])), ## method exact
         CI_upr = map2_dbl(n, n_total, ~ as.numeric(binom.confint( x = .x, n = .y)[5, 6]))) ## method exact

cut_off$No_nestlings <- relevel(as.factor(cut_off$No_nestlings), ref = "1")
cut_off$Habitat <- factor(cut_off$Habitat, levels = c("urban", "rural"))
cut_off$Habitat <- relevel(as.factor(cut_off$Habitat), ref = "urban")

ggplot(cut_off, aes(x = No_nestlings, y = prop, fill = Habitat, ymin = CI_lwr, ymax = CI_upr)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.8) +
  geom_errorbar(position = position_dodge(preserve = "total", width = 0.9), width = 0) +
  scale_fill_grey(labels = c("urban", "rural")) +
  scale_y_continuous("Proportion") +
  scale_x_discrete("Number of nestlings") +
  geom_vline(xintercept = 2.5, linetype = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 0, hjust =  1), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) 

ggsave(filename = "./figures/figure6.pdf", width = 10, height = 7)


##### GLMM
goshawk_nest3 <- table_per_nest[, c("Year", "Location", "No_nestlings", "Laying_begin_day",
                                    "Territory", "Habitat", "Temp_breeding_begin")]
goshawk_nest3$Year <- as.factor(goshawk_nest3$Year)
goshawk_nest3$Location <- as.factor(goshawk_nest3$Location)
goshawk_nest3$Territory <- as.factor(goshawk_nest3$Territory)
goshawk_nest3$Habitat <- as.factor(goshawk_nest3$Habitat)
goshawk_nest3 <- as.data.frame(goshawk_nest3)
goshawk_nest3 <- droplevels(na.omit(goshawk_nest3))

goshawk_nest3$No_nestlings_binary <- goshawk_nest3$No_nestlings > 2

glmm_nestlings <- fitme(No_nestlings_binary ~ Habitat + Laying_begin_day + 
                        Temp_breeding_begin + (1|Location/Territory),
                        family = binomial(link = "logit"), data = goshawk_nest3, 
                        method = "PQL/L")

# effects
glmm_nestlings_no_Habitat <- fitme(No_nestlings_binary ~ Laying_begin_day + 
                                   Temp_breeding_begin + (1|Location/Territory),
                                   family = binomial(link = "logit"), data = goshawk_nest3, 
                                   method = "PQL/L")

glmm_nestlings_no_day <- fitme(No_nestlings_binary ~ Habitat + Temp_breeding_begin + 
                              (1|Location/Territory), family = binomial(link = "logit"), 
                              data = goshawk_nest3, method = "PQL/L")

glmm_nestlings_no_Temp <- fitme(No_nestlings_binary ~ Habitat + Laying_begin_day + 
                                (1|Location/Territory), family = binomial(link = "logit"), 
                                data = goshawk_nest3, method = "PQL/L")

# bootstraps
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
                        re.form = NA, intervals = "predVar") 
pred_day <- data.frame(pred = pred_day_obj, attr(pred_day_obj, "intervals"), day = predictor_day)

predictor_day2 <- seq(min(goshawk_nest3$Laying_begin_day, na.rm = TRUE), max(goshawk_nest3$Laying_begin_day, na.rm = TRUE), length = 30)
pred_day_obj2 <- predict(glmm_nestlings, newdata = expand.grid(Laying_begin_day = predictor_day2,
                                                               Habitat = "rural",
                                                               Temp_breeding_begin = median(goshawk_nest3$Temp_breeding_begin, na.rm = TRUE)),
                         re.form = NA, intervals = "predVar") 
pred_day2 <- data.frame(pred2 = pred_day_obj2, attr(pred_day_obj2, "intervals"), day2 = predictor_day2)

pdf("./figures/figure7.pdf", width = 10, height = 7)
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

# Odd
exp(glmm_nestlings2$fixef["Habitaturban"])
# Habitaturban 
# 2.217902
CI_nestlings <- confint(glmm_nestlings2, "Habitaturban")
exp(CI_nestlings$interval)

#### testing effects########
glmm_nestlings_no_Habitat2 <- fitme(No_nestlings_binary ~ Temp_breeding_begin + 
                                    (1|Location/Territory), family = binomial(link = "logit"),
                                    data = goshawk_nest3, method = "PQL/L")

glmm_nestlings_no_Temp2 <- fitme(No_nestlings_binary ~ Habitat + (1|Location/Territory),
                                 family = binomial(link = "logit"),
                                 data = goshawk_nest3,
                                 method = "PQL/L")

## bootstraps
anova(glmm_nestlings2, glmm_nestlings_no_Habitat2, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_nestlings2, glmm_nestlings_no_Temp2, boot.repl = boot.repl, nb_cores = nb_cores)


###################### Health status of goshawk nestlings ####################################

#### prevalence ####
Tabelle_Statistik_urban_ecology_aktuell <- read_excel("./source_data/goshawk_data_nestlings.xlsx")

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

# subset urban/rural
urban2 <- subset(goshawk_urban, Habitat == "urban")
rural2 <- subset(goshawk_urban, Habitat == "rural")

# Median age
median(urban2$Age)
# 24
median(rural2$Age)
# 23

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
ggsave(filename = "./figures/figure8.pdf", width = 8, height = 8)

## Analysis trichomonas
glmm_tricho_logit <- fitme(Prevalence ~ Habitat + Age + No_nestlings +  Sex + Temp_age + 
                           Year + laying_day +(1|Location/Territory),
                           family = binomial(link = "logit"),
                           data = goshawk_urban,
                           method = "PQL/L")

# odds
exp(glmm_tricho_logit$fixef["Habitaturban"])
# Habitaturban 
# 2.833019
CI_tricho <- confint(glmm_tricho_logit, "Habitaturban")
exp(CI_tricho$interval)

#test effects
glmm_tricho_logit_noHabitat <- fitme(Prevalence ~  Age + No_nestlings + Sex + Temp_age + 
                                     Year + laying_day + (1|Location/Territory),
                                     family = binomial(link = "logit"), 
                                     data = goshawk_urban, method = "PQL/L")

glmm_tricho_logit_nonest <- fitme(Prevalence ~ Habitat + Age + Sex + Temp_age + Year + 
                                  laying_day + (1|Location/Territory),
                                  family = binomial(link = "logit"), 
                                  data = goshawk_urban, method = "PQL/L")

glmm_tricho_logit_noage <- fitme(Prevalence ~ Habitat + No_nestlings + Sex + Temp_age + 
                                 Year + laying_day + (1|Location/Territory),
                                 family = binomial(link = "logit"), data = goshawk_urban, 
                                 method = "PQL/L")

glmm_tricho_logit_noTemp <- fitme(Prevalence ~ Habitat + Age + No_nestlings + Sex + 
                                  Year + laying_day + (1|Location/Territory),
                                  family = binomial(link = "logit"), data = goshawk_urban, 
                                  method = "PQL/L")

glmm_tricho_logit_noSex <- fitme(Prevalence ~ Habitat + Age + No_nestlings + Year + 
                                 laying_day + Temp_age + (1|Location/Territory),
                                 family = binomial(link = "logit"), data = goshawk_urban, 
                                 method = "PQL/L")

glmm_tricho_logit_noYear <- fitme(Prevalence ~ Habitat + Age + No_nestlings + Sex + 
                                  Temp_age + laying_day + (1|Location/Territory),
                                  family = binomial(link = "logit"), data = goshawk_urban, 
                                  method = "PQL/L")

glmm_tricho_logit_noDay <- fitme(Prevalence ~ Habitat + Age + No_nestlings + Sex + 
                                 Temp_age + Year + (1|Location/Territory),
                                 family = binomial(link = "logit"), data = goshawk_urban,
                                 method = "PQL/L")

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

# odds
exp(glmm_clinical$fixef["Habitaturban"])
# Habitaturban 
# 10.89108 
CI_clinical <- confint(glmm_clinical, "Habitaturban")
exp(CI_clinical$interval)

#test effects
glmm_clinical_noHabitat <- fitme(Clinical_binary ~  Age + No_nestlings + Sex + Year + 
                                 laying_day + Temp_age + (1|Location/Territory),
                                 family = binomial(link = "logit"), data = goshawk_urban, 
                                 method = "PQL/L")

glmm_clinical_nonest <- fitme(Clinical_binary ~ Habitat + Age + Sex + Year + laying_day + 
                              Temp_age + (1|Location/Territory),
                              family = binomial(link = "logit"), data = goshawk_urban, 
                              method = "PQL/L")

glmm_clinical_noage <- fitme(Clinical_binary ~ Habitat + No_nestlings + Sex + Year + 
                             laying_day + Temp_age + (1|Location/Territory),
                             family = binomial(link = "logit"), data = goshawk_urban, 
                             method = "PQL/L")

glmm_clinical_noSex <- fitme(Clinical_binary ~ Habitat + Age + No_nestlings + Year + 
                             laying_day + Temp_age + (1|Location/Territory),
                             family = binomial(link = "logit"), data = goshawk_urban, 
                             method = "PQL/L")

glmm_clinical_noTemp <- fitme(Clinical_binary ~ Habitat + Age + No_nestlings + Sex + 
                             Year + laying_day + (1|Location/Territory),
                             family = binomial(link = "logit"), data = goshawk_urban, 
                             method = "PQL/L")

glmm_clinical_noYear <- fitme(Clinical_binary ~ Habitat + Age + No_nestlings + Sex + 
                              laying_day + Temp_age + (1|Location/Territory),
                              family = binomial(link = "logit"), data = goshawk_urban, 
                              method = "PQL/L")

glmm_clinical_noDay <- fitme(Clinical_binary ~ Habitat + Age + No_nestlings + Sex + Year + 
                             Temp_age + (1|Location/Territory),
                             family = binomial(link = "logit"), data = goshawk_urban, 
                             method = "PQL/L")

## bootstraps
anova(glmm_clinical, glmm_clinical_noHabitat, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_clinical, glmm_clinical_nonest, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_clinical, glmm_clinical_noage, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_clinical, glmm_clinical_noSex, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_clinical, glmm_clinical_noTemp, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_clinical, glmm_clinical_noYear, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glmm_clinical, glmm_clinical_noDay, boot.repl = boot.repl, nb_cores = nb_cores)

#### plots predictions clinical signs ######

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
                         re.form = NA, intervals = "predVar") 
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
                         re.form = NA, intervals = "predVar") 
pred_lay4 <- data.frame(pred4 = pred_lay_obj4, attr(pred_lay_obj4, "intervals"), lay4 = predictor_lay4)


pdf("./figures/figure10.pdf", width = 10, height = 7) 
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
                        re.form = NA, intervals = "predVar") 
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
                         re.form = NA, intervals = "predVar") 
pred_age2 <- data.frame(pred2 = pred_age_obj2, attr(pred_age_obj2, "intervals"), age2 = predictor_age2)


pdf("./figures/figure9.pdf", width = 10, height = 7)

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

glm_death <- fitme(Trichomonosiasis ~ location + age + sex,
                   family = binomial(link = "logit"),
                   data = death,
                   method = "PQL/L")

glm_death_GLM <- glm(Trichomonosiasis ~ location + age + sex,
                     family = binomial(link = "logit"),
                     data = death)

# effects
glm_death_no_loc <- fitme(Trichomonosiasis ~ age + sex,
                          family = binomial(link = "logit"),
                          data = death, method = "PQL/L")

glm_death_no_age <- fitme(Trichomonosiasis ~ location + sex,
                          family = binomial(link = "logit"),
                          data = death,
                          method = "PQL/L")

glm_death_no_sex <- fitme(Trichomonosiasis ~ location + age,
                          family = binomial(link = "logit"),
                          data = death,
                          method = "PQL/L")

# odds ratio
exp(glm_death$fixef["locationUrban"]) 
# locationUrban 
# 5.145808
CI_death <- confint(glm_death_GLM, parm = "locationUrban")
exp(CI_death)
# 2.5 %    97.5 % 
# 1.075153 43.764468 


## bootstraps
anova(glm_death_no_loc, glm_death, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glm_death_no_age, glm_death, boot.repl = boot.repl, nb_cores = nb_cores)
anova(glm_death_no_sex, glm_death, boot.repl = boot.repl, nb_cores = nb_cores)

####### window####

glm_death2 <- fitme(Trauma_window ~ location + sex,
                    family = binomial(link = "logit"),
                    data = death,
                    method = "PQL/L")

glm_death2_GLM <- glm(Trauma_window ~ location + sex,
                    family = binomial(link = "logit"),
                    data = death)

# effects
glm_death_no_loc2 <- fitme(Trauma_window ~  sex,
                           family = binomial(link = "logit"),
                           data = death,
                           method = "PQL/L")

glm_death_no_sex2 <- fitme(Trauma_window ~ location,
                           family = binomial(link = "logit"),
                           data = death,
                           method = "PQL/L")

# odds ratio
exp(glm_death2$fixef["locationUrban"]) 
# locationUrban 
# 3.537027
CI_death2loc <- confint(glm_death2_GLM, "locationUrban")
exp(CI_death2loc)
# 2.5 %   97.5 % 
# 1.388356 10.932307 

exp(glm_death2$fixef["sexm"]) 
# sexm 
# 2.311387
CI_death2sex <- confint(glm_death2_GLM, "sexm")
exp(CI_death2sex)
# 2.5 %    97.5 % 
# 1.194481 4.600079 

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


ggplot(data_plot2, aes(y = prop, x = Cause_of_death, fill = location, ymin = CI_lwr, ymax = CI_upr)) + 
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
ggsave(filename = "./figures/figure11.pdf", width = 10, height = 8)




###################### Figures supporting information #####################################################

table_per_nest <- read_excel("./source_data/goshawk_data_nest.xlsx")

goshawk_imper <- subset(table_per_nest, Species == "Goshawk")
goshawk_imper <- as.data.frame(goshawk_imper)

ggplot(goshawk_imper, aes(x = Habitat, y = imperviousness2500)) +
  geom_boxplot() +
  theme_bw() +
  ylab("\n Imperviousness (%) \n")
ggsave(filename = "./figures/figureS1.pdf", width = 8, height = 8)


Rupfungen_per_nest_aktuell <- read_excel("./source_data/rarefaction.xlsx")

urban <- subset(Rupfungen_per_nest_aktuell, Habitat == "Urban")
rural <- subset(Rupfungen_per_nest_aktuell, Habitat == "Rural")

urban1 <- urban[-c(1,2,3)]
rural1 <- rural[-c(1,2,3)]

urban1 <- as.matrix(urban1)

pdf(file = "./figures/figureS2.pdf", width = 7, height = 5)

par(las = 1, mar = c(5,4,1,1))
plot(specaccum(rural1, 
               method = "rarefaction"), 
     ci.type = "polygon", col = "grey", xlim = c(0, 60), ci.lty = 0, xlab = "Territory", ylab = "Number of species")

plot(specaccum(rural1, 
               method = "rarefaction"), lwd = 2,
     ci.type = "polygon", ci = 0, col = "black", add = TRUE)

plot(specaccum(urban1, 
               method = "rarefaction"), 
     add = TRUE, ci.type = "polygon", col = "black", ci.lty = 0)

plot(specaccum(urban1, 
               method = "rarefaction"), lwd = 2, 
     ci.type = "polygon", ci = 0, col = "white", add = TRUE)

legend("bottomright", fill = c("grey", "black"), bty = "n", legend = c("rural", "urban"))

dev.off()



### INFORMATION ON THE SESSION ###
sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] vegan_2.5-6        lattice_0.20-41    permute_0.9-5      doSNOW_1.0.19      snow_0.4-3        
# [6] iterators_1.0.13   foreach_1.5.1      mapdata_2.3.0      maps_3.3.0         ggplot2_3.3.2     
# [11] binom_1.1-1        purrr_0.3.4        car_3.0-10         carData_3.0-4      DHARMa_0.3.3.0    
# [16] readxl_1.3.1       dplyr_1.0.2        MASS_7.3-53        glmmTMB_1.0.2.9000 rstanarm_2.21.1   
# [21] Rcpp_1.0.5         lme4_1.1-25        Matrix_1.2-18      spaMM_3.5.2       
# 
# loaded via a namespace (and not attached):
#   [1] plyr_1.8.6            igraph_1.2.6          TMB_1.7.18            splines_4.0.3         crosstalk_1.1.0.1    
# [6] TH.data_1.0-10        rstantools_2.1.1      inline_0.3.16         digest_0.6.27         htmltools_0.5.0      
# [11] rsconnect_0.8.16      fansi_0.4.1           magrittr_1.5          cluster_2.1.0         openxlsx_4.2.3       
# [16] RcppParallel_5.0.2    matrixStats_0.57.0    xts_0.12.1            sandwich_3.0-0        prettyunits_1.1.1    
# [21] colorspace_2.0-0      haven_2.3.1           callr_3.5.1           crayon_1.3.4          jsonlite_1.7.1       
# [26] Rglpk_0.6-4           survival_3.2-7        zoo_1.8-8             glue_1.4.2            registry_0.5-1       
# [31] gtable_0.3.0          emmeans_1.5.2-1       V8_3.4.0              pkgbuild_1.1.0        rstan_2.21.2         
# [36] abind_1.4-5           scales_1.1.1          mvtnorm_1.1-1         miniUI_0.1.1.1        xtable_1.8-4         
# [41] foreign_0.8-80        proxy_0.4-24          stats4_4.0.3          StanHeaders_2.21.0-6  DT_0.16              
# [46] htmlwidgets_1.5.2     threejs_0.3.3         ellipsis_0.3.1        farver_2.0.3          pkgconfig_2.0.3      
# [51] loo_2.3.1             labeling_0.4.2        tidyselect_1.1.0      rlang_0.4.8           reshape2_1.4.4       
# [56] later_1.1.0.1         munsell_0.5.0         cellranger_1.1.0      tools_4.0.3           cli_2.1.0            
# [61] generics_0.1.0        ROI.plugin.glpk_1.0-0 ggridges_0.5.2        stringr_1.4.0         fastmap_1.0.1        
# [66] yaml_2.2.1            processx_3.4.4        zip_2.1.1             pbapply_1.4-3         nlme_3.1-150         
# [71] mime_0.9              slam_0.1-47           ROI_1.0-0             compiler_4.0.3        bayesplot_1.7.2      
# [76] shinythemes_1.1.2     rstudioapi_0.12       curl_4.3              tibble_3.0.4          statmod_1.4.35       
# [81] stringi_1.5.3         ps_1.4.0              forcats_0.5.0         nloptr_1.2.2.2        markdown_1.1         
# [86] shinyjs_2.0.0         vctrs_0.3.4           pillar_1.4.6          lifecycle_0.2.0       estimability_1.3     
# [91] data.table_1.13.2     httpuv_1.5.4          R6_2.5.0              promises_1.1.1        gridExtra_2.3        
# [96] rio_0.5.16            codetools_0.2-18      boot_1.3-25           colourpicker_1.1.0    gtools_3.8.2         
# [101] assertthat_0.2.1      withr_2.3.0           shinystan_2.5.0       multcomp_1.4-14       mgcv_1.8-33          
# [106] parallel_4.0.3        hms_0.5.3             grid_4.0.3            coda_0.19-4           minqa_1.2.4          
# [111] numDeriv_2016.8-1.1   shiny_1.5.0           base64enc_0.1-3       dygraphs_1.1.1.6-1