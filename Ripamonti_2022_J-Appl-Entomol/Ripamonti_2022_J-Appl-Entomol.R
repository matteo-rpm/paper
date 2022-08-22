#Citation:
#@article{,
#  author =       "Matteo Ripamonti and Luciana Galetto and Federico Maron and Cristina Marzachì and Domenico Bosco",
#  title =        "{Scaphoideus titanus fitness parameters on grapevine varieties with different susceptibility to Flavescence dorée phytoplasma}",
#  journal =      "Journal of Applied Entomology",
#  volume =       "",
#  number =       "",
#  pages =        "",
#  year =         "2022",
#  DOI =          ""
#}

##########Fitness tests for Scaphoideus titanus on different grapevine cultivars#############
library(survival)
library(survminer)
library(tidyverse)
library(cowplot)
library(readxl)
library(patchwork)
library(emmeans)
library(lme4)
library(coxme)
library(DHARMa)
library("performance")

#####DEVELOPMENTAL TIME, cumulative replicates (5)---------
#import dataset PTS
#anadevtRev <- tidyr::as_tibble(rio::import("https://raw.githubusercontent.com/matteo-rpm/projects/master/codes/fitness_db/anadevtR.compl.csv?token=GHSAT0AAAAAABQO3BLMVIAL2JQLTWBB2PJKYPSOEBA"))

#extract Sex and Cultivar columns
anadevtRev <- anadevtRev %>%
  #filter(new.ad == 1) %>% 
  mutate(Sex = case_when(sex == 0 ~ "male",
                         sex == 1 ~ "female",
                         TRUE ~ as.character(c(NA)))) %>% 
  mutate(Cultivar = case_when(treat == "BA84" ~ "Barbera",
                              treat == "BRA20" ~ "Brachetto",
                              treat == "MO190" ~ "Moscato",
                              TRUE ~ "other"))

#import dataset PFER
#anaferR20a <- tidyr::as_tibble(rio::import("https://raw.githubusercontent.com/matteo-rpm/projects/d83fff615c3337f80a7724680237ffeb41996057/codes/fitness_db/anaferR20a.csv?token=GHSAT0AAAAAABQO3BLMJBFJIKKMM2T46NZIYPSNTWQ"))

#anaferR20b <- tidyr::as_tibble(rio::import("https://raw.githubusercontent.com/matteo-rpm/projects/d83fff615c3337f80a7724680237ffeb41996057/codes/fitness_db/anaferR20b.csv?token=GHSAT0AAAAAABQO3BLMJST6OV7X6TP6E7GKYPSNXWQ"))

anaferR20a$manichettaxcv <- as.character(anaferR20a$manichettaxcv)
anaferR20a$moved_to_man_n <- as.character(anaferR20a$moved_to_man_n)
#anaferR20b$`IPSP CODE` <- as.character(anaferR20b$`IPSP CODE`)

#extract infos for columns
anaferR20a <- anaferR20a %>%
  mutate(Sex = case_when(sex == "m" ~ "male",
                         sex == "f" ~ "female",
                         TRUE ~ "other")) %>%
  mutate(Cultivar = case_when(treat == "BA84" ~ "Barbera",
                              treat == "BRA20" ~ "Brachetto",
                              treat == "MO190" ~ "Moscato",
                              TRUE ~ "other"))
#extract infos for columns
anaferR20b <- anaferR20b %>%
  mutate(Sex = case_when(sex == "m" ~ "male",
                         sex == "f" ~ "female",
                         TRUE ~ "other")) %>%
  mutate(Cultivar = case_when(treat == "BA84" ~ "Barbera",
                              treat == "BRA20" ~ "Brachetto",
                              treat == "MO190" ~ "Moscato",
                              TRUE ~ "other"))

#bind rep20A & 20B
anaferR20cum <- bind_rows(anaferR20a, anaferR20b) %>% mutate(dissday = dissected_days)

anaferR20cum$dissday[anaferR20cum$dissected_days == 13] <- 14
anaferR20cum$dissday[anaferR20cum$dissected_days == 24] <- 25

#import dataset PFER21
#anaferR21 <- tidyr::as_tibble(rio::import("https://raw.githubusercontent.com/matteo-rpm/projects/d83fff615c3337f80a7724680237ffeb41996057/codes/fitness_db/anaferR21.csv?token=GHSAT0AAAAAABQO3BLNGUQKL6YA3CP3T6FOYPSNY7A"))

anaferR21$manichettaxcv <- as.character(anaferR21$manichettaxcv)
anaferR21$moved_to_man_n <- as.character(anaferR21$moved_to_man_n)
anaferR21$`eggs_adv.mat` <- as.double(anaferR21$`eggs_adv.mat`)
anaferR21$`eggs_early.mat` <- as.double(anaferR21$`eggs_early.mat`)

anaferR21$dissected_days[anaferR21$dissected_days == 15] <- 14
anaferR21$dissected_days[anaferR21$dissected_days == 24] <- 25
anaferR21$dissected_days[anaferR21$dissected_days == 34] <- 35

#extract infos
anaferR21 <- anaferR21 %>%
  mutate(Sex = case_when(sex == "m" ~ "male",
                         sex == "f" ~ "female",
                         TRUE ~ "other")) %>%
  mutate(Cultivar = case_when(treat == "BA84" ~ "Barbera",
                              treat == "BRA20" ~ "Brachetto",
                              treat == "MO190" ~ "Moscato",
                              TRUE ~ "other"))


#unite PFER datasets
anaferALL <- bind_rows(anaferR20cum, anaferR21)

anaferALL$dissected_days[anaferALL$dissected_days == 13] <- 14
anaferALL$dissected_days[anaferALL$dissected_days == 24] <- 25

anaferALL <- anaferALL %>%
  mutate(rep = case_when(I_instar == "28/05/2020" ~ "Af",
                         I_instar == "18/08/2020" ~ "Bf",
                         I_instar == "20/04/2021" ~ "Cf",
                         TRUE ~ "other"))

anaferALL_sel <- anaferALL %>% 
  dplyr::select(ind, days2emerge, Cultivar, Sex, new.ad, rep) %>%
  rename(days = days2emerge)

#unite replicates in different datasets: PTS + PFER
wkp <- bind_rows(anadevtRev, anaferALL_sel) %>% 
  filter(new.ad != 0) %>%
  select(rep:Cultivar) #%>% 
#mutate(days_sub = days - min(days))

#median and IQR table
wkp_tab <- wkp %>% 
  dplyr::group_by(Cultivar, Sex) %>% 
  dplyr::summarise(n = n(), colquant = quantile(days, c(0.25, 0.75)), quant = c(0.25, 0.75),
                   mean = mean(days), median = median(days), iqr = IQR(days)) %>% 
  tidyr::spread(quant, colquant)

write.csv(wkp_tab, file = "Tab_med-IQR_dev-time.csv", row.names = T)

#dataset analysis with survival package
devt_wkp <- Surv(time = wkp$days, event = wkp$new.ad)

#model fit on Kaplan-Meyer curves, Cultivar+Sex treatment
wkp_fit <- survfit(devt_wkp ~ Cultivar+Sex, data = wkp) #+frailty(rep)
summary(wkp_fit)
surv_summary(wkp_fit)

#Log-Rank test
survdiff(Surv(time = wkp$days, event = wkp$new.ad) ~ Cultivar+Sex, data = wkp)

#Pairwise comparisons
pair_survd <- survminer::pairwise_survdiff(Surv(days, new.ad) ~ Cultivar+Sex, data = wkp, p.adjust.method = "BH")
symnum(pair_survd$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
       symbols = c("****", "***", "**", "*", "+", "NS"),
       abbr.colnames = T, na = "")

# Aalen's Additive Regression Model -> covariate effects not multiplicative (Cox), but additive
aadd <- survival::aareg(Surv(days, new.ad) ~ Cultivar*Sex, data = wkp)
summary(aadd)
summary(aadd)$table #covariate effects, p-val

library(ggfortify)
aaddf <- autoplot(aadd) #check covariates effects, plot

ggsave(filename = "FigS2_wkp_covar-eff_Aalen.tiff", plot = aaddf, dpi = 600, compression = "lzw", width = 10, height = 6.3)

#plot developmental time curves, Kaplan-Meier estimator (non-parametric model, no baseline hazard nor covariates assumptions)
wkp_devplot <- ggsurvplot(wkp_fit, data = wkp, risk.table = "abs_pct",  ylab = "Adult stage reach [%]",
                          xlab = "Time [days after hatch]", fun = "event", linetype = c("strata"),# size = 2, #palette = "grey",#facet.by = "sex_char",
                          xlim= c(19,70), legend.title = "Groups",
                          #log.rank.weights = "n", #pval = TRUE, pval.method = TRUE, pval.method.coord = c(20, 0.25),
                          legend.labs = c("Barbera, S.t. females", "Barbera, S.t. males",
                                          "Brachetto, S.t. females", "Brachetto, S.t. males",
                                          "Moscato, S.t. females", "Moscato, S.t. males"))

wkp_devplot_sel <- wkp_devplot$plot + 
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed", "solid", "dashed")) +
  scale_colour_manual(values = c("red", "red", "green", "green", "blue", "blue"))

wkp_plotdevt_patch <- wkp_devplot_sel / wkp_devplot$table + plot_layout(ncol = 1, heights = c(3, 1)) + 
  theme(axis.text.y.left = element_text(color = rev(c("red", "red", "green", "green", "blue", "blue"))))

ggsave(filename = "Figure1.tiff", plot = wkp_plotdevt_patch, dpi = 600, compression = "lzw", width = 10, height = 6.3)

# Quantification of treatments effects
#cox proportional hazard ratio with fixed effects only + stratification per replicate -> bad model, non-proportional hazards
wkp_forr <- coxph(devt_wkp ~ Cultivar+Sex + strata(rep), data = wkp) #+frailty(rep)
summary(wkp_forr)

cox.zph(wkp_forr) #Test the Proportional Hazards Assumptions -> non-proportional hazards! use another method
ggcoxzph(cox.zph(wkp_forr)) #Test the PH Assumptions, Schoenfeld residuals -> non-proportional hazards! use another method
ggcoxdiagnostics(wkp_forr) #Test PH assumptions, Martingale residuals, bad model

#cox proportional hazard ratio with mixed effects -> bad model, non-proportional hazards
wkp_forrcoxme <- coxme::coxme(devt_wkp ~ Cultivar + Sex + strata(rep) + (1|rep), data = wkp) #+strata() for stratification
summary(wkp_forrcoxme)

cox.zph(wkp_forrcoxme) #Test the Proportional Hazards Assumptions -> non-proportional hazards! use another method
ggcoxzph(cox.zph(wkp_forrcoxme)) #Test the PH Assumptions, Schoenfeld residuals -> non-proportional hazards! use another method

# Additive Cox model - chosen model, Peto's correction for ties + stratification
wkp_addcox <- mgcv::gam(days ~ Cultivar + Sex + strata(rep), data = wkp, family = "cox.ph", weights = new.ad)
ols <- summary(wkp_addcox)
ols
write.csv(ols$p.table, file = "S3_wkp_addcox.csv", row.names = T)

# EMMs pairwise comparisons between groups
mult.wkp_addcox = emmeans::emmeans(wkp_addcox, pairwise ~ Cultivar + Sex, adj="holm")
wkp_GAMmultcomp <- multcomp::cld(mult.wkp_addcox$emmeans, adj="holm", Letters = letters, decreasing = T)
write.csv(wkp_GAMmultcomp, file = "Tab_wkp-GAMmultcomp.csv", row.names = T)

FigS2 <- plot(multcomp::cld(mult.wkp_addcox$emmeans, adj="holm", Letters = letters, decreasing = T))
ggsave(filename = "FigS2_add-cox-mod_effects.tiff", plot = FigS2, dpi = 600, compression = "lzw", width = 10, height = 6.3)


# Nymphs mortality comparison among cultivars and replicates----------

# prepare datasets
fer20a_count <- anaferR20a %>% count(Cultivar, Sex) %>% rename(Af = n)
fer20b_count <- anaferR20b %>% count(Cultivar, Sex) %>% rename(Bf = n)
fer21_count <- anaferR21 %>% count(Cultivar, Sex) %>% rename(Cf = n)
ferALL_count <- left_join(left_join(fer20a_count, fer20b_count), fer21_count)

ferALL_count[is.na(ferALL_count)] <- 0
ferALL_count

# import dead nymphs dataset
#anafer_deadDB <- tidyr::as_tibble(rio::import("https://raw.githubusercontent.com/matteo-rpm/projects/d83fff615c3337f80a7724680237ffeb41996057/codes/fitness_db/anaferR21.csv?token=GHSAT0AAAAAABQO3BLNGUQKL6YA3CP3T6FOYPSNY7A"))

anadevtRev <- anadevtRev %>% mutate(status = case_when(new.ad == 1 ~ "emerged",
                                                       new.ad == 0 ~ "dead",
                                                       TRUE ~ "other"))

anafer_deadDB <- anafer_deadDB %>% mutate(status = case_when(new.ad == 1 ~ "emerged",
                                                             new.ad == 0 ~ "dead",
                                                             TRUE ~ "other"))

# join dead and emerged nymphs databases
fernymphs <- bind_rows(bind_rows(anaferALL_sel %>% mutate(status = case_when(new.ad == 1 ~ "emerged",
                                                                             new.ad == 0 ~ "dead",
                                                                             TRUE ~ "other")), anafer_deadDB), anadevtRev)

# summary of dead vs emerged nymphs per Cultivar
fernymphs_tab <- fernymphs %>% count(Cultivar, Sex, status)
write.csv(fernymphs_tab, file = "Tab_nymph-surv.csv", row.names = T)

#GLMM Bernoulli binomial model on emergence of adults (complementary to death of nymph)

#GLMM, family = "binomial", link = "logit"
glmm_logit <- lme4::glmer(new.ad ~ Cultivar + (1|rep), data = fernymphs, family = binomial(link = "logit"), nAGQ = 11)
summary(glmm_logit)

#GLMM, family = "binomial", link = "probit"
glmm_probit <- lme4::glmer(new.ad ~ Cultivar + (1|rep), data = fernymphs, family = binomial(link = "probit"), nAGQ = 11)
summary(glmm_probit)

#GLMM, family = "binomial", link = "cauchit"
glmm_cauchit <- lme4::glmer(new.ad ~ Cultivar + (1|rep), data = fernymphs, family = binomial(link = "cauchit"), nAGQ = 11)
summary(glmm_cauchit)

# compare models
compare_performance(glmm_logit, glmm_probit, glmm_cauchit, metrics = "all", rank = T)
plot(compare_performance(glmm_logit, glmm_probit, glmm_cauchit, metrics = "all", rank = T))
plot(compare_performance(glmm_logit, glmm_probit, glmm_cauchit)) # Binomial GLMM with Cauchit link fit better
check_distribution(glmm_cauchit) # Bernoulli distribution confirmed
check_model(glmm_cauchit) # check model -> heteroscedasticity present -> differences between cultivars are enormous, no robust SE calculation

dhitri <- simulateResiduals(fittedModel = glmm_cauchit)
plot(dhitri, quantreg = T)
testResiduals(dhitri) #deviation from the ideal curve is small, probably the large amount of observations makes it very sensitive

# EMMs pairwise comparisons between groups
ref.glmm_cauchit <- emmeans::ref_grid(glmm_cauchit)
emm.glmm_cauchit <- emmeans::emmeans(ref.glmm_cauchit, specs = "Cultivar", adj = "holm")
cld.glmm_cauchit <- as_tibble(multcomp::cld(emm.glmm_cauchit, Letters = "abcdef"))
write.csv(cld.glmm_cauchit, file = "Tab_nymph-surv_multcomp.csv", row.names = T)

FigS3 <- plot(multcomp::cld(emm.glmm_cauchit, Letters = "abcdef"))
ggsave(filename = "FigS2_glmm-cauchit_nymphmort.tiff", plot = FigS3, dpi = 600, compression = "lzw", width = 10, height = 6.3)

#report GLMM table with jtools and huxtable
library("huxtable")
library("jtools")
library("broom.mixed")

jtools::export_summs(glmm_logit, glmm_probit, glmm_cauchit, 
                     model.names = c("glmm_logit", "glmm_probit", "glmm_cauchit"),
                     error_format = "({std.error})",
                     to.file = "xlsx", file.name = "FileS2_Tab_GLMMs_nymph-mort.xlsx",
                     exp = F, robust = F)

#####SURVIVAL CURVES, cumulative replicates (5)---------

#dataset import
#anasurvR <- tidyr::as_tibble(rio::import("https://raw.githubusercontent.com/matteo-rpm/projects/d83fff615c3337f80a7724680237ffeb41996057/codes/fitness_db/anasurvR.csv?token=GHSAT0AAAAAABQO3BLNXFCZWZVLEZYFRYHKYPSN2HQ"))

#modify BRA20m outlier -> not dead
anasurvR$day[anasurvR$ind == 33] <- 110 #max value of other treats
anasurvR$death[anasurvR$ind == 33] <- 0 #event "death" not happened
anasurvR %>% filter(ind == 33)

#extract infos
anasurvR <- anasurvR %>%
  mutate(Sex = case_when(sex == 0 ~ "male",
                         sex == 1 ~ "female",
                         TRUE ~ "other")) %>% 
  mutate(Cultivar = case_when(treat == "BA84" ~ "Barbera",
                              treat == "BRA20" ~ "Brachetto",
                              treat == "MO190" ~ "Moscato",
                              TRUE ~ "other"))

#tidy PFER dataset for survival
anaferALLun <- anaferALL %>% unite(event_day, dead_dayspostneosf|dissected_days|days_censored, sep = "", remove = F, na.rm=T)
anaferALLun$event_day <- as.numeric(anaferALLun$event_day)

anaferALLun$event_day[anaferALLun$event_day == 2525] <- 25
anaferALLun$event_day[anaferALLun$event_day == 1414] <- 14
anaferALLun$event_day[anaferALLun$event_day == 1818] <- 18
anaferALLun$event_day[anaferALLun$event_day == 3535] <- 35
anaferALLun$event_day[anaferALLun$event_day == 3435] <- 35

anaferALLun <- anaferALLun %>%
  filter(!is.na(event_day)) %>% # remove "lost" or moved (males to a different cultivar) individuals
  mutate(death = ifelse(is.na(dead_dayspostneosf), 0, 1))
#filter(sex == "f") %>% #`` & ind != 5)

anaferALLun_surv <- anaferALLun %>% 
  dplyr::select(ind, event_day, Cultivar, Sex, death, rep, dead_dayspostneosf)

#unite PLO + PFER datasets for survival
skp <- bind_rows(anasurvR %>%
                   dplyr::select(!c(sex, treat, treat_sex)) %>%
                   rename(sex = sex_char, event_day = day) %>% 
                   mutate(dead_dayspostneosf = event_day),
                 anaferALLun_surv)

#median and IQR table
skp_tab <- skp %>% 
  filter(!is.na(dead_dayspostneosf)) %>% 
  group_by(Cultivar, Sex) %>% 
  summarise(n = n(), colquant = quantile(dead_dayspostneosf, c(0.25, 0.75)), quant = c(0.25, 0.75),
            mean = mean(dead_dayspostneosf), median = median(dead_dayspostneosf), iqr = IQR(dead_dayspostneosf)) %>% 
  tidyr::spread(quant, colquant)

write.csv(skp_tab, file = "Tab_med-IQR_surv.csv", row.names = T)

#dataset analysis with survival package
skp_surv <- Surv(time = skp$event_day, event = skp$death)

#model fit on Kaplan-Meyer curves, treatment = Cultivar+Sex
skp_fit <- survfit(skp_surv ~ Cultivar+Sex, data = skp)
summary(skp_fit)

#Log-Rank test
survdiff(Surv(time = skp$event_day, event = skp$death) ~ Cultivar+Sex, data = skp)

#Pairwise comparisons
pair_survlong <- survminer::pairwise_survdiff(Surv(event_day, death) ~ Cultivar+Sex, data = skp, p.adjust.method = "BH")
symnum(pair_survlong$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
       symbols = c("****", "***", "**", "*", "+", "NS"),
       abbr.colnames = T, na = "")

# Aalen's Additive Regression Model -> covariate effects not multiplicative (Cox), but additive
aass <- survival::aareg(Surv(event_day, death) ~ Cultivar*Sex, data = skp)
summary(aass)
summary(aass)$table #covariate effects, p-val

library(ggfortify)
aassf <- autoplot(aass) #check covariates effects, plot

ggsave(filename = "FigS4_skp_covar-eff_Aalen.tiff", plot = aassf, dpi = 600, compression = "lzw", width = 10, height = 6.3)

#plot survival time curves, Kaplan-Meier estimator (non-parametric model, no baseline hazard nor covariates assumptions)
skp_survplot <- ggsurvplot(skp_fit, data = skp, risk.table = "nrisk_cumcensor", #surv.median.line = "hv",
                           censor.shape = c("X"), linetype = c("strata"), legend.title = "Groups",
                           ylab = "Survival probability [%]", xlab = "Time [days]",
                           legend.labs = c("Barbera, S.t. females", "Barbera, S.t. males",
                                           "Brachetto, S.t. females", "Brachetto, S.t. males",
                                           "Moscato, S.t. females", "Moscato, S.t. males"),
                           xlim= c(0,100))

skp_survplot_sel <- skp_survplot$plot + 
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed", "solid", "dashed")) +
  scale_colour_manual(values = c("red", "red", "green", "green", "blue", "blue"))

skp_survplot_patch <- skp_survplot_sel / skp_survplot$table + plot_layout(ncol = 1, heights = c(3, 1)) +
  theme(axis.text.y.left = element_text(color = rev(c("red", "red", "green", "green", "blue", "blue"))))

ggsave(filename = "Figure2.tiff", plot = skp_survplot_patch, dpi = 600, compression = "lzw", width = 10, height = 6.3)

# Quantification of treatments effects
#cox proportional hazard ratio with fixed effects only + stratification per replicate -> bad model, non-proportional hazards
skp_forr <- coxph(skp_surv ~ Cultivar+Sex +strata(rep), data = skp)
summary(skp_forr)

cox.zph(skp_forr) #Test the Proportional Hazards Assumptions -> non-proportional hazards! use another method
ggcoxzph(cox.zph(skp_forr)) #Test the PH Assumptions, Schoenfeld residuals -> non-proportional hazards! use another method
ggcoxdiagnostics(skp_forr) #Test PH assumptions, Martingale residuals, bad model

#cox proportional hazard ratio with mixed effects -> bad model, non-proportional hazards
skp_forrcoxme <- coxme::coxme(skp_surv ~ Cultivar + Sex + (1|rep), data = skp)
summary(skp_forrcoxme)

cox.zph(skp_forrcoxme) #Test the Proportional Hazards Assumptions -> non-proportional hazards! use another method
ggcoxzph(cox.zph(skp_forrcoxme)) #Test the PH Assumptions, Schoenfeld residuals -> non-proportional hazards! use another method

# Additive Cox model - chosen model, Peto's correction for ties + stratification
skp_addcox <- mgcv::gam(event_day ~ Cultivar + Sex + strata(rep), data = skp, family = "cox.ph", weights = death)
als <- summary(skp_addcox)
als
write.csv(als$p.table, file = "S4_skp_addcox.csv", row.names = T)

# EMMs for the selected model
mult.skp_addcox = emmeans::emmeans(skp_addcox, pairwise ~ Cultivar + Sex, adj="holm")
skp_GAMmultcomp <- multcomp::cld(mult.skp_addcox$emmeans, adj="holm", Letters = letters, decreasing = T)
write.csv(skp_GAMmultcomp, file = "Tab_skp-GAMmultcomp.csv", row.names = T)

FigS4 <- plot(multcomp::cld(mult.skp_addcox$emmeans, adj="holm", Letters = letters, decreasing = T))
ggsave(filename = "FigS4_add-cox-mod_surv-effects.tiff", plot = FigS5, dpi = 600, compression = "lzw", width = 10, height = 6.3)

#Joint developmental time and survival Kaplan-Meier curves---
np_dat <- wkp_devplot$data.survplot
np_dat$time <- np_dat$time - max(np_dat$time)
np_dat$surv <- 1 - np_dat$surv 

df <- rbind(np_dat, skp_survplot$data.survplot)

df <- as_tibble(df) %>% #adding points to start from the zero baseline
  add_row(time = -44, n.risk = 91, n.event = 1, n.censor = 0, surv = 0, std.err = Inf,
          upper = NA, lower = NA, strata = as.factor("Cultivar=Barbera, Sex=female"), Cultivar = as.factor("Barbera"), Sex = as.factor("female")) %>% 
  add_row(time = -47, n.risk = 152, n.event = 1, n.censor = 0, surv = 0, std.err = Inf,
          upper = NA, lower = NA, strata = as.factor("Cultivar=Barbera, Sex=male"), Cultivar = as.factor("Barbera"), Sex = as.factor("male")) %>% 
  add_row(time = -42, n.risk = 67, n.event = 1, n.censor = 0, surv = 0, std.err = Inf,
          upper = NA, lower = NA, strata = as.factor("Cultivar=Brachetto, Sex=female"), Cultivar = as.factor("Brachetto"), Sex = as.factor("female")) %>% 
  add_row(time = -43, n.risk = 100, n.event = 1, n.censor = 0, surv = 0, std.err = Inf,
          upper = NA, lower = NA, strata = as.factor("Cultivar=Brachetto, Sex=male"), Cultivar = as.factor("Brachetto"), Sex = as.factor("male")) %>% 
  add_row(time = -40, n.risk = 41, n.event = 1, n.censor = 0, surv = 0, std.err = Inf,
          upper = NA, lower = NA, strata = as.factor("Cultivar=Moscato, Sex=female"), Cultivar = as.factor("Moscato"), Sex = as.factor("female")) %>% 
  add_row(time = -42, n.risk = 59, n.event = 1, n.censor = 0, surv = 0, std.err = Inf,
          upper = NA, lower = NA, strata = as.factor("Cultivar=Moscato, Sex=male"), Cultivar = as.factor("Moscato"), Sex = as.factor("male"))

devtsurv_joint <- ggplot(df, aes(time, surv, color = Cultivar, linetype = Sex)) + 
  geom_step(size = 0.9) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(y = 'Probability [%]',
       x = 'Time from emergence [days]',
       color = "Cultivar",
       linetype = "Sex")+
  #geom_point(data = df %>% filter(n.censor != 0), aes(time, surv, color = Cultivar))+
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_colour_manual(values = c("red", "green", "blue"))+
  #scale_color_brewer(palette = "Set1")+
  theme_minimal()+
  theme(panel.grid = element_blank(), panel.border = element_rect(fill = NA))

ggsave(filename = "Fig3_devtsurv_joint.tiff", plot = devtsurv_joint, dpi = 600, compression = "lzw", width = 8, height = 4)

#calculate area under developmental time and survival curves
library("pracma")

dfBAf <- df %>% filter(Cultivar == "Barbera") %>% filter(Sex == "female")
AUC_BAf <- trapz(dfBAf$time, dfBAf$surv)

dfBAm <- df %>% filter(Cultivar == "Barbera") %>% filter(Sex == "male")
AUC_BAm <- trapz(dfBAm$time, dfBAm$surv)

dfBRAf <- df %>% filter(Cultivar == "Brachetto") %>% filter(Sex == "female")
AUC_BRAf <- trapz(dfBRAf$time, dfBRAf$surv)

dfBRAm <- df %>% filter(Cultivar == "Brachetto") %>% filter(Sex == "male")
AUC_BRAm <- trapz(dfBRAm$time, dfBRAm$surv)

dfMOf <- df %>% filter(Cultivar == "Moscato") %>% filter(Sex == "female")
AUC_MOf <- trapz(dfMOf$time, dfMOf$surv)

dfMOm <- df %>% filter(Cultivar == "Moscato") %>% filter(Sex == "male")
AUC_MOm <- trapz(dfMOm$time, dfMOm$surv)

areas <- matrix(c("Barbera", "Barbera", "Brachetto", "Brachetto", "Moscato", "Moscato",
                  "female", "male", "female", "male", "female", "male",
                  AUC_BAf, AUC_BAm, AUC_BRAf, AUC_BRAm, AUC_MOf, AUC_MOm), ncol = 3, byrow = F)

colnames(areas) <- c("Cultivar", "Sex", "Area")

write.csv(areas, file = "Tab_areas-curves.csv", row.names = T)

#####PROLIFICACY, cumulative replicates (3)#######

# FERTILITY TEST, eggs + gene exp combined---
#genexpvit <- tidyr::as_tibble(rio::import("https://raw.githubusercontent.com/matteo-rpm/projects/master/codes/fitness_db/PFER2020-2021_genestudycomp_report_Maestro2.2_write.csv?token=GHSAT0AAAAAABQO3BLM2ABC46KVEJ67GY6YYPSUODQ"))

# tidy dataset
genexpvit <- genexpvit %>% filter(Target == "vitellogenin") %>%
  tidyr::separate(Biological_Group_Sample, sep = "dpe_", into = c("Biological_Group", "IPSP_CODE"))

genexpvit$IPSP_CODE <- as.double(genexpvit$IPSP_CODE)

#choose only processed (dissected/extracted) individuals
dob <- anaferALL %>% filter(!is.na(IPSP_CODE))

# rearrange sampling points within +- 1 day
dob$dissected_days[dob$dissected_days == 13] <- 14
dob$dissected_days[dob$dissected_days == 24] <- 25

# unite gene expression db with egg counts db
tmp <- left_join(dob, genexpvit, by = "IPSP_CODE")

tmp <- tmp %>% dplyr::select(IPSP_CODE, sex, treat, dissected_days, eggs_mature, `eggs_adv.mat`, `eggs_early.mat`, Target, `Expression`, Cultivar, Sex, rep) %>% 
  filter(dissected_days != 18) %>% 
  filter(dissected_days != 28) %>% 
  filter(dissected_days != 39)

tmp$`Expression` <- as.numeric(unlist(tmp$`Expression`))
#tmp$`Normalized Expression` <- as.numeric(unlist(tmp$`Normalized Expression`))
#tmp$`Relative Normalized Expression` <- as.numeric(unlist(tmp$`Relative Normalized Expression`))

#number of dissected females, ferALL
Tab_ferALL <- tmp %>% filter(Sex == "female") %>% count(Cultivar, dissected_days) %>%
  pivot_wider(names_from = "dissected_days", values_from = "n")

write.csv(Tab_ferALL, file = "Tab_ferALL.csv", row.names = T)
##check correlations among variables
tmp1 <- tmp %>% dplyr::select(dissected_days, eggs_mature, `eggs_adv.mat`, `eggs_early.mat`,
                              `Expression`) %>%
  filter(!is.na(eggs_mature))

#check correlations -> no correlations
usdm::vifcor(as.data.frame(tmp1), th=5)

#db for ggplot
temp <- tmp %>% 
  filter(!is.na(dissected_days)) %>% 
  filter(dissected_days != 18) %>% 
  filter(Sex != "male") %>% 
  filter(!is.na(Expression)) %>% 
  mutate(diss_time = case_when(dissected_days == "14" ~ "t1_14",
                               dissected_days == "25" ~ "t2_25",
                               dissected_days == "35" ~ "t3_35",
                               TRUE ~ "other"))

is.na(temp$eggs_mature)

# Bi-graph with number of eggs and vitellogenin expression associated---

doubg <- ggplot(temp, aes(x=eggs_mature, y = `Expression`)) +
  geom_point(size = 3, aes(shape = as.factor(dissected_days), color = Cultivar))+
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3"))+
  scale_shape_discrete(name = "Dissection day ")+
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"))+
  #ggpmisc::stat_poly_eq(formula = y ~ s(x, bs = "cs"), 
  #aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
  #parse = TRUE)+
  #RUncommon::stat_smooth_func(geom="text",method="gam",hjust=0,parse=TRUE) +
  theme_bw()+
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0))+
  labs(x = "Number of mature eggs",  y = "Vitellogenin expression")#+
#stat_ellipse(type = "euclid", level = 1, aes(shape = as.factor(dissected_days)))
#ggpmisc::stat_poly_eq(formula =  y ~ s(x, bs = "cs"), 
#            aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#           parse = TRUE) 
#geom_smooth(method = "glm")

ggsave(filename = "FigS7_bi-graph.tiff", plot = doubg, dpi = 600, compression = "lzw", width = 7, height = 6)

# TEST DIFFERENCES in number of mature eggs 
#GLMMs
library(glmmTMB)

#temp %>% filter(Sex != "male") %>% count(Cultivar)

# zero-inflation Poisson GLMM
fit_zipoisson1 <- glmmTMB(eggs_mature ~ Cultivar + as.factor(dissected_days) + (1|rep), data = temp, ziformula = ~1, family = poisson)
fit_zipoisson1
summary(fit_zipoisson1)

# zero-inflation negative binomial GLMM
fit_zinegbin1 <- glmmTMB(eggs_mature ~ Cultivar + as.factor(dissected_days) + (1|rep), 
                         family=list(family="nbinom1",link="log"), data = temp, ziformula = ~1)
fit_zinegbin1
summary(fit_zinegbin1)

# zero-inflation negative binomial GLMM
fit_zinegbin2 <- glmmTMB(eggs_mature ~ Cultivar + as.factor(dissected_days) + (1|rep), 
                         family=list(family="nbinom2", link="log"), data = temp, ziformula = ~1)
fit_zinegbin2
summary(fit_zinegbin2)

# zero-inflation Hurdle GLMM
fit_zihurdle <- glmmTMB(eggs_mature ~ Cultivar + as.factor(dissected_days) + (1|rep), data = temp, ziformula = ~., family = truncated_nbinom1)
fit_zihurdle
summary(fit_zihurdle)

# LMM with nlme
library(nlme)
nlme1 <- nlme::lme(eggs_mature ~ Cultivar+as.factor(dissected_days), random = ~ 1|rep, data = temp)
summary(nlme1)

# LMM with lme4
lmer1 <- lme4::lmer(eggs_mature ~ Cultivar+as.factor(dissected_days) + (1|rep), data = temp)

# compare performances of all the models
compare_performance(fit_zipoisson1, fit_zinegbin1, fit_zinegbin2, fit_zihurdle, nlme1, lmer1, metrics = "all", rank = T)
compare_performance(fit_zihurdle, nlme1, metrics = "all")
plot(compare_performance(fit_zihurdle, nlme1, lmer1))
check_model(fit_zihurdle)
check_model(nlme1) # select nlme1

# EMMs pairwise comparisons between groups
worthy <- emmeans::emmeans(nlme1, ~Cultivar+as.factor(dissected_days))
plot(worthy, horizontal = F)

eggs.mat_multcomp <- as.data.frame(multcomp::cld(worthy, reversed = F, Letters = letters))
emmeans::pwpp(worthy)

newdata <- eggs.mat_multcomp[order("Cultivar", "dissected_days"),]

write.csv(eggs.mat_multcomp, file = "Tab_eggs.mat_multcomp.csv", row.names = T)

#report (G)LMM table with jtools and huxtable
library("huxtable")
library("jtools")

jtools::export_summs(fit_zipoisson1, fit_zinegbin1, fit_zinegbin2, fit_zihurdle, nlme1, lmer1,
                     model.names = c("fit_zipoisson1", "fit_zinegbin1", "fit_zinegbin2", "fit_zihurdle", "nlme1", "lmer1"),
                     error_format = "({std.error})",
                     to.file = "xlsx", file.name = "Tab_GLMMs_eggs-mat.xlsx",
                     exp = F, robust = F)

#facets ggplot for number of mature eggs (eggs_mature) and vitellogenin expression---
#gather columns eggs_mature and Expression to reproduce them side-by-side in ggplot2::geom_boxplot()
tempo <- temp %>% 
  rename(`Number of mature eggs` = eggs_mature) %>% 
  gather(`Number of mature eggs`, `Expression`, key = "variable", value = "number")

#double boxplot with position dodge
#tempoline <- tempo %>% filter(variable == "Expression") #substitute eggs_mature with `Normalized Expression` or viceversa
eggenexp <- ggplot(tempo, aes(x=Cultivar, y = number, fill = variable))+
  geom_boxplot(aes(y=number, fill = variable), position = position_dodge(width = 1), outlier.shape = NA)+
  geom_jitter(aes(y=number, group = variable), width = 0.1)+#, position = position_dodge(width = 1))+
  #geom_smooth(data = tempoline, method = "loess", formula = number ~ as.factor(Cultivar), aes(x=Cultivar, y=number, group = variable))+ #add cultivar line to different dpe
  #geom_line(data = tempoline, aes(x = as.factor(Cultivar), #cultivar
  #                           y = number, #number
  #                          group = variable))+ #variable
  #geom_text(data = as.data.frame(eggs.mat_multcomp), aes(y = 40, label = trimws(.group)))+
  #geom_text(label = c("bc", "de", "e", "ab", "cde", "cde", "a", "bc", "bcd", "", "", "", "", "", "", "", "", ""),
  #          aes(y = 32, x = as.factor(dissected_days)), size = 6)+
  lemon::facet_rep_grid(factor(variable, levels = c("Number of mature eggs", "Expression")) ~ as.factor(dissected_days), labeller = label_value, scales = "free_y")+
  #scale_y_continuous(sec.axis = sec_axis(~./scl, name = "Vitellogenin normalized expression"))+
  #xlab("Dissection day [days after emergence]")+
  xlab("Cultivar")+
  ylab(NULL)+
  #cowplot::theme_cowplot()+
  #cowplot::background_grid()+
  theme_bw()+
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0))+
  scale_fill_manual(name = "Legend", labels = c("Vitellogenin expression", "Number of mature eggs"),
                    values = c("#717171", "grey"))

#add compact letter display for comparisons in number of mature eggs
egg_comp <- eggs.mat_multcomp %>% arrange(dissected_days, Cultivar)

tempo2 <- distinct(tempo, Cultivar, dissected_days, variable) %>%
  arrange(variable, dissected_days, Cultivar)

tempo2$yloc = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 31, 31, 31, 31, 31, 31, 31, 31, 31)

tempo2$label = c( "", "", "", "", "", "", "", "", "", trimws(egg_comp$.group))

eggenexp_cld <- eggenexp + 
  geom_text(data = tempo2, aes(y = yloc, label = label), size = 3.5)

#export graph
ggsave(filename = "Fig3_eggenexp_cld_vLG.tiff", plot = eggenexp_cld, dpi = 600, compression = "lzw", width = 7, height = 6)

##### vitellogenin gene expression ####
# tl;dr -> no differences!
tempexp <- temp %>%
  filter(!is.na(Expression))

tempexp$Expression

# Gamma GLMM
fit_gamma <- glmmTMB(Expression ~ Cultivar + as.factor(dissected_days) + (1|rep), data = tempexp, family = Gamma)
summary(fit_gamma)

# Gaussian GLMM
fit_gaussian <- glmmTMB(Expression ~ Cultivar + as.factor(dissected_days) + (1|rep), data = tempexp, family = gaussian)
summary(fit_gaussian)

# LMM
fit_LMM <- nlme::lme(Expression ~ Cultivar+as.factor(dissected_days), random = ~ 1|rep, data = tempexp)
summary(fit_LMM)

# Hurdle-Gamma GLMM
fit_zigamma <- glmmTMB(Expression ~ Cultivar + as.factor(dissected_days) + (1|rep), data = tempexp, family = ziGamma)
summary(fit_zigamma)

# negative binomial GLMM
fit_glmm_negbin <- lme4::glmer.nb(Expression ~ Cultivar + as.factor(dissected_days) + (1|rep), data = tempexp)
summary(fit_glmm_negbin)

# compare models
compare_performance(fit_gamma, fit_gaussian, fit_zigamma, fit_glmm_negbin, metrics = "all")
compare_performance(fit_gamma, fit_gaussian, metrics = "all", rank = T)
plot(compare_performance(fit_gamma, fit_gaussian, fit_LMM))
check_model(fit_gamma) #choose Gaussian GLMM

# EMMs pairwise comparisons between groups
king <- emmeans::emmeans(fit_gamma, ~Cultivar+as.factor(dissected_days))
plot(king, horizontal = F)

multcomp::cld(king, reversed = F, Letters = letters)
emmeans::pwpp(king)

#report GLMM table with jtools and huxtable
jtools::export_summs(fit_gaussian, fit_LMM, fit_gamma, fit_zigamma, fit_glmm_negbin, 
                     model.names = c("fit_gaussian", "fit_LMM", "fit_gamma", "fit_zigamma", "fit_glmm_negbin"),
                     error_format = "({std.error})",
                     to.file = "xlsx", file.name = "Tab_GLMMs_vitell.xlsx",
                     exp = F, robust = F)

