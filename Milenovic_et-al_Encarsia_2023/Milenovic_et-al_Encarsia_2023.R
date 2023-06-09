library(tidyverse)

#import dataset
paras_surv <- tibble::as_tibble(readxl::read_excel("C:/Users/ripamonti/Documents/Work/Parasitoids/Encarsia/E_formosa_longevity_CC.xlsx"))
#paras_surv <- tibble::as_tibble(readxl::read_excel("E_formosa_longevity_CC.xlsx"))

#tidying
paras_surv <- paras_surv %>% 
  select(!c(notes)) %>% 
  arrange(vial)

paras_surv <- paras_surv %>% 
  group_by(climate, diet, vial) %>%
  mutate(dead_that_day = dead - lag(dead)) %>% 
  mutate(dead_that_day = coalesce(dead_that_day, dead)) %>% 
  mutate(origin = ifelse(grepl('PAR', vial, ignore.case = T), 'LIST_nymphs','Koppert'))

paras_surv %>% group_by(climate, diet) %>% summarise(sum(dead_that_day))

paras_surv %>% filter(dead_that_day < 0)

#uncount dead sum -> one row per individual
paras_surv_4ana <- paras_surv %>%
  uncount(dead_that_day, .remove = T) %>% 
  select(!dead) %>% 
  mutate(death_hapn = 1)

#survival analysis
library(survival)
library(survminer)

paras_surv_4ana <- paras_surv_4ana %>%
  mutate(diet = factor(diet, levels = c("neg_control", "honey")))
paras_surv_4ana <- paras_surv_4ana %>%
  mutate(origin = factor(origin, levels = c("Koppert","LIST_nymphs")))
paras_surv_4ana <- paras_surv_4ana %>%
  mutate(climate = factor(climate, levels = c("future", "present")))


survf <- Surv(time = paras_surv_4ana$day, event = paras_surv_4ana$death_hapn)

#model fit on Kaplan-Meyer curves, climate+diet treatment
survf_fit <- survfit(survf ~ origin+diet+climate, data = paras_surv_4ana) #+frailty(rep) #if needed for stratification
summary(survf_fit)
surv_summary(survf_fit)

survf_fit$strata <- survf_fit$strata  %>% 
  mutate(climate = factor(climate, levels = c("present", "future")))

#preliminar Log-Rank test
survdiff(Surv(time = paras_surv_4ana$day, event = paras_surv_4ana$death_hapn) ~ climate+diet+origin, data = paras_surv_4ana)

# Aalen's Additive Regression Model -> covariate effects not multiplicative (Cox), but additive
aadd <- survival::aareg(Surv(day, death_hapn) ~ origin+diet+climate, data = paras_surv_4ana)
summary(aadd)
summary(aadd)$table #covariate effects, p-val

library(ggfortify)
aaddf <- autoplot(aadd) #check covariates effects, plot


#plot developmental time curves, Kaplan-Meier estimator (non-parametric model, no baseline hazard nor covariates assumptions)
paras_survplot <- ggsurvplot(survf_fit, data = paras_surv_4ana, risk.table = "abs_pct",  ylab = "Survival rate [%]",
                             xlab = "Time [days]", linetype = c("strata"), palette = "Set4",# size = 2,#facet.by = "sex_char",
                             #xlim= c(-1,36),
                             conf.int = T,
                             #conf.int.fill = c("strata"),
                             legend.title = "", surv.median.line = "hv",
                             log.rank.weights = "n", #pval = TRUE, pval.method = TRUE, pval.method.coord = c(20, 0.25),
                              legend.labs = c("Future, no-diet control, F0-commercial nymphs", "Present, no-diet control, F0-commercial nymphs", 
                                              "Future, honey, F0-commercial nymphs", "Present, honey,  F0-commercial nymphs", 
                                              "Future, honey, F1-in-house reared", "Present, honey, F1-in-house reared"), 
                             break.time.by = 5
)

colors = c("black", "#858585", "chocolate3", "darkolivegreen3", "chocolate3", "darkolivegreen3")

paras_survplot_sel <- paras_survplot$plot + 
  #order as in legend.labs above
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "dashed", "dashed")) +
  scale_colour_manual(values = colors) +
  scale_fill_manual(values = colors)

library(patchwork)
paras_survplot_patch <- paras_survplot_sel / paras_survplot$table + plot_layout(ncol = 1, heights = c(3, 1)) + 
  theme(axis.text.y.left = element_text(color = rev(colors)))
paras_survplot_patch
#ggsave(filename = "Figure2.tiff", plot = paras_survplot_patch, dpi = 300, compression = "lzw", units = "px", width = 5700, height = 2700)

# Quantification of treatments effects
#cox proportional hazard ratio with fixed effects only + stratification per replicate -> bad model, non-proportional hazards
paras_forr <- coxph(survf ~ climate+diet+origin + frailty(vial), data = paras_surv_4ana) #+frailty(rep)
summary(paras_forr)

cox.zph(paras_forr) #Test the Proportional Hazards Assumptions -> non-proportional hazards! use another method
ggcoxzph(cox.zph(paras_forr)) #Test the PH Assumptions, Schoenfeld residuals -> non-proportional hazards! use another method
ggcoxdiagnostics(paras_forr) #Test PH assumptions, Martingale residuals, bad model

#cox proportional hazard ratio with mixed effects -> bad model, non-proportional hazards
paras_forrcoxme <- coxme::coxme(survf ~ climate + diet + origin + (1|vial), data = paras_surv_4ana) #+strata(vial) for stratification
summary(paras_forrcoxme)

cox.zph(paras_forrcoxme) #Test the Proportional Hazards Assumptions -> non-proportional hazards! use another method
ggcoxzph(cox.zph(paras_forrcoxme)) #Test the PH Assumptions, Schoenfeld residuals -> non-proportional hazards! use another method

# Additive Cox model - Peto's correction for ties + stratification
paras_surv_4ana$vial <- as.factor(paras_surv_4ana$vial)

paras_addcox <- mgcv::gam(day ~ climate + diet + origin + s(vial, bs = 're'), data = paras_surv_4ana, family = "cox.ph", weights = death_hapn)
mgcv::gam.check(paras_addcox) #check GAM model, ok
summary(paras_addcox) #origin not significant, simplified model below without origin

# Additive Cox model - chosen model, Peto's correction for ties + stratification
paras_addcox_sim <- mgcv::gam(day ~ climate + diet + s(vial, bs = 're'), data = paras_surv_4ana, family = "cox.ph", weights = death_hapn)
mgcv::gam.check(paras_addcox_sim) #check GAM model, ok
summary(paras_addcox_sim)

#EMMs pairwise comparisons between groups
library(emmeans)
ref.paras_addcox <- emmeans::ref_grid(paras_addcox_sim)
mult.paras_addcox <- emmeans::emmeans(ref.paras_addcox, pairwise ~ climate + diet, adj="BH")

paras_GAMmultcomp <- multcomp::cld(mult.paras_addcox, adj="holm", Letters = letters, decreasing = T)
paras_GAMmultcomp

plot(multcomp::cld(mult.paras_addcox, adj="holm", Letters = letters, decreasing = T))
pwpp(mult.paras_addcox)

#median and IQR table
paras_surv_4ana_tab <- paras_surv_4ana %>% 
  dplyr::group_by(climate, diet, origin) %>% 
  dplyr::summarise(n = n(), colquant = quantile(day, c(0.25, 0.75)), quant = c(0.25, 0.75),
                   mean = mean(day), median = median(day), iqr = IQR(day)) %>% 
  tidyr::spread(quant, colquant)

#write.csv(paras_surv_4ana_tab, here::here("Tab_med-IQR_surv.csv"), row.names = T)



##### DEVELOPMENT -----
paras_devt <- tibble::as_tibble(readxl::read_excel("E_formosa_clipcages_survival_col.xlsx"))

'%!in%' <- function(x,y)!('%in%'(x,y))

paras_devt <- paras_devt %>%
  filter(!is.na(emerged_parasitoids)) %>%
  mutate(climate_cond = climate) %>% 
  unite(cage, climate_cond, col = rep_treat) %>% 
  dplyr::filter(rep_treat %!in% c("P5C1_future", "P6C1_future", "P7C1_future", "P7C2_future",
                                  "P7C1_present", "P7C2_present", "P8C1_present"))

#total number of individuals
paras_devt %>% group_by(climate) %>% summarise(sum(emerged_parasitoids))

#tidying
paras_devt <- paras_devt %>% 
  select(!c(survival_test, sampled_EtOH, notes))

#uncount emerged sum -> one row per individual
paras_devt_4ana <- paras_devt %>%
  uncount(emerged_parasitoids, .remove = T) %>% 
  mutate(emergence_hapn = 1)

#survival analysis
library(survival)
library(survminer)

devtf <- Surv(time = paras_devt_4ana$emergence_day, event = paras_devt_4ana$emergence_hapn)

#model fit on Kaplan-Meyer curves, climate treatment
devtf_fit <- survfit(devtf ~ climate, data = paras_devt_4ana) #+frailty(rep) #if needed for stratification
summary(devtf_fit)
surv_summary(devtf_fit)

#check differences with Log-Rank test (non parametric analysis, no assumptions, 1 categorical variable = climate, two groups)
#Log-Rank test -> super significative differences
survdiff(Surv(time = paras_devt_4ana$emergence_day, event = paras_devt_4ana$emergence_hapn) ~ climate, data = paras_devt_4ana)

#median and IQR table
paras_devt_4ana_tab <- paras_devt_4ana %>% 
  dplyr::group_by(climate) %>% 
  dplyr::summarise(n = n(), colquant = quantile(emergence_day, c(0.25, 0.75)), quant = c(0.25, 0.75),
                   mean = mean(emergence_day), median = median(emergence_day), iqr = IQR(emergence_day)) %>% 
  tidyr::spread(quant, colquant)

paras_devt_4ana_tab
#write.csv(paras_devt_4ana_tab, here::here("Tab_med-IQR_dev.csv"), row.names = T)

ggplot(paras_devt_4ana %>% group_by(climate, rep_treat) %>% summarise(total = sum(emergence_hapn)),
       aes(x=climate, y=total , col=climate)) +
  geom_point()

#plot developmental time curves, Kaplan-Meier estimator (non-parametric model, no baseline hazard nor covariates assumptions)
paras_devtplot <- ggsurvplot(devtf_fit, data = paras_devt_4ana, risk.table = "abs_pct",  ylab = "Emergence rate [%]",
                             xlab = "Emergence day [days after oviposition]", #linetype = c("strata"), palette = "Spectral",# size = 2,#facet.by = "sex_char",
                             #xlim = c(15, max(paras_devt_4ana$emergence_day)),
                             fun = "event",
                             conf.int = T,
                             conf.int.style = "ribbon",
                             legend.title = "", surv.median.line = "hv",
                             #log.rank.weights = "n", #pval = TRUE, pval.method = TRUE, pval.method.coord = c(20, 0.25),
                             legend.labs = c("Future", "Present"), break.time.by = 5,

)
paras_devtplot


colors = c("chocolate3", "darkolivegreen3")

paras_devtplot_sel <- paras_devtplot$plot + 
  #order as in legend.labs above
  scale_linetype_manual(values = c("solid", "solid")) +
  scale_colour_manual(values = colors) +
  scale_fill_manual(values = colors)

paras_devtplot_sel


paras_devtplot_patch <- paras_devtplot_sel / paras_devtplot$table + plot_layout(ncol = 1, heights = c(3, 1)) + theme(axis.text.y.left = element_text(color = rev(colors)))
paras_devtplot_patch 

ggsave(filename = "Figure1.tiff", plot = paras_devtplot_patch, dpi = 300, compression = "lzw", units = "px", width = 4000, height = 2700)



#Joint developmental time and survival Kaplan-Meier curves---
np_dat <- paras_devtplot$data.survplot
np_dat$time <- np_dat$time - max(np_dat$time)
np_dat$surv <- 1 - np_dat$surv 

df <- rbind(np_dat %>% mutate(origin = "LIST_nymphs"),
            paras_survplot$data.survplot %>% filter(diet == "honey") %>% select(!diet)) %>% filter(origin != "Koppert")

df <- as_tibble(df) %>% #adding points to start from the zero baseline
  add_row(time = -15, n.risk = 164, n.event = 1, n.censor = 0, surv = 0, std.err = Inf,
          upper = NA, lower = NA, strata = as.factor("climate=future"), climate = as.factor("future"), origin = as.factor("LIST_nymphs")) %>% 
  add_row(time = -8, n.risk = 12, n.event = 1, n.censor = 0, surv = 0, std.err = Inf,
          upper = NA, lower = NA, strata = as.factor("climate=present"), climate = as.factor("present"), origin = as.factor("LIST_nymphs")) #%>% 
  # add_row(time = 0, n.risk = 329, n.event = 1, n.censor = 0, surv = 1, std.err = Inf)
  #         upper = NA, lower = NA, strata = as.factor("climate=future"), climate = as.factor("future"), origin = as.factor("Koppert")) %>% 
  # add_row(time = 0, n.risk = 303, n.event = 1, n.censor = 0, surv = 1, std.err = Inf,
  #         upper = NA, lower = NA, strata = as.factor("climate=present"), climate = as.factor("present"), origin = as.factor("Koppert"))


every_nth <- function(x, nth, empty = TRUE, inverse = FALSE) 
{
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if(empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}

devtsurv_joint <- ggplot(df, aes(time, surv, color = climate)) + 
  geom_step(size = 0.9) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(y = 'Probability [%]',
       x = 'Time from emergence [days]',
       color = "Climate") +
       #linetype = "Nymphs origin")+
  #scale_linetype_manual(values = c(5, 1)) +
  scale_colour_manual(values = c("chocolate3", "darkolivegreen3"))+
  #scale_color_manual(palette = "Set1")+
  theme_minimal()+
  theme(panel.grid = element_blank(), panel.border = element_rect(fill = NA))+
  scale_x_continuous(breaks = seq(-15, 55, 1),
                     labels = every_nth(seq(-15, 55, 1), 10, inverse = TRUE))#+
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 15))

devtsurv_joint
ggsave(filename = "Figure3.tiff", plot = devtsurv_joint, dpi = 300, compression = "lzw", units = "px", width = 3000, height = 2000)

#calculate area under developmental time and survival curves
library("pracma")

df_future_LIST <- df %>% filter(climate == "future") %>% filter(origin == "LIST_nymphs")
AUC_future_LIST <- trapz(df_future_LIST$time, df_future_LIST$surv)

# df_future_Koppert <- df %>% filter(climate == "future") %>% filter(time <= 0 & origin == "LIST_nymphs" | origin == "Koppert")
# AUC_future_Koppert <- trapz(df_future_Koppert$time, df_future_Koppert$surv)

df_present_LIST <- df %>% filter(climate == "present") %>% filter(origin == "LIST_nymphs")
AUC_present_LIST <- trapz(df_present_LIST$time, df_present_LIST$surv)
# 
# df_present_Koppert <- df %>% filter(climate == "present") %>% filter(time <= 0 & origin == "LIST_nymphs" | origin == "Koppert")
# AUC_present_Koppert <- trapz(df_present_Koppert$time, df_present_Koppert$surv)

areas <- matrix(c("Future, LIST nymphs", "Present, LIST nymphs", 
                  AUC_future_LIST, AUC_present_LIST), ncol = 2, byrow = F)

colnames(areas) <- c("Climate", "Area")
areas

write.csv(areas, file = "Tab_areas-curves.csv", row.names = T)



##### PARASITATION RATE -----

paras_rate <- tibble::as_tibble(readxl::read_excel("E_formosa_clipcages_survival.xlsx"))

'%!in%' <- function(x,y)!('%in%'(x,y))

paras_rate <- paras_rate %>%
  dplyr::select(climate:total_emerged_parasitoids) %>% 
  mutate(climate_cond = climate) %>% 
  unite(cage, climate_cond, col = rep_treat) %>% 
  dplyr::filter(rep_treat %!in% c("P5C1_future", "P6C1_future", "P7C1_future", "P7C2_future",
                   "P7C1_present", "P7C2_present", "P8C1_present"))

#total emerged parasitoids, median and IQR table
paras_rate_tab <- paras_rate %>% 
  dplyr::group_by(climate) %>% 
  dplyr::summarise(n = n(), colquant_emerpar = quantile(total_emerged_parasitoids, c(0.25, 0.75)), quant_emerpar = c(0.25, 0.75),
                   mean = mean(total_emerged_parasitoids), median = median(total_emerged_parasitoids), iqr = IQR(total_emerged_parasitoids)) %>% 
                   tidyr::spread(quant_emerpar, colquant_emerpar)

#total number of parasitized nymphs, median and IQR table
paras_nparasn_tab <- paras_rate %>% 
  dplyr::group_by(climate) %>%
  dplyr::summarise(n = n(), colquant_nparasn = quantile(n_paras_4inst_exuviae , c(0.25, 0.75)), quant_nparasn = c(0.25, 0.75),
                  mean = mean(n_paras_4inst_exuviae), median = median(n_paras_4inst_exuviae), iqr = IQR(n_paras_4inst_exuviae)) %>% 
  tidyr::spread(quant_nparasn, colquant_nparasn)

#number of total 4th instar nymphs (available on the leaf to be parasitized), median and IQR table
paras_nnymphs_tab <- paras_rate %>% 
  dplyr::group_by(climate) %>% 
  dplyr::summarise(n = n(), colquant_nnymphs = quantile(n_4inst_exuviae, c(0.25, 0.75)), quant_nnymphs = c(0.25, 0.75),
                   mean = mean(n_4inst_exuviae), median = median(n_4inst_exuviae), iqr = IQR(n_4inst_exuviae)) %>% 
  tidyr::spread(quant_nnymphs, colquant_nnymphs)

#surviving adults at 72h oviposition, median and IQR table
paras_survad_tab <- paras_rate %>% 
  dplyr::group_by(climate) %>% 
  dplyr::summarise(n = n(), colquant_survad = quantile(alive, c(0.25, 0.75)), quant_survad = c(0.25, 0.75),
                   mean = mean(alive), median = median(alive), iqr = IQR(alive)) %>% 
  tidyr::spread(quant_survad, colquant_survad)

bind_rows(paras_survad_tab, paras_nnymphs_tab, paras_nparasn_tab, paras_rate_tab)
paras_summary_tab <- bind_rows(list("Alive adults at 72h oviposition" = paras_survad_tab,
               "Number of 4th instar Bemisia tabaci available nymphs" = paras_nnymphs_tab,
               "Number of visually found parasitized 4th instar nymphs" = paras_nparasn_tab,
               "Total number of emerged parasitoids" = paras_rate_tab), .id = "groups")

write.csv(paras_summary_tab, here::here("Tab_med-IQR_summary-parrate.csv"), row.names = T)

#compare different models for parasitization rate
library(lme4)
library(performance)
parrate_lm <- lm(total_emerged_parasitoids ~ climate  + n_4inst_exuviae + alive, data = paras_rate)
summary(parrate_lm)

parrate_glm <- glm(total_emerged_parasitoids ~ climate  + n_4inst_exuviae + alive, family = poisson, data = paras_rate)
summary(parrate_glm)

parrate_glmm_zipois <- glmmTMB::glmmTMB(total_emerged_parasitoids ~ climate  + n_4inst_exuviae + alive + (1|rep_treat), data = paras_rate, ziformula = ~1, family = poisson)
summary(parrate_glmm_zipois)

parrate_glmm <- lme4::glmer(total_emerged_parasitoids ~ climate  + n_4inst_exuviae + alive + (1|rep_treat), data = paras_rate, family = poisson)
summary(parrate_glmm)

parrate_glmm.nb <- lme4::glmer.nb(total_emerged_parasitoids ~ climate + n_paras_4inst_exuviae + n_4inst_exuviae + alive + (1|rep_treat), data = paras_rate, nAGQ = 0)
summary(parrate_glmm.nb)

parrate_glmm.nb2 <- lme4::glmer.nb(total_emerged_parasitoids ~ climate + n_4inst_exuviae + alive + (1|rep_treat), data = paras_rate)
summary(parrate_glmm.nb2)

compare_performance(parrate_lm, parrate_glm, parrate_glmm_zipois, parrate_glmm, parrate_glmm.nb, parrate_glmm.nb2)
plot(compare_performance(parrate_lm, parrate_glm, parrate_glmm_zipois, parrate_glmm, parrate_glmm.nb, parrate_glmm.nb2))

check_distribution(parrate_glmm.nb)
check_model(parrate_glmm.nb) #model ok, except for overdispersion

#more checks on the negative binomial model
library(DHARMa)
parrate_glmm.nb_simulation_refit <- simulateResiduals(fittedModel = parrate_glmm.nb)
plot(parrate_glmm.nb_simulation_refit)
testDispersion(parrate_glmm.nb_simulation_refit)
testZeroInflation(parrate_glmm.nb_simulation_refit)

#GLMM negative binomial
parrate_glmm.nb <- lme4::glmer.nb(total_emerged_parasitoids ~ climate + n_paras_4inst_exuviae + n_4inst_exuviae + alive + (1|rep_treat), data = paras_rate, nAGQ = 0)
parrate_glmm.nb
summary(parrate_glmm.nb)

#report (G)LMM table with jtools and huxtable
library("huxtable")
library("jtools")

jtools::export_summs(parrate_lm, parrate_glm, parrate_glmm_zipois, parrate_glmm, parrate_glmm.nb, parrate_glmm.nb2,
                     model.names = c("parrate_lm", "parrate_glm", "parrate_glmm_zipois", "parrate_glmm", "parrate_glmm.nb", "parrate_glmm.nb2"),
                     error_format = "({std.error})",
                     to.file = "xlsx", file.name = "Supp_mat_S3.xlsx",
                     exp = F, robust = F)
