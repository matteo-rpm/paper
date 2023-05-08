#parasitoids script
library(tidyverse)

#import dataset
paras_surv <- tibble::as_tibble(read.csv("~/E_eremicus_longevity_CC.csv"))

#total number of individuals
paras_surv %>% distinct(vial, .keep_all = T) %>% summarise(sum(total))

#tidying
paras_surv <- paras_surv %>% 
  select(!c(dead_percent, total))

#add column for number of dead that day, not cumulative
paras_surv$dead_that_day <- c(diff(paras_surv$dead), paras_surv$dead[length(paras_surv$dead)])

paras_surv <- paras_surv %>%
  mutate(dead_that_day=lag(dead_that_day))

paras_surv$dead_that_day[paras_surv$dead_that_day < 0] <- NA
paras_surv$dead_that_day[is.na(paras_surv$dead_that_day)] <- 0

#uncount dead sum -> one row per individual
paras_surv_4ana <- paras_surv %>%
  select(!dead) %>%
  uncount(dead_that_day, .remove = T) %>% 
  mutate(death_hapn = 1)

#survival analysis
library(survival)
library(survminer)

survf <- Surv(time = paras_surv_4ana$day, event = paras_surv_4ana$death_hapn)

#model fit on Kaplan-Meyer curves, climate+diet treatment
survf_fit <- survfit(survf ~ climate+diet, data = paras_surv_4ana) #+frailty(rep) #if needed for stratification
summary(survf_fit)
surv_summary(survf_fit)

#Log-Rank test
survdiff(Surv(time = paras_surv_4ana$day, event = paras_surv_4ana$death_hapn) ~ climate+diet, data = paras_surv_4ana)

#Pairwise comparisons -> super significative differences
pair_survd <- survminer::pairwise_survdiff(Surv(day, death_hapn) ~ climate+diet, data = paras_surv_4ana, p.adjust.method = "BH")
symnum(pair_survd$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
       symbols = c("****", "***", "**", "*", "+", "NS"),
       abbr.colnames = F, na = "")

# Aalen's Additive Regression Model -> covariate effects not multiplicative (Cox), but additive
aadd <- survival::aareg(Surv(day, death_hapn) ~ climate*diet, data = paras_surv_4ana)
summary(aadd)
summary(aadd)$table #covariate effects, p-val

library(ggfortify)
aaddf <- autoplot(aadd) #check covariates effects, plot

#plot developmental time curves, Kaplan-Meier estimator (non-parametric model, no baseline hazard nor covariates assumptions)
paras_survplot <- ggsurvplot(survf_fit, data = paras_surv_4ana, risk.table = "abs_pct",  ylab = "Survival rate [%]",
                          xlab = "Time [days after isolation]", linetype = c("strata"), palette = "Set1",# size = 2,#facet.by = "sex_char",
                          xlim= c(-1,36),
                          conf.int = T,
                          legend.title = "Groups", surv.median.line = "hv",
                          #log.rank.weights = "n", #pval = TRUE, pval.method = TRUE, pval.method.coord = c(20, 0.25),
                          legend.labs = c("Future, honey", "Future, negative control",
                                          "Present, honey", "Present, negative control"), break.time.by = 5
                          )

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

paras_survplot_sel <- paras_survplot$plot + 
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed"))+
  scale_x_continuous(breaks = seq(0, 35, 1),
                     labels = every_nth(seq(0, 35, 1), 5, inverse = TRUE))
  #scale_y_continuous(breaks = seq(0.00, 1.00, 0.1),
   #                  labels = every_nth(seq(0.00, 1.00, 0.1), 0.25, inverse = TRUE))
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 8))
  #scale_colour_manual(values = c("#1b9e77", "#1b9e77", "#d95f02", "#d95f02"))

library(patchwork)
paras_plotdevt_patch <- paras_survplot_sel / paras_survplot$table + plot_layout(ncol = 1, heights = c(3, 1)) #+ 
  #theme(axis.text.y.left = element_text(color = rev(c("#1b9e77", "#1b9e77", "#d95f02", "#d95f02"))))

ggsave(filename = here::here("~", "Work", "Parasitoids", "surv_plot.tiff"), plot = paras_plotdevt_patch, dpi = 320, compression = "lzw", width = 11, height = 6.3)

# Quantification of treatments effects
#cox proportional hazard ratio with fixed effects only + stratification per replicate -> bad model, non-proportional hazards
paras_forr <- coxph(survf ~ climate+diet + frailty(vial), data = paras_surv_4ana) #+frailty(rep)
summary(paras_forr)

cox.zph(paras_forr) #Test the Proportional Hazards Assumptions -> non-proportional hazards! use another method
ggcoxzph(cox.zph(paras_forr)) #Test the PH Assumptions, Schoenfeld residuals -> non-proportional hazards! use another method
ggcoxdiagnostics(paras_forr) #Test PH assumptions, Martingale residuals, bad model

#cox proportional hazard ratio with mixed effects -> bad model, non-proportional hazards
paras_forrcoxme <- coxme::coxme(survf ~ climate + diet + (1|vial), data = paras_surv_4ana) #+strata(vial) for stratification
summary(paras_forrcoxme)

cox.zph(paras_forrcoxme) #Test the Proportional Hazards Assumptions -> non-proportional hazards! use another method
ggcoxzph(cox.zph(paras_forrcoxme)) #Test the PH Assumptions, Schoenfeld residuals -> non-proportional hazards! use another method

# Additive Cox model - chosen model, Peto's correction for ties + stratification
paras_surv_4ana$vial <- as.factor(paras_surv_4ana$vial)
paras_addcox <- mgcv::gam(day ~ climate + diet + s(vial, bs = 're'), data = paras_surv_4ana, family = "cox.ph", weights = death_hapn)
summary(paras_addcox)

#EMMs pairwise comparisons between groups
library(emmeans)
mult.paras_addcox <- emmeans::emmeans(paras_addcox, pairwise ~ climate + diet, adj="BH")
paras_GAMmultcomp <- multcomp::cld(mult.paras_addcox$emmeans, adj="BH", Letters = letters)#, decreasing = T)

plot(multcomp::cld(mult.paras_addcox$emmeans, adj="holm", Letters = letters, decreasing = T))

#median and IQR table
paras_surv_4ana_tab <- paras_surv_4ana %>% 
  dplyr::group_by(climate, diet) %>% 
  dplyr::summarise(n = n(), colquant = quantile(day, c(0.25, 0.75)), quant = c(0.25, 0.75),
                   mean = mean(day), median = median(day), iqr = IQR(day)) %>% 
  tidyr::spread(quant, colquant)

#write.csv(paras_surv_4ana_tab, here::here("~", "Work", "Parasitoids", "Tab_med-IQR_surv.csv"), row.names = T)
