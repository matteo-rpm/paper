#Citation:
@article{,
  author =       "Matteo Ripamonti and Federico Maron and Daniele Cornara and Cristina Marzachì and Alberto Fereres and Domenico Bosco",
  title =        "{Leafhopper feeding behaviour on three grapevine cultivars with different susceptibilities to Flavescence dorée}",
  journal =      "Journal of Insect Physiology",
  volume =       "137",
  number =       "",
  pages =        "",
  year =         "2022",
  DOI =          "https://doi.org/10.1016/j.jinsphys.2022.104366"
}

#Rwaves ANALYSIS------
#load packages
devtools::install_github("mchiapello/Rwaves")
library("Rwaves")

#use with raw data as ".ANA" files----
#set working directory, containing all the ".ANA" marking files
#setwd("~/folder/subfolder")
#import marked file ".ANA"
#x <- readData(".", "ANA") #use when raw data as ".ANA" files

#import raw data from GitHub-----
#as table resulting from Rwaves::readData() function
x <- tidyr::as_tibble(rio::import("https://raw.githubusercontent.com/matteo-rpm/papers/main/Ripamonti_2021_Journal-of-Insect-Physiology/Ripamonti_2021_INSPHY_rawdata.csv"))

#run analysis
y <- rwaves(x)

#load packages
library("dplyr")
library(ggplot2)
library("ggpubr")

#DATASET MODIFICATIONS-------
#rename Variables
y <- y %>% 
  rename(n_np = f1_1,
         s_np = f2_1,
         s_2np = f3_1,
         n_Pr = f14,
         s_Pr = f24,
         s_C = f2_2,
         n_G = f1_7,
         s_G = f2_7,
         n_E1 = f1_4,
         s_E1 = f2_4,
         n_E2 = f89_5,
         n_sE2 = f90_5,
         s_E2 = f91,
         mean_E2 = f92_5,
         s_longestE2 = f93_5,
         s_E = f96_345,
         s_notE = f98,
         t_LE1.Z = f107,
         t_1sE2.exp = f109,
         t_1E2.exp = f112,
         percprobtime_E1 = f115_4,
         percprobtime_E2 = f119,
         percprobtime_E = f119E,
         percprobtime_C = f115_2,
         percprobtime_G = f115_7,
         E2index = f95,
         mean_fr_Ninterrup = f200,
         percNinterrup_E2 = f201,
         n_Ninterrupt_E2 = f202,
         s_npto1stprobe = f210,
         t_1st_sE2 = f150
  )

#modify variables' values when phloem phase does not occur
y$mean_E2[y$mean_E2 == 0] <- NA 
y$t_1st_sE2[y$t_1st_sE2 == 0] <- 28800 
y$t_1E2.exp[y$t_1E2.exp == 0] <- 28800 
y$t_1sE2.exp[y$t_1sE2.exp == 0] <- 28800 

#extract column "treatment" or "cultivar" from marked file name
ynew <- y %>%
  mutate(Cultivar = ifelse(grepl('BA84', File, ignore.case = T), 'Barbera',
                           ifelse(grepl('BRA20', File, ignore.case = T), 'Brachetto',
                                  ifelse(grepl('MO190', File, ignore.case = T), 'Moscato', 'marking_error')
                           )))

#extract column "sex" from marked file name
ynew2 <- ynew %>%
  mutate(Sex = ifelse(grepl('male', File, ignore.case = T), 'male',
                      ifelse(grepl('fem', File, ignore.case = T), 'female', 'sexmarking_error')
  ))

#extract db with hours-variables
hourvar <- ynew2 %>% select(File, Cultivar, Sex, f191_1:f191_28799.9, s_E2)

#remove hours-variables
ynew2 <- ynew2 %>% select(!f191_1:f191_28799.9)

#remove empty variables
ynew2 <- ynew2 %>% select(!c(n_E1, s_E1, t_LE1.Z, percprobtime_E1, s_E, percprobtime_E))

#FIGURE S4 & S5--------
#exclude non-phloem recordings
ynew2_phloem <- ynew2[!(ynew2$s_E2 == "0"), ]

#new dummy dataset
ynew3_phloem <-
  ynew2_phloem %>% 
  tidyr::pivot_longer(cols = n_np:s_npto1stprobe,names_to = "Variable", values_to = "values") %>% 
  arrange(Variable)

#order in boxplot facets
ynew3_phloem$Variable_f = factor(ynew3_phloem$Variable, levels=c("n_np", "s_np", "s_2np", "s_npto1stprobe", "n_Pr", "s_Pr", "s_C", "n_G",
                                                                 "s_G", "n_E2", "n_sE2", "s_E2", "mean_E2", "s_longestE2", "s_notE", 
                                                                 "t_1E2.exp", "t_1sE2.exp", "t_1st_sE2", "percprobtime_E2",
                                                                 "percprobtime_C", "percprobtime_G", "E2index",
                                                                 "mean_fr_Ninterrup", "percNinterrup_E2", "n_Ninterrupt_E2"))

#boxplot for the dummy dataset
sp_phloem <- ggplot(ynew3_phloem, aes(x=Cultivar, y=values, fill=Sex)) + geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=Sex))

#small boxplots for different variables x cultivar x sex
StEPG_phloem <- sp_phloem + facet_wrap( ~ Variable_f, scales = "free_y") +
  theme_minimal()+
  theme(panel.border = element_rect(fill="transparent"), axis.text.x = element_text(angle = 45, vjust = 0.7))+
  scale_fill_manual(values = c("#717171","#aaaaaa"))

StEPG_phloem

ggsave(filename = "Figure_S4.tiff", plot = StEPG_phloem, dpi = 320, compression = "lzw", width = 12, height = 10)

#exclude phloem recordings
ynew2_exclphlo <- ynew2[(ynew2$s_E2 == "0"), ]

#select non-phloem variables
ynew2_exclphlo <- ynew2_exclphlo %>% 
  select(File, Cultivar, Sex, percprobtime_C, percprobtime_G, n_G, n_np, n_Pr, s_C, s_G, s_notE, s_np, s_Pr, s_2np, s_npto1stprobe)

#new dummy dataset
ynew3_exclphlo <-
  ynew2_exclphlo %>% 
  select(File, Cultivar, Sex, percprobtime_C, percprobtime_G, n_G, n_np, n_Pr, s_C, s_G, s_notE, s_np, s_Pr, s_2np, s_npto1stprobe) %>% 
  tidyr::pivot_longer(cols = percprobtime_C:s_npto1stprobe,
                      names_to = "Variable", values_to = "values") %>% 
  arrange(Variable)

#order in boxplot facets
ynew3_exclphlo$Variable_f = factor(ynew3_exclphlo$Variable,
                                   levels=c("n_np", "s_np", "s_2np", "s_npto1stprobe", "n_Pr", "s_Pr", "s_C", "n_G",
                                            "s_G", "s_notE", "percprobtime_C", "percprobtime_G"))

#boxplot for the dummy dataset
sp_exclphlo <- ggplot(ynew3_exclphlo, aes(x=Cultivar, y=values, fill=Sex)) + geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=Sex))

#define comparisons
my_comp <- list(c("Barbera", "Moscato"), c("Barbera", "Brachetto"), c("Brachetto", "Moscato"))

#small boxplots for different variables x cultivar x sex
StEPG_exclphloem <- sp_exclphlo + facet_wrap( ~ Variable_f, scales = "free_y") +
  stat_compare_means(aes(group=Cultivar), comparisons=my_comp , label = "p.format",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("****", "***", "**", "*", "")))+
  theme_minimal()+
  theme(panel.border = element_rect(fill="transparent"), axis.text.x = element_text(angle = 45, vjust = 0.7))+
  scale_fill_manual(values = c("#717171","#aaaaaa"))

StEPG_exclphloem

ggsave(filename = "Figure_S5.tiff", plot = StEPG_exclphloem, dpi = 320, compression = "lzw", width = 12, height = 12)

##FIGURE 1-------------
#FIGURE 1A--------
#operations on dataset
ta <- hourvar %>%
  mutate(f190 = 1) %>% 
  select(File, f190, f191_3600:f191_28799.9, Cultivar, Sex)

#new dummy dataset
library(tidyr)
tat <- ta %>% 
  pivot_longer(
    cols = c(f190:f191_28799.9), 
    names_to = "hour", 
    values_to = "waveform",
    values_drop_na = T)

db <- tat %>% 
  group_by(Cultivar, hour, waveform) %>% 
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>% 
  mutate(perchour = n/tot*100)

library(tibble)
# fake db to avoid errors in temporal development graph
#(i.e. areas with empty spaces when one waveform has value = 0 for a specific hour)

#fake db_0 for all recordings
db_0 <- matrix(c("Barbera", 0, 2, 0, 31, 0, "Barbera", 0, 5, 0, 31, 0, "Barbera", 0, 7, 0, 31, 0,
                 "Brachetto", 0, 2, 0, 32, 0, "Brachetto", 0, 5, 0, 32, 0, "Brachetto", 0, 7, 0, 32, 0,
                 "Brachetto", 8, 7, 0, 32, 0,
                 "Moscato", 0, 2, 0, 37, 0, "Moscato", 0, 5, 0, 37, 0, "Moscato", 0, 7, 0, 37, 0), nrow = 10, byrow = T)

#operations on fake dataset
colnames(db_0) <- c("Cultivar", "hour", "waveform", "n", "tot", "perchour")
db_0 <- as_tibble(db_0)
db_0$hour <- as.numeric(db_0$hour)
db_0$n <- as.numeric(db_0$n)
db_0$tot <- as.numeric(db_0$tot)
db_0$perchour <- as.numeric(db_0$perchour)

#make hours column explict
db$hour[db$hour == "f190"] <- "0"
db$hour[db$hour == "f191_3600"] <- "1"
db$hour[db$hour == "f191_7200"] <- "2"
db$hour[db$hour == "f191_10800"] <- "3"
db$hour[db$hour == "f191_14400"] <- "4"
db$hour[db$hour == "f191_18000"] <- "5"
db$hour[db$hour == "f191_21600"] <- "6"
db$hour[db$hour == "f191_25200"] <- "7"
db$hour[db$hour == "f191_28799.9"] <- "8"

db$hour <- as.numeric(db$hour)
db$waveform <- as.character(db$waveform)

#join real and fake dataset
dbgraphs <- bind_rows(db, db_0)

#order waveforms
dbgraphs$waveform <- factor(dbgraphs$waveform , levels=c("5", "7", "2", "1"))

#single figures
#Barbera, all recordings
ggbaA <- ggplot(dbgraphs %>% filter(Cultivar=="Barbera"), aes(x=hour, y=perchour, fill=waveform))+
  geom_area(alpha=1 , size=1, colour="black", outline.type = "full", position = position_stack(reverse = F))+
  labs(title = "Barbera", tag = "A")+
  scale_fill_grey(name = "Waveform",
                  breaks=c("5", "7", "2", "1"),
                  labels=c("Phloem phase", "Xylem phase", "Pathway phase", "non-probing"),
                  start = 0)+
  theme_minimal()+
  ylab("Percentage of individuals [%]")+
  xlab(NULL)+
  theme(panel.grid = element_blank(), 
        axis.title.x=element_blank(), axis.text.x=element_blank())

#Brachetto, all recordings
ggbraA <- ggplot(dbgraphs %>% filter(Cultivar=="Brachetto"), aes(x=hour, y=perchour, fill=waveform))+
  geom_area(alpha=1 , size=1, colour="black", outline.type = "full", position = position_stack(reverse = F))+
  labs(title = "Brachetto")+
  scale_fill_grey(name = "Waveform",
                  breaks=c("5", "7", "2", "1"),
                  labels=c("Phloem phase", "Xylem phase", "Pathway phase", "non-probing"),
                  start = 0)+
  theme_minimal()+
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.title.x=element_blank(), axis.text.x=element_blank())

#Moscato, all recordings
ggmoA <- ggplot(dbgraphs %>% filter(Cultivar=="Moscato"), aes(x=hour, y=perchour, fill=waveform))+
  geom_area(alpha=1 , size=1, colour="black", outline.type = "full", position = position_stack(reverse = F))+
  labs(title = "Moscato")+
  scale_fill_grey(name = "Waveform",
                  breaks=c("5", "7", "2", "1"),
                  labels=c("Phloem phase", "Xylem phase", "Pathway phase", "non-probing"),
                  start = 0)+
  theme_minimal()+
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.title.x=element_blank(), axis.text.x=element_blank())

#Figure1B-----
#operations on dataset
taB <- hourvar %>%
  mutate(f190 = 1) %>% 
  filter(s_E2 != 0) %>% 
  select(File, f190, f191_3600:f191_28799.9, Cultivar, Sex)

#new dummy dataset
library(tidyr)
tatB <- taB %>% 
  pivot_longer(
    cols = c(f190:f191_28799.9), 
    names_to = "hour", 
    values_to = "waveform",
    values_drop_na = T)

dbB <- tatB %>% 
  group_by(Cultivar, hour, waveform) %>% 
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>% 
  mutate(perchour = n/tot*100)

library(tibble)
# fake db to avoid errors in temporal development graph
#(i.e. areas with empty spaces when one waveform has value = 0 for a specific hour)

#fake dbB_0 for phloem recordings only
dbB_0B <- matrix(c("Barbera", 0, 2, 0, 31, 0, "Barbera", 0, 7, 0, 31, 0, "Barbera", 0, 5, 0, 31, 0, "Barbera", 2, 7, 0, 31, 0,
                 "Brachetto", 0, 2, 0, 32, 0, "Brachetto", 0, 7, 0, 32, 0, "Brachetto", 0, 5, 0, 32, 0, "Brachetto", 1, 5, 0, 32, 0, "Brachetto", 8, 7, 0, 32, 0,
                 "Moscato", 0, 2, 0, 37, 0, "Moscato", 0, 7, 0, 37, 0, "Moscato", 0, 5, 0, 37, 0, "Moscato", 0, 2, 0, 37, 0, "Moscato", 4, 7, 0, 37, 0, "Moscato", 7, 7, 0, 37, 0), nrow = 15, byrow = T)

#operations on fake dataset
colnames(dbB_0B) <- c("Cultivar", "hour", "waveform", "n", "tot", "perchour")
dbB_0B <- as_tibble(dbB_0B)
dbB_0B$hour <- as.numeric(dbB_0B$hour)
dbB_0B$n <- as.numeric(dbB_0B$n)
dbB_0B$tot <- as.numeric(dbB_0B$tot)
dbB_0B$perchour <- as.numeric(dbB_0B$perchour)

#make hours column explicit
dbB$hour[dbB$hour == "f190"] <- "0"
dbB$hour[dbB$hour == "f191_3600"] <- "1"
dbB$hour[dbB$hour == "f191_7200"] <- "2"
dbB$hour[dbB$hour == "f191_10800"] <- "3"
dbB$hour[dbB$hour == "f191_14400"] <- "4"
dbB$hour[dbB$hour == "f191_18000"] <- "5"
dbB$hour[dbB$hour == "f191_21600"] <- "6"
dbB$hour[dbB$hour == "f191_25200"] <- "7"
dbB$hour[dbB$hour == "f191_28799.9"] <- "8"

dbB$hour <- as.numeric(dbB$hour)
dbB$waveform <- as.character(dbB$waveform)

#join real and fake dataset
dbBgraphs <- bind_rows(dbB, dbB_0B)

#order waveforms
dbBgraphs$waveform <- factor(dbBgraphs$waveform , levels=c("5", "7", "2", "1"))

#single figures
#Barbera, phloem recordings
ggbaB <- ggplot(dbBgraphs %>% filter(Cultivar=="Barbera"), aes(x=hour, y=perchour, fill=waveform))+
  geom_area(alpha=1 , size=1, colour="black", outline.type = "full", position = position_stack(reverse = F))+
  labs(tag = "B")+
  scale_fill_grey(name = "Waveform",
                  breaks=c("5", "7", "2", "1"),
                  labels=c("Phloem phase", "Xylem phase", "Pathway phase", "non-probing"),
                  start = 0)+
  theme_minimal()+
  theme(panel.grid = element_blank(), axis.title.x = element_blank())+
  ylab("Percentage of individuals [%]")

#Brachetto, phloem recordings
ggbraB <- ggplot(dbBgraphs %>% filter(Cultivar=="Brachetto"), aes(x=hour, y=perchour, fill=waveform))+
  geom_area(alpha=1 , size=1, colour="black", outline.type = "full", position = position_stack(reverse = F))+
  scale_fill_grey(name = "Waveform",
                  breaks=c("5", "7", "2", "1"),
                  labels=c("Phloem phase", "Xylem phase", "Pathway phase", "non-probing"),
                  start = 0)+
  theme_minimal()+
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())+
  xlab("Time [hour]")

#Moscato, phloem recordings
ggmoB <- ggplot(dbBgraphs %>% filter(Cultivar=="Moscato"), aes(x=hour, y=perchour, fill=waveform))+
  geom_area(alpha=1 , size=1, colour="black", outline.type = "full", position = position_stack(reverse = F))+
  scale_fill_grey(name = "Waveform",
                  breaks=c("5", "7", "2", "1"),
                  labels=c("Phloem phase", "Xylem phase", "Pathway phase", "non-probing"),
                  start = 0)+
  theme_minimal()+
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.title.x = element_blank())

#Join Figure 1A and Figure 1B
library(patchwork)
gghour <- ((ggbaA + ggbraA + ggmoA)/(ggbaB + ggbraB + ggmoB))+
  plot_layout(guides="collect")

gghour

ggsave(filename = "Figure1.tiff", plot = gghour, dpi = 320, compression = "lzw", width = 10, height = 7.25)

###UNIVARIATE ANALYSES-----------
##GLMs------------

dabob <- ynew2_phloem #OR ynew2 for Supplementary File S2
library(emmeans)
library(stringr)
library(car)
library(hnp)

#n_np
#glm - quasipoisson distribution
gqp.n_np <- glm(n_np ~ Cultivar, family = quasipoisson, data = dabob)
summary(gqp.n_np) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(n_np ~ Cultivar, family = quasipoisson, data = dabob) #test homogeneity of variance
hnp(gqp.n_np, main = "gqp.n_np") #fits perfectly
ref.gqp.n_np <- emmeans::ref_grid(gqp.n_np)
emm.gqp.n_np <- emmeans::emmeans(ref.gqp.n_np, specs = "Cultivar")
cld.gqp.n_np <- as_tibble(multcomp::cld(emm.gqp.n_np, Letters = "abcdef"))
cld.gqp.n_np <- rbind(cld.gqp.n_np, cld.gqp.n_np %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.n_np = tolower(.group))

#s_np
#glm - gamma distribution
gqp.s_np <- glm(s_np ~ Cultivar, family = Gamma, data = dabob)
summary(gqp.s_np) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(s_np ~ Cultivar, family = Gamma, data = dabob) #test homogeneity of variance
hnp(gqp.s_np, main = "gqp.s_np", how.many.out = T) #fits perfectly
ref.gqp.s_np <- emmeans::ref_grid(gqp.s_np)
emm.gqp.s_np <- emmeans::emmeans(ref.gqp.s_np, specs = "Cultivar")
cld.gqp.s_np <- as_tibble(multcomp::cld(emm.gqp.s_np, Letters = "abcdef"))
cld.gqp.s_np <- rbind(cld.gqp.s_np, cld.gqp.s_np %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.s_np = tolower(.group))

#s_2np
#glm - inverse.gaussian distribution
gqp.s_2np <- glm(s_2np ~ Cultivar, family = inverse.gaussian, data = dabob) #family = Gamma out of goodness-of-fit, inverse-gaussian fit precisely
summary(gqp.s_2np) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(s_2np ~ Cultivar, family = inverse.gaussian, data = dabob) #test homogeneity of variance
hnp(gqp.s_2np, main = "gqp.s_2np") #fits perfectly
ref.gqp.s_2np <- emmeans::ref_grid(gqp.s_2np)
emm.gqp.s_2np <- emmeans::emmeans(ref.gqp.s_2np, specs = "Cultivar")
cld.gqp.s_2np <- as_tibble(multcomp::cld(emm.gqp.s_2np, Letters = "abcdef"))
cld.gqp.s_2np <- rbind(cld.gqp.s_2np, cld.gqp.s_2np %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.s_2np = tolower(.group))

#n_Pr
#glm - quasipoisson distribution
gqp.n_Pr <- glm(n_Pr ~ Cultivar, family = quasipoisson, data = dabob)
summary(gqp.n_Pr) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(n_Pr ~ Cultivar, family = quasipoisson, data = dabob) #test homogeneity of variance
hnp(gqp.n_Pr, main = "gqp.n_Pr") #fits perfectly
ref.gqp.n_Pr <- emmeans::ref_grid(gqp.n_Pr)
emm.gqp.n_Pr <- emmeans::emmeans(ref.gqp.n_Pr, specs = "Cultivar")
cld.gqp.n_Pr <- as_tibble(multcomp::cld(emm.gqp.n_Pr, Letters = "abcdef")) 
cld.gqp.n_Pr <- rbind(cld.gqp.n_Pr, cld.gqp.n_Pr %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.n_Pr = tolower(.group))

#s_Pr
#glm - gamma distribution
gqp.s_Pr <- glm(s_Pr ~ Cultivar, family = Gamma, data = dabob)
summary(gqp.s_Pr) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(s_Pr ~ Cultivar, family = Gamma, data = dabob) #test homogeneity of variance
hnp(gqp.s_Pr, main = "gqp.s_Pr") #fits perfectly
ref.gqp.s_Pr <- emmeans::ref_grid(gqp.s_Pr)
emm.gqp.s_Pr <- emmeans::emmeans(ref.gqp.s_Pr, specs = "Cultivar")
cld.gqp.s_Pr <- as_tibble(multcomp::cld(emm.gqp.s_Pr, Letters = "abcdef"))
cld.gqp.s_Pr <- rbind(cld.gqp.s_Pr, cld.gqp.s_Pr %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.s_Pr = tolower(.group))

#s_C
#glm - gamma distribution
gqp.s_C <- glm(s_C ~ Cultivar * Sex, family = Gamma, data = dabob)
summary(gqp.s_C) # DIFFERENCEs for Sex (0.08) AND Cultivar*Sex (BA84female vs MO190male 0.0563), keep Cultivar * Sex
leveneTest(s_C ~ Cultivar * Sex, family = Gamma, data = dabob) #test homogeneity of variance
hnp(gqp.s_C, main = "gqp.s_C") #fits perfectly
ref.gqp.s_C <- emmeans::ref_grid(gqp.s_C)
emm.gqp.s_C <- emmeans::emmeans(ref.gqp.s_C, specs = c("Cultivar", "Sex"))
cld.gqp.s_C <- as_tibble(multcomp::cld(emm.gqp.s_C, Letters = "abcdef")) %>% 
  mutate(group.s_C = .group)

#n_G
#glm - negative-binomial distribution
gqp.n_G <- MASS::glm.nb(n_G ~ Cultivar * Sex, data = dabob[-11,]) #removed 1 outlier to better fit the model
#gqp.n_G <- glm(n_G ~ Cultivar*Sex, family = quasipoisson, data = dabob)
summary(gqp.n_G) # keep Cultivar * Sex
leveneTest(n_G ~ Cultivar*Sex, family = negative.binomial, data = dabob) #test homogeneity of variance
hnp(gqp.n_G, main = "gqp.n_G", how.many.out = T, sim = 999) #fits perfectly
ref.gqp.n_G <- emmeans::ref_grid(gqp.n_G)
emm.gqp.n_G <- emmeans::emmeans(ref.gqp.n_G, specs = c("Cultivar", "Sex"))
cld.gqp.n_G <- as_tibble(multcomp::cld(emm.gqp.n_G, Letters = "abcdef")) %>% 
  mutate(group.n_G = .group)

#s_G
#glm - gamma distribution
gqp.s_G <- glm(s_G ~ Cultivar, family = Gamma, data = dabob)
summary(gqp.s_G) # No differences for Sex NOR interaction (per poco), keep cultivar
leveneTest(s_G ~ Cultivar, family = Gamma, data = dabob) #test homogeneity of variance
hnp(gqp.s_G, main = "gqp.s_G") #fits perfectly
ref.gqp.s_G <- emmeans::ref_grid(gqp.s_G)
emm.gqp.s_G <- emmeans::emmeans(ref.gqp.s_G, specs = c("Cultivar"))
cld.gqp.s_G <- as_tibble(multcomp::cld(emm.gqp.s_G, Letters = "abcdef")) 
cld.gqp.s_G <- rbind(cld.gqp.s_G, cld.gqp.s_G %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.s_G = tolower(.group))

#n_E2
#glm - quasipoisson distribution
gqp.n_E2 <- glm(n_E2 ~ Cultivar, family = quasipoisson, data = dabob)
summary(gqp.n_E2) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(n_E2 ~ Cultivar, family = quasipoisson, data = dabob) #test homogeneity of variance *, maxvar/minvar = ~6
hnp(gqp.n_E2, main = "gqp.n_E2") #fits perfectly
ref.gqp.n_E2 <- emmeans::ref_grid(gqp.n_E2, vcov. = sandwich::vcovHC)
emm.gqp.n_E2 <- emmeans::emmeans(ref.gqp.n_E2, specs = c("Cultivar"))
cld.gqp.n_E2 <- as_tibble(multcomp::cld(emm.gqp.n_E2, Letters = "abcdef"))
cld.gqp.n_E2 <- rbind(cld.gqp.n_E2, cld.gqp.n_E2 %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.n_E2 = tolower(.group))

#n_sE2
#glm - quasipoisson distribution
gqp.n_sE2 <- glm(n_sE2 ~ Cultivar, family = quasipoisson, data = dabob)
summary(gqp.n_sE2) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(n_sE2 ~ Cultivar, family = quasipoisson, data = dabob) #test homogeneity of variance
hnp(gqp.n_sE2, main = "gqp.n_sE2") #fits perfectly
ref.gqp.n_sE2 <- emmeans::ref_grid(gqp.n_sE2)
emm.gqp.n_sE2 <- emmeans::emmeans(ref.gqp.n_sE2, specs = c("Cultivar"))
cld.gqp.n_sE2 <- as_tibble(multcomp::cld(emm.gqp.n_sE2, Letters = "abcdef"))
cld.gqp.n_sE2 <- rbind(cld.gqp.n_sE2, cld.gqp.n_sE2 %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.n_sE2 = tolower(.group))

#s_E2
#glm - gamma distribution
gqp.s_E2 <- glm(s_E2 ~ Cultivar, family = Gamma, data = dabob) #family = Gamma IF dabob is ynew2_phloem, gaussian if ynew2
summary(gqp.s_E2) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(s_E2 ~ Cultivar, family = Gamma, data = dabob) #test homogeneity of variance * , maxvar/minvar = 2.39
var(dabob %>% filter(Cultivar == "Barbera") %>% dplyr::select(s_E2))/var(dabob %>% filter(Cultivar == "Brachetto") %>% dplyr::select(s_E2))
#var(dabob %>% filter(Cultivar == "Moscato") %>% dplyr::select(s_E2)) #calculate variances between groups: if maxvar < 4 times minvar no need for heteroscedastically robust SE
hnp(gqp.s_E2, main = "gqp.s_E2") #fits perfectly
ref.gqp.s_E2 <- emmeans::ref_grid(gqp.s_E2)
emm.gqp.s_E2 <- emmeans::emmeans(ref.gqp.s_E2, specs = c("Cultivar"))
cld.gqp.s_E2 <- as_tibble(multcomp::cld(emm.gqp.s_E2, Letters = "abcdef"))
cld.gqp.s_E2 <- rbind(cld.gqp.s_E2, cld.gqp.s_E2 %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.s_E2 = tolower(.group))

#mean_E2
#glm - inverse.gaussian distribution
gqp.mean_E2 <- glm(mean_E2 ~ Cultivar, family = inverse.gaussian, data = dabob)
summary(gqp.mean_E2) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(mean_E2 ~ Cultivar, family = inverse.gaussian, data = dabob) #test homogeneity of variance
hnp(gqp.mean_E2, main = "gqp.mean_E2") #fits perfectly
ref.gqp.mean_E2 <- emmeans::ref_grid(gqp.mean_E2)
emm.gqp.mean_E2 <- emmeans::emmeans(ref.gqp.mean_E2, specs = c("Cultivar"))
cld.gqp.mean_E2 <- as_tibble(multcomp::cld(emm.gqp.mean_E2, Letters = "abcdef"))
cld.gqp.mean_E2 <- rbind(cld.gqp.mean_E2, cld.gqp.mean_E2 %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.mean_E2 = tolower(.group))

#s_longestE2
#glm - gamma distribution
gqp.s_longestE2 <- glm(s_longestE2 ~ Cultivar, family = Gamma, data = dabob[-42,])
#gqp.s_longestE2 <- update(gqp.s_longestE2, subset = -42) #exclude 1 outlier
summary(gqp.s_longestE2) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(s_longestE2 ~ Cultivar, family = Gamma, data = dabob[-42,]) #test homogeneity of variance ** , maxvar/minvar = 5.08
hnp(gqp.s_longestE2, main = "gqp.s_longestE2", how.many.out = T) #fits well (~7% obs out of the envelope)
ref.gqp.s_longestE2 <- emmeans::ref_grid(gqp.s_longestE2)#, vcov. = sandwich::vcovHC) outlier in MO190 exclude. Since hnp is good enough and variance between groups is ~5, exclude robust SE calc
emm.gqp.s_longestE2 <- emmeans::emmeans(ref.gqp.s_longestE2, specs = c("Cultivar"))
cld.gqp.s_longestE2 <- as_tibble(multcomp::cld(emm.gqp.s_longestE2, Letters = "abcdef"))
cld.gqp.s_longestE2 <- rbind(cld.gqp.s_longestE2, cld.gqp.s_longestE2 %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.s_longestE2 = tolower(.group))

#s_notE
#glm - gamma distribution
gqp.s_notE <- glm(s_notE ~ Cultivar * Sex, family = Gamma, data = dabob)
summary(gqp.s_notE) # NO differences for Cultivar*Sex, keep Cultivar + Sex
leveneTest(s_notE ~ Cultivar * Sex, family = Gamma, data = dabob) #test homogeneity of variance
hnp(gqp.s_notE, main = "gqp.s_notE") #fits perfectly
ref.gqp.s_notE <- emmeans::ref_grid(gqp.s_notE)
emm.gqp.s_notE <- emmeans::emmeans(ref.gqp.s_notE, specs = c("Cultivar", "Sex"))
cld.gqp.s_notE <- as_tibble(multcomp::cld(emm.gqp.s_notE, Letters = "abcdef")) %>% 
  mutate(group.s_notE = .group)

#t_1E2.exp
#glm - gamma distribution
gqp.t_1E2.exp <- glm(t_1E2.exp ~ Cultivar, family = Gamma, data = dabob)
summary(gqp.t_1E2.exp) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(t_1E2.exp ~ Cultivar, family = Gamma, data = dabob) #test homogeneity of variance * , maxvar/minvar = 6.12
hnp(gqp.t_1E2.exp, main = "gqp.t_1E2.exp") #fits perfectly
ref.gqp.t_1E2.exp <- emmeans::ref_grid(gqp.t_1E2.exp, vcov. = sandwich::vcovHC)
emm.gqp.t_1E2.exp <- emmeans::emmeans(ref.gqp.t_1E2.exp, specs = c("Cultivar"))
cld.gqp.t_1E2.exp <- as_tibble(multcomp::cld(emm.gqp.t_1E2.exp, Letters = "abcdef")) 
cld.gqp.t_1E2.exp <- rbind(cld.gqp.t_1E2.exp, cld.gqp.t_1E2.exp %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.t_1E2.exp = tolower(.group))

#t_1sE2.exp
#glm - inverse.gaussian distribution
gqp.t_1sE2.exp <- glm(t_1sE2.exp ~ Cultivar, family = inverse.gaussian, data = dabob)
summary(gqp.t_1sE2.exp) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(t_1sE2.exp ~ Cultivar, family = inverse.gaussian, data = dabob) #test homogeneity of variance ** , maxvar/minvar = 4.49
hnp(gqp.t_1sE2.exp, main = "gqp.t_1sE2.exp") #fits perfectly
ref.gqp.t_1sE2.exp <- emmeans::ref_grid(gqp.t_1sE2.exp, vcov. = sandwich::vcovHC)
emm.gqp.t_1sE2.exp <- emmeans::emmeans(ref.gqp.t_1sE2.exp, "Cultivar")
cld.gqp.t_1sE2.exp <- as_tibble(multcomp::cld(emm.gqp.t_1sE2.exp, Letters = "abcdef"))
cld.gqp.t_1sE2.exp <- rbind(cld.gqp.t_1sE2.exp, cld.gqp.t_1sE2.exp %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.t_1sE2.exp = tolower(.group))

#percprobtime_E2
#glm - beta-regression

y.transf.betareg <- function(y){
  n.obs <- sum(!is.na(y))
  (y * (n.obs - 1) + 0.5) / n.obs
}

library(betareg)
dabob <- dabob %>% mutate(percprobtime_E2 = percprobtime_E2/100)
gbr.percprobtime_E2 <- betareg(y.transf.betareg(percprobtime_E2) ~ Cultivar, data = dabob)
summary(gbr.percprobtime_E2) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar 
leveneTest(y.transf.betareg(percprobtime_E2) ~ Cultivar, data = dabob) #test homogeneity of variance
plot(gbr.percprobtime_E2, main = "gbr.percprobtime_E2", which = 5, type = "deviance", nsim = 99, level = .95) #fits perfectly
ref.percprobtime_E2 <- emmeans::ref_grid(gbr.percprobtime_E2)
emm.percprobtime_E2 <- emmeans::emmeans(ref.percprobtime_E2, "Cultivar")
cld.gqp.percprobtime_E2 <- as_tibble(multcomp::cld(emm.percprobtime_E2, Letters = "abcdef")) 
cld.gqp.percprobtime_E2 <- rbind(cld.gqp.percprobtime_E2, cld.gqp.percprobtime_E2 %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.percprobtime_E2 = tolower(.group))

#percprobtime_C
#glm - beta-regression
dabob <- dabob %>% mutate(percprobtime_C = percprobtime_C/100)
gbr.percprobtime_C <- betareg(percprobtime_C ~ Cultivar*Sex, data = dabob[-42,])
#gbr.percprobtime_C <- update(gbr.percprobtime_C, subset = -42) #excluding an outlier to better fit the model
leveneTest(y.transf.betareg(percprobtime_C) ~ Cultivar * Sex, data = dabob[-42,]) #test homogeneity of variance
summary(gbr.percprobtime_C) # keep all vars
plot(gbr.percprobtime_C, main = "gbr.percprobtime_C", which = 5, type = "deviance", nsim = 999, level = .95) # fits well (~5% points out of the envelope)
ref.percprobtime_C <- emmeans::ref_grid(gbr.percprobtime_C)
emm.percprobtime_C <- emmeans::emmeans(ref.percprobtime_C, specs = c("Cultivar", "Sex"))
cld.gqp.percprobtime_C <- as_tibble(multcomp::cld(emm.percprobtime_C, Letters = "abcdef")) %>% 
  mutate(group.percprobtime_C = .group)

#percprobtime_G
#glm - beta-regression
dabob <- dabob %>% mutate(percprobtime_G = percprobtime_G/100)
gbr.percprobtime_G <- betareg(percprobtime_G ~ Cultivar * Sex, data = dabob)
leveneTest(y.transf.betareg(percprobtime_G) ~ Cultivar * Sex, data = dabob) #test homogeneity of variance
summary(gbr.percprobtime_G) # keep all vars
plot(gbr.percprobtime_G, main = "gbr.percprobtime_G", which = 5, type = "deviance", nsim = 99, level = .95) # fits perfectly
ref.percprobtime_G <- emmeans::ref_grid(gbr.percprobtime_G)
emm.percprobtime_G <- emmeans::emmeans(ref.percprobtime_G, specs = c("Cultivar", "Sex"))
cld.gqp.percprobtime_G <- as_tibble(multcomp::cld(emm.percprobtime_G, Letters = "abcdef")) %>% 
  mutate(group.percprobtime_G = .group)


#E2index
#glm - beta-regression

#dabob$E2index[is.na(dabob$E2index)] <- median(dabob$E2index, na.rm=TRUE) #use if NAs are present (like in ynew2 || all recs)
dabob <- dabob %>% mutate(E2index = E2index/100)
gbr.E2index <- betareg(y.transf.betareg(E2index) ~ Cultivar, data = dabob)
summary(gbr.E2index) # NO differences for Sex NOR interaction
leveneTest(y.transf.betareg(E2index) ~ Cultivar, data = dabob) #test homogeneity of variance
shapiro.test(gbr.E2index$residuals) # test normality of residuals
plot(gbr.E2index, main = "gbr.E2index", which = 5, type = "deviance", nsim = 99, level = .95) # fits perfectly
ref.E2index <- emmeans::ref_grid(gbr.E2index)
emm.E2index <- emmeans::emmeans(ref.E2index, specs = c("Cultivar"))
cld.gqp.E2index <- as_tibble(multcomp::cld(emm.E2index, Letters = "abcdef"))
cld.gqp.E2index <- rbind(cld.gqp.E2index, cld.gqp.E2index %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.E2index = tolower(.group))


#mean_fr_Ninterrup
#glm - gamma distribution
gqp.mean_fr_Ninterrup <- glm(mean_fr_Ninterrup ~ Cultivar*Sex, family = Gamma, data = dabob)
summary(gqp.mean_fr_Ninterrup) # Sex and Cultivar appears differents
leveneTest(mean_fr_Ninterrup ~ Cultivar*Sex, family = Gamma, data = dabob) #test homogeneity of variance .
hnp(gqp.mean_fr_Ninterrup, main = "gqp.mean_fr_Ninterrup") # fits perfectly
ref.mean_fr_Ninterrup <- emmeans::ref_grid(gqp.mean_fr_Ninterrup)
emm.mean_fr_Ninterrup <- emmeans::emmeans(ref.mean_fr_Ninterrup, specs = c("Cultivar", "Sex"))
cld.gqp.mean_fr_Ninterrup <- as_tibble(multcomp::cld(emm.mean_fr_Ninterrup, Letters = "abcdef")) %>% 
  mutate(group.mean_fr_Ninterrup = .group)

#percNinterrup_E2
#glm - betareg distribution
dabob <- dabob %>% mutate(percNinterrup_E2 = percNinterrup_E2/100)
gqp.percNinterrup_E2 <- betareg(y.transf.betareg(percNinterrup_E2) ~ Cultivar, data = dabob) 
summary(gqp.percNinterrup_E2) # NO differences for Sex NOR interaction
leveneTest(y.transf.betareg(percNinterrup_E2) ~ Cultivar, data = dabob) #test homogeneity of variance * , maxvar/minvar = 33.25
plot(gqp.percNinterrup_E2, main = "gqp.percNinterrup_E2", which = 5, type = "deviance", nsim = 999, level = .95)  #fits well (~6% obs out of the envelope) 
ref.percNinterrup_E2 <- emmeans::ref_grid(gqp.percNinterrup_E2, vcov. = sandwich::vcovHC)
emm.percNinterrup_E2 <- emmeans::emmeans(ref.percNinterrup_E2, specs = c("Cultivar"))
cld.gqp.percNinterrup_E2 <- as_tibble(multcomp::cld(emm.percNinterrup_E2, Letters = "abcdef"))
cld.gqp.percNinterrup_E2 <- rbind(cld.gqp.percNinterrup_E2, cld.gqp.percNinterrup_E2 %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.percNinterrup_E2 = tolower(.group))

#n_Ninterrupt_E2
#glm - quasipoisson distribution
gqp.n_Ninterrupt_E2 <- glm(n_Ninterrupt_E2 ~ Cultivar, family = quasipoisson, data = dabob)
summary(gqp.n_Ninterrupt_E2) # NO differences for Sex NOR interaction
leveneTest(n_Ninterrupt_E2 ~ Cultivar, family = quasipoisson, data = dabob) #test homogeneity of variance
hnp(gqp.n_Ninterrupt_E2, main = "gqp.n_Ninterrupt_E2") #fits perfectly
ref.n_Ninterrupt_E2 <- emmeans::ref_grid(gqp.n_Ninterrupt_E2)
emm.n_Ninterrupt_E2 <- emmeans::emmeans(ref.n_Ninterrupt_E2, specs = c("Cultivar"))
cld.gqp.n_Ninterrupt_E2 <- as_tibble(multcomp::cld(emm.n_Ninterrupt_E2, Letters = "abcdef"))
cld.gqp.n_Ninterrupt_E2 <- rbind(cld.gqp.n_Ninterrupt_E2, cld.gqp.n_Ninterrupt_E2 %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.n_Ninterrupt_E2 = tolower(.group))

#s_npto1stprobe
#glm - inverse.gaussian distribution
gqp.s_npto1stprobe <- glm(s_npto1stprobe ~ Cultivar, family = inverse.gaussian, data = dabob)
summary(gqp.s_npto1stprobe) # NO differences
leveneTest(s_npto1stprobe ~ Cultivar, family = inverse.gaussian, data = dabob) #test homogeneity of variance
hnp(gqp.s_npto1stprobe, main = "gqp.s_npto1stprobe", how.many.out = T) # fits perfectly
ref.s_npto1stprobe <- emmeans::ref_grid(gqp.s_npto1stprobe)
emm.s_npto1stprobe <- emmeans::emmeans(ref.s_npto1stprobe, specs = c("Cultivar"))
cld.gqp.s_npto1stprobe <- as_tibble(multcomp::cld(emm.s_npto1stprobe, Letters = "abcdef")) 
cld.gqp.s_npto1stprobe <- rbind(cld.gqp.s_npto1stprobe, cld.gqp.s_npto1stprobe %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.s_npto1stprobe = tolower(.group))

#t_1st_sE2
#glm - inverse.gaussian distribution
gqp.t_1st_sE2 <- glm(t_1st_sE2 ~ Cultivar, family = inverse.gaussian, data = dabob)
summary(gqp.t_1st_sE2) # NO differences for Sex NOR interaction
leveneTest(t_1st_sE2 ~ Cultivar, family = inverse.gaussian, data = dabob) #test homogeneity of variance ** , maxvar/minvar = 4.68
hnp(gqp.t_1st_sE2, main = "gqp.t_1st_sE2", how.many.out = T) # fits perfectly
ref.t_1st_sE2 <- emmeans::ref_grid(gqp.t_1st_sE2, vcov. = sandwich::vcovHC)
emm.t_1st_sE2 <- emmeans::emmeans(ref.t_1st_sE2, specs = c("Cultivar"))
cld.gqp.t_1st_sE2 <- as_tibble(multcomp::cld(emm.t_1st_sE2, Letters = "abcdef"))
cld.gqp.t_1st_sE2 <- rbind(cld.gqp.t_1st_sE2, cld.gqp.t_1st_sE2 %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.t_1st_sE2 = tolower(.group))

#SUPPLEMENTARY FILE S3a-b----------
#report GLMs details in a supplementary table
library(huxtable)
library(jtools)

#GLMs
tables2a1 <- jtools::export_summs(gqp.n_np, gqp.s_np, gqp.n_Pr, gqp.s_Pr, gqp.s_C,
                                  gqp.s_G, gqp.n_sE2, gqp.s_E2, 
                                  gqp.s_notE, gbr.percprobtime_E2,
                                  gbr.percprobtime_G, gbr.E2index, gqp.mean_fr_Ninterrup,
                                  gqp.n_Ninterrupt_E2, 
                                  model.names = c("quasiP.n_np", "gamma.s_np", "quasiP.n_Pr", "gamma.s_Pr", "gamma.s_C",
                                                  "gamma.s_G", "quasiP.n_sE2",  "gamma.s_E2",
                                                  "gamma.s_notE", "betareg.percprobtime_E2",
                                                  "betareg.percprobtime_G", "betareg.E2index", "gamma.mean_fr_Ninterrup",
                                                  "quasiP.n_Ninterrupt_E2"),
                                  error_format = "({std.error})",
                                  to.file = "xlsx", file.name = "Table_S2a1_glms.xlsx",
                                  exp = F, robust = F)

#inverse.gaussian GLMs -- if not working, use 'summ' command and export coefs manually
jtools::export_summs(gqp.s_2np, gqp.mean_E2, gqp.s_npto1stprobe,
                     model.names = c("inv.gau.s_2np", "inv.gau.mean_E2", "inv.gau.s_npto1stprobe"),
                     error_format = "({std.error})",
                     to.file = "xlsx", file.name = "Table_S2a12_glms.xlsx")

jtools::summ(gqp.s_2np)
jtools::summ(gqp.mean_E2)
jtools::summ(gqp.s_npto1stprobe)

#GLMs with 1 outlier removed (N inferior to other GLMs -> gives error if exported together)
tables2a2 <- jtools::export_summs(gqp.n_G, gqp.s_longestE2, gbr.percprobtime_C, 
                                  model.names = c("neg.bin.n_G", "gamma.s_longestE2","betareg.percprobtime_C"),
                                  error_format = "({std.error})",
                                  to.file = "xlsx", file.name = "Table_S2a2_glms.xlsx",
                                  exp = F, robust = F)

#GLMs with heteroskedasticity robust SE
tables2b <- jtools::export_summs(gqp.n_E2, gqp.t_1E2.exp, gqp.t_1sE2.exp,
                                 gqp.percNinterrup_E2, gqp.t_1st_sE2,
                                 model.names = c("quasiP.n_E2", "gamma.t_1E2.exp", "inv.gau.t_1sE2.exp",
                                                 "betareg.percNinterrup_E2", "inv.gau.t_1st_sE2"),
                                 error_format = "({std.error})",
                                 to.file = "xlsx", file.name = "Table_S2b_heterosked-rob.xlsx",
                                 exp = F, robust = TRUE)

detach("package:hnp", unload = TRUE)
detach("package:MASS", unload = TRUE)

#TABLES 3 to 6------------
#preliminary operations: calculation of median and SE for every variable
library(dplyr)
dub <- ynew2_phloem

#new tibble: column "median.n_np"
a <-
  dub %>% 
  group_by(Cultivar, Sex, n_np) %>% 
  summarise(n = n()) %>% 
  mutate(median.n_np = median(n_np)) %>% 
  mutate(se.n_np = sd(n_np)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.n_np, se.n_np)

#new tibble: column "median.s_np"
b <-
  dub %>% 
  group_by(Cultivar, Sex, s_np) %>% 
  summarise(n = n()) %>% 
  mutate(median.s_np = median(s_np)) %>% 
  mutate(se.s_np = sd(s_np)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.s_np, se.s_np)

#new tibble: column "median.s_2np"
c <-
  dub %>% 
  group_by(Cultivar, Sex, s_2np) %>% 
  summarise(n = n()) %>% 
  mutate(median.s_2np = median(s_2np))%>% 
  mutate(se.s_2np = sd(s_2np)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.s_2np, se.s_2np)

#new tibble: column "median.n_Pr"
d <-
  dub %>% 
  group_by(Cultivar, Sex, n_Pr) %>% 
  summarise(n = n()) %>% 
  mutate(median.n_Pr = median(n_Pr)) %>% 
  mutate(se.n_Pr = sd(n_Pr)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.n_Pr, se.n_Pr)

#new tibble: column "median.s_Pr"
e <-
  dub %>% 
  group_by(Cultivar, Sex, s_Pr) %>% 
  summarise(n = n()) %>% 
  mutate(median.s_Pr = median(s_Pr)) %>% 
  mutate(se.s_Pr = sd(s_Pr)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.s_Pr, se.s_Pr)

#new tibble: column "median.s_C"
f <-
  dub %>% 
  group_by(Cultivar, Sex, s_C) %>% 
  summarise(n = n()) %>% 
  mutate(median.s_C = median(s_C)) %>% 
  mutate(se.s_C = sd(s_C)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.s_C, se.s_C)

#new tibble: column "median.n_G"
g <-
  dub %>% 
  group_by(Cultivar, Sex, n_G) %>% 
  summarise(n = n()) %>% 
  mutate(median.n_G = median(n_G)) %>% 
  mutate(se.n_G = sd(n_G)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.n_G, se.n_G)

#new tibble: column "median.s_G"
h <-
  dub %>% 
  group_by(Cultivar, Sex, s_G) %>% 
  summarise(n = n()) %>% 
  mutate(median.s_G = median(s_G)) %>% 
  mutate(se.s_G = sd(s_G)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.s_G, se.s_G)

#new tibble: column "median.n_E2"
i <-
  dub %>% 
  group_by(Cultivar, Sex, n_E2) %>% 
  summarise(n = n()) %>% 
  mutate(median.n_E2 = median(n_E2)) %>% 
  mutate(se.n_E2 = sd(n_E2)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.n_E2, se.n_E2)

#new tibble: column "median.n_sE2"
j <-
  dub %>% 
  group_by(Cultivar, Sex, n_sE2) %>% 
  summarise(n = n()) %>% 
  mutate(median.n_sE2 = median(n_sE2)) %>% 
  mutate(se.n_sE2 = sd(n_sE2)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.n_sE2, se.n_sE2)

#new tibble: column "median.s_E2"
k <-
  dub %>% 
  group_by(Cultivar, Sex, s_E2) %>% 
  summarise(n = n()) %>% 
  mutate(median.s_E2 = median(s_E2)) %>% 
  mutate(se.s_E2 = sd(s_E2)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.s_E2, se.s_E2)

#new tibble: column "median.mean_E2"
l <-
  dub %>% 
  group_by(Cultivar, Sex, mean_E2) %>% 
  summarise(n = n()) %>% 
  mutate(median.mean_E2 = median(mean_E2)) %>% 
  mutate(se.mean_E2 = sd(mean_E2)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.mean_E2, se.mean_E2)

#new tibble: column "median.s_longestE2"
m <-
  dub %>% 
  group_by(Cultivar, Sex, s_longestE2) %>% 
  summarise(n = n()) %>% 
  mutate(median.s_longestE2 = median(s_longestE2)) %>% 
  mutate(se.s_longestE2 = sd(s_longestE2)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.s_longestE2, se.s_longestE2)

#new tibble: column "median.s_notE"
n <-
  dub %>% 
  group_by(Cultivar, Sex, s_notE) %>% 
  summarise(n = n()) %>% 
  mutate(median.s_notE = median(s_notE)) %>% 
  mutate(se.s_notE = sd(s_notE)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.s_notE, se.s_notE)

#new tibble: column "median.t_1sE2.exp"
o <-
  dub %>% 
  group_by(Cultivar, Sex, t_1sE2.exp) %>% 
  summarise(n = n()) %>% 
  mutate(median.t_1sE2.exp = median(t_1sE2.exp)) %>% 
  mutate(se.t_1sE2.exp = sd(t_1sE2.exp)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.t_1sE2.exp, se.t_1sE2.exp)

#new tibble: column "median.t_1E2.exp"
p <-
  dub %>% 
  group_by(Cultivar, Sex, t_1E2.exp) %>% 
  summarise(n = n()) %>% 
  mutate(median.t_1E2.exp = median(t_1E2.exp)) %>% 
  mutate(se.t_1E2.exp = sd(t_1E2.exp)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.t_1E2.exp, se.t_1E2.exp)

#new tibble: column "median.percprobtime_E2"
q <-
  dub %>% 
  group_by(Cultivar, Sex, percprobtime_E2) %>% 
  summarise(n = n()) %>% 
  mutate(median.percprobtime_E2 = median(percprobtime_E2)) %>% 
  mutate(se.percprobtime_E2 = sd(percprobtime_E2)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.percprobtime_E2, se.percprobtime_E2)

#new tibble: column "median.percprobtime_C"
r <-
  dub %>% 
  group_by(Cultivar, Sex, percprobtime_C) %>% 
  summarise(n = n()) %>% 
  mutate(median.percprobtime_C = median(percprobtime_C)) %>% 
  mutate(se.percprobtime_C = sd(percprobtime_C)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.percprobtime_C, se.percprobtime_C)

#new tibble: column "median.percprobtime_G"
s <-
  dub %>% 
  group_by(Cultivar, Sex, percprobtime_G) %>% 
  summarise(n = n()) %>% 
  mutate(median.percprobtime_G = median(percprobtime_G)) %>% 
  mutate(se.percprobtime_G = sd(percprobtime_G)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.percprobtime_G, se.percprobtime_G)

#new tibble: column "median.E2index"
dub$E2index[is.na(dub$E2index)] <- with(dub, ave(E2index, Cultivar,
                                                 FUN = function(x) median(x, na.rm = TRUE)))[is.na(dub$E2index)]
t <-
  dub %>% 
  group_by(Cultivar, Sex, E2index) %>% 
  summarise(n = n()) %>% 
  mutate(median.E2index = median(E2index)) %>% 
  mutate(se.E2index = sd(E2index)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.E2index, se.E2index)

#new tibble: column "median.mean_fr_Ninterrup"
u <-
  dub %>% 
  group_by(Cultivar, Sex, mean_fr_Ninterrup) %>% 
  summarise(n = n()) %>% 
  mutate(median.mean_fr_Ninterrup = median(mean_fr_Ninterrup)*1000) %>% 
  mutate(se.mean_fr_Ninterrup = sd(mean_fr_Ninterrup)/sqrt(sum(n))*1000) %>%
  distinct(Cultivar, Sex, median.mean_fr_Ninterrup, se.mean_fr_Ninterrup)

#new tibble: column "median.percNinterrup_E2"
v <-
  dub %>% 
  group_by(Cultivar, Sex, percNinterrup_E2) %>% 
  summarise(n = n()) %>% 
  mutate(median.percNinterrup_E2 = median(percNinterrup_E2)) %>% 
  mutate(se.percNinterrup_E2 = sd(percNinterrup_E2)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.percNinterrup_E2, se.percNinterrup_E2)

#new tibble: column "median.n_Ninterrupt_E2"
w <-
  dub %>% 
  group_by(Cultivar, Sex, n_Ninterrupt_E2) %>% 
  summarise(n = n()) %>% 
  mutate(median.n_Ninterrupt_E2 = median(n_Ninterrupt_E2)) %>% 
  mutate(se.n_Ninterrupt_E2 = sd(n_Ninterrupt_E2)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.n_Ninterrupt_E2, se.n_Ninterrupt_E2)

#new tibble: column "median.s_npto1stprobe"
at <-
  dub %>% 
  group_by(Cultivar, Sex, s_npto1stprobe) %>% 
  summarise(n = n()) %>% 
  mutate(median.s_npto1stprobe = median(s_npto1stprobe)) %>% 
  mutate(se.s_npto1stprobe = sd(s_npto1stprobe)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.s_npto1stprobe, se.s_npto1stprobe)

#new tibble: column "median.t_1st_sE2"
bt <-
  dub %>% 
  group_by(Cultivar, Sex, t_1st_sE2) %>% 
  summarise(n = n()) %>% 
  mutate(median.t_1st_sE2 = median(t_1st_sE2)) %>% 
  mutate(se.t_1st_sE2 = sd(t_1st_sE2)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.t_1st_sE2, se.t_1st_sE2)

#report median, SE, GLMs post-hoc in Tables 3 to 6

#join all variables dfs
dfr <- list(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,at,bt,
            cld.gqp.n_np, cld.gqp.s_np, cld.gqp.s_2np, cld.gqp.n_Pr, cld.gqp.s_Pr, cld.gqp.s_C, cld.gqp.n_G,
            cld.gqp.s_G, cld.gqp.n_E2, cld.gqp.n_sE2, cld.gqp.n_Ninterrupt_E2, cld.gqp.s_E2,
            cld.gqp.mean_E2, cld.gqp.s_longestE2, cld.gqp.s_notE, cld.gqp.t_1E2.exp, cld.gqp.t_1sE2.exp,
            cld.gqp.percprobtime_E2, cld.gqp.percprobtime_C, cld.gqp.percprobtime_G, cld.gqp.E2index,
            cld.gqp.mean_fr_Ninterrup, cld.gqp.percNinterrup_E2, cld.gqp.s_npto1stprobe, cld.gqp.t_1st_sE2)

#intermediate table
median_gentab <- plyr::join_all(dfr, by = c("Cultivar", "Sex"), type='full') %>%
  mutate(across(where(is.numeric), round, 1))

#general table
#unite median +- SE
library(tidyr)

median_gentab <- median_gentab %>% 
  unite(medianse.n_np, median.n_np|se.n_np, sep = " ± ", remove = T) %>%
  unite("Number of non-probing periods *", medianse.n_np|group.n_np, sep = " ", remove = T) %>% 
  unite(medianse.s_np, median.s_np|se.s_np, sep = " ± ", remove = T) %>%
  unite("Total duration of non-probing periods [min] *", medianse.s_np|group.s_np, sep = " ", remove = T) %>%
  unite(medianse.s_2np, median.s_2np|se.s_2np, sep = " ± ", remove = T) %>%
  unite("Duration of the 2nd non-probing period [s] *", medianse.s_2np|group.s_2np, sep = " ", remove = T) %>%
  unite(medianse.n_Pr, median.n_Pr|se.n_Pr, sep = " ± ", remove = T) %>% 
  unite("Number of probes *", medianse.n_Pr|group.n_Pr, sep = " ", remove = T) %>% 
  unite(medianse.s_Pr, median.s_Pr|se.s_Pr, sep = " ± ", remove = T) %>% 
  unite("Total probing time [min] *", medianse.s_Pr|group.s_Pr, sep = " ", remove = T) %>% 
  unite(medianse.s_C, median.s_C|se.s_C, sep = " ± ", remove = T) %>% 
  unite("Total duration of pathway phase [min] **", medianse.s_C|group.s_C, sep = " ", remove = T) %>% 
  unite(medianse.n_G, median.n_G|se.n_G, sep = " ± ", remove = T) %>% 
  unite("Number of active ingestion phases **", medianse.n_G|group.n_G, sep = " ", remove = T) %>% 
  unite(medianse.s_G, median.s_G|se.s_G, sep = " ± ", remove = T) %>% 
  unite("Total duration of active ingestion [min] *", medianse.s_G|group.s_G, sep = " ", remove = T) %>% 
  unite(medianse.n_E2, median.n_E2|se.n_E2, sep = " ± ", remove = T) %>% 
  unite("Number of phloem ingestions *", medianse.n_E2|group.n_E2, sep = " ", remove = T) %>% 
  unite(medianse.n_sE2, median.n_sE2|se.n_sE2, sep = " ± ", remove = T) %>% 
  unite("Number of sustained (> 10 min) phloem ingestion *", medianse.n_sE2|group.n_sE2, sep = " ", remove = T) %>% 
  unite(medianse.s_E2, median.s_E2|se.s_E2, sep = " ± ", remove = T) %>% 
  unite("Total duration of phloem ingestions [min] *", medianse.s_E2|group.s_E2, sep = " ", remove = T) %>% 
  unite(medianse.mean_E2, median.mean_E2|se.mean_E2, sep = " ± ", remove = T) %>% 
  unite("Mean duration of a single event of phloem ingestion [min] *", medianse.mean_E2|group.mean_E2, sep = " ", remove = T) %>% 
  unite(medianse.s_longestE2, median.s_longestE2|se.s_longestE2, sep = " ± ", remove = T) %>% 
  unite("Duration of the longest phloem ingestion [min] *", medianse.s_longestE2|group.s_longestE2, sep = " ", remove = T) %>% 
  unite(medianse.s_notE, median.s_notE|se.s_notE, sep = " ± ", remove = T) %>% 
  unite("Total duration of non-phloematic phases [min] **", medianse.s_notE|group.s_notE, sep = " ", remove = T) %>% 
  unite(medianse.t_1sE2.exp, median.t_1sE2.exp|se.t_1sE2.exp, sep = " ± ", remove = T) %>% 
  unite("Time from 1st probe to 1st sustained (> 10 min) phloem ingestion [min] *", medianse.t_1sE2.exp|group.t_1sE2.exp, sep = " ", remove = T) %>% 
  unite(medianse.t_1E2.exp, median.t_1E2.exp|se.t_1E2.exp, sep = " ± ", remove = T) %>% 
  unite("Time from 1st probe to 1st phloem ingestion [min] *", medianse.t_1E2.exp|group.t_1E2.exp, sep = " ", remove = T) %>% 
  unite(medianse.percprobtime_E2, median.percprobtime_E2|se.percprobtime_E2, sep = " ± ", remove = T) %>% 
  unite("Percentage of probing time spent in phloem ingestion [%] *", medianse.percprobtime_E2|group.percprobtime_E2, sep = " ", remove = T) %>% 
  unite(medianse.percprobtime_C, median.percprobtime_C|se.percprobtime_C, sep = " ± ", remove = T) %>% 
  unite("Percentage of probing time spent in pathway-phase [%] **", medianse.percprobtime_C|group.percprobtime_C, sep = " ", remove = T) %>% 
  unite(medianse.percprobtime_G, median.percprobtime_G|se.percprobtime_G, sep = " ± ", remove = T) %>% 
  unite("Percentage of probing time spent in active ingestion [%] **", medianse.percprobtime_G|group.percprobtime_G, sep = " ", remove = T) %>% 
  unite(medianse.E2index, median.E2index|se.E2index, sep = " ± ", remove = T) %>% 
  unite("Potential E2 index [%] *", medianse.E2index|group.E2index, sep = " ", remove = T) %>% 
  unite(medianse.mean_fr_Ninterrup, median.mean_fr_Ninterrup|se.mean_fr_Ninterrup, sep = " ± ", remove = T) %>% 
  unite("Mean frequency of interruptions during phloem phase [mHz] **", medianse.mean_fr_Ninterrup|group.mean_fr_Ninterrup, sep = " ", remove = T) %>% 
  unite(medianse.percNinterrup_E2, median.percNinterrup_E2|se.percNinterrup_E2, sep = " ± ", remove = T) %>% 
  unite("Percentage of time spent in interruption during phloem phase [%] *", medianse.percNinterrup_E2|group.percNinterrup_E2, sep = " ", remove = T) %>% 
  unite(medianse.n_Ninterrupt_E2, median.n_Ninterrupt_E2|se.n_Ninterrupt_E2, sep = " ± ", remove = T) %>% 
  unite("Number of interruptions during phloem ingestion *", medianse.n_Ninterrupt_E2|group.n_Ninterrupt_E2, sep = " ", remove = T) %>% 
  unite(medianse.s_npto1stprobe, median.s_npto1stprobe|se.s_npto1stprobe, sep = " ± ", remove = T) %>% 
  unite("Time from 1st np to 1st probe [s] *", medianse.s_npto1stprobe|group.s_npto1stprobe, sep = " ", remove = T) %>% 
  unite(medianse.t_1st_sE2, median.t_1st_sE2|se.t_1st_sE2, sep = " ± ", remove = T) %>% 
  unite("Time of 1st sustained phloem phase [min] *", medianse.t_1st_sE2|group.t_1st_sE2, sep = " ", remove = T) %>% 
  dplyr::select(!c(emmean:response))

#add numbers of recs per cvs and sex
n_tab <- dub %>% 
  group_by(Cultivar, Sex) %>% 
  summarise(n = n())

median_gentab_n <- full_join(n_tab, median_gentab, by = c("Cultivar", "Sex"))

#Table 3
#transposed table (reviewer2 suggestion)
treno <- as_tibble(t(median_gentab_n), rownames = "row_names")                                               
write.csv(treno, file = "uniqueTab.csv", row.names = F)
                                                 
#SUPPLEMENTARY FILE S2----------
dacomp <- ynew2
library(emmeans)
library(stringr)
library(car)
library(hnp)

#n_np
#glm - quasipoisson distribution
gqcomp.n_np <- MASS::glm.nb(n_np ~ Cultivar, data = dacomp)
summary(gqcomp.n_np) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(n_np ~ Cultivar, family = quasipoisson, data = dacomp) #test homogeneity of variance
hnp(gqcomp.n_np, main = "gqcomp.n_np", how.many.out = T, sim = 999, conf = 0.99)
ref.gqcomp.n_np <- emmeans::ref_grid(gqcomp.n_np)
emm.gqcomp.n_np <- emmeans::emmeans(ref.gqcomp.n_np, specs = "Cultivar")
cld.gqcomp.n_np <- as_tibble(multcomp::cld(emm.gqcomp.n_np, Letters = "abcdef"))
cld.gqcomp.n_np <- rbind(cld.gqcomp.n_np, cld.gqcomp.n_np %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.n_np = tolower(.group))

#s_np
#glm - gamma distribution
gqcomp.s_np <- glm(s_np ~ Cultivar, family = gaussian, data = dacomp)
summary(gqcomp.s_np) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(s_np ~ Cultivar, family = Gamma, data = dacomp) #test homogeneity of variance
hnp(gqcomp.s_np, main = "gqcomp.s_np", how.many.out = T, sim = 999, conf = 0.99) #fits poorly
ref.gqcomp.s_np <- emmeans::ref_grid(gqcomp.s_np)
emm.gqcomp.s_np <- emmeans::emmeans(ref.gqcomp.s_np, specs = "Cultivar")
cld.gqcomp.s_np <- as_tibble(multcomp::cld(emm.gqcomp.s_np, Letters = "abcdef"))
cld.gqcomp.s_np <- rbind(cld.gqcomp.s_np, cld.gqcomp.s_np %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.s_np = tolower(.group))

#s_2np
#glm - inverse.gaussian distribution
gqcomp.s_2np <- glm(s_2np ~ Cultivar, family = inverse.gaussian, data = dacomp[-5,]) #family = Gamma out of goodness-of-fit, inverse-gaussian fit precisely
summary(gqcomp.s_2np) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(s_2np ~ Cultivar, family = inverse.gaussian, data = dacomp) #test homogeneity of variance
hnp(gqcomp.s_2np, main = "gqcomp.s_2np", how.many.out = T, sim = 999, conf = 0.95) #fits perfectly
ref.gqcomp.s_2np <- emmeans::ref_grid(gqcomp.s_2np)
emm.gqcomp.s_2np <- emmeans::emmeans(ref.gqcomp.s_2np, specs = "Cultivar")
cld.gqcomp.s_2np <- as_tibble(multcomp::cld(emm.gqcomp.s_2np, Letters = "abcdef"))
cld.gqcomp.s_2np <- rbind(cld.gqcomp.s_2np, cld.gqcomp.s_2np %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.s_2np = tolower(.group))

#n_Pr
#glm - quasipoisson distribution
gqcomp.n_Pr <- MASS::glm.nb(n_Pr ~ Cultivar, data = dacomp)
summary(gqcomp.n_Pr) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(n_Pr ~ Cultivar, family = quasipoisson, data = dacomp) #test homogeneity of variance
hnp(gqcomp.n_Pr, main = "gqcomp.n_Pr", how.many.out = T, sim = 999, conf = 0.99) #fits poorly
ref.gqcomp.n_Pr <- emmeans::ref_grid(gqcomp.n_Pr)
emm.gqcomp.n_Pr <- emmeans::emmeans(ref.gqcomp.n_Pr, specs = "Cultivar")
cld.gqcomp.n_Pr <- as_tibble(multcomp::cld(emm.gqcomp.n_Pr, Letters = "abcdef")) 
cld.gqcomp.n_Pr <- rbind(cld.gqcomp.n_Pr, cld.gqcomp.n_Pr %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.n_Pr = tolower(.group))

#s_Pr
#glm - gamma distribution
gqcomp.s_Pr <- glm(s_Pr ~ Cultivar, family = inverse.gaussian, data = dacomp)
summary(gqcomp.s_Pr) # NO differences for Sex NOR Cultivar*Sex, keep Cultivar only
leveneTest(s_Pr ~ Cultivar, family = Gamma, data = dacomp) #test homogeneity of variance
hnp(gqcomp.s_Pr, main = "gqcomp.s_Pr", how.many.out = T, sim = 999, conf = 0.95) #fits perfectly
ref.gqcomp.s_Pr <- emmeans::ref_grid(gqcomp.s_Pr)
emm.gqcomp.s_Pr <- emmeans::emmeans(ref.gqcomp.s_Pr, specs = "Cultivar")
cld.gqcomp.s_Pr <- as_tibble(multcomp::cld(emm.gqcomp.s_Pr, Letters = "abcdef"))
cld.gqcomp.s_Pr <- rbind(cld.gqcomp.s_Pr, cld.gqcomp.s_Pr %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.s_Pr = tolower(.group))

#s_C
#glm - gamma distribution
gqcomp.s_C <- glm(s_C ~ Cultivar, family = inverse.gaussian, data = dacomp)
summary(gqcomp.s_C) # keep Cultivar
leveneTest(s_C ~ Cultivar, family = Gamma, data = dacomp) #test homogeneity of variance
hnp(gqcomp.s_C, main = "gqcomp.s_C", how.many.out = T, sim = 999, conf = 0.95) #fits perfectly
ref.gqcomp.s_C <- emmeans::ref_grid(gqcomp.s_C)
emm.gqcomp.s_C <- emmeans::emmeans(ref.gqcomp.s_C, specs = c("Cultivar"))
cld.gqcomp.s_C <- as_tibble(multcomp::cld(emm.gqcomp.s_C, Letters = "abcdef"))
cld.gqcomp.s_C <- rbind(cld.gqcomp.s_C, cld.gqcomp.s_C %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.s_C = tolower(.group))

#n_G
#glm - negative-binomial distribution
gqcomp.n_G <- MASS::glm.nb(n_G ~ Cultivar, data = dacomp) #
#gqcomp.n_G <- glm(n_G ~ Cultivar*Sex, family = quasipoisson, data = dacomp)
summary(gqcomp.n_G) # keep Cultivar 
leveneTest(n_G ~ Cultivar*Sex, family = negative.binomial, data = dacomp) #test homogeneity of variance
hnp(gqcomp.n_G, main = "gqcomp.n_G", how.many.out = T, sim = 999, conf = 0.95) #fits perfectly
ref.gqcomp.n_G <- emmeans::ref_grid(gqcomp.n_G)
emm.gqcomp.n_G <- emmeans::emmeans(ref.gqcomp.n_G, specs = c("Cultivar"))
cld.gqcomp.n_G <- as_tibble(multcomp::cld(emm.gqcomp.n_G, Letters = "abcdef")) 
cld.gqcomp.n_G <- rbind(cld.gqcomp.n_G, cld.gqcomp.n_G %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.n_G = tolower(.group))

#s_G
#glm - gamma distribution
gqcomp.s_G <- glm(s_G ~ Cultivar, family = Gamma, data = dacomp)
summary(gqcomp.s_G) # No differences for Sex NOR interaction, keep cultivar
leveneTest(s_G ~ Cultivar, family = Gamma, data = dacomp) #test homogeneity of variance
hnp(gqcomp.s_G, main = "gqcomp.s_G", how.many.out = T, sim = 999, conf = 0.95) #fits perfectly
ref.gqcomp.s_G <- emmeans::ref_grid(gqcomp.s_G)
emm.gqcomp.s_G <- emmeans::emmeans(ref.gqcomp.s_G, specs = c("Cultivar"))
cld.gqcomp.s_G <- as_tibble(multcomp::cld(emm.gqcomp.s_G, Letters = "abcdef")) 
cld.gqcomp.s_G <- rbind(cld.gqcomp.s_G, cld.gqcomp.s_G %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.s_G = tolower(.group))

#s_npto1stprobe
#glm - inverse.gaussian distribution
gqcomp.s_npto1stprobe <- glm(s_npto1stprobe ~ Cultivar, family = inverse.gaussian, data = dacomp)
summary(gqcomp.s_npto1stprobe) # NO differences
leveneTest(s_npto1stprobe ~ Cultivar, family = inverse.gaussian, data = dacomp) #test homogeneity of variance
tmp <- hnp(gqcomp.s_npto1stprobe, main = "gqcomp.s_npto1stprobe", how.many.out = T, sim = 999, conf = 0.95, paint.out = T) # fits perfectly
ref.s_npto1stprobe <- emmeans::ref_grid(gqcomp.s_npto1stprobe)
emm.s_npto1stprobe <- emmeans::emmeans(ref.s_npto1stprobe, specs = c("Cultivar"))
cld.gqcomp.s_npto1stprobe <- as_tibble(multcomp::cld(emm.s_npto1stprobe, Letters = "abcdef")) 
cld.gqcomp.s_npto1stprobe <- rbind(cld.gqcomp.s_npto1stprobe, cld.gqcomp.s_npto1stprobe %>% mutate(Cultivar=Cultivar, .group=toupper(.group))) %>% 
  mutate(Sex = ifelse (.group == tolower(.group), "female", "male")) %>% 
  mutate(group.s_npto1stprobe = tolower(.group))

detach("package:hnp", unload = TRUE)
detach("package:MASS", unload = TRUE)

#Table S2, dataset calculations
library(dplyr)

#table all recordings
dubcomp <- ynew2

#new tibble: column "median.n_np"
acomp <-
  dubcomp %>% 
  group_by(Cultivar, Sex, n_np) %>% 
  summarise(n = n()) %>% 
  mutate(median.n_np = median(n_np)) %>% 
  mutate(se.n_np = sd(n_np)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.n_np, se.n_np)

#new tibble: column "median.s_np"
bcomp <-
  dubcomp %>% 
  group_by(Cultivar, Sex, s_np) %>% 
  summarise(n = n()) %>% 
  mutate(median.s_np = median(s_np)/60) %>% 
  mutate(se.s_np = sd(s_np)/sqrt(sum(n))/60) %>%
  distinct(Cultivar, Sex, median.s_np, se.s_np)

#new tibble: column "median.s_2np"
ccomp <-
  dubcomp %>% 
  group_by(Cultivar, Sex, s_2np) %>% 
  summarise(n = n()) %>% 
  mutate(median.s_2np = median(s_2np))%>% 
  mutate(se.s_2np = sd(s_2np)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.s_2np, se.s_2np)

#new tibble: column "median.n_Pr"
dcomp <-
  dubcomp %>% 
  group_by(Cultivar, Sex, n_Pr) %>% 
  summarise(n = n()) %>% 
  mutate(median.n_Pr = median(n_Pr)) %>% 
  mutate(se.n_Pr = sd(n_Pr)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.n_Pr, se.n_Pr)

#new tibble: column "median.s_Pr"
ecomp <-
  dubcomp %>% 
  group_by(Cultivar, Sex, s_Pr) %>% 
  summarise(n = n()) %>% 
  mutate(median.s_Pr = median(s_Pr)/60) %>% 
  mutate(se.s_Pr = sd(s_Pr)/sqrt(sum(n))/60) %>%
  distinct(Cultivar, Sex, median.s_Pr, se.s_Pr)

#new tibble: column "median.s_C"
fcomp <-
  dubcomp %>% 
  group_by(Cultivar, Sex, s_C) %>% 
  summarise(n = n()) %>% 
  mutate(median.s_C = median(s_C)/60) %>% 
  mutate(se.s_C = sd(s_C)/sqrt(sum(n))/60) %>%
  distinct(Cultivar, Sex, median.s_C, se.s_C)

#new tibble: column "median.n_G"
gcomp <-
  dubcomp %>% 
  group_by(Cultivar, Sex, n_G) %>% 
  summarise(n = n()) %>% 
  mutate(median.n_G = median(n_G)) %>% 
  mutate(se.n_G = sd(n_G)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.n_G, se.n_G)

#new tibble: column "median.s_G"
hcomp <-
  dubcomp %>% 
  group_by(Cultivar, Sex, s_G) %>% 
  summarise(n = n()) %>% 
  mutate(median.s_G = median(s_G)/60) %>% 
  mutate(se.s_G = sd(s_G)/sqrt(sum(n))/60) %>%
  distinct(Cultivar, Sex, median.s_G, se.s_G)

#new tibble: column "median.s_npto1stprobe"
atcomp <-
  dubcomp %>% 
  group_by(Cultivar, Sex, s_npto1stprobe) %>% 
  summarise(n = n()) %>% 
  mutate(median.s_npto1stprobe = median(s_npto1stprobe)) %>% 
  mutate(se.s_npto1stprobe = sd(s_npto1stprobe)/sqrt(sum(n))) %>%
  distinct(Cultivar, Sex, median.s_npto1stprobe, se.s_npto1stprobe)

#join datasets
dfrcomp <- list(acomp,bcomp,ccomp,dcomp,ecomp,fcomp,gcomp,hcomp,atcomp,
            cld.gqcomp.n_np, cld.gqcomp.s_np, cld.gqcomp.s_2np, cld.gqcomp.n_Pr, cld.gqcomp.s_Pr, cld.gqcomp.s_C, cld.gqcomp.n_G,
            cld.gqcomp.s_G, cld.gqcomp.s_npto1stprobe)

#intermediate table
median_gentabcomp <- plyr::join_all(dfrcomp, by = c("Cultivar", "Sex"), type='full') %>%
  mutate(across(where(is.numeric), round, 1))

library(tidyr)

#general table
#unite median +- SE
median_gentabcomp <- median_gentabcomp %>% 
  unite(medianse.n_np, median.n_np|se.n_np, sep = " ? ", remove = T) %>%
  unite("Number of non-probing periods *", medianse.n_np|group.n_np, sep = " ", remove = T) %>% 
  unite(medianse.s_np, median.s_np|se.s_np, sep = " ? ", remove = T) %>%
  unite("Total duration of non-probing periods [min] *", medianse.s_np|group.s_np, sep = " ", remove = T) %>%
  unite(medianse.s_2np, median.s_2np|se.s_2np, sep = " ? ", remove = T) %>%
  unite("Duration of the 2nd non-probing period [s] *", medianse.s_2np|group.s_2np, sep = " ", remove = T) %>%
  unite(medianse.n_Pr, median.n_Pr|se.n_Pr, sep = " ? ", remove = T) %>% 
  unite("Number of probes *", medianse.n_Pr|group.n_Pr, sep = " ", remove = T) %>% 
  unite(medianse.s_Pr, median.s_Pr|se.s_Pr, sep = " ? ", remove = T) %>% 
  unite("Total probing time [min] *", medianse.s_Pr|group.s_Pr, sep = " ", remove = T) %>% 
  unite(medianse.s_C, median.s_C|se.s_C, sep = " ? ", remove = T) %>% 
  unite("Total duration of pathway phase [min] *", medianse.s_C|group.s_C, sep = " ", remove = T) %>% 
  unite(medianse.n_G, median.n_G|se.n_G, sep = " ? ", remove = T) %>% 
  unite("Number of active ingestion phases *", medianse.n_G|group.n_G, sep = " ", remove = T) %>% 
  unite(medianse.s_G, median.s_G|se.s_G, sep = " ? ", remove = T) %>% 
  unite("Total duration of active ingestion [min] *", medianse.s_G|group.s_G, sep = " ", remove = T) %>% 
  unite(medianse.s_npto1stprobe, median.s_npto1stprobe|se.s_npto1stprobe, sep = " ? ", remove = T) %>% 
  unite("Time from 1st np to 1st probe [s] *", medianse.s_npto1stprobe|group.s_npto1stprobe, sep = " ", remove = T) %>% 
  dplyr::select(!c(emmean:.group))

#add numbers of recs per cvs and sex
n_tabcomp <- dubcomp %>% 
  group_by(Cultivar, Sex) %>% 
  summarise(n = n())

#complete table
median_gentab_ncomp <- full_join(n_tabcomp, median_gentabcomp, by = c("Cultivar", "Sex"))

#Table S2
#table: non phloem variables
write.csv(median_gentab_n %>% 
            dplyr::select(Cultivar, Sex, n, c(1:12)), file = "TableS2.csv", row.names = F)


###MULTIVARIATE ANALYSIS---------
#FIGURE 2------------
#CCA

#add cultivar-related median to specifics columns
dab <- ynew2_phloem

dab$E2index[is.na(dab$E2index)] <- with(dab, ave(E2index, Cultivar,
                                                 FUN = function(x) median(x, na.rm = TRUE)))[is.na(dab$E2index)]

#exclude collinear variables
vif <- usdm::vifcor(as.data.frame(dab[,3:27]), th=0.95)
dabex <- usdm::exclude(dab, vif)

dabex <- dplyr::bind_cols(dabex, dab %>% dplyr::select(File, Cultivar, Sex))

#standardise dataset values
library(vegan)
library(ggordiplots)
ccazz_stand <- decostand(dabex[,1:20], na.rm = T, #always control column numbers
                         method = "hellinger")

#CCA
ccazz <- cca(ccazz_stand ~ Cultivar*Sex, na.action = na.exclude,
             data = dabex)

#raw plot
ccazz_ggordi <- gg_ordiplot(ccazz, groups = dabex$Cultivar, plot = T, kind = "se", conf = 0.99)

#define observations
ccazz_scores <- as.data.frame(scores(ccazz, choices = c(1,2), display = "species"))

library(tibble)
ccazz_scores <- rownames_to_column(ccazz_scores, 
                                   var = "species")

#define environmental variables vectors
ccazz_envi_scores <- envfit(ccazz ~ Cultivar*Sex, 
                            data = dabex)
ccazz_envi_scores_vector <- as.data.frame(ccazz_envi_scores$factors$centroids)
ccazz_envi_scores_vector <- rownames_to_column(ccazz_envi_scores_vector, 
                                               var = "envi")

#define groups
ord.ccazz <- ccazz_ggordi$df_ord
ord.ccazz$Sex <- dabex$Sex
colnames(ord.ccazz) <- c("x", "y", "Cultivar", "Sex")

#Figure 2
library(ggrepel)
ccazz_plot <- ggplot() + 
  geom_point(data = ord.ccazz, 
             aes(x = x, y = y, colour = Cultivar, shape = Sex), 
             size = 3.5, alpha = 0.9, stroke = 0) + 
  geom_point(data = ccazz_ggordi$df_spiders, 
             aes(x = cntr.x, y = cntr.y, colour = Group), 
             size = 4.5, show.legend = FALSE) + 
  geom_path(data = ccazz_ggordi$df_ellipse,   
            aes(x = x, y = y, colour = Group),  
            size = 1.1) + 
  geom_segment(data = ccazz_scores,
               aes(x=0, y=0, xend = CCA1, yend = CCA2), alpha = .2)+
  geom_text_repel(data = ccazz_scores,
                  aes(x = CCA1, y = CCA2, 
                      label = species), arrow = arrow (angle = 10, length=unit(4, "mm"), ends = "last", type = "closed"),
                  alpha = .6) +
  labs(x = "CCA 1 (28.7%)",        y = "CCA 2 (1.9%)", #control gg_ordiplot every time! ccazz_ggordi
       colour = "Cultivar") + 
  scale_colour_manual(values = c("#5e3c99", "#fdb863", "#e66101"))+
  theme_minimal()+
  theme(panel.grid = element_blank(), panel.border = element_rect(fill="transparent"))+
  geom_hline(yintercept = 0, lty = 2, alpha=0.3) +
  geom_vline(xintercept = 0, lty = 2, alpha=0.3)+
  guides(color = guide_legend(title = "Cultivar"))

ccazz_plot

ggsave(filename = "Figure2.tiff", plot = ccazz_plot, dpi = 320, compression = "lzw", width = 10, height = 8)

#TABLE 7-------
#perMANOVA
permanova <- adonis2(ccazz_stand ~ Cultivar*Sex, data = dabex, method='bray', permutations = 9999)

write.csv(permanova, file = "Table7.csv", row.names = T)
