##############################
## surivorship analysis
##
## Matt Brachmann (PhDMattyB)
##
## 16.12.2024
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')

library(tidyverse)
library(ggridges)
library(viridis)
library(patchwork)

theme_set(theme_bw())

data = read_csv('convertedeusurvdata.csv') %>% 
  na.omit()%>% 
  mutate(group_treatment = paste(group, 
                                 treatment, sep = "_"))
data$treatment = as.factor(data$treatment)


data_clean = data %>% 
  mutate(.data = data,
         group2 = as.factor(case_when(
           group == 'AshnW' ~ 'ASHN - Warm',
           group == 'AshnC' ~ 'ASHN - Cold',
           group == 'SKRW' ~ 'SKR - Warm', 
           group == 'SKRC' ~ 'SKR - Cold', 
           group == 'MyvatW' ~ 'MYV - Warm', 
           group == 'MyvatC' ~ 'MYV - Cold', 
           group == 'GTS' ~ 'GTS - Warm', 
           group == 'CAUSEWAY' ~ 'GAR - Cold'))) %>% 
  mutate(.data = data, 
         treatment2 = as.factor(case_when(
           treatment == '1212' ~ '12 @ 12', 
           treatment == '1218' ~ '12 @ 18', 
           treatment == '1812' ~ '18 @ 12', 
           treatment == '1818' ~ '18 @ 18'
         )))

data_clean$group2 = factor(data_clean$group2, 
                           levels = c('ASHN - Warm', 
                                      'ASHN - Cold', 
                                      'MYV - Warm', 
                                      'MYV - Cold', 
                                      'SKR - Warm', 
                                      'SKR - Cold', 
                                      'GTS - Warm', 
                                      'GAR - Cold'))
  

# whole experiment survival -----------------------------------------------


# whole_exp_aov = aov(adjustedSurvWholeExp ~ group_treatment, 
#                     data = data)
# summary(whole_exp_aov)

# whole_exp_aov2 = aov(adjustedSurvWholeExp ~ treatment*ecotype*pair, 
#                     data = data)
# summary(whole_exp_aov2)

# Anova(whole_exp_aov)
# Anova(whole_exp_aov2)

# TukeyHSD(whole_exp_aov)


  survive_whole_exp_plot = data_clean %>% 
    ggplot(aes(x = adjustedSurvWholeExp, 
             y = group2, 
             fill = group2))+
  facet_grid(~treatment2)+
  geom_density_ridges()+
  scale_fill_viridis(discrete = T)+
  labs(x = 'Survival over whole experiment', 
       title = 'Treatment')+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 14), 
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(size = 12, 
                                  face = 'bold'),
        plot.title = element_text(hjust = 0.5, 
                                  face = 'bold', 
                                  size = 14),
        legend.position = 'none')


  
  
  
  
  ggsave('~/Parsons_Postdoc/Written_things/Stickle_genomic_paper/Fat_Phenotype.tiff', 
         plot = lipid_plot, 
         dpi = 'retina', 
         units = 'cm', 
         width = 15, 
         height = 10)
  
  
# 12@12 treatments --------------------------------------------------------

T1212 = data %>% 
  filter(treatment == '1212')

whole_exp_1212_aov = aov(adjustedSurvWholeExp ~ group, 
                         data = T1212)

summary(whole_exp_1212_aov)

TukeyHSD(whole_exp_1212_aov)


# 12@18 treatments --------------------------------------------------------

T1218 = data %>% 
  filter(treatment == '1218')

whole_exp_1218_aov = aov(adjustedSurvWholeExp ~ group, 
                         data = T1218)

summary(whole_exp_1218_aov)


# 18@12 treatments --------------------------------------------------------

T1812 = data %>% 
  filter(treatment == '1812')

whole_exp_1812_aov = aov(adjustedSurvWholeExp ~ group, 
                         data = T1812)

summary(whole_exp_1812_aov)


# 18@18 treatments --------------------------------------------------------

T1818 = data %>% 
  filter(treatment == '1818')

whole_exp_1818_aov = aov(adjustedSurvWholeExp ~ group, 
                         data = T1818)

summary(whole_exp_1818_aov)

TukeyHSD(whole_exp_1818_aov)



# EU survival -------------------------------------------------------------

EU_data = data %>% 
  rename(EU_survival = longdata...8)


eu_survive_aov = aov(Eusurvivalafter6weeks ~ EU_survival*treatment,
                    data = EU_data)
summary(eu_survive_aov)

car::Anova(eu_survive_aov, 
           type = 'II')

TukeyHSD(eu_survive_aov)



# 6 week survival ---------------------------------------------------------


