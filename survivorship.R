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
  


# log rank survival -------------------------------------------------------

fit <- survfit(Surv(group, treatment) ~ early6weeksurvivalrate, 
               data = data_clean)
print(fit)
summary(fit)$table

surv_diff <- survdiff(Surv(early6weeksurvivalrate) ~ group + treatment, 
                      data = data_clean)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw())

# Family details ----------------------------------------------------------

data_clean %>% 
  group_by(group_treatment) %>% View()


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


surv_diff <- survdiff(Surv(adjustedSurvWholeExp) ~  group2+treatment, 
                         data = data_clean)


  survive_whole_exp_plot = data_clean %>% 
    ggplot(aes(x = adjustedSurvWholeExp, 
             y = group2, 
             fill = group2))+
  facet_grid(~treatment2)+
  geom_density_ridges()+
  scale_fill_viridis(discrete = T)+
  labs(x = 'Survival over experiment')+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 14), 
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(size = 12, 
                                  face = 'bold'),
        # plot.title = element_text(hjust = 0.5, 
        #                           face = 'bold', 
        #                           size = 14),
        legend.position = 'none')


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

EU_data = data_clean %>% 
  rename(EU_survival = longdata...8) %>% 
  mutate(.data = data, 
         EU_survival_grouping = as.factor(case_when(
           EU_survival == 'eu1surv' ~ 'Experimental unit 1', 
           EU_survival == 'eu2surv' ~ 'Experimental unit 2', 
           EU_survival == 'eu3surv' ~ 'Experimental unit 3', 
           EU_survival == 'eu4surv' ~ 'Experimental unit 4', 
           EU_survival == 'eu5surv' ~ 'Experimental unit 5'
         )))


EU_surv_diff <- survdiff(Surv(Eusurvivalafter6weeks) ~ EU_survival_grouping + treatment, 
                      data = EU_data)


eu_survive_aov = aov(Eusurvivalafter6weeks ~ EU_survival_grouping*treatment,
                    data = EU_data)
summary(eu_survive_aov)

car::Anova(eu_survive_aov, 
           type = 'II')

TukeyHSD(eu_survive_aov)


EU_survival_plot = EU_data %>% 
  ggplot(aes(x = Eusurvivalafter6weeks, 
             y = EU_survival_grouping, 
             fill = EU_survival_grouping))+
  facet_grid(~treatment2)+
  geom_density_ridges()+
  scale_fill_viridis(discrete = T)+
  labs(x = 'Survival per experimental unit')+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 14), 
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(size = 12, 
                                  face = 'bold'),
        legend.position = 'none')



# Combine and save the plot -----------------------------------------------
combo_survial_plot = EU_survival_plot/survive_whole_exp_plot


ggsave('Survival_plots.tiff', 
       plot = combo_survial_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 15)



# six week survival -------------------------------------------------------


eu_survive_aov = aov(Eusurvivalafter6weeks ~ EU_survival_grouping*treatment,
                     data = EU_data)
summary(eu_survive_aov)

car::Anova(eu_survive_aov, 
           type = 'II')

TukeyHSD(eu_survive_aov)


EU_survival_plot = EU_data %>% 
  ggplot(aes(x = Eusurvivalafter6weeks, 
             y = EU_survival_grouping, 
             fill = EU_survival_grouping))+
  facet_grid(~treatment2)+
  geom_density_ridges()+
  scale_fill_viridis(discrete = T)+
  labs(x = 'Survival per experimental unit')+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 14), 
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(size = 12, 
                                  face = 'bold'),
        legend.position = 'none')




# Investigate GTS shit survival -------------------------------------------

EU_data %>% 
  filter(group == 'GTS') %>% 
  ggplot(aes(x = adjustedSurvWholeExp, 
             # y = treatment2, 
             fill = EU_survival_grouping))+
  geom_histogram()+
  facet_grid(~treatment2)+
  # geom_density_ridges()+
  scale_fill_viridis(discrete = T)+
  labs(x = 'Survival per experimental unit')+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 14), 
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(size = 12, 
                                  face = 'bold'))

