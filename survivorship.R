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

theme_set(theme_bw())

data = read_csv('convertedeusurvdata.csv') %>% 
  na.omit()%>% 
  mutate(group_treatment = paste(group, 
                                 treatment, sep = "_"))

data %>% 
  ggplot(aes(x = Eusurvivalafter6weeks))+
  geom_histogram()

data %>% 
  ggplot(aes(x = early6weeksurvivalrate))+
  geom_histogram()

data %>% 
  ggplot(aes(x = adjustedSurvWholeExp))+
  geom_histogram()

data$treatment = as.factor(data$treatment)

whole_exp_aov = aov(adjustedSurvWholeExp ~ group_treatment, 
                    data = data)
summary(whole_exp_aov)

# whole_exp_aov2 = aov(adjustedSurvWholeExp ~ treatment*ecotype*pair, 
#                     data = data)
# summary(whole_exp_aov2)

Anova(whole_exp_aov)
# Anova(whole_exp_aov2)

TukeyHSD(whole_exp_aov)


data %>% 
  ggplot(aes(x = adjustedSurvWholeExp, 
             y = group, 
             fill = group))+
  facet_grid(~treatment)+
  geom_density_ridges()
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


