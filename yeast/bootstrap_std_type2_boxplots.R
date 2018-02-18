rm(list = ls())

pacman::p_load(tidyverse, forcats)

####get the boxplots
rank_df <- read_csv('bootstrap_std_type2_ranks.csv') %>% 
  gather(Method, Estimates) %>% 
  mutate(type = 'Rank')

card_df <- read_csv('bootstrap_std_type2_cards.csv') %>% 
  gather(Method, Estimates) %>% 
  mutate(type = 'Number of selected predictors')

comb_df <- bind_rows(rank_df, card_df) %>%
  filter(Method != 'FSCV') %>%
  mutate(type = fct_relevel(type, 'Rank', 'Number of selected predictors'),
         Method = fct_relevel(Method, 'PSCV', 'CV')) %>%
  mutate(Method = factor(Method, labels = c('SCV', 'CV')))


ggplot(comb_df, aes(x = Method, y = Estimates)) +
  geom_boxplot(fatten = 3.5) +
  facet_wrap(~ type, scales = 'free_y') +
  labs(x = '', y = '') +
  theme_bw()
####

####get the median DFs
rank_df2 <- read_csv('bootstrap_std_type2_ranks.csv') %>% 
  gather(Method, Estimates) %>% 
  dplyr::rename(rank = Estimates)

card_df2 <- read_csv('bootstrap_std_type2_cards.csv') %>% 
  gather(Method, Estimates) %>% 
  dplyr::rename(card = Estimates)

comb_df2 <- bind_cols(rank_df2, card_df2) %>%
  dplyr::select(-Method1) %>% 
  group_by(Method) %>% 
  mutate(DF = (card + 18 - rank) * rank) %>% 
  summarise(median = median(DF))

comb_df2
####
