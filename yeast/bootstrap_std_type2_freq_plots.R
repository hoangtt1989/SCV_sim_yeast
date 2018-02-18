rm(list = ls())

pacman::p_load(tidyverse, forcats, ggrepel, readr, readxl, gridExtra, ggpubr, directlabels)


TF_df <- readRDS('TF_names.rds') %>% 
  dplyr::rename(idx = pos)

sheet_names <- c('PSCV', 'FSCV', 'CV')

brks <- seq(1, 0, by = -.01)
#######TF21 frequencies
file_name <- 'bootstrap_std_type2_TF21_freqs.xls'

TF21_freq_df <- data_frame(Method = sheet_names) %>% 
  mutate(data = map(Method, function(x) {read_xls(file_name, x) %>% 
      mutate(TF = TF_df$TF[idx])}))

TF21_freq_curves <- TF21_freq_df %>% 
  mutate(curve_data = map(data, function(x) {num_TFs <- map_dbl(brks, function(y) {sum(x$freq / 200 >= y)})
  pct <- brks
  data_frame(num_TFs = num_TFs, pct = pct)}))

TF21_curves_df <- TF21_freq_curves %>% 
  unnest(curve_data) %>% 
  # filter(Method != 'PSCV') %>%
  filter(Method != 'FSCV') %>%
  # mutate(Method = fct_relevel(Method, 'FSCV', 'CV')) %>% 
  mutate(Method = fct_relevel(Method, 'PSCV', 'CV')) %>%
  mutate(Method = factor(Method, labels = c('SCV', 'CV'))) %>% 
  mutate(label_pos = pct == .75) %>% 
  # mutate(label_pos2 = case_when(
  #   label_pos & Method == 'FSCV' ~ 'FSCV',
  #   label_pos & Method == 'CV' ~ 'CV',
  #   !label_pos ~ ''
  # ))
mutate(label_pos2 = case_when(
  label_pos & Method == 'SCV' ~ 'SCV',
  label_pos & Method == 'CV' ~ 'CV',
  !label_pos ~ ''
))



TF21_freq_df2 <- TF21_freq_df %>% 
  unnest(data) %>% 
  mutate(pct = freq / 200 * 100) %>% 
  filter(Method != 'FSCV')
  # filter(Method != 'PSCV')
  
TF21_freq_df3 <- TF21_freq_df2 %>% 
  dplyr::select(-idx, -freq) %>% 
  spread(Method, pct)

freq_brks <- seq(0, 100, by = 25)

leave_text <- c('MCM1', 'FKH2', 'MBP1', 'ACE2', 'SWI6', 'STE12', 'SWI4', 'NDD1', 'SWI5')

TF21_freq_df4 <- TF21_freq_df3 %>% 
  mutate(leave_text = !(TF %in% leave_text)) %>% 
  mutate()

TF21_ident <- ggplot(TF21_freq_df4, aes(x = CV, y = PSCV, label = TF, alpha = leave_text)) +
  geom_point(size = 1) +
  geom_text_repel(size = 2.5, segment.size = .25) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_alpha_discrete(name = '', labels = '', breaks = '', range = c(.6, 1)) +
  scale_x_continuous(breaks = freq_brks, limits = c(0, 100)) +
  scale_y_continuous(breaks = freq_brks, limits = c(0, 100)) +
  labs(y = 'SCV') +
  theme_bw()


m_num <- .15

TF21_curves <- ggplot(TF21_curves_df %>% filter(pct >= .5), aes(x = pct * 100, y = num_TFs / 21 * 100, linetype = Method)) +
  geom_line() +
  geom_text(aes(x = pct * 100, y = num_TFs / 21 * 100, label = label_pos2), nudge_y = 2) +
  scale_x_reverse() +
  guides(linetype = F) +
  labs(x = 'Cut-off frequency (%)', y = 'Percent of TFs') +
  theme_bw() +
  theme(plot.margin = margin(m_num, m_num, m_num, 1.3, 'cm'))

ggarrange(TF21_ident, TF21_curves, ncol = 2)


###############colored version

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

TF21_ident <- ggplot(TF21_freq_df4, aes(x = CV, y = PSCV, label = TF, alpha = leave_text)) +
  # geom_point(size = 1, aes(color = leave_text)) +
  geom_point(size = 1, color = 'red1') +
  geom_text_repel(size = 2.5, segment.size = .25, color = 'forestgreen', force = 1) +
  # geom_text_repel(size = 2.5, segment.size = .25, aes(color = leave_text)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = cbPalette[3]) +
  scale_alpha_manual(values = c(.6, 1)) +
  # scale_fill_manual(values = )
  # scale_color_manual(values = c('forestgreen', 'red1')) +
  # scale_color_manual(values = c(cbPalette[4], cbPalette[8])) +
  scale_x_continuous(breaks = freq_brks, limits = c(0, 100)) +
  scale_y_continuous(breaks = freq_brks, limits = c(0, 100)) +
  guides(color = F, alpha = F) +
  labs(y = 'SCV') +
  theme_bw()

TF21_ident
# grid.arrange(TF21_ident, TF21_curves, ncol = 2)

m_num <- .15

TF21_curves <- ggplot(TF21_curves_df %>% filter(pct >= .5), aes(x = pct * 100, y = num_TFs / 21 * 100, color = Method)) +
  geom_line() +
  geom_text(aes(x = pct * 100, y = num_TFs / 21 * 100, label = label_pos2), nudge_y = 2) +
  scale_x_reverse() +
  guides(linetype = F, color = F) +
  labs(x = 'Cut-off frequency (%)', y = 'Percent of TFs') +
  theme_bw() +
  theme(plot.margin = margin(m_num, m_num, m_num, 1.3, 'cm'))

ggarrange(TF21_ident, TF21_curves, ncol = 2)
#################
