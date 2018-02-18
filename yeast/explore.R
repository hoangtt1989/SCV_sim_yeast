rm(list = ls())

pacman::p_load(tidyverse, spls, stringr)


data(yeast)


yeast$y
yeast$x


write_csv(as.data.frame(yeast$x), 'X_mat.csv', col_names = F)
write_csv(as.data.frame(yeast$y), 'Y_mat.csv', col_names = F)


###the 21 experimentally verified TFs
exp_TF <- c('ACE2', 'SWI4', 'SWI5', 'SWI6', 'MBP1', 'STB1', 'FKH1', 'FKH2',
            'NDD1', 'MCM1', 'ABF1', 'BAS1', 'CBF1', 'GCN4', 'GCR1', 'GCR2',
            'LEU3', 'MET31', 'REB1', 'SKN7', 'STE12')

simp_names <- str_split(colnames(yeast$x), '_', simplify = T)[, 1]

X_names <- yeast$x
colnames(X_names) <- simp_names

not_exp <- setdiff(colnames(X_names), exp_TF)

X_reorder <- X_names[, c(exp_TF, not_exp)]

TF_df <- data_frame(TF = str_split(colnames(yeast$x), '_', simplify = T)[, 1], pos = 1:length(colnames(yeast$x))) %>% 
  mutate(in_21 = TF %in% exp_TF) %>% 
  mutate(in_co_pairs = TF %in% unique(co_pairs))


saveRDS(TF_df, 'TF_names.rds')
