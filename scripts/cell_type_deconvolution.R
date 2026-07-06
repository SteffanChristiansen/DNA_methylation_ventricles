# loading packages --------------------------------------------------------
packages <- c("tidyverse", "here", "EpiSCORE")
invisible(lapply(packages, require, character.only = TRUE)) 

# loading data ------------------------------------------------------------
beta_values <- readRDS("results/processed/beta_values_combat_2024_06_27.rds")
sample_sheet <- readRDS("data/sample_sheet_2024_10_23.rds")
selected_colors <- readRDS("data/selected_colors_2023_12_19.rds")
sample_sheet <- sample_sheet %>% 
  mutate(id = factor(id))

# running episcore --------------------------------------------------------
promotor_data_matrix_EPICv2 <- constAvBetaTSS(
  beta.m = beta_values,
  type = "EPICv2")

est_ventricles <- wRPC(
  promotor_data_matrix_EPICv2,
  ref=biosino_heart_ref_mat,
  useW=TRUE,
  wth=0.4,
  maxit=200 
)

est_ventricles_df <- data.frame(est_ventricles$estF) %>% 
  rownames_to_column(., var = "id_pos") %>% 
  left_join(., sample_sheet)

est_ventricles_collapsed <- est_ventricles_df %>% 
  mutate(CM_SMC = CM + SMC) %>% 
  select(id, origin, EC, FB, MP, CM_SMC) 

# difference in CM and SMC cell type composition --------------------------
model <- aov(CM_SMC ~ origin + Error(id/origin), data = est_ventricles_collapsed)
model
mean(est_ventricles_collapsed$CM_SMC) # Equal to "Grand Mean"

est_ventricles_collapsed %>% 
  group_by(origin) %>% 
  reframe(mean = mean(CM_SMC))

summary(model)