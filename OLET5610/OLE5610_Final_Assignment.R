# Load libraries
library(tidyverse)
library(purrr)
library(theft)
library(factoextra)
library(FactoMineR)
library(MVN)
library(rstatix)
library(tidymodels)
library(randomForest)
library(vip)
library(patchwork)
library(cowplot)
theme_set(theme_cowplot())

################################################################################
# Data prep
################################################################################

# Define data path
data_dir <- "D:/Virtual_Machines/Shared_Folder/github/PhD_courses/OLET5610/data/"

# clinical info data
clinical_info_data <- read.csv(paste0(data_dir, "ADNIMERGE.csv")) %>%
  dplyr::select(RID, DX_bl)

# tau PET data
tau_PET_data_wide <- read.csv(paste0(data_dir, "UCBERKELEYAV1451_PVC_04_29_22.csv")) %>%
  mutate(EXAMDATE = as.Date(EXAMDATE)) %>%
  group_by(RID) %>%
  filter(EXAMDATE == min(EXAMDATE, na.rm=T)) %>%
  ungroup() %>%
  left_join(., clinical_info_data) %>%
  filter(!is.na(DX_bl), DX_bl != "") %>%
  dplyr::select(RID, LEFT_MIDDLEFR_SUVR:CTX_RH_TRANSVERSETEMPORAL_VOLUME, DX_bl) %>%
  pivot_longer(cols = c(LEFT_MIDDLEFR_SUVR:CTX_RH_TRANSVERSETEMPORAL_VOLUME),
               names_to = "Variable",
               values_to = "Values") %>%
  mutate(Measurement = ifelse(str_detect(Variable, "SUVR"),
                              "SUVR", "Volume")) %>%
  mutate(Brain_Region = str_replace_all(Variable, "_SUVR|_VOLUME", ""),
         .keep = "unused") %>%
  filter(Measurement == "SUVR") %>%
  mutate(Brain_Region = stringr::str_to_lower(Brain_Region)) %>%
  distinct() %>%
  filter(Values > 0) %>%
  pivot_wider(id_cols=c(RID, DX_bl), names_from=Brain_Region, values_from=Values) %>%
  mutate(DX_bl = factor(DX_bl, levels = c("CN", "SMC", "EMCI", "LMCI", "AD")))

# Visualize distribution of classes
tau_PET_data_wide %>%
  group_by(DX_bl) %>%
  count() %>%
  ungroup() %>%
  ggplot(data = ., mapping = aes(x=DX_bl, y=n)) +
  geom_bar(stat="identity", aes(fill=DX_bl)) +
  geom_text(aes(label = n),vjust = -0.3) +
  scale_fill_manual(values = viridis::viridis(5, direction=1)) +
  ylab("# Participants") +
  xlab("Baseline Cognitive Diagnosis") +
  ggtitle("Number of Participants per Diagnostic Group") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5))
ggsave("plots/Figure1_Number_Participants_per_DX.png",
       width=6, height=4, units="in", dpi=300)


# Scale data
zscore <- function(values) {
  zscored_values <- (values - mean(values, na.rm=T))/sd(values, na.rm=T)
  return(zscored_values)
}

# Create numerical matrix
tau_num_mat <- tau_PET_data_wide %>%
  dplyr::select(-RID, -DX_bl) %>%
  as.matrix()

# Apply z-scoring
tau_num_mat_z <- apply(tau_num_mat, 2, zscore)

# Check for multivariate normality assumptions
mvn(tau_num_mat_z, mvnTest = "mardia")

# Check for homogeneity of variance-covariance matrices
box_m(tau_num_mat_z, as.character(tau_PET_data_wide$DX_bl))


################################################################################
# Method 1: RF
################################################################################

# Prep data
tau_num_z_with_dx <- as.data.frame(cbind(DX_bl = as.character(tau_PET_data_wide$DX_bl),
                                             tau_num_mat_z)) %>%
  mutate_at(vars(-("DX_bl")), as.numeric)

# Split data into 80% training, 20% testing
set.seed(127)
split_tau_PET_data <- initial_split(tau_num_z_with_dx, 
                                    prop = 0.8)

# Preprocess training data
recipe_tau_PET <-
  training(split_tau_PET_data) %>% # Use the training data
  recipe(DX_bl ~ .) %>% # Define the RF formula
  step_corr(all_predictors()) %>% # Remove correlated predictors > 0.9
  prep() # Preprocess

# Extract preprocessed training data
train_tau_PET <- juice(recipe_tau_PET)

# Apply defined preprocessing to test data
test_tau_PET <-
  recipe_tau_PET %>%
  bake(testing(split_tau_PET_data))

rforest_tau_PET <-
  rand_forest(trees = 100, mode = "classification") %>%
  set_engine("randomForest") %>%
  fit(factor(DX_bl) ~ ., data = train_tau_PET)

# Apply RF classifier to test data
tau_PET_predicted <- predict(rforest_tau_PET, test_tau_PET) 

# Model validation
validate <-
  test_tau_PET %>%
  bind_cols(tau_PET_predicted) %>% # add tau_PET_predicted into test data frame
  metrics(truth = DX_bl, # then perform comparison of original data in truth
          estimate = .pred_class) # against predicted data from model
validate

vi_model(rforest_tau_PET) %>%
  arrange(desc(Importance)) %>%
  slice_head(n=10) %>%
  mutate(Variable = fct_reorder(Variable, Importance, .desc=F)) %>%
  ggplot(data = ., mapping = aes(x = Variable, 
                                 y = Importance,
                                 fill = Importance)) +
  geom_bar(stat = "identity") +
  ggtitle("Variable Importance in\nRF Approach #1") +
  xlab("Brain Region") +
  scale_fill_viridis_b() +
  coord_flip() +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5))

ggsave("plots/Figure3_Approach1_VIP.png",
       width=5, height=6, units="in", dpi=300)

################################################################################
# Method 2: PCA + RF
################################################################################

# Correlation plot
corr_mat <- tau_PET_data_wide %>%
  dplyr::select(-RID, -DX_bl) %>%
  cor(., method = "pearson", use = "complete.obs")

ggplot(reshape2::melt(corr_mat), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8) +
  scale_fill_gradient2(low="blue", mid="white", high="red") +
  theme_minimal() +
  coord_equal() +
  ggtitle("Pearson Correlation among\nPredictor Variables") +
  labs(x = "", y = "",fill="Pearson's\nR") +
  theme(axis.text.x=element_text(size=5, angle=45, vjust=1, hjust=1, 
                                 margin=margin(-3,0,0,0)),
        axis.text.y=element_text(size=6, margin=margin(0,-3,0,0)),
        plot.title = element_text(hjust=0.5, face="bold", size=14),
        panel.grid.major=element_blank()) 

ggsave("plots/Figure2_Tau_PET_SUVR_corrplot.png",
       width=6, height=5, units="in", dpi=300)


tau_PCA <- stats::prcomp(tau_num_mat_z,
                         center = F,
                         scale. = F)

# Scree plot + eigenvalue plot
as.data.frame(get_eig(tau_PCA)) %>%
  dplyr::select(eigenvalue, cumulative.variance.percent) %>%
  dplyr::rename("Cumulative Variance" = "cumulative.variance.percent") %>%
  mutate(num_PC = row_number()) %>%
  mutate(kaiser_PC = max(num_PC[eigenvalue > 1]),
         cumvar_PC = min(num_PC[`Cumulative Variance` >= 80])) %>%
  pivot_longer(cols = c(eigenvalue,
                        `Cumulative Variance`)) %>%
  mutate(PC_to_plot = ifelse(name=="eigenvalue", 
                               kaiser_PC,
                               cumvar_PC),
         line_to_plot = ifelse(name=="eigenvalue", 
                               1,
                               80)) %>%
  dplyr::select(-kaiser_PC, -cumvar_PC) %>%
  ggplot(data = ., mapping = aes(x=num_PC, y=value)) +
  geom_bar(stat = "identity", fill = "lightsteelblue") +
  geom_hline(aes(yintercept = line_to_plot),
             linetype = 2) +
  geom_vline(aes(xintercept = PC_to_plot),
             linetype = 1) +
  facet_grid(name ~ ., scales="free", switch = "both") +
  ggtitle("Visualization of PCA Metrics by # PCs") +
  xlab("Number of PCs") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5),
        axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle=0))
ggsave("plots/Figure4_PCA_CumVar_Eigen.png",
       width=6, height=3.5, units="in", dpi=300)

# num features = 73, num observations = 690

# PC1 vs. PC2
pc12 <- fviz_pca_ind(tau_PCA,
                     addEllipses = TRUE,
                     geom = "point",
                     col.ind = tau_PET_data_wide$DX_bl) +
  labs(x = "PC1", y = "PC2", title = "") +
  theme(legend.position = "none",
        plot.title = element_blank())

# PC1 vs. PC3
pc13 <- fviz_pca_ind(tau_PCA,
                     addEllipses = TRUE,
                     axes = c(1,3),
                     geom = "point",
                     col.ind = tau_PET_data_wide$DX_bl) +
  labs(x = "PC1", y = "PC3", title = "", color = "Group", 
       fill = "Group",
       shape = "Group") +
  theme(plot.title = element_blank())

pc12 + pc13 +
  plot_annotation(title = "ADNI Subjects by diagnostic group") &
  theme(plot.title = element_text(hjust=0.5))
ggsave("plots/Figure5_PC_Scores_vs_Group_First3.png",
       width=6, height=3, units="in", dpi=300)

# Plot loadings for all brain regions onto first 9 PCs

# Extract the first 9 PCs to be used as input for random forest classifier
PC_scores_for_SVM <- as.data.frame(cbind(DX_bl = as.character(tau_PET_data_wide$DX_bl),
                                         tau_PCA$x[, 1:9]))  %>%
  mutate_at(vars(-("DX_bl")), as.numeric)

# Split data into 80% training, 20% testing
set.seed(127)
split_tau_PET_data_PC <- initial_split(PC_scores_for_SVM, 
                                    prop = 0.8)

# Preprocess training data
recipe_tau_PET_PC <-
  training(split_tau_PET_data_PC) %>% # Use the training data
  recipe(DX_bl ~ .) %>% # Define the RF formula
  step_corr(all_predictors()) %>% # Remove correlated predictors > 0.9
  prep() # Preprocess

# Extract preprocessed training data
train_tau_PET_PC <- juice(recipe_tau_PET_PC)

# Apply defined preprocessing to test data
test_tau_PET_PC <-
  recipe_tau_PET_PC %>%
  bake(testing(split_tau_PET_data_PC))

# Construct RF classifier using training data
rforest_tau_PET_PC <-
  rand_forest(trees = 100, mode = "classification") %>%
  set_engine("randomForest") %>%
  fit(DX_bl ~ ., data = train_tau_PET_PC)

# Apply RF classifier to test data
tau_PET_predicted_PC <- predict(rforest_tau_PET_PC, test_tau_PET_PC) 

# Model validation
validate_PC <-
  test_tau_PET_PC %>%
  bind_cols(tau_PET_predicted_PC) %>% # add tau_PET_predicted_PC into test data frame
  metrics(truth = DX_bl, # then perform comparison of original data in truth
          estimate = .pred_class) # against predicted data from model
validate_PC


vi_model(rforest_tau_PET_PC) %>%
  arrange(desc(Importance)) %>%
  mutate(Variable = fct_reorder(Variable, Importance, .desc=F)) %>%
  ggplot(data = ., mapping = aes(x = Variable, 
                                 y = Importance,
                                 fill = Importance)) +
  geom_bar(stat = "identity") +
  ggtitle("Variable Importance in\nPCA+RF Approach #2") +
  xlab("Brain Region") +
  scale_fill_viridis_b() +
  coord_flip() +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5))

ggsave("plots/Figure6_Approach2_VIP.png",
       width=3, height=6, units="in", dpi=300)
