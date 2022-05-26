################################################################################

# Code for linear modeling and figure generation for OLET5608 May 2022 Assignment
# Annie G. Bryant
# The University of Sydney
# 26 May 2022

################################################################################

# Load packages
library(tidyverse)
library(ggfortify)
library(lmtest)
library(cowplot)
library(fastDummies)
library(patchwork)
library(knitr)
library(kableExtra)
library(car)
library(gridExtra)
library(lmtest)
library(sandwich)
library(MASS)

# Set plot theme to cowplot
theme_set(theme_cowplot())

# Image output folder
plot_path <- "../Report_Images/"


################################################################################
# Data cleaning
################################################################################

# We can start by reading in the data, provided in a CSV format. Here's the first four values and description of each feature in the dataset:
data_file <- "../Project_Dataset/ADNIMERGE.csv"
description_file <- read.csv("../Project_Dataset/Variable_Descriptions.csv")
ADNI_data <- read.csv(data_file, na.strings = c("NA", "")) %>%
  
  # Select a subset of columns of interest
  dplyr::select(RID, EXAMDATE, AGE, PTGENDER, PTEDUCAT, PTRACCAT, PTMARRY, APOE4, FDG, ABETA, PTAU, MMSE, CDRSB, ADAS11, ADAS13, Ventricles, Hippocampus, WholeBrain, Entorhinal, Fusiform, MidTemp, ICV) %>%
  
  # Convert exam date to date format
  mutate(EXAMDATE = as.Date(EXAMDATE)) %>%
  
  # Convert numeric variables to be treated as numeric
  mutate_at(c("ABETA", "PTAU"),
            function(x) as.numeric(x)) %>%
  group_by(RID) %>%
  
  # Drop na values for these columns and only keep
  # participants where marital status is known
  filter(!is.na(MMSE), !is.na(PTAU), !is.na(Entorhinal),
         !is.na(Hippocampus), !is.na(ABETA), 
         PTMARRY != "Unknown", PTRACCAT != "Unknown") %>%
  
  # Pick the first visit date available per subject
  filter(EXAMDATE == min(EXAMDATE, na.rm=T)) %>%
  dplyr::select(-EXAMDATE) %>%
  drop_na() %>%
  ungroup() %>%
  # Set APOE4 to boolean and subject ID to character
  mutate(APOE4 = ifelse(APOE4>0, TRUE, FALSE),
         RID = as.character(RID))

head_data <- head(ADNI_data, 4) %>%
  t()

# Table 1
as.data.frame(head_data) %>%
  unite("First 4 Values", V1:V4, sep=", ") %>%
  rownames_to_column(var="Variable") %>%
  left_join(., description_file)


# I manually encoded `PTGENDER`, `PTRACCAT`, `PTMARRY`, and `APOE4` as factors, since they are classes rather than continuously distributed numerical variables. I also discretized PTEDUCAT into a new variable (PTHIGHERED) two groups: (1) 16 or fewer years and (2) greater than 16 years, since this is not a truly continuous variable.
# Manually encode factors
ADNI_data <- ADNI_data %>%
  mutate(PTGENDER = factor(PTGENDER, levels=unique(PTGENDER)),
         PTRACCAT = factor(PTRACCAT, levels=unique(PTRACCAT)),
         PTMARRY = factor(PTMARRY, levels=unique(PTMARRY)),
         APOE4 = factor(APOE4, levels = c(FALSE, TRUE)))

# Discretize PTHIGHERED
ADNI_data$PTHIGHERED <- ifelse(ADNI_data$PTEDUCAT > 16, TRUE, FALSE)
ADNI_data <- ADNI_data %>% 
  dplyr::select(-PTEDUCAT) %>%
  mutate(PTHIGHERED = factor(PTHIGHERED, levels = c(FALSE, TRUE)))

# Figure 1 -- distribution of variables in ADNI dataset

# Function to emulate the ggplot default color palette for a given number of features
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Function to create either a histogram or bar chart for each feature
# Depending on whether it's numerical or categorical, respectively
plot_func <- function(df, plot_list) {
  num_features <- length(colnames(df))
  plot_palette <- gg_color_hue(num_features)
  for (i in 1:length(colnames(df))) {
    feature <- sort(colnames(df))[i]
    if (is.numeric(df %>% pull(feature))) {
      feature_p <- df %>%
        dplyr::select(all_of(feature)) %>%
        dplyr::mutate(row_id = dplyr::row_number()) %>%
        pivot_longer(cols=c(-row_id),
                     names_to = "Variable",
                     values_to = "Raw_Value") %>%
        ggplot(data=., mapping=aes(x=Raw_Value)) +
        geom_histogram(fill = plot_palette[i]) 
    } else {
      feature_p <- df %>%
        dplyr::select(all_of(feature)) %>%
        dplyr::mutate(row_id = dplyr::row_number()) %>%
        pivot_longer(cols=c(-row_id),
                     names_to = "Variable",
                     values_to = "Value") %>%
        group_by(Variable, Value) %>%
        count() %>%
        ggplot(data=., mapping=aes(x = Value, y = n)) +
        geom_bar(stat = "identity", fill = plot_palette[i])
    }
    feature_p <- feature_p  +
      facet_wrap(Variable ~ ., scales="free") +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.4, size=7),
            axis.text.y = element_text(size=7),
            plot.title = element_text(hjust=0.5, size=9),
            axis.title = element_blank(),
            strip.text = element_text(size=8),
            legend.position = "none")
    plot_list <- rlist::list.append(plot_list, feature_p)
  }
  
  return(plot_list)
}

# Generate list of plots, one per ADNI dataset feature (other than RID)
plot_list <- plot_func(ADNI_data %>% dplyr::select(-RID),
                       list())

# Draw the plots
patchwork::wrap_plots(plot_list, ncol=5) + plot_annotation(
  title = 'Distribution of Features in ADNI Dataset'
) & theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"),
          plot.title = element_text(size=12))
grid::grid.draw(grid::textGrob("Number of Subjects", x = 0.01, rot = 90))
grid::grid.draw(grid::textGrob("Feature Value", y = 0.01))

# Save to a PNG
ggsave(paste0(plot_path, "Figure1_Distribution_of_Numerical_Features_in_ADNI_Dataset.png"),
       width=6, height=5, units="in", dpi=300)

# Convert factors to dummy binary variables
ADNI_data_dummy <- dummy_cols(ADNI_data, 
                              select_columns = c("PTGENDER", "PTRACCAT", "PTMARRY", "PTHIGHERED", "APOE4"),
                              remove_first_dummy = TRUE,
                              remove_selected_columns = TRUE)

# Function to calculate z-score per feature
z_score <- function(variable_values) {
  z_scores <- (variable_values - mean(variable_values, na.rm=T))/sd(variable_values, na.rm=T)
  return(z_scores)
}

# Extract qualitative features + outcome variable
ADNI_qual <- ADNI_data %>%
  dplyr::select(RID, PTAU, PTGENDER, PTRACCAT, PTMARRY, APOE4, PTHIGHERED)

# Apply z-score function to quantitative features
ADNI_z_quant <- ADNI_data %>%
  dplyr::select(-PTAU) %>%
  dplyr::select(AGE, FDG:ICV) %>%
  dplyr::select_if(is.numeric) %>%
  # Create arbitrary ID to enable pivoting
  dplyr::mutate(row_id = dplyr::row_number()) %>%
  pivot_longer(cols=c(-row_id),
               names_to = "Variable",
               values_to = "Raw_Value") %>%
  # Group by variable
  group_by(Variable) %>%
  mutate(Z_Value = z_score(Raw_Value)) %>%
  pivot_wider(id_cols = row_id, names_from=Variable, values_from=Z_Value) %>%
  dplyr::select(-row_id) %>%
  cbind(., RID=ADNI_data$RID)

# Merge z-scored data back with qualitative features
ADNI_data_z <- left_join(ADNI_z_quant, ADNI_qual)

################################################################################
# Building models and testing assumptions
################################################################################


# Assumption 1: Linear relationship to outcome variable

# Figure 2
ADNI_z_quant %>%
  cbind(., dplyr::select(ADNI_data, c("PTAU"))) %>%
  pivot_longer(cols = c(-RID, -PTAU),
               names_to = "Predictor_Variable",
               values_to = "Value") %>%
  ggplot(data=., mapping=aes(x=Value, y=log(PTAU))) +
  geom_point(aes(color = Predictor_Variable), alpha = 0.5) +
  facet_wrap(Predictor_Variable ~ ., scales="free", ncol=5) +
  theme(legend.position = "none") +
  ylab("log-PTAU in CSF") +
  xlab("Predictor Value (Z-Scored)") +
  geom_smooth(method = "lm", color = "black")

ggsave(paste0(plot_path, "Figure2_Linear_Relationship_to_Outcome_Variable.png"),
       width=8, height=4.5, units="in", dpi=300)

# Construct OLS
model1 <- lm(PTAU ~ .,
             data = dplyr::select(ADNI_data_z,
                                  c(-RID)))

# Correlation plot among predictor variables
corr_mat <- ADNI_data_dummy %>%
  dplyr::select(-RID, -PTAU) %>%
  cor(., method = "pearson", use = "complete.obs")

ggplot(reshape2::melt(corr_mat), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8) +
  scale_fill_gradient2(low="blue", mid="white", high="red") +
  theme_minimal() +
  coord_equal() +
  ggtitle("Pearson Correlation among Predictor Variables") +
  labs(x = "", y = "",fill="Pearson's\nR") +
  theme(axis.text.x=element_text(size=7, angle=45, vjust=1, hjust=1, 
                                 margin=margin(-3,0,0,0)),
        axis.text.y=element_text(size=8, margin=margin(0,-3,0,0)),
        plot.title = element_text(hjust=0.5, face="bold", size=14),
        panel.grid.major=element_blank()) 

ggsave(paste0(plot_path, "Figure3_ADNI_Feature_Corrplot.png"),
       width=6, height=5, units="in", dpi=300)


# Multicollinearity -- calculating variance inflation factor (VIF)
cat("\nVariance Inflation Factor (VIF) values above 5:\n")
vif(model1)[,1][vif(model1)[,1] > 5]


# Construct OLS without multicollinearity
model2 <- lm(PTAU ~ .,
             data = dplyr::select(ADNI_data_z,
                                  c(-RID, -ADAS11, -ICV)))

# Assumption 2: Normal distribution of residuals

# Figure 4: Q-Q plot for Model 2
qq <- autoplot(model2, which = 2)
residual_histogram <- data.frame(Residual = model2$residuals) %>%
  ggplot(data=., mapping=aes(x=Residual)) +
  geom_histogram() +
  ggtitle("Distribution of OLS\nModel Residuals")
png(paste0(plot_path, "Figure4_Normal_Distribution_of_OLS_Residuals.png"),
    width=7, height=3.5, units="in", res=300)
qq + residual_histogram
dev.off()

# Figure 5: Box-Cox plot for Model 2
png(paste0(plot_path, "Figure5_OLS_BoxCox_Plot.png"),
    width=5, height=3.5, units="in", res=300)
boxcox(model2, plotit=T)
dev.off()

# Figure 6: Possible PTAU transforms
ADNI_data %>%
  dplyr::select(PTAU) %>%
  dplyr::rename("Original_Value" = "PTAU") %>%
  mutate(Log_Value = log(Original_Value),
         Sqrt_Value = sqrt(Original_Value),
         Inv_Value = 1/Original_Value) %>%
  pivot_longer(cols=c(Log_Value, Sqrt_Value, Original_Value, Inv_Value)) %>%
  mutate(name = factor(name, levels = c("Original_Value",
                                        "Log_Value",
                                        "Sqrt_Value",
                                        "Inv_Value"))) %>%
  ggplot(data=., mapping=aes(x=value)) +
  ggtitle("Distribution of PTAU with various power transforms") +
  geom_histogram(aes(fill=name)) +
  facet_grid(. ~ name, scales="free") +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none",
        axis.text.x = element_text(angle=90, vjust=0.4, hjust=1)) +
  xlab("Value") +
  ylab("Number of Subjects")
ggsave(paste0(plot_path, "Figure6_PTAU_Transformations.png"),
       width=7, height=3, units="in", dpi=300)


# Construct Model 3
model3 <- lm(log(PTAU) ~ ., 
             data = dplyr::select(ADNI_data_z,
                                  c(-RID, -ADAS11, -ICV)))

# Figure 7: Q-Q plot for Model 3
qq <- autoplot(model3, which = 2)
residual_histogram <- data.frame(Residual = model3$residuals) %>%
  ggplot(data=., mapping=aes(x=Residual)) +
  geom_histogram() +
  ggtitle("Distribution of Residuals\nfor Model 3")
png(paste0(plot_path, "Figure7_Normal_Distribution_of_OLS_Residuals_logPTAU_Model3.png"),
    width=7, height=3.5, units="in", res=300)
qq + residual_histogram
dev.off()

# Assumption 3: Constant variance of residuals

# Figure 8: Residual variance for models 2 and 3
png(paste0(plot_path, "Figure8_Constant_Residual_Variance_autoplot_models_2_3.png"),
    width=7, height=6, units="in", res=300)
model2_autoplot <- autoplot(model2, which = c(1,3))
model3_autoplot <- autoplot(model3, which = c(1,3))

grid.arrange(model2_autoplot@plots[[1]], model2_autoplot@plots[[2]],
             model3_autoplot@plots[[1]], model3_autoplot@plots[[2]], nrow = 2)
dev.off()

# Run Breusch-Pagan tests
# Model 2
bptest(model2)

# Model 3
bptest(model3)

# Model 4 = robust standard errors for model 3
model4 <- coeftest(model3, vcov = vcovHC(model3, "HC3"), save=T)


################################################################################
# Results + model selection
################################################################################

# Stepwise regression
model5 <- stepAIC(model3, direction = "both", 
                  trace = FALSE)
summary(model5)


# Model summary
# Stepwise regression model
model5_df <- data.frame(model5_coefs = round(model5$coefficients, 3),
                        model5_SE = round(summary(model5)$coefficients[,2], 3),
                        model5_pvals = summary(model5)$coefficients[,2], 3) %>%
  rownames_to_column(var="Term") %>%
  filter(Term != "(Intercept)") %>%
  dplyr::select(-X3)

# R2 and adjusted R2
R2_df <- data.frame(Term = c("R2", "Adjusted R2"),
                    Model3 = c(summary(model3)$r.squared, summary(model3)$adj.r.squared),
                    Model4 = c(summary(model3)$r.squared, summary(model3)$adj.r.squared),
                    Model5 = c(summary(model5)$r.squared, summary(model5)$adj.r.squared))

data.frame(model3_coefs = round(model3$coefficients, 3),
           model3_SE = round(summary(model3)$coefficients[,2], 3),
           model3_pvals = summary(model3)$coefficients[,4],
           model4_coefs = round(model4[,1], 3),
           model4_SE = round(model4[,2], 3),
           model4_pvals = model4[,4]) %>%
  rownames_to_column(var="Term") %>%
  filter(Term != "(Intercept)") %>%
  left_join(., model5_df) %>%
  arrange(model3_pvals) %>%
  mutate(model3_Stars = case_when(model3_pvals < 0.0001 ~ "***",
                                  model3_pvals < 0.001 ~ "**",
                                  model3_pvals < 0.05 ~ "*",
                                  T ~ ""),
         model4_Stars = case_when(model4_pvals < 0.0001 ~ "***",
                                  model4_pvals < 0.001 ~ "**",
                                  model4_pvals < 0.05 ~ "*",
                                  T ~ ""),
         model5_Stars = case_when(model5_pvals < 0.0001 ~ "***",
                                  model5_pvals < 0.001 ~ "**",
                                  model5_pvals < 0.05 ~ "*",
                                  T ~ "")) %>%
  dplyr::select(-model3_pvals, -model4_pvals, -model5_pvals) %>%
  mutate("Model3" = paste0(model3_coefs, model3_Stars, " (", model3_SE, ")"),
         "Model4" = paste0(model4_coefs, model4_Stars, " (", model4_SE, ")"),
         "Model5" = ifelse(is.na(model5_coefs), "", (paste0(model5_coefs, model5_Stars, " (", model5_SE, ")"))),
         .keep = "unused") %>%
  plyr::rbind.fill(., R2_df)