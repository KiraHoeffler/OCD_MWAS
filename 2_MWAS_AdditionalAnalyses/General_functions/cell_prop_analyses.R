

#' Extract variables needed for the cell proportions
#'
#' @param samplesheet needs to contain the columns Percent_improvement, Sex, Age_Collection_years, smoking_CpG_beta, Epi, Fib, B,
#' NK, CD4T, Mono, Neutro, Time_point, Resp_nr
#' @param ctrl_PCA PCA from the positive control probes with PC1, PC2, PC3 as columns
#' @param SNP_PCA PCA from the ancestry information with PC1, PC2, PC3 as columns
#'
#' @return df with variables as columns, samples as rows
#' @export
extract_variables_cellprop <- function(samplesheet, ctrl_PCA, SNP_PCA){
  variables_df <- data.frame(
    Response = as.numeric(samplesheet$Percent_improvement),
    Sex  = as.factor(samplesheet$Sex),
    Age = scale(as.numeric(samplesheet$Age)),
    smoking = scale(as.numeric(samplesheet$smoking_CpG)),
    Epi = scale(as.numeric(samplesheet$Epi)),
    Fib = scale(as.numeric(samplesheet$Fib)),
    B = scale(as.numeric(samplesheet$B)),
    NK = scale(as.numeric(samplesheet$NK)),
    CD4T = scale(as.numeric(samplesheet$CD4T)),
    Mono = scale(as.numeric(samplesheet$Mono)),
    Neutro = scale(as.numeric(samplesheet$Neutro)),
    ctrl_PC1 = ctrl_PCA$PC1,
    ctrl_PC2 = ctrl_PCA$PC2,
    ctrl_PC3 = ctrl_PCA$PC3,
    ctrl_PC5 = ctrl_PCA$PC5,
    ctrl_PC6 = ctrl_PCA$PC6,
    ctrl_PC7 = ctrl_PCA$PC7,
    ctrl_PC8 = ctrl_PCA$PC8,
    ctrl_PC9 = ctrl_PCA$PC9,
    ctrl_PC10 = ctrl_PCA$PC10,
    ctrl_PC11 = ctrl_PCA$PC11,
    ctrl_PC12 = ctrl_PCA$PC12,
    ctrl_PC13 = ctrl_PCA$PC13,
    ctrl_PC14 = ctrl_PCA$PC14,
    ctrl_PC15 = ctrl_PCA$PC15,
    snps_PC1 = SNP_PCA$PC1,
    snps_PC2 = SNP_PCA$PC2,
    snps_PC3 = SNP_PCA$PC3,
    snps_PC4 = SNP_PCA$PC4,
    snps_PC5 = SNP_PCA$PC5,
    Time_point = as.factor(samplesheet$Time_point),
    Resp_nr = as.factor(samplesheet$Resp_nr)
  )
  return(variables_df)
}







################################################################################

#' Run cross-sectional model for cell proportion analysis
#'
#' @param variables_df see function above (should be already filtered on the time point)
#' @param time_point time point to label
#' @param cell_type_vector all cell types you want to include in the analysis (should be column names in variables_df)
#'
#' @return results_df with statistics on Response, Sex, Age, and smoking
#' @export
cross_sectional_cellprop <- function(variables_df, time_point, cell_type_vector){
  
  # combine cell types in list
  cell_type_list <- list()
  for (cellprop in cell_type_vector){
    cell_type_list[[cellprop]] <- variables_df[, cellprop]
  }
  
  # create empty df for results
  comb_results <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 18))
  colnames(comb_results) <- c("Time_point", "Cell_type", 
                              "Estimate_Response", "Std_Error_Response", "t_value_Response", "p_value_Response",
                              "Estimate_Sex", "Std_Error_Sex", "t_value_Sex", "p_value_Sex",
                              "Estimate_Age", "Std_Error_Age", "t_value_Age", "p_value_Age", 
                              "Estimate_smoking", "Std_Error_smoking", "t_value_smoking", "p_value_smoking")
  
  # run model and add results
  for (i in c(1:length(cell_type_vector))){
    cell_type <- cell_type_list[[i]]
    
    variables_df_copy <- variables_df
    variables_df_copy$cell_type <- cell_type
    
    model <- lm(cell_type ~ Response + Sex + Age + smoking + 
                  ctrl_PC1 + ctrl_PC2 + ctrl_PC3 + ctrl_PC5 + ctrl_PC6 + ctrl_PC7 + ctrl_PC8 + ctrl_PC9 + ctrl_PC10 +
                  ctrl_PC11 + ctrl_PC12 + ctrl_PC13 + ctrl_PC14 + ctrl_PC15 +
                  snps_PC1 + snps_PC2 + snps_PC3 + snps_PC4 + snps_PC5, data = variables_df_copy)
    model_summary <- summary(model)
    rm(variables_df_copy)
    
    results <- data.frame(Time_point = time_point, 
                          Cell_type = cell_type_vector[i],
                          Estimate_Response = model_summary$coefficients["Response", "Estimate"],
                          Std_Error_Response = model_summary$coefficients["Response", "Std. Error"],
                          t_value_Response = model_summary$coefficients["Response", "t value"],
                          p_value_Response = model_summary$coefficients["Response", "Pr(>|t|)"],
                          Estimate_Sex = model_summary$coefficients["SexM", "Estimate"],
                          Std_Error_Sex = model_summary$coefficients["SexM", "Std. Error"],
                          t_value_Sex = model_summary$coefficients["SexM", "t value"],
                          p_value_Sex = model_summary$coefficients["SexM", "Pr(>|t|)"],
                          Estimate_Age = model_summary$coefficients["Age", "Estimate"],
                          Std_Error_Age = model_summary$coefficients["Age", "Std. Error"],
                          t_value_Age = model_summary$coefficients["Age", "t value"],
                          p_value_Age = model_summary$coefficients["Age", "Pr(>|t|)"],
                          Estimate_smoking =  model_summary$coefficients["smoking", "Estimate"],
                          Std_Error_smoking = model_summary$coefficients["smoking", "Std. Error"],
                          t_value_smoking = model_summary$coefficients["smoking", "t value"],
                          p_value_smoking = model_summary$coefficients["smoking", "Pr(>|t|)"])
    
    comb_results <- rbind(comb_results, results)
    
  }
  
  return(comb_results)
}









################################################################################

#' Run "Independent-of-Time" model for cell proportion analysis
#'
#' @param variables_df see function above (should be already filtered on the time point)
#' @param cell_type_vector all cell types you want to include in the analysis (should be column names in variables_df)
#'
#' @return results_df with results on Response, Time point, Sex, Age, and smoking
#' @export
IndepOfTime_cellprop <- function(variables_df, cell_type_vector){
  
  # combine cell types in list
  cell_type_list <- list()
  for (cellprop in cell_type_vector){
    cell_type_list[[cellprop]] <- variables_df[, cellprop]
  }
  
  colname_vector <- c("Cell_type", "AIC", 
                      "Estimate_Resp", "Std_Error_Resp", "t_value_Resp", "p_value_Resp", 
                      "Estimate_D4", "Std_Error_D4", "t_value_D4", "p_value_D4", 
                      "Estimate_M3", "Std_Error_M3", "t_value_M3", "p_value_M3", 
                      "Estimate_Sex", "Std_Error_Sex", "t_value_Sex", "p_value_Sex", 
                      "Estimate_Age", "Std_Error_Age", "t_value_Age", "p_value_Age", 
                      "Estimate_smoking", "Std_Error_smoking", "t_value_smoking", "p_value_smoking")
  
  # create empty df for results
  comb_results <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 26))
  colnames(comb_results) <- colname_vector
  
  # run model and add results
  for (i in c(1:length(cell_type_vector))){
    cell_type <- cell_type_list[[i]]
    
    variables_df_copy <- variables_df
    variables_df_copy$cell_type <- cell_type
    
    
    model <- lmerTest::lmer(cell_type ~ Response + Time_point + Sex + Age + smoking + 
                              ctrl_PC1 + ctrl_PC2 + ctrl_PC3 + ctrl_PC5 + ctrl_PC6 + ctrl_PC7 + ctrl_PC8 + ctrl_PC9 + ctrl_PC10 +
                              ctrl_PC11 + ctrl_PC12 + ctrl_PC13 + ctrl_PC14 + ctrl_PC15 +
                              snps_PC1 +  snps_PC2 + snps_PC3  + snps_PC4 + snps_PC5 + (1 | Resp_nr), data = variables_df_copy)
    model_summary <- summary(model)
    rm(variables_df_copy)
    
    results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 26))
    colnames(results) <- colname_vector
    
    
    results[1,1] <- cell_type_vector[i]
    results[1,2] <- AIC(model)
    
    results[1,3] <- model_summary$coefficients["Response", "Estimate"]
    results[1,4] <- model_summary$coefficients["Response", "Std. Error"]
    results[1,5] <- model_summary$coefficients["Response", "t value"]
    results[1,6] <- model_summary$coefficients["Response", "Pr(>|t|)"]
    
    results[1,7] <- model_summary$coefficients["Time_pointD4", "Estimate"]
    results[1,8] <- model_summary$coefficients["Time_pointD4", "Std. Error"]
    results[1,9] <- model_summary$coefficients["Time_pointD4", "t value"]
    results[1,10] <- model_summary$coefficients["Time_pointD4", "Pr(>|t|)"]
    
    results[1,11] <- model_summary$coefficients["Time_pointM3", "Estimate"]
    results[1,12] <- model_summary$coefficients["Time_pointM3", "Std. Error"]
    results[1,13] <- model_summary$coefficients["Time_pointM3", "t value"]
    results[1,14] <- model_summary$coefficients["Time_pointM3", "Pr(>|t|)"]
    
    results[1,15] <- model_summary$coefficients["SexM", "Estimate"]
    results[1,16] <- model_summary$coefficients["SexM", "Std. Error"]
    results[1,17] <- model_summary$coefficients["SexM", "t value"]
    results[1,18] <- model_summary$coefficients["SexM", "Pr(>|t|)"]
    
    results[1,19] <- model_summary$coefficients["Age", "Estimate"]
    results[1,20] <- model_summary$coefficients["Age", "Std. Error"]
    results[1,21] <- model_summary$coefficients["Age", "t value"]
    results[1,22] <- model_summary$coefficients["Age", "Pr(>|t|)"]
    
    results[1,23] <- model_summary$coefficients["smoking", "Estimate"]
    results[1,24] <- model_summary$coefficients["smoking", "Std. Error"]
    results[1,25] <- model_summary$coefficients["smoking", "t value"]
    results[1,26] <- model_summary$coefficients["smoking", "Pr(>|t|)"]
    
    
    comb_results <- rbind(comb_results, results)
    
  }
  
  return(comb_results)
}










################################################################################


#' Run "Independent-of-Time" model for cell proportion analysis
#'
#' @param variables_df 
#' @param cell_type_vector 
#'
#' @return 
#' @export

#' Run longitudinal model for cell proportion analysis
#'
#' @param variables_df see function above (should be already filtered on the time point)
#' @param cell_type_vector all cell types you want to include in the analysis (should be column names in variables_df)
#'
#' @return results_df with results on Response x Time point
#' @export
Longitudinal_cellprop <- function(variables_df, cell_type_vector){
  
  # combine cell types in list
  cell_type_list <- list()
  for (cellprop in cell_type_vector){
    cell_type_list[[cellprop]] <- variables_df[, cellprop]
  }
  
  colname_vector <- c("Cell_type",
                      "Estimate_RespD4", "Std_Error_RespD4", "t_value_RespD4", "p_value_RespD4", 
                      "Estimate_RespM3", "Std_Error_RespM3", "t_value_RespM3", "p_value_RespM3")
  
  # create empty df for results
  comb_results <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 9))
  colnames(comb_results) <- colname_vector
  
  # run model and add results
  for (i in c(1:length(cell_type_vector))){
    cell_type <- cell_type_list[[i]]
    
    variables_df_copy <- variables_df
    variables_df_copy$cell_type <- cell_type
    
    model <- lm(cell_type ~ Response * Time_point + Sex + Age + smoking + 
                  ctrl_PC1 + ctrl_PC2 + ctrl_PC3 + ctrl_PC5 + ctrl_PC6 + ctrl_PC7 + ctrl_PC8 + ctrl_PC9 + ctrl_PC10 +
                  ctrl_PC11 + ctrl_PC12 + ctrl_PC13 + ctrl_PC14 + ctrl_PC15 +
                  snps_PC1 +  snps_PC2 + snps_PC3  + snps_PC4 + snps_PC5, data = variables_df_copy)
    clustered_covariance <- vcovCL(model, cluster = ~variables_df_copy$Resp_nr)
    model_summary <- coeftest(model, clustered_covariance)
    
    results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 9))
    colnames(results) <- colname_vector
    
    results[1, 1] <- cell_type_vector[i]
    results[1, 2] <- model_summary["Response:Time_pointD4", "Estimate"]
    results[1, 3] <- model_summary["Response:Time_pointD4", "Std. Error"]
    results[1, 4] <- model_summary["Response:Time_pointD4", "t value"]
    results[1, 5] <- model_summary["Response:Time_pointD4", "Pr(>|t|)"]
    results[1, 6] <- model_summary["Response:Time_pointM3", "Estimate"]
    results[1, 7] <- model_summary["Response:Time_pointM3", "Std. Error"]
    results[1, 8] <- model_summary["Response:Time_pointM3", "t value"]
    results[1, 9] <- model_summary["Response:Time_pointM3", "Pr(>|t|)"]
    
    comb_results <- rbind(comb_results, results)
    
  }
  
  return(comb_results)
}









################################################################################

#' Categorical response plot for cell type proportions
#'
#' @param cell_type cell type how the column is called in the samplesheet (e.g. "Mono")
#' @param cell_type_label written label for plot (e.g. "Monocytes")
#' @param samplesheet samplesheet for which plot should be generated (incl. columns Time_point, the respsective cell type proportion and Response_status)
#' @param avg_cell_type df with two columns: Cell_type and average -> shows the average of cell type proportions for the controls
#'
#' @return plot
#' @export
#'
#' @examples
categ_cell_plot <- function(cell_type, cell_type_label, samplesheet, avg_cell_type){
  cell_types_df <- samplesheet[, c("Time_point", cell_type, "Response_status")]
  colnames(cell_types_df) <- c("Time_point", "cell_type", "Response_status")
  
  df_avg <- cell_types_df %>%
    group_by(Time_point, Response_status) %>%
    summarise(
      avg = mean(cell_type, na.rm = TRUE),
      se = sd(cell_type, na.rm = TRUE) / sqrt(n()))
  
  Ctrl_value <- avg_cell_type$average[which(avg_cell_type$Cell_type == cell_type)]
  
  if (mean(df_avg$avg) > 0.05){
    added_nr <- 0.003
  }
  
  if (mean(df_avg$avg) <= 0.05){
    added_nr <- 0.0003
  }
  
  CellTypes_plot <- ggplot(df_avg, aes(x = Time_point, y = avg, color = Response_status, group = Response_status)) +
    geom_point(size = 3) +
    geom_line() +
    geom_hline(yintercept = Ctrl_value, color = ctrl_color, linetype = "dashed") +
    annotate("text", x=3.08, y = Ctrl_value + added_nr, label = "CTRL", color = ctrl_color, hjust = 0, fontface = "bold") +
    geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width = 0.2) +
    th + th_transparent +
    ggtitle(cell_type_label) +
    labs(
      x = "Time point",
      y = "Estimated cell type proportion"
    ) + 
    theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(), plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=c("responder"= resp_color, "non-responder" = non_resp_color))
  
  return(CellTypes_plot)
}








################################################################################

#' Correlation plot for cell type proportions
#'
#' @param cell_type cell type how the column is called in the samplesheet (e.g. "Mono")
#' @param cell_type_label written label for plot (e.g. "Monocytes")
#' @param samplesheet samplesheet for which plot should be generated (incl. columns Time_point, the respsective cell type proportion and Response_status)
#' @param TP2 time point of second time point ("D4" or "M3")
#'
#' @return plot
#' @export
cell_corr_plot <- function(cell_type, cell_type_label, samplesheet, TP2){
  TP1 <- "D1"
  cell_types_df <- samplesheet[, c("Resp_nr", "Time_point", cell_type, "Percent_improvement")]
  colnames(cell_types_df) <- c("Resp_nr", "Time_point", "cell_type", "Percent_improvement")
  
  cell_types_df_comb <- unique(samplesheet[, c("Resp_nr", "Percent_improvement")])
  cell_types_df_comb$change <- NA
  
  for (i in c(1:nrow(cell_types_df_comb))){
    Resp_nr <- cell_types_df_comb$Resp_nr[i]
    samples_indiv <- cell_types_df[which(cell_types_df$Resp_nr == Resp_nr), ]
    
    if (TP1 %in% samples_indiv$Time_point & TP2 %in% samples_indiv$Time_point){
      TP1_sample <- samples_indiv$cell_type[which(samples_indiv$Time_point == TP1)]
      TP2_sample <- samples_indiv$cell_type[which(samples_indiv$Time_point == TP2)]
      change <- TP2_sample - TP1_sample
      cell_types_df_comb$change[i] <- change
    }
  }
  
  cell_types_df_comb <- na.omit(cell_types_df_comb)
  
  if (TP2 == "D4"){
    yaxistitle <- expression(bold(paste(Delta, "estimated cell proportion (Post - Pre)")))
  }
  
  if (TP2 == "M3"){
    yaxistitle <- expression(bold(paste(Delta, "estimated cell proportion (Follow-up - Pre)")))
  }
  
  # CORRELATION
  Pearson_cor <- cor(cell_types_df_comb$Percent_improvement, cell_types_df_comb$change)
  corr <- round(Pearson_cor,3)
  
  
  Corr_plot <- ggplot(cell_types_df_comb, aes(x = Percent_improvement, y = change)) +
    geom_point(color = "#084081", alpha = 0.15) + 
    geom_smooth(color = "#084081", method = "lm", linewidth=1.2, alpha = 0.35) +
    scale_color_identity() +
    ggtitle(cell_type_label) +
    xlab("% clinical improvement") + ylab(yaxistitle) +
    th + th_transparent + 
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate("text", x=0, y = min(cell_types_df_comb$change), label = paste0("r = ", corr), color = "black")
  
  return(Corr_plot)
}

