# This script includes all functions that were created and used during the analytical process of the meta-analysis. 
# Each function is also listed in the scripts were they are used.

# -------------------------------------------
### Load packages ###

load_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

packages_needed <- c(
  "readxl", "writexl", "tidyverse", "grid", "gridExtra", "ggplot2", "patchwork", "gt", "webshot2", "gridGraphics",
  "ggpattern", "knitr", "esc", "metafor", "meta", "dmetar", "DeLuciatoR", "RColorBrewer", "ggsignif", "broom", "cowplot"
)

# -------------------------------------------
### Delete empty rows ###

clean_empty_rows <- function(df) {
  df <- df[!apply(is.na(df) | df == "", 1, all), ] 
}

# -------------------------------------------
### Run meta-analysis in first script ###

run_and_annotate <- function(data, metabolite_name, results_dir) {
  subset_data <- dataset %>% filter(Neurometabolite_Type == metabolite_name) 
  
  if(nrow(subset_data) < 3) {
    message(paste("Skipping", metabolite_name, "- need at least 3 studies"))
    return(NULL)
  }
  
  meta_result <- metagen(
    TE = SDM, seTE = `SDM Std Err`, data = subset_data, studlab = Study,
    common = FALSE, random = TRUE, method.tau = "DL", 
    method.random.ci = "classic", prediction = TRUE, sm = "SMD"
  )
  
  saveRDS(meta_result, file = file.path(results_dir, paste0("meta_", metabolite_name, ".rds")))
  return(subset_data)
}

# -------------------------------------------
### Run baseline meta-analysis in second script ###

run_baseline_meta_pipeline <- function(results_dir, metab_colors = c("Glu" = "#C76B72", "Glx" = "#E6B65E", "Overall" = "#05668d")) {
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  dataset_all <- bind_rows(
    if(exists("dataset_glu")) dataset_glu %>% mutate(Neurometabolite_Type = "Glu") else NULL,
    if(exists("dataset_glx")) dataset_glx %>% mutate(Neurometabolite_Type = "Glx") else NULL
  )
  
  overall_meta <- metagen(TE = SDM, seTE = `SDM Std Err`, data = dataset_all, studlab = Study,
                          common = FALSE, random = TRUE, method.tau = "DL", sm = "SMD")

  sink(file.path(results_dir, "Overall_Metabolite_Stats_Full.txt"))
  cat("--- OVERALL META-ANALYSIS (ALL METABOLITES COMBINED) ---\n")
  print(summary(overall_meta))
  sink()
  
  metabolites <- c("Glu", "Glx")
  meta_results_table <- data.frame()
  
  for (m in metabolites) {
    df_sub <- dataset_all %>% filter(Neurometabolite_Type == m)
    if (nrow(df_sub) >= 3) {
      m_meta <- metagen(TE = SDM, seTE = `SDM Std Err`, data = df_sub, studlab = Study,
                        common = FALSE, random = TRUE, method.tau = "DL", sm = "SMD")
      meta_results_table <- rbind(meta_results_table, data.frame(
        Metabolite = m, k = m_meta$k, Effect_Size = as.numeric(m_meta$TE.random),
        Lower = as.numeric(m_meta$lower.random), Upper = as.numeric(m_meta$upper.random),
        P_Value = as.numeric(m_meta$pval.random), Q_Stat = as.numeric(m_meta$Q),
        Q_P_Value = as.numeric(m_meta$pval.Q), I2 = as.numeric(m_meta$I2) * 100, stringsAsFactors = FALSE
      ))
    }
  }
  
    overall_row <- data.frame(
    Metabolite = "Overall", k = overall_meta$k, Effect_Size = as.numeric(overall_meta$TE.random),
    Lower = as.numeric(overall_meta$lower.random), Upper = as.numeric(overall_meta$upper.random),
    P_Value = as.numeric(overall_meta$pval.random), Q_Stat = as.numeric(overall_meta$Q),
    Q_P_Value = as.numeric(overall_meta$pval.Q), I2 = as.numeric(overall_meta$I2) * 100, stringsAsFactors = FALSE
  )
  
  meta_results_table <- bind_rows(overall_row, meta_results_table) %>%
    mutate(
      Stars = case_when(P_Value < 0.001 ~ "***", P_Value < 0.01 ~ "**", P_Value < 0.05 ~ "*", TRUE ~ ""),
      Metabolite = factor(Metabolite, levels = c("Overall", "Glu", "Glx"))
    )
  
  metab_levels <- unique(as.character(dataset_all$Neurometabolite_Type))
  pairwise_results <- data.frame()
  
  if (length(metab_levels) >= 2) {
    level_pairs <- utils::combn(metab_levels, 2, simplify = FALSE)
    for (pair in level_pairs) {
      p_subset <- dataset_all %>% filter(Neurometabolite_Type %in% pair)
      try({
        comp_meta_base <- metagen(TE = SDM, 
                                  seTE = `SDM Std Err`, 
                                  data = p_subset, 
                                  common = FALSE, 
                                  random = TRUE, 
                                  method.tau = "DL", 
                                  sm = "SMD")
        comp_meta <- update(comp_meta_base, subgroup = p_subset$Neurometabolite_Type)
        pairwise_results <- rbind(pairwise_results, data.frame(
          Group_1 = pair[1], 
          Group_2 = pair[2], 
          Q_Between = comp_meta$Q.b.random, 
          P_Value = comp_meta$pval.Q.b.random, 
          stringsAsFactors = FALSE))
      }, silent = TRUE)
    }
  }
  
  if (nrow(pairwise_results) > 0) {
    pairwise_results$FDR_p <- p.adjust(pairwise_results$P_Value, method = "fdr")
  }
  
  fig_main_results <- ggplot(meta_results_table %>% filter(Metabolite != "Overall"), 
                             aes(x = Metabolite, y = Effect_Size, fill = Metabolite)) + 
    geom_col(width = 0.6, color = "black", show.legend = FALSE) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15, linewidth = 0.8) +
    geom_hline(yintercept = 0.0, linetype = "solid", color = "grey30") +
    geom_text(aes(label = Stars, y = Upper + 0.05), size = 6, fontface = "bold", vjust = 0) +
    geom_text(aes(label = paste0("k=", k), y = Lower - 0.1), size = 4, fontface = "italic") +
    scale_y_continuous(name = "Pooled Effect Size (SMD)") +
    scale_fill_manual(values = metab_colors) +
    labs(title = "Summary of Metabolite Changes") +
    theme_linedraw() + 
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(), axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5), axis.text = element_text(color = "black", size = 12), 
      axis.title.x = element_blank(), axis.text.x  = element_text(angle = 0, hjust = 0.5) 
    )
  
  ggsave(filename = file.path(results_dir, "Figure_Main_Meta_Summary.pdf"), plot = fig_main_results, width = 7, height = 5)
  writexl::write_xlsx(list(Main_Results = meta_results_table, Pairwise_Comparisons = pairwise_results), file.path(results_dir, "Table_Main_Meta_Results.xlsx"))
  
  if (nrow(pairwise_results) > 0) {
    post_hoc_table <- pairwise_results %>%
      select(Group_1, Group_2, Q_Between, P_Value, FDR_p) %>%
      gt() %>%
      tab_header(title = "Post-Hoc Pairwise Comparisons Between Metabolites") %>%
      cols_label(Group_1 = "Metabolite", Group_2 = "Comparison", Q_Between = "Q-Statistic", P_Value = "p-value", FDR_p = "FDR Adjusted p") %>%
      fmt_number(columns = c(Q_Between), decimals = 2) %>%
      fmt_number(columns = c(P_Value, FDR_p), decimals = 3) %>%
      tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels()) %>%
      tab_style(style = cell_fill(color = "#f8f9fa"), locations = cells_body(rows = FDR_p < 0.05)) %>%
      tab_style(style = cell_text(color = "#A11332", weight = "bold"), locations = cells_body(columns = FDR_p, rows = FDR_p < 0.05)) %>%
      tab_options(table.border.top.color = "black", column_labels.border.bottom.color = "black", table.width = pct(80))
    
    gtsave(data = post_hoc_table, filename = file.path(results_dir, "Table_PostHoc_Final.pdf"))
  }
  
  total_sample_size <- dataset_all %>%
    distinct(Study_ID, .keep_all = TRUE) %>%
    summarise(Total_Cases = sum(Case_Sample_size, na.rm = TRUE),
              Total_Controls = sum(Control_Sample_size, na.rm = TRUE),
              Grand_Total = Total_Cases + Total_Controls)
  
  print(total_sample_size)
  
  return(list(summary_table = meta_results_table, pairwise = pairwise_results))
}

# -------------------------------------------
### Categorical moderator analyses in second script ###

run_categorical_moderator_pipeline <- function(subgroup_var, 
                                               subgroup_name, 
                                               results_dir, 
                                               x_axis_order = NULL) {
  
  safe_name <- gsub("[^a-zA-Z0-9_.-]", "_", subgroup_name)
  output_path <- file.path(results_dir, "Moderator_Analyses", safe_name)
  if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
  
  pipeline_theme <- theme_bw(base_size = 18) + 
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(color = "black", size = 12), 
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 16)
    )
  
  df_overall <- dataset_all %>%
    filter(!is.na(!!sym(subgroup_var)), !is.na(SDM)) %>%
    group_by(!!sym(subgroup_var)) %>%
    filter(n() >= 3) %>%
    ungroup()
  
  overall_output <- NULL
  
  if (nrow(df_overall) > 0) {
    meta_overall <- metagen(
      TE = SDM, seTE = `SDM Std Err`, data = df_overall, studlab = Study,
      common = FALSE, random = TRUE, method.tau = "DL", sm = "SMD",
      subgroup = df_overall[[subgroup_var]]
    )
    
    n_counts_overall <- df_overall %>%
      group_by(!!sym(subgroup_var)) %>%
      summarise(k = n(), N_studies = n_distinct(Study_ID), .groups = 'drop') %>%
      rename(subgroup = !!sym(subgroup_var))
    
    summary_overall <- data.frame(
      subgroup = names(meta_overall$k.w),
      estimate = as.numeric(meta_overall$TE.random.w),      
      p_value = as.numeric(meta_overall$pval.random.w),     
      conf.low = as.numeric(meta_overall$lower.random.w),   
      conf.high = as.numeric(meta_overall$upper.random.w),  
      stringsAsFactors = FALSE
    ) %>%
      left_join(n_counts_overall, by = "subgroup") %>%
      mutate(
        significance = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**", p_value < 0.05 ~ "*", TRUE ~ "ns"),
        plot_label = paste0(subgroup, "\n(k=", k, ")")
      )
    
    p_overall <- ggplot(summary_overall, aes(x = plot_label, y = estimate)) +
      geom_col(fill = "#5F6DAE", colour = "black", width = 0.6) + 
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
      geom_hline(yintercept = 0, linetype = "solid", color = "grey30") + 
      geom_text(aes(y = conf.high + 0.05, label = ifelse(significance == "ns", "", significance)), size = 6) +
      labs(title = paste("Overall Effect of", subgroup_name), x = NULL, y = "Pooled SMD") +
      pipeline_theme
    
    ggsave(file.path(output_path, "Overall_Subgroup_Plot.png"), p_overall, width = 8, height = 6)
    
    overall_output <- list(model = meta_overall, table_data = summary_overall)
  }
  
  metabolite_datasets <- list(
    Glu = if(exists("dataset_glu")) dataset_glu else NULL,
    Glx = if(exists("dataset_glx")) dataset_glx else NULL
  )
  metabolite_datasets <- metabolite_datasets[!sapply(metabolite_datasets, is.null)]
  
  summary_list <- list()
  all_pairwise_results <- data.frame()
  metab_colors <- c("Glu" = "#C76B72", "Glx" = "#E6B65E")
  
  for (metab in names(metabolite_datasets)) {
    df_current <- metabolite_datasets[[metab]] %>% filter(!is.na(!!sym(subgroup_var)))
    if (nrow(df_current) < 3) next
    
    meta_sub <- metagen(TE = SDM, seTE = `SDM Std Err`, data = df_current, studlab = Study,
                        common = FALSE, random = TRUE, method.tau = "DL", sm = "SMD")
    meta_subgroup <- update(meta_sub, subgroup = df_current[[subgroup_var]], subgroup.name = subgroup_name)
    
    sink(file.path(output_path, paste0("meta_results_", metab, ".txt")))
    cat(paste0("Analysis for Metabolite: ", metab, "\n"))
    print(summary(meta_subgroup))
    sink()
    
    n_counts_sub <- df_current %>%
      group_by(!!sym(subgroup_var)) %>%
      summarise(
        k = n(), N_studies = n_distinct(Study_ID),
        Total_N_Case = sum(distinct(pick(Study_ID, Case_Sample_size))$Case_Sample_size, na.rm = TRUE),
        Total_N_Control = sum(distinct(pick(Study_ID, Control_Sample_size))$Control_Sample_size, na.rm = TRUE),
        .groups = 'drop'
      ) %>% rename(subgroup = !!sym(subgroup_var))
    
    if (!is.null(meta_subgroup$bylevs)) {
      summary_df <- data.frame(
        Neurometabolite_Type = metab,
        subgroup = as.character(meta_subgroup$bylevs),
        Effect_Size = as.numeric(meta_subgroup$TE.random.w),
        Lower = as.numeric(meta_subgroup$lower.random.w),
        Upper = as.numeric(meta_subgroup$upper.random.w),
        Pvalue = as.numeric(meta_subgroup$pval.random.w),
        Q_stat = as.numeric(meta_subgroup$Q.w),    
        I2 = as.numeric(meta_subgroup$I2.w) * 100,
        stringsAsFactors = FALSE
      ) %>%
        left_join(n_counts_sub, by = "subgroup") %>%
        mutate(
          plot_label = paste0(subgroup, "\n(k=", k, ", N=", N_studies, ")\n", "C/H: ", Total_N_Case, "/", Total_N_Control),
          Significance_Stars = case_when(Pvalue < 0.001 ~ "***", Pvalue < 0.01  ~ "**", Pvalue < 0.05  ~ "*", TRUE ~ "")
        )
      
      summary_list[[metab]] <- summary_df
      
      p_metab <- ggplot(summary_df, aes(x = plot_label, y = Effect_Size, fill = Neurometabolite_Type)) +
        geom_col(colour = "black", width = 0.6, show.legend = FALSE) +
        geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.1) +
        geom_hline(yintercept = 0, linetype = "solid", color = "grey30") + 
        geom_text(aes(y = Upper + 0.05, label = Significance_Stars), size = 6) +
        scale_fill_manual(values = metab_colors) + 
        labs(title = paste("Average Effect Size:", metab), subtitle = paste("Moderator:", subgroup_name), x = NULL, y = "Effect Size (SMD)") +
        pipeline_theme
      
      ggsave(file.path(output_path, paste0(metab, "_Subgroup_Plot.png")), p_metab, width = 7, height = 6)
      
      mod_levels <- as.character(meta_subgroup$bylevs)
      if (length(mod_levels) >= 2) {
        level_pairs <- utils::combn(mod_levels, 2, simplify = FALSE)
        for (pair in level_pairs) {
          pair_df <- df_current[df_current[[subgroup_var]] %in% pair, ]
          try({
            pair_base <- metagen(TE = SDM, seTE = `SDM Std Err`, data = pair_df, common = FALSE, random = TRUE, method.tau = "DL", sm = "SMD")
            pair_meta <- update(pair_base, subgroup = pair_df[[subgroup_var]])
            all_pairwise_results <- rbind(all_pairwise_results, data.frame(
              Metabolite = metab, Group_1 = pair[1], Group_2 = pair[2],
              Q_Test = pair_meta$Q.b.random, P_Value = pair_meta$pval.Q.b.random, stringsAsFactors = FALSE
            ))
          }, silent = TRUE)
        }
      }
    }
  }
  
  excel_out_list <- list()
  
  if (length(summary_list) > 0) { 
    combined_summary <- bind_rows(summary_list) %>%
      mutate(Neurometabolite_Type = factor(Neurometabolite_Type, levels = c("Glu", "Glx"))) %>%
      arrange(Neurometabolite_Type)
    
    excel_out_list$Subgroup_Stats <- combined_summary
    
    if (length(summary_list) > 1) {
      comp_plot <- ggplot(combined_summary, aes(x = plot_label, y = Effect_Size, fill = Neurometabolite_Type)) +
        geom_col(position = position_dodge(width = 0.8), colour = "black", width = 0.7) +
        geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.8), width = 0.1) +
        geom_hline(yintercept = 0) +
        scale_fill_manual(values = metab_colors, name = "Metabolite") +
        geom_text(aes(y = Upper + 0.1, label = Significance_Stars), position = position_dodge(width = 0.8), size = 6, fontface = "bold") +
        labs(title = paste("Metabolite Comparison:", subgroup_name), x = NULL, y = "Effect Size (SMD)") +
        pipeline_theme + theme(legend.position = "right")
      
      ggsave(file.path(output_path, "Comparison_Across_Metabolites_Plot.png"), comp_plot, width = 12, height = 7)
    }
  }
  
  if (nrow(all_pairwise_results) > 0) {
    all_pairwise_results$FDR_p <- p.adjust(all_pairwise_results$P_Value, method = "fdr")
    excel_out_list$Pairwise_Comparisons <- all_pairwise_results
  }
  
  if (length(excel_out_list) > 0) {
    write_xlsx(excel_out_list, path = file.path(output_path, paste0(safe_name, "_Results_Summary.xlsx")))
  }
  
  return(list(overall = overall_output, metabolites = summary_list)) 
}

# -------------------------------------------
### Run continuous meta-regressions in the second script ###

run_continuous_regression_pipeline <- function(df, moderator, moderator_name, results_dir) {
  safe_mod_name <- gsub("[^a-zA-Z0-9_.-]", "_", gsub("#", "Number", gsub("%", "Percent", moderator_name)))
  df_clean <- df %>%
    filter(!is.na(!!sym(moderator))) %>%
    mutate(temp_mod = as.numeric(!!sym(moderator))) %>%
    filter(!is.na(temp_mod))
  
  if (nrow(df_clean) < 3) {
    message(paste("Skipping Regressions for", moderator, ": Too few studies."))
    return(NULL)
  }
  
  pub_theme <- theme_linedraw() + 
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_blank(), panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 13),
      legend.background = element_rect(fill = "white"), legend.position   = "bottom"
    )
  
  plot_list <- list()
  
  overall_folder <- file.path(results_dir, "Overall_Regressions", safe_mod_name)
  if (!dir.exists(overall_folder)) dir.create(overall_folder, recursive = TRUE)
  
  reg_overall <- rma(yi = SDM, sei = `SDM Std Err`, mods = ~ temp_mod, data = df_clean, method = "DL")
  p_val_ov <- reg_overall$pval[2]
  p_label_ov <- if (p_val_ov < 0.001) "p < 0.001" else paste0("p = ", round(p_val_ov, 3))
  
  results_txt <- paste0("Analysis: OVERALL (All Metabolites Combined)\nModerator: ", moderator_name, 
                        "\nk (obs): ", nrow(df_clean), "\nN (studies): ", n_distinct(df_clean$Study_ID),
                        "\nSlope: ", round(reg_overall$b[2], 4), "\nP-value: ", p_label_ov)
  writeLines(results_txt, file.path(overall_folder, paste0("Overall_Summary_", safe_mod_name, ".txt")))
  
  p_overall <- ggplot(df_clean, aes(x = temp_mod, y = SDM)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey80") +
    geom_point(aes(size = 1/`SDM Std Err`), shape = 21, fill = "#26547c", alpha = 0.5) +
    geom_smooth(method = "lm", color = "black", fill = "grey90", linewidth = 1) +
    labs(title = paste("Overall Analysis vs", moderator_name),
         subtitle = paste0("Slope = ", round(reg_overall$b[2], 3), " | ", p_label_ov),
         x = moderator_name, y = "Effect Size (SMD)", size = "Precision (1/SE)") +
    pub_theme
  
  ggsave(file.path(overall_folder, paste0("Overall_Plot_", safe_mod_name, ".png")), p_overall, width = 8, height = 6, dpi = 300)
  plot_list[["Overall"]] <- p_overall
  
  split_folder <- file.path(results_dir, "Regression_Glu_vs_Glx", safe_mod_name)
  if (!dir.exists(split_folder)) dir.create(split_folder, recursive = TRUE)
  
  metab_colors <- c("Glu" = "#C76B72", "Glx" = "#E6B65E")
  
  for (metab in c("Glu", "Glx")) {
    metab_df <- df_clean %>% filter(Neurometabolite_Type == metab)
    if (nrow(metab_df) <= 2) next
    
    reg_indiv <- rma(yi = SDM, sei = `SDM Std Err`, mods = ~ temp_mod, data = metab_df, method = "DL")
    variable_p <- if (reg_indiv$pval[2] < 0.001) "p < 0.001" else paste0("p = ", round(reg_indiv$pval[2], 3))
    
    sink(file.path(split_folder, paste0(metab, "_Meta_Regression.txt")))
    cat(paste0("--- ", metab, " META-REGRESSION ---\nk = ", nrow(metab_df), "\n"))
    print(summary(reg_indiv))
    sink()
    
    p_indiv <- ggplot(metab_df, aes(x = temp_mod, y = SDM)) +
      geom_hline(yintercept = 0, linetype = "solid", color = "grey80") +
      geom_point(shape = 21, size = 3, fill = metab_colors[metab], color = "black") +
      geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey90", alpha = 0.4) +
      labs(title = paste0(metab, " Meta-Regression: ", moderator_name),
           subtitle = paste0("Slope=", round(reg_indiv$b[2], 3), " | ", variable_p),
           x = moderator_name, y = "Effect Size (SMD)") +
      pub_theme
    
    ggsave(file.path(split_folder, paste0(metab, "_Plot.png")), p_indiv, width = 7, height = 5)
    plot_list[[metab]] <- p_indiv
  }
  
  if (length(unique(df_clean$Neurometabolite_Type)) >= 2) {
    reg_int <- rma(yi = SDM, sei = `SDM Std Err`, mods = ~ temp_mod * Neurometabolite_Type, data = df_clean, method = "DL")
    write.csv(as.data.frame(coef(summary(reg_int))), file.path(split_folder, "Interaction_Comparison_Results.csv"))
    
    comp_plot <- ggplot(df_clean, aes(x = temp_mod, y = SDM, color = Neurometabolite_Type, fill = Neurometabolite_Type)) +
      geom_hline(yintercept = 0, linetype = "solid", color = "grey80") +
      geom_point(aes(size = 1/`SDM Std Err`), shape = 21, alpha = 0.5, color = "black") +
      geom_smooth(method = "lm", alpha = 0.1, fullrange = TRUE) +
      scale_color_manual(values = metab_colors) + scale_fill_manual(values = metab_colors) +
      labs(title = paste("Comparison Across Metabolites:", moderator_name), subtitle = "Interaction Model Line Fit",
           x = moderator_name, y = "Effect Size (SMD)", size = "Precision (1/SE)") +
      pub_theme
    
    ggsave(file.path(split_folder, "Comparison_Interaction_Plot.png"), comp_plot, width = 8, height = 6)
    plot_list[["Comparison"]] <- comp_plot
  }
  
  return(plot_list)
}

# -------------------------------------------
### Comparison between metabolite statistics in the second script ###

run_metabolite_comparison_stats <- function(subgroup_var) {
  
  combined_data <- bind_rows(
    if(exists("dataset_glu")) dataset_glu else NULL,
    if(exists("dataset_glx")) dataset_glx else NULL
  )
  
  combined_data <- combined_data %>% filter(!is.na(!!sym(subgroup_var)))
  diff_stats <- data.frame()
  subtypes <- unique(combined_data[[subgroup_var]])
  
  for (s in subtypes) {
    sub_df <- combined_data %>% filter(!!sym(subgroup_var) == s)
    metab_count <- length(unique(sub_df$Neurometabolite_Type))
    if(nrow(sub_df) > 3 & metab_count > 1) {
      try({
        res <- rma(yi = SDM, sei = `SDM Std Err`, 
                   mods = ~ Neurometabolite_Type, 
                   data = sub_df, method = "DL")
        diff_stats <- rbind(diff_stats, data.frame(
          subgroup = s,
          Q_Value = res$QM,
          P_Value = res$QMp,
          stringsAsFactors = FALSE
        ))
      }, silent = TRUE)
    }
  }
  
  return(diff_stats)
}

# -------------------------------------------
### Define statistical significance ###

get_stars <- function(p) {
  if (length(p) == 0 || is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

get_stars_comparisons <- function(p) {
  if (length(p) == 0 || is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

# -------------------------------------------
### Pairwise comparison calculation ###

calc_pairwise_p <- function(g1, g2, model) {
  te1 <- model$TE.random.w[g1]
  te2 <- model$TE.random.w[g2]
  se1 <- model$seTE.random.w[g1]
  se2 <- model$seTE.random.w[g2]
  z <- (te1 - te2) / sqrt(se1^2 + se2^2)
  p <- 2 * (1 - pnorm(abs(z)))
  return(p)
}
