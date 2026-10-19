rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(parallel)
  library(readxl)
  library(openxlsx)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("FATAL ERROR: No task index provided from SLURM!")
}

TARGET_INDEX <- as.numeric(args[1])
cat(">>> Starting Monte Carlo Simulation for Global Task Index =", TARGET_INDEX, "<<<\n")

BASE_DIR <- "/N/slate/renzime/R/final/all_in_one"
setwd(BASE_DIR)
FUNCTION_FILE <- "/N/slate/renzime/R/final/functions.R"
TRUTH_TABLE_PATH <- "/N/slate/renzime/R/final/Truth_Table.xlsx"
TRUTH_TABLE_SHEET <- "tb"

SIMULATION_TIMES <- 10000L
SAMPLE_SIZE_LIST <- c(10L, 15L, 20L, 25L, 30L, 40L, 50L, 60L, 80L, 100L, 130L, 160L, 200L)
REQUESTED_CORES <- 16L
USE_LOGICAL_CORES <- FALSE
LEAVE_FREE_CORES <- 0L

BETA_INT <- 20.5
BETA_BASE <- 0.4
BETA_P2 <- -2
WRITE_XLSX <- TRUE
SAVE_RDS <- TRUE

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

source(FUNCTION_FILE)

make_design_A <- function() {
  tibble::tribble(
    ~seq, ~BASELINE, ~ACTIVE, ~PERIOD2, ~CARRYOVER, ~SUBG,
    1, 1, 0, 0, 0, 1,
    1, 0, 1, 0, 0, 1,
    1, 0, 0, 1, 1, 1,
    1, 1, 0, 0, 0, 0,
    1, 0, 1, 0, 0, 0,
    1, 0, 0, 1, 1, 0
  )
}

make_design_B <- function() {
  tibble::tribble(
    ~seq, ~BASELINE, ~ACTIVE, ~PERIOD2, ~CARRYOVER, ~SUBG,
    2, 1, 0, 0, 0, 1,
    2, 0, 0, 0, 0, 1,
    2, 0, 1, 1, 0, 1,
    2, 1, 0, 0, 0, 0,
    2, 0, 0, 0, 0, 0,
    2, 0, 1, 1, 0, 0
  )
}

make_cov_dict <- function() {
  tibble::tribble(
    ~Matrix_Index, ~Matrix_Name, ~e11, ~e12, ~e13, ~e21, ~e22, ~e23, ~e31, ~e32, ~e33,
    1, "CS_3_0.4", 3.0, 1.2, 1.2, 1.2, 3.0, 1.2, 1.2, 1.2, 3.0,
    2, "CS_6_0.4", 6.0, 2.4, 2.4, 2.4, 6.0, 2.4, 2.4, 2.4, 6.0,
    3, "CS_9_0.4", 9.0, 3.6, 3.6, 3.6, 9.0, 3.6, 3.6, 3.6, 9.0,
    4, "CS_12_0.4", 12.0, 4.8, 4.8, 4.8, 12.0, 4.8, 4.8, 4.8, 12.0,
    5, "CS_15_0.4", 15.0, 6.0, 6.0, 6.0, 15.0, 6.0, 6.0, 6.0, 15.0,
    6, "CS_18_0.4", 18.0, 7.2, 7.2, 7.2, 18.0, 7.2, 7.2, 7.2, 18.0,
    7, "CS_20_0.4", 20.0, 8.0, 8.0, 8.0, 20.0, 8.0, 8.0, 8.0, 20.0,
    8, "CS_24_0.4", 24.0, 9.6, 9.6, 9.6, 24.0, 9.6, 9.6, 9.6, 24.0,
    9, "CS_3_0.6", 3.0, 1.8, 1.8, 1.8, 3.0, 1.8, 1.8, 1.8, 3.0,
    10, "CS_6_0.6", 6.0, 3.6, 3.6, 3.6, 6.0, 3.6, 3.6, 3.6, 6.0,
    11, "CS_9_0.6", 9.0, 5.4, 5.4, 5.4, 9.0, 5.4, 5.4, 5.4, 9.0,
    12, "CS_12_0.6", 12.0, 7.2, 7.2, 7.2, 12.0, 7.2, 7.2, 7.2, 12.0,
    13, "CS_15_0.6", 15.0, 9.0, 9.0, 9.0, 15.0, 9.0, 9.0, 9.0, 15.0,
    14, "CS_18_0.6", 18.0, 10.8, 10.8, 10.8, 18.0, 10.8, 10.8, 10.8, 18.0,
    15, "CS_20_0.6", 20.0, 12.0, 12.0, 12.0, 20.0, 12.0, 12.0, 12.0, 20.0,
    16, "CS_24_0.6", 24.0, 14.4, 14.4, 14.4, 24.0, 14.4, 14.4, 14.4, 24.0,
    17, "UN_UN_1", 3.0, 0.8*sqrt(3.0*6.0), 0.1*sqrt(3.0*9.0), 0.8*sqrt(3.0*6.0), 6.0, 0.2*sqrt(6.0*9.0), 0.1*sqrt(3.0*9.0), 0.2*sqrt(6.0*9.0), 9.0,
    18, "UN_UN_2", 6.0, 0.8*sqrt(6.0*15.0), 0.1*sqrt(6.0*9.0), 0.8*sqrt(6.0*15.0), 15.0, 0.2*sqrt(15.0*9.0), 0.1*sqrt(6.0*9.0), 0.2*sqrt(15.0*9.0), 9.0,
    19, "UN_UN_3", 3.0, 0.8*sqrt(3.0*15.0), 0.1*sqrt(3.0*3.0), 0.8*sqrt(3.0*15.0), 15.0, 0.2*sqrt(15.0*3.0), 0.1*sqrt(3.0*3.0), 0.2*sqrt(15.0*3.0), 3.0,
    20, "UN_RW_4", 8.0, 0.5*sqrt(8.0*16.0), 0.2*sqrt(8.0*24.0), 0.5*sqrt(8.0*16.0), 16.0, 0.5*sqrt(16.0*24.0), 0.2*sqrt(8.0*24.0), 0.5*sqrt(16.0*24.0), 24.0,
    21, "UN_RW_5", 9.0, 0.5*sqrt(9.0*18.0), 0.2*sqrt(9.0*27.0), 0.5*sqrt(9.0*18.0), 18.0, 0.5*sqrt(18.0*27.0), 0.2*sqrt(9.0*27.0), 0.5*sqrt(18.0*27.0), 27.0,
    22, "UN_RW_6", 10.0, 0.5*sqrt(10.0*20.0), 0.2*sqrt(10.0*30.0), 0.5*sqrt(10.0*20.0), 20.0, 0.5*sqrt(20.0*30.0), 0.2*sqrt(10.0*30.0), 0.5*sqrt(20.0*30.0), 30.0,
    23, "UN_RW_7", 8.0, 0.7*sqrt(8.0*18.0), 0.3*sqrt(8.0*28.0), 0.7*sqrt(8.0*18.0), 18.0, 0.7*sqrt(18.0*28.0), 0.3*sqrt(8.0*28.0), 0.7*sqrt(18.0*28.0), 28.0,
    24, "UN_RW_8", 9.0, 0.7*sqrt(9.0*22.0), 0.3*sqrt(9.0*28.0), 0.7*sqrt(9.0*22.0), 22.0, 0.7*sqrt(22.0*28.0), 0.3*sqrt(9.0*28.0), 0.7*sqrt(22.0*28.0), 28.0,
    25, "UN_RW_9", 10.0, 0.7*sqrt(10.0*18.0), 0.3*sqrt(10.0*28.0), 0.7*sqrt(10.0*18.0), 18.0, 0.6*sqrt(18.0*28.0), 0.3*sqrt(10.0*28.0), 0.6*sqrt(18.0*28.0), 28.0
  )
}

make_pattern_dict <- function() {
  tidyr::crossing(
    Seq_Name = c("1_1", "8_7", "9_7"),
    Subg_Name = c("1_1", "3_5", "5_3")
  ) %>%
    dplyr::mutate(
      Pattern_Index = dplyr::row_number(),
      Seq_Ratio_1 = as.numeric(sub("([0-9]+)_.*", "\\1", Seq_Name)),
      Seq_Ratio_2 = as.numeric(sub(".*_([0-9]+)", "\\1", Seq_Name)),
      Subg_Ratio_1 = as.numeric(sub("([0-9]+)_.*", "\\1", Subg_Name)),
      Subg_Ratio_2 = as.numeric(sub(".*_([0-9]+)", "\\1", Subg_Name))
    )
}

make_cov_sas <- function(e11, e12, e13, e21, e22, e23, e31, e32, e33) {
  tibble::tibble(
    `_TYPE_` = c("MEAN", "COV", "COV", "COV"),
    `_NAME_` = c(".", "e1", "e2", "e3"),
    e1 = c(0, e11, e21, e31),
    e2 = c(0, e12, e22, e32),
    e3 = c(0, e13, e23, e33)
  )
}

extract_sigma_matrix <- function(tempdata) {
  if (is.matrix(tempdata)) return(tempdata)
  sigma <- tempdata %>%
    dplyr::filter(.data[["_TYPE_"]] == "COV") %>%
    dplyr::select(e1, e2, e3) %>%
    as.matrix()
  storage.mode(sigma) <- "double"
  sigma
}

repeat_df <- function(df, times) {
  if (times <= 0) return(df[0, , drop = FALSE])
  dplyr::bind_rows(replicate(times, df, simplify = FALSE))
}

get_sas_suffix <- function(dataset_source) {
  parts <- strsplit(as.character(dataset_source), "_", fixed = TRUE)[[1]]
  if (length(parts) >= 2L) parts[2L] else as.character(dataset_source)
}

alloc <- function(N, ratio) {
  x <- N * ratio / sum(ratio)
  n <- floor(x)
  k <- N - sum(n)
  if (k > 0) {
    id <- order(x - n, decreasing = TRUE)[seq_len(k)]
    n[id] <- n[id] + 1L
  }
  n
}

build_sas_like_dataset <- function(samplesize, tempdata, seed, seq_ratio, subg_ratio) {
  A <- make_design_A()
  B <- make_design_B()

  ns <- alloc(samplesize, seq_ratio)

  nA <- alloc(ns[1], subg_ratio)
  nB <- alloc(ns[2], subg_ratio)

  n <- c(nA, nB)

  X <- dplyr::bind_rows(
    repeat_df(A[A$SUBG == 0, ], n[1]),
    repeat_df(A[A$SUBG == 1, ], n[2]),
    repeat_df(B[B$SUBG == 0, ], n[3]),
    repeat_df(B[B$SUBG == 1, ], n[4])
  )
  
  sigma_mat <- extract_sigma_matrix(tempdata)
  set.seed(seed)
  eps <- as.data.frame(MASS::mvrnorm(n = samplesize, mu = c(0, 0, 0), Sigma = sigma_mat))
  names(eps) <- c("e1", "e2", "e3")
  
  err_long <- eps %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = c("e1", "e2", "e3"), names_to = "eps_name", values_to = "epsilon") %>%
    dplyr::arrange(row_id, factor(eps_name, levels = c("e1", "e2", "e3"))) %>%
    dplyr::select(epsilon)
  
  n_use <- min(nrow(X), nrow(err_long))
  X <- X[seq_len(n_use), , drop = FALSE]
  X$epsilon <- err_long$epsilon[seq_len(n_use)]
  X
}

prepare_model_data <- function(dat) {
  dat %>%
    dplyr::group_by(SUBJID) %>%
    dplyr::mutate(Visit = factor(dplyr::row_number(), levels = 1:3)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(SUBJID = factor(SUBJID))
}

add_all_outcomes <- function(dat, truth_table, beta_int, beta_base, beta_p2) {
  for (j in seq_len(nrow(truth_table))) {
    ds_name <- truth_table$Dataset_Source[j]
    out_name <- paste0("outcome_", get_sas_suffix(ds_name))
    
    dat[[out_name]] <- beta_int +
      beta_base * dat$BASELINE +
      beta_p2 * dat$PERIOD2 +
      truth_table$SUBGT_2_0[[j]] * dat$ACTIVE +
      truth_table$SUBGC_2_0[[j]] * dat$CARRYOVER +
      0 * dat$SUBG +
      truth_table$interaction_trt[[j]] * dat$ACTIVE * dat$SUBG +
      truth_table$interaction_cry[[j]] * dat$CARRYOVER * dat$SUBG +
      dat$epsilon
  }
  dat
}

build_task_grid <- function(sample_size_list) {
  tibble::tibble(
    Sample_Size = as.integer(sample_size_list),
    Task_ID = seq_along(sample_size_list)
  )
}

sanitize_sheet_name <- function(x) {
  x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
  substr(x, 1, 31)
}

bind_and_summarize_results <- function(task_results, truth_table) {
  model_names <- get_model_names()
  final_summarized <- setNames(vector("list", length(model_names)), model_names)
  
  ordered_ds <- unique(truth_table$Dataset_Source)
  
  raw_all <- setNames(vector("list", length(model_names)), model_names)
  for (nm in model_names) raw_all[[nm]] <- tibble::tibble()
  for (res in task_results) {
    for (nm in model_names) {
      if (!is.null(res[[nm]]) && nrow(res[[nm]]) > 0) {
        raw_all[[nm]] <- dplyr::bind_rows(raw_all[[nm]], res[[nm]])
      }
    }
  }
  
  truth_ref <- truth_table %>% 
    dplyr::select(label, Dataset_Source, SUBGT_2_0, SUBGC_2_0, interaction_trt, interaction_cry) %>% 
    dplyr::distinct()
  
  for (nm in model_names) {
    if (nrow(raw_all[[nm]]) > 0) {
      final_summarized[[nm]] <- raw_all[[nm]] %>%
        dplyr::left_join(truth_ref, by = "Dataset_Source") %>%
        dplyr::mutate(
          Target = dplyr::case_when(
            Effect %in% c("ACTIVE", "CARRYOVER", "Trt Effect (SUBG=0)", "Cry Effect (SUBG=0)", "Trt Effect (Naive)", "Cry Effect (Naive)") ~ "Primary",
            Effect %in% c("ACTIVE:SUBG", "CARRYOVER:SUBG", "Trt*SUBG Interaction", "Cry*SUBG Interaction") ~ "Interaction",
            TRUE ~ "Unknown"
          ),
          Truth = dplyr::case_when(
            Effect %in% c("ACTIVE", "Trt Effect (SUBG=0)", "Trt Effect (Naive)") ~ SUBGT_2_0,
            Effect %in% c("CARRYOVER", "Cry Effect (SUBG=0)", "Cry Effect (Naive)") ~ SUBGC_2_0,
            Effect %in% c("ACTIVE:SUBG", "Trt*SUBG Interaction") ~ interaction_trt,
            Effect %in% c("CARRYOVER:SUBG", "Cry*SUBG Interaction") ~ interaction_cry,
            TRUE ~ NA_real_
          ),
          Bias_raw = mean - Truth,
          Covered = ifelse(!is.na(Truth) & !is.na(CI_LOWERBOUND) & !is.na(CI_HIGHERBOUND) & 
                             CI_LOWERBOUND <= Truth & CI_HIGHERBOUND >= Truth, 1, 0),
          Signif = ifelse(!is.na(p_value) & p_value < 0.05, 1, 0),
          
          Dataset_Source = factor(Dataset_Source, levels = ordered_ds),
          Source_Group = dplyr::case_when(
            grepl("^Mixed_S", Dataset_Source) ~ 1,
            grepl("^Mixed_N", Dataset_Source) ~ 2,
            grepl("^Mixed_T", Dataset_Source) ~ 3,
            TRUE ~ 4
          )
        ) %>%
        dplyr::group_by(Matrix_Name, Seq_Name, Subg_Name, label, Dataset_Source, Source_Group, Sample_Size, Effect, Target, Truth) %>%
        dplyr::summarise(
          Mean_Est = mean(mean, na.rm = TRUE),
          SD_Est   = sd(mean, na.rm = TRUE),
          Bias     = mean(Bias_raw, na.rm = TRUE),
          Coverage = mean(Covered, na.rm = TRUE),
          Power    = mean(Signif, na.rm = TRUE),
          .groups  = "drop"
        ) %>%
        dplyr::arrange(Matrix_Name, Seq_Name, Subg_Name, Source_Group, Dataset_Source, Sample_Size, Effect) %>%
        dplyr::mutate(Dataset_Source = as.character(Dataset_Source)) %>%
        dplyr::select(Seq_Name, Subg_Name, label, Sample_Size, Dataset_Source, Effect, Target, Mean_Est, SD_Est, Bias, Coverage, Power, Truth, Matrix_Name)
    }
  }
  
  final_summarized
}

write_model_workbook <- function(df, file_path) {
  wb <- openxlsx::createWorkbook()
  
  if (nrow(df) == 0) {
    openxlsx::addWorksheet(wb, "EMPTY")
    openxlsx::writeData(wb, "EMPTY", data.frame(Message = "No rows"))
    openxlsx::saveWorkbook(wb, file = file_path, overwrite = TRUE)
    return(invisible(file_path))
  }
  
  matrix_names <- unique(df$Matrix_Name)
  total_sheets <- length(matrix_names)
  
  start_time <- Sys.time()
  message(sprintf("\n>>> Starting: %s", basename(file_path)))
  message(sprintf(">>> Total Sheet...", total_sheets))
  
  for (i in seq_along(matrix_names)) {
    mn <- matrix_names[i]
    sheet_name <- sanitize_sheet_name(mn)
    openxlsx::addWorksheet(wb, sheet_name)
    
    sheet_df <- df %>%
      dplyr::filter(Matrix_Name == mn) %>%
      dplyr::select(-Matrix_Name)
    
    openxlsx::writeData(wb, sheet_name, sheet_df)
    
    current_time <- Sys.time()
    elapsed_secs <- as.numeric(difftime(current_time, start_time, units = "secs"))
    avg_secs_per_sheet <- elapsed_secs / i
    remaining_secs <- avg_secs_per_sheet * (total_sheets - i)
    eta_time <- current_time + remaining_secs
    
    msg <- sprintf(
      "  [%d/%d] finish Sheet: %-10s | Current time: %s | ETA: %s",
      i, total_sheets, mn,
      format(current_time, "%H:%M:%S"),
      ifelse(i == total_sheets, "Finish!", format(eta_time, "%H:%M:%S"))
    )
    message(msg)
  }
  
  message(">>> Saving Works...")
  openxlsx::saveWorkbook(wb, file = file_path, overwrite = TRUE)
  message(sprintf(">>> Finish output: %s\n", basename(file_path)))
  
  invisible(file_path)
}

write_all_model_workbooks <- function(final_all, output_dir, file_suffix) {
  for (nm in names(final_all)) {
    if (!is.null(final_all[[nm]]) && nrow(final_all[[nm]]) > 0) {
      out_file <- file.path(output_dir, paste0("Result_", nm, "_", file_suffix, ".xlsx"))
      write_model_workbook(final_all[[nm]], out_file)
    }
  }
  invisible(TRUE)
}

run_one_task <- function(cov_row, sample_size, sim_times, truth_table, beta_int, beta_base, beta_p2, seq_ratio, subg_ratio, seq_name, subg_name) {
  current_cov <- make_cov_sas(
    cov_row$e11[[1]], cov_row$e12[[1]], cov_row$e13[[1]],
    cov_row$e21[[1]], cov_row$e22[[1]], cov_row$e23[[1]],
    cov_row$e31[[1]], cov_row$e32[[1]], cov_row$e33[[1]]
  )
  model_names <- get_model_names()
  task_out <- setNames(vector("list", length(model_names)), model_names)
  
  for (s in seq_len(sim_times)) {
    dat <- build_sas_like_dataset(
      samplesize = sample_size,
      tempdata = current_cov,
      seed = 123456L + s,
      seq_ratio = seq_ratio,
      subg_ratio = subg_ratio
    )
    
    stopifnot(nrow(dat) == sample_size * 3L)
    
    dat <- dat %>%
      dplyr::mutate(
        SUBJID = rep(seq_len(sample_size), each = 3L)
      )
    
    dat <- add_all_outcomes(dat, truth_table, beta_int, beta_base, beta_p2)
    dat_model <- prepare_model_data(dat)
    
    for (j in seq_len(nrow(truth_table))) {
      ds_name <- truth_table$Dataset_Source[j]
      response <- paste0("outcome_", get_sas_suffix(ds_name))
      fit_list <- fit_all_models(dat_model, response)
      
      for (nm in names(fit_list)) {
        fit_obj <- fit_list[[nm]]
        if (is.null(fit_obj)) next
        
        tmp <- extract_estimate_rows(fit_obj, nm) %>%
          dplyr::mutate(
            Dataset_Source = ds_name,
            Matrix_Name = cov_row$Matrix_Name[[1]],
            Sim_Run = s,
            Sample_Size = sample_size,
            Seq_Name = seq_name,
            Subg_Name = subg_name
          )
        task_out[[nm]] <- dplyr::bind_rows(task_out[[nm]], tmp)
      }
    }
  }
  task_out
}

run_task_index <- function(task_index, task_grid, cov_row, pattern_row, truth_table, sim_times, beta_int=20.5, beta_base=0.4, beta_p2=-2) {
  sample_size <- task_grid$Sample_Size[task_index]
  seq_ratio <- c(pattern_row$Seq_Ratio_1, pattern_row$Seq_Ratio_2)
  subg_ratio <- c(pattern_row$Subg_Ratio_1, pattern_row$Subg_Ratio_2)
  seq_name <- pattern_row$Seq_Name
  subg_name <- pattern_row$Subg_Name
  
  run_one_task(cov_row, sample_size, sim_times, truth_table, beta_int, beta_base, beta_p2, seq_ratio, subg_ratio, seq_name, subg_name)
}

make_truth_table <- function(path, sheet = 1) { 
  out <- readxl::read_xlsx(path = path, sheet = sheet) 
  required_cols <- c("label", "Dataset_Source", "SUBGT_2_0", "SUBGC_2_0", "interaction_trt", "interaction_cry") 
   if (!all(required_cols %in% names(out))){ 
     stop("Truth_Table is missing required columns. Please check: label, Dataset_Source, SUBGT_2_0, SUBGC_2_0, interaction_trt, interaction_cry")
   } 
  out 
}

truth_table <- make_truth_table(TRUTH_TABLE_PATH, TRUTH_TABLE_SHEET)

master_grid <- tidyr::crossing(
  Matrix_Index = 1:25,
  Pattern_Index = 1:9
)

if (TARGET_INDEX < 1 || TARGET_INDEX > nrow(master_grid)) {
  stop("FATAL ERROR: TARGET_INDEX is out of bounds (must be 1 to 225).")
}

current_matrix_idx <- master_grid$Matrix_Index[TARGET_INDEX]
current_pattern_idx <- master_grid$Pattern_Index[TARGET_INDEX]

full_cov_dict <- make_cov_dict()
current_cov_row <- full_cov_dict %>% dplyr::filter(Matrix_Index == current_matrix_idx)

full_pattern_dict <- make_pattern_dict()
current_pattern_row <- full_pattern_dict %>% dplyr::filter(Pattern_Index == current_pattern_idx)

CURRENT_MATRIX_NAME <- current_cov_row$Matrix_Name[1]
CURRENT_SEQ_NAME <- current_pattern_row$Seq_Name[1]
CURRENT_SUBG_NAME <- current_pattern_row$Subg_Name[1]

FILE_SUFFIX <- paste0("Sim_", CURRENT_MATRIX_NAME)
CURRENT_OUT_DIR <- file.path(BASE_DIR, CURRENT_SEQ_NAME, CURRENT_SUBG_NAME)
dir.create(CURRENT_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

task_grid <- build_task_grid(SAMPLE_SIZE_LIST)

available_cores <- parallel::detectCores(logical = USE_LOGICAL_CORES)
if (is.na(available_cores)) available_cores <- 1L
usable_cores <- max(1L, as.integer(available_cores) - as.integer(LEAVE_FREE_CORES))
N_CORES <- min(as.integer(REQUESTED_CORES), usable_cores, nrow(task_grid))
N_CORES <- max(1L, N_CORES)

main_objects_to_export <- c(
  "truth_table", "current_cov_row", "current_pattern_row", "task_grid", 
  "SIMULATION_TIMES", "alloc", "BETA_INT", "BETA_BASE", "BETA_P2", 
  "make_design_A", "make_design_B", "make_cov_sas", "extract_sigma_matrix", 
  "repeat_df", "build_sas_like_dataset", "get_sas_suffix", "prepare_model_data",
  "add_all_outcomes", "run_one_task", "run_task_index"
)

if (N_CORES > 1L) {
  cl <- parallel::makePSOCKcluster(N_CORES)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, varlist = c("BASE_DIR", "FUNCTION_FILE"), envir = environment())
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(dplyr)
      library(tidyr)
      library(MASS)
      library(mmrm)
    })
    source(FUNCTION_FILE)
    NULL
  })
  parallel::clusterExport(cl, varlist = main_objects_to_export, envir = environment())
  
  task_results <- parallel::parLapplyLB(cl, seq_len(nrow(task_grid)), function(task_index) {
    run_task_index(task_index, task_grid, current_cov_row, current_pattern_row, truth_table, SIMULATION_TIMES, BETA_INT, BETA_BASE, BETA_P2)
  })
} else {
  task_results <- lapply(seq_len(nrow(task_grid)), function(task_index) {
    run_task_index(task_index, task_grid, current_cov_row, current_pattern_row, truth_table, SIMULATION_TIMES, BETA_INT, BETA_BASE, BETA_P2)
  })
}

final_all <- bind_and_summarize_results(task_results, truth_table)

sub_final <- setNames(vector("list", length(names(final_all))), names(final_all))
for (nm in names(final_all)) {
  if (!is.null(final_all[[nm]]) && nrow(final_all[[nm]]) > 0) {
    sub_final[[nm]] <- final_all[[nm]] %>% 
      dplyr::select(-Seq_Name, -Subg_Name)
  }
}

if (SAVE_RDS) {
  output_filename <- paste0("ALL_Final_", FILE_SUFFIX, ".rds")
  saveRDS(sub_final, file = file.path(CURRENT_OUT_DIR, output_filename))
}
if (WRITE_XLSX) {
  write_all_model_workbooks(sub_final, CURRENT_OUT_DIR, FILE_SUFFIX)
}
