suppressPackageStartupMessages({
  library(dplyr)
  library(mmrm)
})

get_model_names <- function() {
  c("CS", "UN", "ANOVA", "CS_Naive", "UN_Naive", "ANOVA_Naive")
}

safe_fit <- function(expr) {
  tryCatch(expr, error = function(e) NULL)
}

#fit models CS, UN, ANOVA and ingnore subgroup CS, UN,ANOVA
fit_cs <- function(dat, response) {
  mmrm::mmrm(
    formula = as.formula(paste0(
      response,
      " ~ BASELINE + PERIOD2 + ACTIVE + SUBG + CARRYOVER + ",
      "ACTIVE:SUBG + CARRYOVER:SUBG + cs(Visit | SUBJID)"
    )),
    data = dat,
    reml = TRUE,
    control = mmrm::mmrm_control(
      method = "Kenward-Roger",
      vcov = "Kenward-Roger-Linear",
      n_cores = 1L
    )
  )
}

fit_un <- function(dat, response) {
  mmrm::mmrm(
    formula = as.formula(paste0(
      response,
      " ~ BASELINE + PERIOD2 + ACTIVE + SUBG + CARRYOVER + ACTIVE:SUBG + CARRYOVER:SUBG + us(Visit | SUBJID)"
    )),
    data = dat,
    reml = TRUE,
    control = mmrm::mmrm_control(
      method = "Kenward-Roger",
      vcov = "Kenward-Roger-Linear",
      n_cores = 1L
    )
  )
}

fit_anova <- function(dat, response) {
  stats::lm(
    as.formula(paste0(
      response,
      " ~ BASELINE + PERIOD2 + ACTIVE + SUBG + CARRYOVER + ACTIVE:SUBG + CARRYOVER:SUBG"
    )),
    data = dat
  )
}

fit_cs_naive <- function(dat, response) {
  mmrm::mmrm(
    formula = as.formula(paste0(
      response,
      " ~ BASELINE + PERIOD2 + ACTIVE + CARRYOVER + cs(Visit | SUBJID)"
    )),
    data = dat,
    reml = TRUE,
    control = mmrm::mmrm_control(
      method = "Kenward-Roger",
      vcov = "Kenward-Roger-Linear",
      n_cores = 1L
    )
  )
}

fit_un_naive <- function(dat, response) {
  mmrm::mmrm(
    formula = as.formula(paste0(
      response,
      " ~ BASELINE + PERIOD2 + ACTIVE + CARRYOVER + us(Visit | SUBJID)"
    )),
    data = dat,
    reml = TRUE,
    control = mmrm::mmrm_control(
      method = "Kenward-Roger",
      vcov = "Kenward-Roger-Linear",
      n_cores = 1L
    )
  )
}

fit_anova_naive <- function(dat, response) {
  stats::lm(
    as.formula(paste0(
      response,
      " ~ BASELINE + PERIOD2 + ACTIVE + CARRYOVER"
    )),
    data = dat
  )
}

fit_all_models <- function(dat, response) {
  list(
    CS = safe_fit(fit_cs(dat, response)),
    UN = safe_fit(fit_un(dat, response)),
    ANOVA = safe_fit(fit_anova(dat, response)),
    CS_Naive = safe_fit(fit_cs_naive(dat, response)),
    UN_Naive = safe_fit(fit_un_naive(dat, response)),
    ANOVA_Naive = safe_fit(fit_anova_naive(dat, response))
  )
}

get_coef_table_lme <- function(fit) {
  tt <- as.data.frame(summary(fit)$tTable)
  tt$Label <- rownames(tt)
  rownames(tt) <- NULL
  
  if (any(is.na(tt$DF)) || any(tt$DF <= 0)) {
    tt$DF <- pmax(tt$DF, 1)
  }
  
  tt %>%
    dplyr::transmute(
      Estimate = .data$Value,
      StdErr = .data$Std.Error,
      DF = .data$DF,
      tValue = .data$`t-value`,
      Probt = .data$`p-value`,
      Label = .data$Label
    )
}

get_coef_table_mmrm <- function(fit) {
  tt <- as.data.frame(summary(fit)$coefficients)
  tt$Label <- rownames(tt)
  rownames(tt) <- NULL

  tt %>%
    dplyr::transmute(
      Estimate = .data$Estimate,
      StdErr   = .data$`Std. Error`,
      DF       = .data$df,
      tValue   = .data$`t value`,
      Probt    = .data$`Pr(>|t|)`,
      Label    = .data$Label
    )
}

get_coef_table_lm <- function(fit) {
  tt <- as.data.frame(summary(fit)$coefficients)
  tt$Label <- rownames(tt)
  rownames(tt) <- NULL

  tt$DF <- stats::df.residual(fit)
  
  tt %>%
    dplyr::mutate(DF = stats::df.residual(fit)) %>%
    dplyr::transmute(
      Estimate = .data$Estimate,
      StdErr   = .data$`Std. Error`,
      DF       = .data$DF,
      tValue   = .data$`t value`,
      Probt    = .data$`Pr(>|t|)`,
      Label    = .data$Label
    )
}

get_coef_table <- function(fit) {
  if (inherits(fit, "mmrm")) {
    return(get_coef_table_mmrm(fit))
  }
  if (inherits(fit, "lme")) {
    return(get_coef_table_lme(fit))
  }
  if (inherits(fit, "lm")) {
    return(get_coef_table_lm(fit))
  }
  stop("Unsupported model class.")
}


extract_estimate_rows <- function(fit, model_name) {
  coef_tab <- get_coef_table(fit)
  
  if (model_name %in% c("CS", "UN", "ANOVA")) {
    coef_tab <- coef_tab %>%
      dplyr::filter(
        Label == "ACTIVE" | 
          Label == "CARRYOVER" | 
          (grepl("ACTIVE", Label) & grepl("SUBG", Label)) | 
          (grepl("CARRYOVER", Label) & grepl("SUBG", Label))
      ) %>%
      dplyr::mutate(
        Effect = dplyr::case_when(
          Label == "ACTIVE" ~ "Trt Effect (SUBG=0)",
          Label == "CARRYOVER" ~ "Cry Effect (SUBG=0)",
          grepl("ACTIVE", Label) & grepl("SUBG", Label) ~ "Trt*SUBG Interaction",
          grepl("CARRYOVER", Label) & grepl("SUBG", Label) ~ "Cry*SUBG Interaction",
          TRUE ~ Label
        )
      )
  } else {
    coef_tab <- coef_tab %>%
      dplyr::filter(Label %in% c("ACTIVE", "CARRYOVER")) %>%
      dplyr::mutate(
        Effect = dplyr::case_when(
          Label == "ACTIVE" ~ "Trt Effect (Naive)",
          Label == "CARRYOVER" ~ "Cry Effect (Naive)",
          TRUE ~ Label
        )
      )
  }
  
  coef_tab %>%
    dplyr::mutate(
      mean = .data$Estimate,
      se = .data$StdErr,
      CI_LOWERBOUND = .data$Estimate - stats::qt(0.975, df = pmax(.data$DF, 1)) * .data$StdErr,
      CI_HIGHERBOUND = .data$Estimate + stats::qt(0.975, df = pmax(.data$DF, 1)) * .data$StdErr,
      p_value = .data$Probt
    ) %>%
    dplyr::select(Effect, mean, se, CI_LOWERBOUND, CI_HIGHERBOUND, p_value)
}


