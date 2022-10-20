grouped_mean <- function(data, groupVars, var) {
  data$temp <- data[ , var]

  data_temp <-
    data %>%
    select(.data$Attributes, .data$Levels, .data$temp) %>%
    group_by(.data$Attributes, .data$Levels) %>%
    summarise(temp = mean(.data$temp),
              .groups = "drop")

  colnames(data_temp)[colnames(data_temp) == "temp"] <- var
  return(data_temp)
}


#' DCE samples
#'
#' A simulated sample of DCE data
#'
#' @format ## `dce_sample`
#' A data frame with 14 rows and 8 columns:
#' \describe{
#'   \item{Samples}{sample number}
#'   \item{Attributes}{attributes of each treatment}
#'   \item{Levels}{levels of each attribute}
#'   \item{no trt, trt1, trt2, trt3, trt4}{treatments}
#' }
#' @source simulated data
"dce_sample"

#' DCE samples
#'
#' A simulated sample of patient data for a DCE analysis
#'
#' @format ## `patient_sample`
#' A data frame with 210 rows and 9 columns:
#' \describe{
#'   \item{id}{patient id}
#'   \item{response}{data outcomes}
#'   \item{X1, X2}{covariates with their original values}
#'   \item{X1_1. X1_2, X2_1, X2_2, X2_3}{covariates converted to their corresponding dummy variables}
#' }
#' @source simulated data
"patient_sample"


#' Generate bootstrap samples
#'
#' @description
#' Generate bootstrap samples based on the probabilistic distribution of attributes of different treatments provided in the dca_samples data
#'
#' @param dce_samples samples of probabilistic distribution for the different treatments' attributes. The dataset needs to have the following columns:
#'                    Samples (sample number), Attributes (treatments' attributes), Levels (levels of different attributes), and one column per treatment. See dca_sample data in the package as an example of this dataset.
#' @param boot_n number of required bootstrap samples
#' @param na_zero logical (defaul is FALSE). If true, replace all the NA's in the dataset by zero
#' @param seed set seed for reproducibility purposes
#'
#' @return two datasets: 1. generated bootstrap samples, and 2. aggregated probabilistic distribution of attributes of each treatment
#' @export
#'
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom dplyr summarise
#' @importFrom rlang .data
#'
#' @examples
#' dce_boot_obj <- dce_bootstrap(dce_samples = dce_sample, boot_n = 100, seed = 1)
#' head(dce_boot_obj$boot_data)
#' head(dce_boot_obj$boot_prob)
dce_bootstrap <- function(dce_samples,
                          boot_n,
                          na_zero = FALSE,
                          seed = NULL) {

  if (! is.null(seed)) set.seed(seed = seed)

  if (na_zero) {
    dce_samples[is.na(dce_samples)] <- 0
  }

  trt_list <- colnames(dce_samples)[! colnames(dce_samples) %in% c("Attributes", "Levels", "Samples")]

  ## aggregate samples
  dce_agg <- dce_samples[dce_samples$Samples == 1 , colnames(dce_samples) != "Samples"]
  for (col in trt_list) {
    dce_agg_temp <-
      grouped_mean(data = dce_samples, groupVars = c("Attributes", "Levels"), var = col)

    dce_agg %<>%
      select(! col) %>%
      left_join(dce_agg_temp, by = c("Attributes", "Levels"))

    rm(dce_agg_temp)
  }

  dce_boot <- data.frame(Sample = rep(seq(1 : boot_n),
                                      times = length(unique(dce_agg$Attributes))),
                         Attributes = rep(unique(dce_agg$Attributes), each = boot_n))

  for (col in trt_list) {
    dce_boot_temp <-
      dce_boot %>%
      select(.data$Sample, .data$Attributes) %>%
      group_by(.data$Sample, .data$Attributes) %>%
      mutate(temp = sample(x = dce_agg$Levels[dce_agg$Attributes == .data$Attributes],
                           size = 1, replace = T,
                           prob = dce_agg[dce_agg$Attributes == .data$Attributes, col]))

    dce_boot %<>%
      left_join(dce_boot_temp, by = c("Sample", "Attributes"))
    colnames(dce_boot)[colnames(dce_boot) == "temp"] <- col

    rm(dce_boot_temp)
  }

  dce_boot %<>%
    arrange(.data$Sample, .data$Attributes)


  dce_probs <-
    dce_agg %>%
    select(.data$Attributes, .data$Levels) %>%
    group_by(.data$Attributes, .data$Levels) %>%
    summarise(temp_sim = mean(dce_boot[dce_boot$Attributes == .data$Attributes, trt_list[1]] == .data$Levels) * 100,
              .groups = 'drop') %>%
    left_join(dce_agg[ , c("Attributes", "Levels", trt_list[1])], by = c("Attributes", "Levels"))
  colnames(dce_probs)[colnames(dce_probs) %in% c(trt_list[1], "temp_sim")] <- c(paste0(trt_list[1], "_obs"),
                                                                                paste0(trt_list[1], "_sim"))

  for (col in trt_list[-1]) {
    dce_probs_temp <-
      dce_agg %>%
      select(.data$Attributes, .data$Levels) %>%
      group_by(.data$Attributes, .data$Levels, .drop = TRUE) %>%
      summarise(temp_sim = mean(dce_boot[dce_boot$Attributes == .data$Attributes, col] == .data$Levels) * 100,
                .groups = 'drop') %>%
      left_join(dce_agg[ , c("Attributes", "Levels", col)], by = c("Attributes", "Levels"))
    colnames(dce_probs_temp)[colnames(dce_probs_temp) %in% c(col, "temp_sim")] <- c(paste0(col, "_obs"),
                                                                                    paste0(col, "_sim"))

    dce_probs %<>%
      left_join(dce_probs_temp, by = c("Attributes", "Levels"))

    rm(dce_probs_temp)
  }

  return(list(boot_data = dce_boot, boot_prob = dce_probs))
}


# original_data = DCE_org_data
# dce_sim_data = dce_boot_obj$boot_data
# covariate_list = c("eaccess2", "eaccess3", "eseizure2",
#                    "eseizure3", "eseizure4", "emajorrisk2",
#                    "emajorrisk3", "emajorrisk4", "eminorrisk2",
#                    "eminorrisk3", "eminorrisk4", "eevidence2",
#                    "eevidence3", "ecost2", "ecost3", "ecost4")
# response = "choice"
# ID = "newid"
# trt_list = list(RNS = "eaccess3 + eseizure3 +
#                        emajorrisk3 + eminorrisk3 +
#                        eevidence2 + ecost3",
#                 MRgLITT = "eaccess2 + eseizure3 +
#                            emajorrisk3 + eminorrisk3 +
#                            eevidence2 - ecost2 -
#                            ecost3 - ecost4")


#' Generate the bootstrap samples amd fit a conditional logit model
#'
#' @description
#' Generate the bootstrap samples of attributes of different treatments, fit a conditional logit model on the patient data, and estimate effect of each covariates on the outcomes. The function provides SE of the effect estimates that are based on two uncertanities: 1. uncertanity of using sample of patients (not population), and 2. uncertanity in attributes of different treatments.
#'
#' @param original_data patient data, which need to include the folowing columns: id (patient id), response (a binary outcome), and covariates. See patient_sample data provided in the package for an example of this argument
#' @param dce_samples samples of probabilistic distribution for the different treatments' attributes. The dataset needs to have the following columns:
#'                    Samples (sample number), Attributes (treatments' attributes), Levels (levels of different attributes), and one column per treatment. See dca_sample data in the package as an example of this dataset.
#' @param dce_sim_data generated bootstrape samples of the attributes of different treatments (to estimate the second source of uncertainity in estimated SE's). These data can either be provided by user (e.g., by using the dce_bootstrap function provided in this package), or can be left null (the default value) and then this function will generate a dataset automatically
#' @param covariate_list list of covariates' name to be used in the conditional logit model
#' @param response response column name
#' @param ID id column name
#' @param coef_sim_n number of generated samples from the estimated sampling distribution of covariates' effect (to estimate the first source of uncertainity in estimated SE's)
#' @param trt_list list of treatment names
#' @param reference_trt name of reference treatment in estimating other treatments' uptake. If left null (default), no treatment will be used as the reference treatment.
#' @param boot_n number of required bootstrap samples
#' @param na_zero logical (defaul is FALSE). If true, replace all the NA's in the dataset by zero
#' @param seed set seed for reproducibility purposes
#'
#' @return a summary of the fitted conditional logit model, the estimated utility of each treatment, and the estimated uptake of each treatment relative to the refrence treatment
#' @export
#'
#'
#' @importFrom stats as.formula confint
#' @importFrom survival clogit strata coxph Surv
#' @importFrom MASS mvrnorm
#' @importFrom survey svycontrast
#' @importFrom dplyr select left_join group_by mutate arrange
#' @importFrom rlang .data
#'
#'
#' @examples
#' dce_CLogit_model <-
#'   dce_CLogitFit(original_data = patient_sample,
#'                 dce_samples = dce_sample,
#'                 covariate_list = c("X1_1", "X1_2", "X2_1", "X2_2", "X2_3"),
#'                 response = "response",
#'                 ID = "id",
#'                 coef_sim_n = 100,
#'                 trt_list = list(trt1 = "X1_1 + X1_2 + X2_1",
#'                                 trt2 = "X1_1 + X2_1 + X2_2 + X2_3"),
#'                 boot_n = 100, seed = 1)
#' dce_CLogit_model$model_summary
#' dce_CLogit_model$coefficient_estimates
#' dce_CLogit_model$reference_TRT_utility
#' dce_CLogit_model$`trt1 utility`
#' dce_CLogit_model$`trt1 uptake`
#'
dce_CLogitFit <- function(original_data,
                          dce_samples,
                          dce_sim_data = NULL,
                          covariate_list,
                          response,
                          ID,
                          coef_sim_n = 1000,
                          trt_list,
                          reference_trt = NULL,
                          boot_n = NULL,
                          na_zero = FALSE,
                          seed = NULL) {

  if (is.null(dce_sim_data)) {
    dce_boot_obj <- dce_bootstrap(dce_samples = dce_samples, boot_n = boot_n, seed = seed)
    dce_sim_data <- dce_boot_obj$boot_data
  }

  res <- list()

  model_formula <- as.formula(paste0(response, " ~ ",
                                     paste(covariate_list, collapse = " + "),
                                     paste0(" + strata(", ID, ")")))

  CLogitMod <- clogit(model_formula, data = original_data)

  res$model_summary <- summary(CLogitMod)
  res$coefficient_estimates <- cbind(summary(CLogitMod)$coef,
                                     confint(CLogitMod))


  ## for bootstrap simulations
  coef_sim <- mvrnorm(n = coef_sim_n, mu = CLogitMod$coefficients, Sigma = CLogitMod$var)

  ## no treatment utility
  if (is.null(reference_trt)) refTRT_call <- paste("-", paste(covariate_list, collapse = " - "))
  else refTRT_call <- reference_trt
  refTRT_util <- svycontrast(CLogitMod, parse(text = refTRT_call))
  res$reference_TRT_utility <- cbind(refTRT_util, confint(refTRT_util))

  ## treatments
  # i <- 1
  for (i in 1 : length(trt_list)) {
    # trt
    trt_temp <- names(trt_list)[i]
    # utility
    utility_temp <- svycontrast(CLogitMod, parse(text = trt_list[i]))
    res$temp <- cbind(utility_temp, confint(utility_temp))
    names(res)[which(names(res) == "temp")] <- paste0(trt_temp, " utility")
    # uptake
    if (as.character(parse(text = trt_list[i])) == as.character(parse(text = refTRT_call))) {
      res$temp <- c(1, NA, NA)
    }
    else {
      uptake_temp <- svycontrast(CLogitMod,
                                 parse(text = paste("exp(", trt_list[i], ") / (exp(",
                                                    trt_list[i], ") + exp(",
                                                    refTRT_call, "))")))
      res$temp <- cbind(uptake_temp, confint(uptake_temp))
    }
    names(res)[which(names(res) == "temp")] <- paste0(trt_temp, " uptake")
  }
  return(res)
}


