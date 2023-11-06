library(tidyverse)
library(magrittr)
library(useful)
library(scam)
library(mixtools)
library(scales)
library(ggthemes)
library(dr4pl)
library(drc)

#----
# LOAD THE RAW DATA
#----

# compound_metadata <- data.table::fread("data/PRISM_Oncology_Reference_23Q4_Compound_List.csv") 
analyte_meta = data.table::fread("data/RELEASE FILES/PRISM_Oncology_Reference_23Q4_Analyte_Meta.csv")
inst_meta <- data.table::fread("data/RELEASE FILES/PRISM_Oncology_Reference_23Q4_Inst_Meta.csv")
LMFI.long <- data.table::fread("data/RELEASE FILES/PRISM_Oncology_Reference_23Q4_LMFI.csv") 


# -----
# FILTER WELLS AND ANALYTES WITH LOW BEAD-COUNTS
# -----

# Analytes with less than 10 beads are filtered. Furthermore, if more than 25% of the analytes are filtered out the rest of the well is omitted.
LMFI.long <- LMFI.long %>%  
  dplyr::left_join(analyte_meta ) %>% 
  dplyr::left_join(inst_meta) %>% 
  dplyr::filter(is.finite(LMFI), is.finite(COUNT), !is.na(ccle_name))%>% 
  dplyr::group_by(prism_replicate, pert_well) %>%
  dplyr::filter(quantile(COUNT, probs = 0.25, na.rm = T) >= 10) %>% 
  dplyr::ungroup() %>%
  dplyr::filter(COUNT >= 10) 

# -----
# FILTER WELLS THAT CONTROL BARCODES DOESN'T COVER AT LEAST 25% OF THE RANGE OF CELL LINE BARCODES OR HAS SPEARMAN CORRELATION LESS THAN 0.5 WITH THE MEDIAN OF CONTROL BARCODES ACROSS NEGATIVE CONTROL WELLS.
# ----

CB.quantiles <- LMFI.long %>% 
  dplyr::filter(!is.na(ccle_name)) %>% 
  dplyr::group_by(prism_replicate, pert_well) %>%
  dplyr::arrange(LMFI) %>% 
  dplyr::mutate(quant = (1:n()) / n()) %>% 
  dplyr::filter(pool_id == "CTLBC") %>% 
  dplyr::distinct(prism_replicate, pert_well, pert_type,
                  ccle_name, LMFI, quant) %>% 
  dplyr::group_by(prism_replicate, pert_well) %>% 
  dplyr::mutate(range = max(quant, na.rm = T) - min(quant, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(range >= 0.25) 



CB.quantiles %<>% 
  dplyr::filter(pert_type == "ctl_vehicle") %>% 
  dplyr::group_by(prism_replicate, pert_well) %>% 
  dplyr::mutate(m = median(LMFI, na.rm = T)) %>%
  dplyr::group_by(prism_replicate) %>% 
  dplyr::mutate(mm = median(m, na.rm = T)) %>% 
  dplyr::group_by(prism_replicate, ccle_name) %>% 
  dplyr::summarise(m.LMFI = median(LMFI - m + mm, na.rm = T)) %>% 
  dplyr::ungroup() %>%
  dplyr::left_join(CB.quantiles) 


LMFI.long <- CB.quantiles %>% 
  dplyr::group_by(prism_replicate, pert_well, pert_type) %>%
  dplyr::summarise(concordance = cor(m.LMFI, LMFI, use = "p", method = "spearman")) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(concordance >= 0.5) %>%
  dplyr::distinct(prism_replicate, pert_well) %>%
  dplyr::left_join(LMFI.long) 

rm(CB.quantiles)

# ----
# COMPUTE THE REFERENCE VALUES FOR CONTROL BARCODES
# ----


REF = LMFI.long %>%  
  dplyr::filter(pool_id == "CTLBC") %>% 
  dplyr::filter(pert_type == "ctl_vehicle") %>% 
  dplyr::group_by(prism_replicate, profile_id) %>% 
  dplyr::mutate(m = median(LMFI, na.rm = T)) %>% 
  dplyr::group_by(prism_replicate) %>%  
  dplyr::mutate(m = m - median(m, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(LMFI = LMFI - m) %>% 
  dplyr::group_by(prism_replicate, analyte_id) %>%
  dplyr::summarise(LMFI.ref = median(LMFI, na.rm = T)) %>% 
  dplyr::ungroup()


# ----
# NORMALIZE USING SPLINE FITS
# ----

LMFI.long %<>% dplyr::left_join(REF) 

profiles = LMFI.long$profile_id %>% unique

LMFI.normalized.long = list()
ix = 1
for (profile in profiles) {
  print(ix)
  temp = LMFI.long %>%
    dplyr::filter(profile_id == profile) 
  
  
  g <- tryCatch(
    scam(
      LMFI.ref ~ s(LMFI, bs = "micx", k = 4),
      data = temp %>% dplyr::filter(is.finite(LMFI.ref)),
      optimizer = "nlm.fd"
    ),
    error = function(e)
      NULL
  )
 
  if (!is.null(g)) {
    p = predict(g, temp, se.fit = T)
    temp$LMFI.normalized = p$fit
    temp$LMFI.normalized.se = p$se.fit
  }

  LMFI.normalized.long[[ix]] = temp
  ix = ix + 1
}

LMFI.normalized.long %<>% dplyr::bind_rows()

rm(LMFI.long, REF, temp, p, g, ix, profile, profiles)

# ----
# QC TABLE
# ----

QC = LMFI.normalized.long %>%
  dplyr::filter(pert_type %in% c("ctl_vehicle", "trt_poscon")) %>%
  dplyr::group_by(compound_plate, prism_replicate, cellset, analyte_id) %>%
  dplyr::filter(is.finite(LMFI.normalized)) %>% 
  dplyr::summarise(
    error_rate =  min(PRROC::roc.curve(scores.class0 = LMFI.normalized,
                                       weights.class0 = pert_type == "ctl_vehicle",
                                       curve = TRUE)$curve[,1] + 1 -
                        PRROC::roc.curve(scores.class0 = LMFI.normalized,
                                         weights.class0 = pert_type == "ctl_vehicle",
                                         curve = TRUE )$curve[,2])/2,
    NC.median = median(LMFI.normalized[pert_type == "ctl_vehicle"], na.rm = T),
    NC.mad = mad(LMFI.normalized[pert_type == "ctl_vehicle"], na.rm = T),
    PC.median = median(LMFI.normalized[pert_type == "trt_poscon"], na.rm = T),
    PC.mad = mad( LMFI.normalized[pert_type == "trt_poscon"], na.rm = T)) %>%
  dplyr::mutate(DR = NC.median - PC.median,
                SSMD = DR / sqrt(NC.mad ^ 2 + PC.mad ^ 2)) %>%
  dplyr::mutate(PASS = (error_rate <= 0.05) & (DR > 2)) %>%
  dplyr::distinct()


# ----
# COMPUTE LOG-FOLD-CHANGES
# ----

LFC<- LMFI.normalized.long %>%  
  dplyr::left_join(QC %>% dplyr::select(analyte_id, prism_replicate,
                                        NC.median, PC.median, NC.mad, PC.mad, error_rate, PASS)) %>%
  dplyr::mutate(LFC = LMFI.normalized - NC.median) %>%
  dplyr::select(profile_id, analyte_id, cellset, screen, LFC, PASS) 



# ----
# FLAG AND FILTER OUTLIER WELLS
# ----

POOL.MEDIAN.VIABILITIES <- LFC %>% 
  dplyr::left_join(inst_meta) %>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::filter(pert_type == "trt_cp", PASS) %>% 
  dplyr::group_by(compound_plate, pool_id, pert_mfc_id, pert_dose, ccle_name) %>% 
  dplyr::mutate(LFC.median = median(LFC, na.rm = T)) %>% 
  dplyr::group_by(compound_plate, prism_replicate, pert_well, pool_id, pert_mfc_id, pert_dose) %>%
  dplyr::summarise(median.Viab = median(pmin(2,2^LFC), na.rm = T),
                   median.Ref.Viab = median(pmin(2, 2^LFC.median), na.rm = T)) %>%
  dplyr::ungroup()

dose_indices <- inst_meta %>% 
  dplyr::distinct(compound_plate, prism_replicate, pert_mfc_id, pert_dose) %>% 
  dplyr::group_by(compound_plate, prism_replicate, pert_mfc_id) %>% 
  dplyr::arrange(pert_dose) %>% 
  dplyr::mutate(dose_ix = 1:n()) %>%
  dplyr::ungroup()

filtered_pools <- dose_indices %>% 
  dplyr::left_join(dose_indices %>% 
                     dplyr::rename(pert_dose.prev = pert_dose) %>% 
                     dplyr::mutate(dose_ix = dose_ix + 1)) %>%   
  dplyr::left_join(dose_indices %>% 
                     dplyr::rename(pert_dose.next = pert_dose) %>% 
                     dplyr::mutate(dose_ix = dose_ix - 1)) %>%   
  dplyr::select(-dose_ix) %>% 
  dplyr::left_join(POOL.MEDIAN.VIABILITIES) %>% 
  dplyr::left_join(POOL.MEDIAN.VIABILITIES %>% 
                     dplyr::rename(pert_dose.prev = pert_dose) %>%
                     dplyr::group_by(prism_replicate, pool_id, pert_mfc_id, pert_dose.prev) %>%
                     dplyr::summarise(median.Viab.prev = median(median.Viab, na.rm = T)) %>%
                     dplyr::ungroup()) %>% 
  dplyr::left_join(POOL.MEDIAN.VIABILITIES %>% 
                     dplyr::rename(pert_dose.next = pert_dose) %>%
                     dplyr::group_by(prism_replicate, pool_id, pert_mfc_id, pert_dose.next) %>%
                     dplyr::summarise(median.Viab.next = median(median.Viab, na.rm = T)) %>%
                     dplyr::ungroup()) %>% 
  dplyr::mutate(outlier.flag = abs(median.Viab - median.Ref.Viab) > 0.5,
                mon.flag = pmin(median.Viab - median.Viab.prev, median.Viab.next - median.Viab, na.rm = T) > 0.5) %>%
  dplyr::filter(outlier.flag | mon.flag) %>% 
  dplyr::distinct(compound_plate, prism_replicate, pert_well, pool_id)


filtered_wells <- LFC %>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::left_join(inst_meta) %>% 
  dplyr::filter(pert_type == "trt_cp", pool_id != "CTLBC") %>% 
  dplyr::distinct(compound_plate, prism_replicate, pert_well, pool_id) %>% 
  dplyr::left_join(filtered_pools %>% dplyr::mutate(flag = TRUE)) %>% 
  dplyr::mutate(flag = ifelse( is.na(flag), 0, 1)) %>% 
  dplyr::group_by(prism_replicate, pert_well) %>% 
  dplyr::summarise(flagged.fraction = mean(flag, na.rm = T)) %>% 
  dplyr::filter(flagged.fraction > 0.5) %>%
  dplyr::ungroup() 
  

LFC.FILTERED <- LFC %>% 
  dplyr::left_join(inst_meta) %>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::anti_join(filtered_wells) %>% 
  dplyr::anti_join(filtered_pools) %>% 
  dplyr::filter(PASS) %>% 
  dplyr::distinct(profile_id, analyte_id, cellset, screen, LFC)
  
rm(dose_indices, POOL.MEDIAN.VIABILITIES)

# -----
# APPLY COMBAT FOR THE POOL EFFECTS 
# -----

apply_combat <- function(Y) {
  # create "cond" column to be used as "batches"
  df <- Y %>%
    dplyr::distinct(profile_id, pool_id, cellset, analyte_id, LFC) %>%
    tidyr::unite(cond, pool_id, profile_id, cellset, sep = "::") %>%
    dplyr::filter(is.finite(LFC))
  
  # calculate means and sd's of each condition
  batch <- df$cond
  m <- rbind(df$LFC,
             rnorm(length(df$LFC),
                   mean =  mean(df$LFC, na.rm = TRUE),
                   sd = sd(df$LFC, na.rm = TRUE)))
  
  # use ComBat to align means and sd's of conditions
  combat <- sva::ComBat(dat = m, batch = batch) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(analyte_id = df$analyte_id, cond = df$cond) %>%
    dplyr::rename(LFC_cb = V1) %>%
    dplyr::mutate(pool_id = stringr::word(cond, 1, sep = stringr::fixed("::")),
                  profile_id = stringr::word(cond, 2, sep = stringr::fixed("::")),
                  cellset = stringr::word(cond, 3, sep = stringr::fixed("::"))) %>%
    dplyr::select(-cond, -V2)
  
  combat_corrected <- Y %>%
    dplyr::left_join(combat, by = c("profile_id", "analyte_id", "pool_id", "cellset")) %>%
    .$LFC_cb
  
  return(combat_corrected)
}

LFC.FILTERED %<>%  
  dplyr::left_join(analyte_meta %>% dplyr::distinct(pool_id, analyte_id, cellset, screen)) %>% 
  dplyr::left_join(inst_meta %>% dplyr::distinct(profile_id, compound_plate, pert_type, pert_mfc_id, pert_dose)) %>%
  dplyr::filter(is.finite(LFC), pool_id != "CTLBC", !is.na(pool_id), pert_type == "trt_cp") %>% 
  tidyr::unite(col = "condition", pert_mfc_id, pert_dose, compound_plate, screen, sep = "::", remove = F) %>% 
  split(.$condition) %>% 
  purrr::map_dfr(~dplyr::mutate(.x, LFC_cb = apply_combat(.))) %>%
  dplyr::select(-condition)


# -----
# COLLAPSE THE REPLICATES
# -----


LFC.collapsed <- LFC.FILTERED %>% 
  dplyr::filter(is.finite(LFC_cb)) %>%
  dplyr::group_by(analyte_id, pool_id, cellset, screen, compound_plate, pert_mfc_id, pert_dose) %>% 
  dplyr::filter(n() > 1) %>% 
  dplyr::summarise(LFC_cb = median(LFC_cb, na.rm = T)) %>%
  dplyr::ungroup() 

# ----
# FIT DOSE RESPONSE CURVES
# ----


compute_auc <- function(LL, UL, Inflection, Slope, md, MD) {
  f1 = function(x) pmax(pmin((UL + (LL - UL)/(1 + (2^x/Inflection)^Slope)), 1, na.rm = T), 0, na.rm = T)
  return(tryCatch(integrate(f1, log2(md), log2(MD))$value/(log2(MD/md)),
                  error = function(e) {print(e); NA}))
}
compute_log_ic50 <- function(LL, UL, Inflection, Slope, md, MD) {
  if((LL >= 0.5) | (UL <= 0.5)) {
    return(NA)
  } else {
    f1 = function(x) (UL + (LL - UL)/(1 + (2^x/Inflection)^Slope)- 0.5)
    return(tryCatch(uniroot(f1, c(log2(md), log2(MD)))$root,
                    error = function(x) NA))
  }
}
compute_MSE_MAD <- function(LFC,  UL, LL,  Slope, Inflection, FC_column="FC", dose_column="dose") {
  #Slope = -Slope
  mse_compute <- LFC %>% 
    dplyr::filter(is.finite(.[[FC_column]]),is.finite(.[[dose_column]]) ) %>% ## in case there are some na values accidentally passed.
    dplyr::mutate(FC.pred = UL  + (LL -UL )/(1 + (.[[dose_column]]/Inflection)^Slope) ) %>% 
    dplyr::mutate(squared_deviation = (.[[FC_column]]-FC.pred)^2, abs_deviation = abs(.[[FC_column]]-FC.pred)) %>%
    dplyr::summarise(mse = mean(squared_deviation,na.rm=T), mad= median(abs_deviation,na.rm=T))
  return (mse_compute)
}

get_best_fit <- function(FC, dose, UL_low=0.8, UL_up=1.01, slope_decreasing=TRUE){
  ## get best fit among different attempts at fitting, and see if this fit works sufficiently well to be reported.
  LFC_filtered = tibble(pert_dose = dose, FC = FC) %>%
    tidyr::drop_na()
  
  var_data <- LFC_filtered$FC %>% var()
  riemann_AUC <- pmin(LFC_filtered$FC,1) %>% mean() ## mean fold-change after rounding FC to 1.
  
  all_fits.df <- fit_4param_drc(LFC_filtered, "pert_dose",  var_data, 
                                UL_low, UL_up, slope_decreasing)
  
  res.df  <- data.frame(successful_fit=FALSE, AUC_Riemann = riemann_AUC) ## default return value if fit is unsuccessful
  
  
  if (nrow(all_fits.df)>0){
    all_fits.df %<>%
      dplyr::filter(!is.na(frac_var_explained))
    } ## remove entries with NA in variance explained 
  
  if (nrow(all_fits.df)>0){
    res.df <- all_fits.df %>%
      slice_max(frac_var_explained, n = 1, with_ties = FALSE) %>%  ## return best fit. if tied, return first of the ties
      dplyr::mutate(successful_fit = TRUE, AUC_Riemann = as.numeric(riemann_AUC) ) ## fit has to be at least as good as predicting just the mean of the data to be called successful
  }
  
  return (res.df)
}

fit_4param_drc <- function(LFC_filtered, dose_var,  var_data, UL_low=0.8, UL_up=1.01, slope_decreasing=TRUE) {
  #fits a number of alternate models  to the DRC and passes the results to the calling function (which chooses the best fit.)
  
  # UL low is the lowerbound of UL we pass to the optimizer and UL_up is the upper bound of UL that we pass to the optimizer
  # fomat of output will be:-
  # results.df <- data.frame("fit_name"=character(),"Lower_Limit"=double(),
  #                          "Upper_Limit"=double(), 
  #                          "Slope"=double(),
  #                          "Inflection"=double(), 
  #                          "MSE"=double(), "MAD" =double(),
  #                          "frac_var_explained"=double(),
  #                          "Input_Parameters"=character())

  
  results.df <- list(); ix = 1
  
  slope_bound <- ifelse(slope_decreasing, 1e-5, Inf)  # bound the slopes by default unless passed another option
  
  dr4pl_initL <- tryCatch(dr4pl(as.formula(paste("FC ~ ", dose_var)), data = LFC_filtered,
                                init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.005),
                                method.init = "logistic",
                                lowerl = c(UL_low, -Inf, -Inf, 0),
                                upperl = c(UL_up, Inf, slope_bound, 1.01)),
                          error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )
  if (dr4pl_initL$convergence==FALSE){
    if (!is.null(dr4pl_initL$dr4pl.robust)) { 
      dr4pl_initL <- dr4pl_initL$dr4pl.robust
    }
  }
  
  if (dr4pl_initL$convergence){
    mse_df <- compute_MSE_MAD(LFC_filtered, dr4pl_initL$parameters[[1]], dr4pl_initL$parameters[[4]],
                              dr4pl_initL$parameters[[3]], dr4pl_initL$parameters [[2]],
                              "FC", dose_var)
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_initL", 
                                Lower_Limit = as.numeric(dr4pl_initL$parameters [[4]]),
                                Upper_Limit = as.numeric(dr4pl_initL$parameters [[1]]),
                                Slope = as.numeric(dr4pl_initL$parameters [[3]]),
                                Inflection = as.numeric(dr4pl_initL$parameters [[2]]),
                                MSE = mse_df$mse, MAD = mse_df$mad, frac_var_explained = 1-mse_df$mse/var_data,
                                Input_Parameters = "constrained|init_logistic")
    ix = ix + 1 
    
  }
  
  
  dr4pl_initM_optNM <- tryCatch(dr4pl(as.formula(paste("FC ~ ", dose_var)), data = LFC_filtered,
                                      init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.01),
                                      method.init = "Mead",
                                      lowerl = c(UL_low, -Inf, -Inf, 0),
                                      upperl = c(UL_up, Inf, slope_bound, 1.01),
                                      method.optim="Nelder-Mead"),
                                error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )
  if (dr4pl_initM_optNM$convergence==FALSE){
    if (!is.null(dr4pl_initM_optNM$dr4pl.robust)) { 
      dr4pl_initM_optNM <- dr4pl_initM_optNM$dr4pl.robust
    }
  }
  
  if (dr4pl_initM_optNM$convergence){
    mse_df <- compute_MSE_MAD(LFC_filtered, dr4pl_initM_optNM$parameters[[1]], dr4pl_initM_optNM$parameters[[4]],
                              dr4pl_initM_optNM$parameters[[3]], dr4pl_initM_optNM$parameters [[2]],
                              "FC", dose_var)
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_initM_optNM", 
                                Lower_Limit = as.numeric(dr4pl_initM_optNM$parameters [[4]]),
                                Upper_Limit = as.numeric(dr4pl_initM_optNM$parameters [[1]]),
                                Slope = as.numeric(dr4pl_initM_optNM$parameters [[3]]),
                                Inflection = as.numeric(dr4pl_initM_optNM$parameters [[2]]),
                                MSE = mse_df$mse, MAD = mse_df$mad, frac_var_explained = 1-mse_df$mse/var_data,
                                Input_Parameters = "constrained|init_Mead|optim_Nelder-Mead")
    ix = ix + 1 
  }
  
  
  dr4pl_initL_optB <- tryCatch(dr4pl(as.formula(paste("FC ~ ", dose_var)), data = LFC_filtered,
                                     init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.01),
                                     method.init = "logistic",
                                     lowerl = c(UL_low, -Inf, -Inf, 0),
                                     upperl = c(UL_up, Inf, slope_bound, 1.01),
                                     method.optim="BFGS"),
                               error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )
  if (dr4pl_initL_optB$convergence==FALSE){
    if (!is.null(dr4pl_initL_optB$dr4pl.robust)) { 
      dr4pl_initL_optB <- dr4pl_initL_optB$dr4pl.robust
    }
  }
  
  if (dr4pl_initL_optB$convergence){
    mse_df <- compute_MSE_MAD(LFC_filtered, dr4pl_initL_optB$parameters[[1]], dr4pl_initL_optB$parameters[[4]],
                              dr4pl_initL_optB$parameters[[3]], dr4pl_initL_optB$parameters [[2]],
                              "FC", dose_var)
    
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_optB", 
                                Lower_Limit = as.numeric(dr4pl_initL_optB$parameters [[4]]),
                                Upper_Limit = as.numeric(dr4pl_initL_optB$parameters [[1]]),
                                Slope = as.numeric(dr4pl_initL_optB$parameters [[3]]),
                                Inflection = as.numeric(dr4pl_initL_optB$parameters [[2]]),
                                MSE = mse_df$mse, MAD = mse_df$mad, frac_var_explained = 1-mse_df$mse/var_data,
                                Input_Parameters = "constrained|init_logistic|optim_BFGS")
    ix = ix + 1 
    
  }
  
  dr4pl_lossHuber<- tryCatch(dr4pl(as.formula(paste("FC ~ ", dose_var)), data = LFC_filtered,
                                   init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.01),
                                   method.robust="Huber",
                                   lowerl = c(UL_low, -Inf, -Inf, 0),
                                   upperl = c(UL_up, Inf, slope_bound, 1.01)),
                             error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )
  if (dr4pl_lossHuber$convergence==FALSE){
    if (!is.null(dr4pl_lossHuber$dr4pl.robust)) { 
      dr4pl_lossHuber <- dr4pl_lossHuber$dr4pl.robust
    }
  }
  
  
  if (dr4pl_lossHuber$convergence){
    mse_df <- compute_MSE_MAD(LFC_filtered, dr4pl_lossHuber$parameters[[1]], dr4pl_lossHuber$parameters[[4]],
                              dr4pl_lossHuber$parameters[[3]], dr4pl_lossHuber$parameters [[2]],
                              "FC", dose_var)
    
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_lossHuber", 
                                Lower_Limit = as.numeric(dr4pl_lossHuber$parameters [[4]]),
                                Upper_Limit = as.numeric(dr4pl_lossHuber$parameters [[1]]),
                                Slope = as.numeric(dr4pl_lossHuber$parameters [[3]]),
                                Inflection = as.numeric(dr4pl_lossHuber$parameters [[2]]),
                                MSE = mse_df$mse, MAD = mse_df$mad, frac_var_explained = 1-mse_df$mse/var_data,
                                Input_Parameters = "constrained|loss_Huber")
    ix = ix + 1 
    
  }
  
  ### add in original default drc into pipeline just to compare. ######
  dr4pl_unconstrained <- tryCatch(dr4pl(dose = LFC_filtered[[dose_var]], response = LFC_filtered$FC,
                                        init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.3),
                                        method.init = "logistic",
                                        lowerl = c(0.99, -Inf, -Inf, 0),
                                        upperl = c(1.01, Inf, Inf, 1.01)),
                                  error = function(e) {print(e); return(NA)})
  
  # if it fits and doesn't converge grab robust fit
  if (!all(is.na(dr4pl_unconstrained))) {
    if (!dr4pl_unconstrained$convergence) {
      dr4pl_unconstrained <- dr4pl_unconstrained$dr4pl.robust 
    }
  }
  # get parameters
  param <- tryCatch(dr4pl_unconstrained$parameters, error = function(e) return(NA))
  if (!all(is.na(param))){
    if(as.numeric(dr4pl_unconstrained$parameters [[3]])<slope_bound){ ### while slope bound is not passed to this last optimizer, we do not accept a solution not within the bound
      mse_df <- compute_MSE_MAD(LFC_filtered, dr4pl_unconstrained$parameters[[1]], dr4pl_unconstrained$parameters[[4]],
                                dr4pl_unconstrained$parameters[[3]], dr4pl_unconstrained$parameters [[2]],
                                "FC", dose_var)
      results.df[[ix]] <- tibble( fit_name = "original", 
                                  Lower_Limit = as.numeric(dr4pl_unconstrained$parameters [[4]]),
                                  Upper_Limit = as.numeric(dr4pl_unconstrained$parameters [[1]]),
                                  Slope = as.numeric(dr4pl_unconstrained$parameters [[3]]),
                                  Inflection = as.numeric(dr4pl_unconstrained$parameters [[2]]),
                                  MSE = mse_df$mse, MAD = mse_df$mad, frac_var_explained = 1-mse_df$mse/var_data,
                                  Input_Parameters = "unconstrained_optim_dr4pl")
      ix = ix + 1 
    }
  }
  
  ### two additional manual ways of initializing the fit to cover all bases.
  dr4pl_initMan_optB <- tryCatch(dr4pl(as.formula(paste("FC ~ ", dose_var)), data = LFC_filtered,
                                       init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_2 = 1.1*min(LFC_filtered[[dose_var]]),
                                                                      theta_3= -.5, theta_4 = 0.01),
                                       lowerl = c(UL_low, -Inf, -Inf, 0),
                                       upperl = c(UL_up, Inf, slope_bound, 1.01),
                                       method.optim="BFGS"),
                                 error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )
  if (dr4pl_initMan_optB$convergence==FALSE){
    if (!is.null(dr4pl_initMan_optB$dr4pl.robust)) { 
      dr4pl_initMan_optB <- dr4pl_initMan_optB$dr4pl.robust
    }
  }
  
  
  if (dr4pl_initMan_optB$convergence){
    mse_df <- compute_MSE_MAD(LFC_filtered, dr4pl_initMan_optB$parameters[[1]], dr4pl_initMan_optB$parameters[[4]],
                              dr4pl_initMan_optB$parameters[[3]], dr4pl_initMan_optB$parameters [[2]],
                              "FC", dose_var)
    
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_initMan_optB", 
                                Lower_Limit = as.numeric(dr4pl_initMan_optB$parameters [[4]]),
                                Upper_Limit = as.numeric(dr4pl_initMan_optB$parameters [[1]]),
                                Slope = as.numeric(dr4pl_initMan_optB$parameters [[3]]),
                                Inflection = as.numeric(dr4pl_initMan_optB$parameters [[2]]),
                                MSE = mse_df$mse, MAD = mse_df$mad, frac_var_explained = 1-mse_df$mse/var_data,
                                Input_Parameters = "constrained|init_Manual|optim_BFGS")
    ix = ix + 1 
  }
  
  dr4pl_initMan_optNM <- tryCatch(dr4pl(as.formula(paste("FC ~ ", dose_var)), data = LFC_filtered,
                                        init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_2 = 8*min(LFC_filtered[[dose_var]]),
                                                                       theta_3= -3, theta_4 = 0.01),
                                        lowerl = c(UL_low, -Inf, -Inf, 0),
                                        upperl = c(UL_up, Inf, slope_bound, 1.01),
                                        method.optim="Nelder-Mead"),
                                  error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )
  if (dr4pl_initMan_optNM$convergence==FALSE){
    if (!is.null(dr4pl_initMan_optNM$dr4pl.robust)) { 
      dr4pl_initMan_optNM <- dr4pl_initMan_optNM$dr4pl.robust
    }
  }
  
  if (dr4pl_initMan_optNM$convergence){
    mse_df <- compute_MSE_MAD(LFC_filtered, dr4pl_initMan_optNM$parameters[[1]], dr4pl_initMan_optNM$parameters[[4]],
                              dr4pl_initMan_optNM$parameters[[3]], dr4pl_initMan_optNM$parameters [[2]],
                              "FC", dose_var)
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_initMan_optNM", 
                                Lower_Limit = as.numeric(dr4pl_initMan_optNM$parameters [[4]]),
                                Upper_Limit = as.numeric(dr4pl_initMan_optNM$parameters [[1]]),
                                Slope = as.numeric(dr4pl_initMan_optNM$parameters [[3]]),
                                Inflection = as.numeric(dr4pl_initMan_optNM$parameters [[2]]),
                                MSE = mse_df$mse, MAD = mse_df$mad, frac_var_explained = 1-mse_df$mse/var_data,
                                Input_Parameters = "constrained|init_Manual|optim_Nelder-Mead")
    ix = ix + 1 
  }
  

  drc_model <-  tryCatch(drc::drm(as.formula(paste("FC ~ ", dose_var)), data=LFC_filtered,
                                  fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                  lowerl = c(-slope_bound,0.0, UL_low, -Inf),upperl = c(Inf,1.01,UL_up, Inf)),
                         error = function(e)
                         {return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))})
  # "slope" in drc package is -ve of slope in dr4pl package
  
  
  if (drc_model$fit$convergence){
    
    mse_df <- compute_MSE_MAD(LFC_filtered, as.numeric(drc_model$coefficients [[3]]), as.numeric(drc_model$coefficients [[2]]),
                              -as.numeric(drc_model$coefficients [[1]]), as.numeric(drc_model$coefficients [[4]]),
                              "FC", dose_var)
    # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.
    
    results.df[[ix]] <- tibble( fit_name = "drc_drm",
                                Lower_Limit = as.numeric(drc_model$coefficients [[2]]),
                                Upper_Limit = as.numeric(drc_model$coefficients [[3]]),
                                Slope = -as.numeric(drc_model$coefficients [[1]]),
                                Inflection = as.numeric(drc_model$coefficients [[4]]),
                                MSE = mse_df$mse, MAD = mse_df$mad, frac_var_explained = 1-mse_df$mse/var_data,
                                Input_Parameters = "constrained-drc")
    ix = ix + 1
  }
  
  
  return (dplyr::bind_rows(results.df))
}


# !!! THIS TAKES A LONG TIME
DRC <-  LFC.FILTERED %>%
 dplyr::filter(is.finite(LFC_cb)) %>%
 dplyr::group_by(pert_mfc_id, pert_dose, compound_plate, analyte_id, cellset) %>%
 dplyr::filter(n() > 1) %>%
 dplyr::group_by(pert_mfc_id, compound_plate, analyte_id, cellset) %>%
 dplyr::summarise(get_best_fit(pmin(2^LFC_cb,2), pert_dose)) %>%
 dplyr::ungroup() 




DRC <- inst_meta %>% 
  dplyr::group_by(screen, compound_plate, pert_mfc_id) %>% 
  dplyr::summarise(md = min(pert_dose, na.rm = T),
                   MD = max(pert_dose, na.rm = T)) %>%
  dplyr::ungroup() %>% 
  dplyr::inner_join(DRC) %>% 
  dplyr::filter(successful_fit) %>% 
  dplyr::mutate(Lower_Limit = as.numeric(Lower_Limit), 
                Upper_Limit = as.numeric(Upper_Limit), 
                Inflection = as.numeric(Inflection), 
                Slope = as.numeric(Slope)) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(AUC = compute_auc(Lower_Limit, Upper_Limit, Inflection, Slope, md, MD),
                log2.IC50 = compute_log_ic50(Lower_Limit, Upper_Limit, Inflection, Slope, md, MD)) 

# Sanity check
DRC %>% 
  dplyr::ungroup() %>% 
  dplyr::sample_n(100) %>% 
  dplyr::distinct(AUC_Riemann, AUC, log2.IC50) %>% 
  psych::pairs.panels()


# -----
# WRITE PORTAL FILES
# -----

M <- LFC.FILTERED %>% 
  dplyr::filter(is.finite(LFC_cb)) %>%
  dplyr::group_by(pert_mfc_id, pert_dose, compound_plate, analyte_id, cellset) %>%
  dplyr::filter(n() > 1) %>% 
  dplyr::ungroup() %>%
  dplyr::filter(is.finite(LFC_cb)) %>%
  dplyr::left_join(analyte_meta) %>%
  dplyr::left_join(inst_meta) %>% 
  tidyr::unite(cn, pert_mfc_id, pert_dose, compound_plate, replicate, sep = "::") %>% 
  reshape2::acast(cn ~ depmap_id, value.var = "LFC_cb",
                  fun.aggregate = median)


Model.annotations <- tibble(cell_line_name = colnames(M),
                            index = 0:(dim(M)[2]-1)) 



Condition.annotations <- tibble(cn = rownames(M), index = 0:(dim(M)[1]-1)) %>%
  dplyr::mutate(pert_mfc_id = word(cn, sep = fixed("::")),
                pert_dose = as.numeric(word(cn, 2, sep = fixed("::"))),
                compound_plate = word(cn, 3, sep = fixed("::")),
                replicate = word(cn, 4, sep = fixed("::"))) %>% 
  dplyr::left_join(inst_meta) %>% 
  dplyr::distinct(index, pert_mfc_id, pert_dose, replicate, compound_plate) %>%
  #dplyr::left_join(compound_list, by = c("pert_mfc_id" = "IDs")) %>%
  dplyr::rename(compound_name = pert_mfc_id,
                dose = pert_dose) %>%
  dplyr::mutate(compound_name = paste0("BRD:", compound_name)) %>% 
  dplyr::mutate(masked = "FALSE") %>% 
  dplyr::group_by(compound_name, dose) %>% 
  dplyr::arrange(compound_plate, replicate) %>% dplyr::mutate(replicate = 1:n()) %>% 
  dplyr::ungroup() %>% dplyr::select(-compound_plate) %>% 
  dplyr::distinct() 


colnames(M) <- 0:(dim(M)[2]-1)
rownames(M) <- 0:(dim(M)[1]-1)


DRC %>% 
  dplyr::left_join(analyte_meta) %>%  
  dplyr::rename(cell_line_name = depmap_id,
                compound_name = pert_mfc_id,
                ec50 = Inflection,
                lower_asymptote = Lower_Limit,
                upper_asymptote = Upper_Limit,
                slope = Slope,
  ) %>% 
  dplyr::mutate(compound_name = paste0("BRD:", compound_name)) %>% 
  dplyr::group_by(cell_line_name, compound_name, cellset) %>%  
  dplyr::top_n(1, frac_var_explained) %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct(cell_line_name, compound_name, ec50, lower_asymptote, upper_asymptote, slope) %>% 
  write_csv("data/PORTAL FILES/Response_curves.csv")



DRC %>% 
  dplyr::left_join(analyte_meta) %>%  
  dplyr::filter(is.finite(AUC)) %>% 
  dplyr::mutate(compound_name = paste0("BRD:", pert_mfc_id)) %>% 
  reshape2::acast(compound_name ~ depmap_id,
                  value.var = "AUC", fun.aggregate = median) %>%
  write.csv("data/PORTAL FILES/AUC_matrix.csv")


DRC %>% 
  dplyr::left_join(analyte_meta) %>%  
  dplyr::filter(is.finite(log2.IC50)) %>% 
  dplyr::mutate(compound_name = paste0("BRD:", pert_mfc_id)) %>% 
  reshape2::acast(compound_name ~ depmap_id,
                  value.var = "log2.IC50", fun.aggregate = median) %>%
  write.csv("data/PORTAL FILES/log2_IC50_matrix.csv")

(2^M) %>% 
  write.csv("data/PORTAL FILES/Viability_matrix.csv")

Condition.annotations %>% 
  write_csv("data/PORTAL FILES/Condition_annotations.csv")

Model.annotations %>%  
  write_csv("data/PORTAL FILES/Model_annotations.csv")


# -----
# WRITE THE RELEASE FILES
# -----


LFC.FILTERED %>% 
  dplyr::group_by(pert_mfc_id, pert_dose, compound_plate, analyte_id, cellset) %>% 
  dplyr::mutate(n.PASS = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(profile_id, analyte_id, LFC, LFC_cb, n.PASS, screen, compound_plate, cellset) %>%
  dplyr::distinct() %>% 
  write_csv("data/RELEASE FILES/PRISM_Oncology_Reference_23Q4_LFC.csv")

LFC.collapsed %>% 
  dplyr::filter(is.finite(LFC_cb)) %>% 
  dplyr::left_join(analyte_meta) %>% 
  tidyr::unite(cn, pert_mfc_id, pert_dose, compound_plate, screen, sep = "::") %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "LFC_cb", fun.aggregate = median) %>% 
  write.csv("data/RELEASE FILES/PRISM_Oncology_Reference_23Q4_LFC_Matrix.csv")

DRC %>% 
  dplyr::select(pert_mfc_id, compound_plate, screen, md, MD, analyte_id, cellset,
                fit_name, Lower_Limit, Upper_Limit, Slope, Inflection,
                MSE, MAD, frac_var_explained, Input_Parameters, AUC, log2.IC50, AUC_Riemann) %>%
  dplyr::distinct() %>% 
  write_csv("data/RELEASE FILES/PRISM_Oncology_Reference_23Q4_Dose_Response_Parameters.csv")

QC %>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::select(screen, compound_plate, prism_replicate, cellset, analyte_id, pool_id, ccle_name, depmap_id,
                error_rate, NC.median, NC.mad, PC.median, PC.mad, DR, SSMD, PASS) %>%
  dplyr::distinct() %>%
  write_csv("data/RELEASE FILES/PRISM_Oncology_Reference_23Q4_QC_Table.csv")


failed_pools <- LFC %>%
  dplyr::left_join(inst_meta) %>%
  dplyr::left_join(analyte_meta) %>% 
  dplyr::distinct(screen, compound_plate, prism_replicate, pool_id, pert_well) 

dplyr::inner_join(failed_pools, filtered_wells) %>%
  dplyr::select(-flagged.fraction) %>% 
  dplyr::bind_rows(dplyr::inner_join(failed_pools, filtered_pools)) %>% 
  distinct() %>% 
  write_csv("data/RELEASE FILES/PRISM_Oncology_Reference_23Q4_Filtered_Pools.csv")



# ----------

