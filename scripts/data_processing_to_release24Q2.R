library(tidyverse)
library(magrittr)
library(useful)
library(scam)
library(mixtools)
library(scales)
library(ggthemes)
library(parallel)
library(dr4pl)
library(drc)




#----
# LOAD THE RAW DATA
#----

compound_list <- data.table::fread("data/24Q2/PRISMOncologyReferenceCompoundList.csv")
analyte_meta = data.table::fread("data/24Q2/PRISMOncologyReferenceAnalyteMeta2.csv")
inst_meta <- data.table::fread("data/24Q2/PRISMOncologyReferenceInstMeta.csv")
LMFI.long <- data.table::fread("data/24Q2/PRISMOncologyReferenceLMFI.csv") 




# -----
# FILTER WELLS AND ANALYTES WITH LOW BEAD-COUNTS
# -----

# Analytes with less than 10 beads are filtered. Furthermore, if more than 25% of the analytes are filtered out the rest of the well is omitted.
LMFI.long <- LMFI.long %>%  
  dplyr::left_join(analyte_meta) %>% 
  dplyr::left_join(inst_meta) %>% 
  dplyr::filter(is.finite(LMFI), is.finite(COUNT), !is.na(depmap_id))%>% 
  dplyr::group_by(prism_replicate, pert_well) %>%
  dplyr::filter(quantile(COUNT, probs = 0.25, na.rm = T) >= 10) %>% 
  dplyr::ungroup() %>%
  dplyr::filter(COUNT >= 10) 



# -----
# IMPUTE MISSING NOISE FLOOR ESTIMATES WITH THE LOWEST READ-OUT IN EACH WELL
# -----

# Sanity check
LMFI.long %>% 
  dplyr::group_by(CompoundPlate, cellset, replicate, pert_type,
                  prism_replicate, pert_well, NF.mean, NF.sd) %>% 
  dplyr::summarize(NF.emp = min(LMFI, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>%
  ggplot(aes(x = NF.emp, y = NF.mean)) +
  geom_point(aes(color = pert_type), size = .25) +
  geom_density2d() +
  geom_abline(color = "red") +
  facet_grid(cellset ~ CompoundPlate)

LMFI.long %>% 
  dplyr::distinct(CompoundPlate, screen, replicate, pert_type, 
                  NF.mean, NF.sd, pert_well, prism_replicate, cellset) %>% 
  ggplot() +
  geom_density(aes(x = NF.mean, color = CompoundPlate,
                   group = prism_replicate)) +
  facet_grid(screen ~ cellset) 

LMFI.long %>%
  dplyr::distinct(cellset, CompoundPlate, pert_type, pert_well,
                  NF.mean, NF.sd) %>%
  ggplot() +
  geom_boxplot(aes(x = CompoundPlate, y = NF.sd, color = pert_type)) +
  facet_grid(. ~ cellset )


# Imputation
LMFI.long <- LMFI.long %>% 
  dplyr::group_by(prism_replicate, pert_well) %>% 
  dplyr::mutate(NF.mean = ifelse(is.na(NF.mean), min(LMFI, na.rm = TRUE), NF.mean),
                NF.sd = ifelse(is.na(NF.mean), 0.25, NF.sd)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(LMFI = log2(2^LMFI - 2^NF.mean)) %>% 
  dplyr::mutate(LMFI = ifelse(is.finite(LMFI), LMFI, NF.sd))


# -----
# FILTER WELLS THAT CONTROL BARCODES DOESN'T COVER AT LEAST 25% OF THE RANGE OF CELL LINE BARCODES OR HAS SPEARMAN CORRELATION LESS THAN 0.5 WITH THE MEDIAN OF CONTROL BARCODES ACROSS NEGATIVE CONTROL WELLS.
# ----

CB.quantiles <- LMFI.long %>% 
  dplyr::group_by(prism_replicate, pert_well) %>%
  dplyr::arrange(LMFI) %>% 
  dplyr::mutate(quant = (1:n()) / n()) %>% 
  dplyr::filter(pool_id == "CTLBC") %>% 
  dplyr::distinct(prism_replicate, pert_well, pert_type,
                  depmap_id, LMFI, quant) %>% 
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
  dplyr::group_by(prism_replicate, depmap_id) %>% 
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
  dplyr::group_by(screen, CompoundPlate, prism_replicate, cellset, analyte_id) %>%
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
  dplyr::distinct() %>%
  dplyr::ungroup()


QC %>% 
  dplyr::distinct(CompoundPlate, prism_replicate, cellset, analyte_id, PASS) %>%  
  dplyr::group_by(CompoundPlate, cellset, analyte_id) %>% 
  dplyr::summarise(n.PASS = sum(PASS, na.rm = T)) %>% 
  ggplot() +
  geom_bar(aes(x = CompoundPlate, fill = as.factor(n.PASS))) +
  facet_wrap(cellset ~ ., ncol = 2) +
  labs(fill = "Passing Rep.") +
  coord_flip()
  
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
  dplyr::group_by(CompoundPlate, pool_id, SampleID, pert_dose, depmap_id) %>% 
  dplyr::mutate(LFC.median = median(LFC, na.rm = T)) %>% 
  dplyr::group_by(CompoundPlate, prism_replicate, pert_well, pool_id, SampleID, pert_dose) %>%
  dplyr::summarise(median.Viab = median(pmin(2,2^LFC), na.rm = T),
                   median.Ref.Viab = median(pmin(2, 2^LFC.median), na.rm = T)) %>%
  dplyr::ungroup()

dose_indices <- inst_meta %>% 
  dplyr::distinct(CompoundPlate, prism_replicate, SampleID, pert_dose) %>% 
  dplyr::group_by(CompoundPlate, prism_replicate, SampleID) %>% 
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
                     dplyr::group_by(prism_replicate, pool_id, SampleID, pert_dose.prev) %>%
                     dplyr::summarise(median.Viab.prev = median(median.Viab, na.rm = T)) %>%
                     dplyr::ungroup()) %>% 
  dplyr::left_join(POOL.MEDIAN.VIABILITIES %>% 
                     dplyr::rename(pert_dose.next = pert_dose) %>%
                     dplyr::group_by(prism_replicate, pool_id, SampleID, pert_dose.next) %>%
                     dplyr::summarise(median.Viab.next = median(median.Viab, na.rm = T)) %>%
                     dplyr::ungroup()) %>% 
  dplyr::mutate(outlier.flag = abs(median.Viab - median.Ref.Viab) > 0.5,
                mon.flag = pmin(median.Viab - median.Viab.prev, median.Viab.next - median.Viab, na.rm = T) > 0.5) %>%
  dplyr::filter(outlier.flag | mon.flag) %>% 
  dplyr::distinct(CompoundPlate, prism_replicate, pert_well, pool_id)


filtered_wells <- LFC %>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::left_join(inst_meta) %>% 
  dplyr::filter(pert_type == "trt_cp", pool_id != "CTLBC") %>% 
  dplyr::distinct(CompoundPlate, prism_replicate, pert_well, pool_id) %>% 
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


# Lost profiles
inst_meta %>% 
  dplyr::filter(pert_type == "trt_cp") %>% 
  dplyr::distinct(profile_id, CompoundPlate, SampleID, pert_dose, prism_replicate) %>% 
  dplyr::anti_join(LFC.FILTERED) %>% 
  dplyr::count(CompoundPlate, SampleID, pert_dose) %>% 
  dplyr::filter(n > 1) %>% 
  dplyr::left_join(compound_list) %>%
  head

# -----
# APPLY COMBAT FOR THE POOL EFFECTS 
# -----

apply_combat <- function(Y) {
  # create "cond" column to be used as "batches"
  df <- Y %>%
    dplyr::distinct(profile_id, pool_id, cellset, analyte_id, LFC) %>%
    tidyr::unite(cond, pool_id, profile_id, cellset, sep = "::") %>%
    dplyr::filter(is.finite(LFC))
  
  set.seed(25) # ! setting the seed 
  
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

LFC.FILTERED <- LFC.FILTERED %>%  
  dplyr::left_join(analyte_meta %>% dplyr::distinct(pool_id, analyte_id, cellset, screen)) %>% 
  dplyr::left_join(inst_meta %>% dplyr::distinct(profile_id, CompoundPlate, pert_type, SampleID, pert_dose)) %>%
  dplyr::filter(is.finite(LFC), pool_id != "CTLBC", !is.na(pool_id), pert_type == "trt_cp") %>% 
  tidyr::unite(col = "condition", SampleID, pert_dose, CompoundPlate, screen, sep = "::", remove = F) %>% 
  split(.$condition) %>% 
  purrr::map_dfr(~dplyr::mutate(.x, LFC_cb = apply_combat(.))) %>%
  dplyr::select(-condition)


# -----
# COLLAPSE THE REPLICATES
# -----

LFC.collapsed <- LFC.FILTERED %>% 
  dplyr::filter(is.finite(LFC_cb)) %>%
  dplyr::group_by(analyte_id, pool_id, cellset, screen, CompoundPlate, SampleID, pert_dose) %>% 
  dplyr::filter(n() > 1) %>% 
  dplyr::summarise(LFC_cb = median(LFC_cb, na.rm = T)) %>%
  dplyr::ungroup() 

# ----
# FIT DOSE RESPONSE CURVES
# ----

# Auxilary functions
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
compute_MSE_MAD <- function(FC, dose,  UL, LL,  Slope, Inflection) {
  FC.pred = UL  + (LL -UL )/(1 + (dose/Inflection)^Slope)
  residuals = FC - FC.pred
  return(list(mse = mean(residuals^2), mad = median(abs(residuals))))
}
get_best_fit <- function(FC, dose, UL_low=0.8, UL_up=1.01, slope_decreasing=TRUE) {
  require(dr4pl)
  require(drc)
  require(tidyverse)
  require(magrittr)
  
  # Fits a number of alternate models  to the DRC and chooses the best fit.
  
  # UL low is the lowerbound of UL we pass to the optimizer and UL_up is the upper bound of UL that we pass to the optimizer
  # fomat of output will be:-
  # results.df <- data.frame("fit_name"=character(),"Lower_Limit"=double(),
  #                          "Upper_Limit"=double(),
  #                          "Slope"=double(),
  #                          "Inflection"=double(),
  #                          "MSE"=double(), "MAD" =double(),
  #                          "frac_var_explained"=double())
  
  
  
  riemann_AUC <- mean(pmin(1,FC)) ## mean fold-change after rounding FC to 1.
  var_data = var(FC)
  slope_bound <- ifelse(slope_decreasing, 1e-5, Inf)  # bound the slopes by default unless passed another option
  
  
  results.df <- list(); ix = 1
  
  
  # FIT 1 ---
  drc_model <-  tryCatch(drc::drm(FC ~ dose, data= data.frame(FC = FC, dose = dose),
                                  fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                  lowerl = c(-slope_bound,0.0, UL_low, -Inf),upperl = c(Inf,1.01,UL_up, Inf)), # ??
                         error = function(e)
                         {return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))})
  # "slope" in drc package is -ve of slope in dr4pl package
  
  
  if (drc_model$fit$convergence){
    mse_mad <- compute_MSE_MAD(FC, dose, as.numeric(drc_model$coefficients[[3]]), as.numeric(drc_model$coefficients[[2]]),
                               -as.numeric(drc_model$coefficients[[1]]), as.numeric(drc_model$coefficients[[4]]))
    # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.
    
    results.df[[ix]] <- tibble( fit_name = "drc_drm_constrained",
                                Lower_Limit = as.numeric(drc_model$coefficients[[2]]),
                                Upper_Limit = as.numeric(drc_model$coefficients[[3]]),
                                Slope = -as.numeric(drc_model$coefficients[[1]]),
                                Inflection = as.numeric(drc_model$coefficients[[4]]),
                                MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }
  
  # FIT 2 ---
  drc_model <-  tryCatch(drc::drm(FC ~ dose, data= data.frame(FC = FC, dose = dose),
                                  fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                  lowerl = c(-Inf,0.0, UL_low, -Inf),upperl = c(Inf,1.01,UL_up, Inf)), # ??
                         error = function(e)
                         {return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))})
  # "slope" in drc package is -ve of slope in dr4pl package
  
  
  if (drc_model$fit$convergence){
    mse_mad <- compute_MSE_MAD(FC, dose, as.numeric(drc_model$coefficients[[3]]), as.numeric(drc_model$coefficients[[2]]),
                               -as.numeric(drc_model$coefficients[[1]]), as.numeric(drc_model$coefficients[[4]]))
    # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.
    
    results.df[[ix]] <- tibble( fit_name = "drc_drm_unconstrained",
                                Lower_Limit = as.numeric(drc_model$coefficients[[2]]),
                                Upper_Limit = as.numeric(drc_model$coefficients[[3]]),
                                Slope = -as.numeric(drc_model$coefficients[[1]]),
                                Inflection = as.numeric(drc_model$coefficients[[4]]),
                                MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }
  
  # FIT 3 ---
  dr4pl_initMan_optNM <- tryCatch(dr4pl(dose, FC,
                                        init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_2 = 8*min(dose),
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
    mse_mad <- compute_MSE_MAD(FC, dose, dr4pl_initMan_optNM$parameters[[1]], dr4pl_initMan_optNM$parameters[[4]],
                               dr4pl_initMan_optNM$parameters[[3]], dr4pl_initMan_optNM$parameters[[2]])
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_initMan_constrained_optNM",
                                Lower_Limit = as.numeric(dr4pl_initMan_optNM$parameters[[4]]),
                                Upper_Limit = as.numeric(dr4pl_initMan_optNM$parameters[[1]]),
                                Slope = as.numeric(dr4pl_initMan_optNM$parameters[[3]]),
                                Inflection = as.numeric(dr4pl_initMan_optNM$parameters[[2]]),
                                MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }
  
  # FIT 4 ---
  dr4pl_unconstrained <- tryCatch(dr4pl(dose, FC,
                                        init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.3),
                                        method.init = "logistic",
                                        lowerl = c(0.99, -Inf, -Inf, 0),
                                        upperl = c(1.01, Inf, Inf, 1.01)),
                                  error = function(e) {print(e); return(NA)})
  
  if (!all(is.na(dr4pl_unconstrained))) {
    if (!dr4pl_unconstrained$convergence) {
      dr4pl_unconstrained <- dr4pl_unconstrained$dr4pl.robust
    }
  }
  
  
  param <- tryCatch(dr4pl_unconstrained$parameters, error = function(e) return(NA))
  if (!all(is.na(param))){
    if(as.numeric(dr4pl_unconstrained$parameters[[3]])<slope_bound){ ### while slope bound is not passed to this last optimizer, we do not accept a solution not within the bound
      mse_mad <- compute_MSE_MAD(FC, dose, dr4pl_unconstrained$parameters[[1]], dr4pl_unconstrained$parameters[[4]],
                                 dr4pl_unconstrained$parameters[[3]], dr4pl_unconstrained$parameters[[2]])
      results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_unconstrained",
                                  Lower_Limit = as.numeric(dr4pl_unconstrained$parameters[[4]]),
                                  Upper_Limit = as.numeric(dr4pl_unconstrained$parameters[[1]]),
                                  Slope = as.numeric(dr4pl_unconstrained$parameters[[3]]),
                                  Inflection = as.numeric(dr4pl_unconstrained$parameters[[2]]),
                                  MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
  }
  
  # FIT 5 ---
  dr4pl_initL <- tryCatch(dr4pl(dose, FC,
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
    mse_mad <- compute_MSE_MAD(FC,dose, dr4pl_initL$parameters[[1]], dr4pl_initL$parameters[[4]],
                               dr4pl_initL$parameters[[3]], dr4pl_initL$parameters[[2]])
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_constrained",
                                Lower_Limit = as.numeric(dr4pl_initL$parameters[[4]]),
                                Upper_Limit = as.numeric(dr4pl_initL$parameters[[1]]),
                                Slope = as.numeric(dr4pl_initL$parameters[[3]]),
                                Inflection = as.numeric(dr4pl_initL$parameters[[2]]),
                                MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }
  
  # Choose the best fit among the successful fits---
  results.df <- dplyr::bind_rows(results.df)
  
  # if (nrow(results.df)>0){
  #   results.df <-  dplyr::filter(results.df, frac_var_explained > 0)
  # }
  
  if (nrow(results.df)>0){
    results.df <- results.df %>%
      dplyr::arrange(desc(frac_var_explained)) %>%
      head(1) %>%
      dplyr::mutate(successful_fit = TRUE, AUC_Riemann = as.numeric(riemann_AUC) )
  }else{
    results.df  <- data.frame(successful_fit=FALSE, AUC_Riemann = riemann_AUC)
  }
  
  return (results.df)
}

## Sequential implementation (note this may take several hours to run, please see the alternative below) ----

# DRC <-  LFC.FILTERED %>%
#  dplyr::filter(is.finite(LFC_cb)) %>%
#  dplyr::group_by(SampleID, pert_dose, CompoundPlate, analyte_id, cellset) %>%
#  dplyr::filter(n() > 1) %>%
#  dplyr::group_by(SampleID, CompoundPlate, analyte_id, cellset) %>%
#  dplyr::summarise(get_best_fit(pmin(2^LFC_cb,2), pert_dose)) %>%
#  dplyr::ungroup() 



# Parallel implementation for for curve fitting ----

LFC.parallel <- LFC.FILTERED %>%
  dplyr::filter(is.finite(LFC_cb)) %>%
  dplyr::group_by(SampleID, pert_dose, CompoundPlate, analyte_id, cellset) %>%
  dplyr::filter(n() > 1) %>% 
  dplyr::ungroup() %>% 
  tidyr::unite(px, SampleID, CompoundPlate, cellset, remove = FALSE) 

LFC.parallel <- lapply(unique(LFC.parallel$px), function(x) dplyr::filter(LFC.parallel, px == x))

# Encapsulating get_best_fit and compute_MSE_MAD functions for parallelization
f <- function(df) { 
  require(tidyverse)
  compute_MSE_MAD <- function(FC, dose,  UL, LL,  Slope, Inflection) {
    FC.pred = UL  + (LL -UL )/(1 + (dose/Inflection)^Slope)
    residuals = FC - FC.pred
    return(list(mse = mean(residuals^2), mad = median(abs(residuals))))
  }
  get_best_fit <- function(FC, dose, UL_low=0.8, UL_up=1.01, slope_decreasing=TRUE) {
    require(dr4pl)
    require(drc)
    require(tidyverse)
    require(magrittr)
    
    # Fits a number of alternate models  to the DRC and chooses the best fit.
    
    # UL low is the lowerbound of UL we pass to the optimizer and UL_up is the upper bound of UL that we pass to the optimizer
    # fomat of output will be:-
    # results.df <- data.frame("fit_name"=character(),"Lower_Limit"=double(),
    #                          "Upper_Limit"=double(), 
    #                          "Slope"=double(),
    #                          "Inflection"=double(), 
    #                          "MSE"=double(), "MAD" =double(),
    #                          "frac_var_explained"=double())
    
    
    
    riemann_AUC <- mean(pmin(1,FC)) ## mean fold-change after rounding FC to 1.
    var_data = var(FC)
    slope_bound <- ifelse(slope_decreasing, 1e-5, Inf)  # bound the slopes by default unless passed another option
    
    
    results.df <- list(); ix = 1
    
    
    # FIT 1 ---------
    drc_model <-  tryCatch(drc::drm(FC ~ dose, data= data.frame(FC = FC, dose = dose),
                                    fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                    lowerl = c(-slope_bound,0.0, UL_low, -Inf),upperl = c(Inf,1.01,UL_up, Inf)), # ?? 
                           error = function(e)
                           {return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))})
    # "slope" in drc package is -ve of slope in dr4pl package
    
    
    if (drc_model$fit$convergence){
      mse_mad <- compute_MSE_MAD(FC, dose, as.numeric(drc_model$coefficients[[3]]), as.numeric(drc_model$coefficients[[2]]),
                                 -as.numeric(drc_model$coefficients[[1]]), as.numeric(drc_model$coefficients[[4]]))
      # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.
      
      results.df[[ix]] <- tibble( fit_name = "drc_drm_constrained",
                                  Lower_Limit = as.numeric(drc_model$coefficients[[2]]),
                                  Upper_Limit = as.numeric(drc_model$coefficients[[3]]),
                                  Slope = -as.numeric(drc_model$coefficients[[1]]),
                                  Inflection = as.numeric(drc_model$coefficients[[4]]),
                                  MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
    
    # FIT 2 ---------
    drc_model <-  tryCatch(drc::drm(FC ~ dose, data= data.frame(FC = FC, dose = dose),
                                    fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                    lowerl = c(-Inf,0.0, UL_low, -Inf),upperl = c(Inf,1.01,UL_up, Inf)), # ?? 
                           error = function(e)
                           {return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))})
    # "slope" in drc package is -ve of slope in dr4pl package
    
    
    if (drc_model$fit$convergence){
      mse_mad <- compute_MSE_MAD(FC, dose, as.numeric(drc_model$coefficients[[3]]), as.numeric(drc_model$coefficients[[2]]),
                                 -as.numeric(drc_model$coefficients[[1]]), as.numeric(drc_model$coefficients[[4]]))
      # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.
      
      results.df[[ix]] <- tibble( fit_name = "drc_drm_unconstrained",
                                  Lower_Limit = as.numeric(drc_model$coefficients[[2]]),
                                  Upper_Limit = as.numeric(drc_model$coefficients[[3]]),
                                  Slope = -as.numeric(drc_model$coefficients[[1]]),
                                  Inflection = as.numeric(drc_model$coefficients[[4]]),
                                  MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
    
    # FIT 3 ------
    dr4pl_initMan_optNM <- tryCatch(dr4pl(dose, FC,
                                          init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_2 = 8*min(dose),
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
      mse_mad <- compute_MSE_MAD(FC, dose, dr4pl_initMan_optNM$parameters[[1]], dr4pl_initMan_optNM$parameters[[4]],
                                 dr4pl_initMan_optNM$parameters[[3]], dr4pl_initMan_optNM$parameters[[2]])
      
      results.df[[ix]] <- tibble( fit_name = "dr4pl_initMan_constrained_optNM",
                                  Lower_Limit = as.numeric(dr4pl_initMan_optNM$parameters[[4]]),
                                  Upper_Limit = as.numeric(dr4pl_initMan_optNM$parameters[[1]]),
                                  Slope = as.numeric(dr4pl_initMan_optNM$parameters[[3]]),
                                  Inflection = as.numeric(dr4pl_initMan_optNM$parameters[[2]]),
                                  MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
    
    # FIT 4 -----
    dr4pl_unconstrained <- tryCatch(dr4pl(dose, FC,
                                          init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.3),
                                          method.init = "logistic",
                                          lowerl = c(0.99, -Inf, -Inf, 0),
                                          upperl = c(1.01, Inf, Inf, 1.01)),
                                    error = function(e) {print(e); return(NA)})
    
    if (!all(is.na(dr4pl_unconstrained))) {
      if (!dr4pl_unconstrained$convergence) {
        dr4pl_unconstrained <- dr4pl_unconstrained$dr4pl.robust
      }
    }
    
    
    param <- tryCatch(dr4pl_unconstrained$parameters, error = function(e) return(NA))
    if (!all(is.na(param))){
      if(as.numeric(dr4pl_unconstrained$parameters[[3]])<slope_bound){ ### while slope bound is not passed to this last optimizer, we do not accept a solution not within the bound
        mse_mad <- compute_MSE_MAD(FC, dose, dr4pl_unconstrained$parameters[[1]], dr4pl_unconstrained$parameters[[4]],
                                   dr4pl_unconstrained$parameters[[3]], dr4pl_unconstrained$parameters[[2]])
        results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_unconstrained",
                                    Lower_Limit = as.numeric(dr4pl_unconstrained$parameters[[4]]),
                                    Upper_Limit = as.numeric(dr4pl_unconstrained$parameters[[1]]),
                                    Slope = as.numeric(dr4pl_unconstrained$parameters[[3]]),
                                    Inflection = as.numeric(dr4pl_unconstrained$parameters[[2]]),
                                    MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
        ix = ix + 1
      }
    }
    
    # FIT 5 ----
    dr4pl_initL <- tryCatch(dr4pl(dose, FC,
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
      mse_mad <- compute_MSE_MAD(FC,dose, dr4pl_initL$parameters[[1]], dr4pl_initL$parameters[[4]],
                                 dr4pl_initL$parameters[[3]], dr4pl_initL$parameters[[2]])
      
      results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_constrained",
                                  Lower_Limit = as.numeric(dr4pl_initL$parameters[[4]]),
                                  Upper_Limit = as.numeric(dr4pl_initL$parameters[[1]]),
                                  Slope = as.numeric(dr4pl_initL$parameters[[3]]),
                                  Inflection = as.numeric(dr4pl_initL$parameters[[2]]),
                                  MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
    
    # Choose the best fit among the successful fits------
    results.df <- dplyr::bind_rows(results.df) 
    
    # if (nrow(results.df)>0){
    #   results.df <-  dplyr::filter(results.df, frac_var_explained > 0)
    # }
    
    if (nrow(results.df)>0){
      results.df <- results.df %>%
        dplyr::arrange(desc(frac_var_explained)) %>% 
        head(1) %>% 
        dplyr::mutate(successful_fit = TRUE, AUC_Riemann = as.numeric(riemann_AUC) ) 
    }else{
      results.df  <- data.frame(successful_fit=FALSE, AUC_Riemann = riemann_AUC) 
    }
    
    return (results.df)
  }
  
  df %>%
    dplyr::group_by(SampleID, CompoundPlate, analyte_id, cellset) %>%
    dplyr::summarise(get_best_fit(pmin(2^LFC_cb,2), pert_dose)) %>%
    dplyr::ungroup() 
}

# Create a cluster
cl <- makeCluster(detectCores() - 1)
# Fit the curves
DRC <- parLapply(cl, LFC.parallel, f)
# Stop the cluster
stopCluster(cl)
rm(LFC.parallel, f)
DRC %<>% dplyr::bind_rows()


DRC <- inst_meta %>% 
  dplyr::group_by(screen, CompoundPlate, SampleID) %>% 
  dplyr::summarise(md = min(pert_dose, na.rm = T),
                   MD = max(pert_dose, na.rm = T)) %>%
  dplyr::ungroup() %>% 
  dplyr::inner_join(DRC) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(AUC = ifelse(successful_fit, compute_auc(Lower_Limit, Upper_Limit, Inflection, Slope, md, MD), NA),
                log2.IC50 = ifelse(successful_fit, compute_log_ic50(Lower_Limit, Upper_Limit, Inflection, Slope, md, MD), NA)) 



# Sanity check
DRC %>% 
  dplyr::ungroup() %>% 
  dplyr::sample_n(1000) %>% 
  dplyr::distinct(AUC_Riemann, AUC, log2.IC50) %>% 
  psych::pairs.panels()

DRC %>% dplyr::mutate(FF = is.finite(AUC)) %>% dplyr::count(FF)
  head

# -----
# WRITE THE RELEASE FILES
# -----


LFC.FILTERED %>% 
  dplyr::group_by(SampleID, pert_dose, CompoundPlate, analyte_id, cellset) %>% 
  dplyr::mutate(n.PASS = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(profile_id, analyte_id, LFC, LFC_cb, n.PASS, screen, CompoundPlate, cellset) %>%
  dplyr::mutate(LFC = as.numeric(LFC)) %>% 
  dplyr::distinct() %>% 
  write_csv("data/24Q2/PRISMOncologyReferenceLFC.csv")

LFC.collapsed %>% 
  dplyr::filter(is.finite(LFC_cb)) %>% 
  dplyr::left_join(analyte_meta) %>% 
  tidyr::unite(cn, SampleID, pert_dose, CompoundPlate, screen, sep = "::") %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "LFC_cb", fun.aggregate = median) %>% 
  write.csv("data/24Q2/PRISMOncologyReferenceLFCMatrix.csv")

DRC %>% 
  dplyr::select(SampleID, CompoundPlate, screen, md, MD, analyte_id, cellset,
                fit_name, Lower_Limit, Upper_Limit, Slope, Inflection,
                MSE, MAD, frac_var_explained, AUC, log2.IC50, AUC_Riemann, successful_fit) %>%
  dplyr::distinct() %>% 
  write_csv("data/24Q2/PRISMOncologyReferenceDoseResponseParameters.csv")

QC %>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::select(screen, CompoundPlate, prism_replicate, cellset, analyte_id, pool_id, depmap_id,
                error_rate, NC.median, NC.mad, PC.median, PC.mad, DR, SSMD, PASS) %>% 
  dplyr::distinct() %>% 
  write_csv("data/24Q2/PRISMOncologyReferenceQCTable.csv")


failed_pools <- LFC %>%
  dplyr::left_join(inst_meta) %>%
  dplyr::left_join(analyte_meta) %>% 
  dplyr::distinct(screen, CompoundPlate, prism_replicate, pool_id, pert_well) 

dplyr::inner_join(failed_pools, filtered_wells) %>%
  dplyr::select(-flagged.fraction) %>% 
  dplyr::bind_rows(dplyr::inner_join(failed_pools, filtered_pools)) %>% 
  distinct() %>% 
  write_csv("data/24Q2/PRISMOncologyReferenceFilteredPools.csv")

# -----
# WRITE PORTAL FILES
# -----


M <- LFC.FILTERED %>% 
  dplyr::filter(is.finite(LFC_cb)) %>%
  dplyr::group_by(SampleID, pert_dose, CompoundPlate, analyte_id, cellset) %>%
  dplyr::filter(n() > 1) %>% 
  dplyr::ungroup() %>%
  dplyr::filter(is.finite(LFC_cb)) %>%
  dplyr::left_join(analyte_meta) %>%
  dplyr::left_join(inst_meta) %>% 
  tidyr::unite(cn, SampleID, pert_dose, CompoundPlate, replicate, pert_well, sep = "::") %>% 
  dplyr::group_by(cn, depmap_id) %>% 
  dplyr::top_n(1, cellset) %>% # This line assumes PR500 and PR300P as cell set !
  dplyr::ungroup() %>%
  dplyr::distinct(cn, depmap_id, LFC_cb) %>% 
  reshape2::acast(cn ~ depmap_id, value.var = "LFC_cb")


Condition.annotations <- tibble(cn = rownames(M), Label = 0:(dim(M)[1]-1)) %>%
  dplyr::mutate(SampleID = word(cn, sep = fixed("::")),
                pert_dose = round(as.numeric(word(cn, 2, sep = fixed("::"))),5),
                CompoundPlate = word(cn, 3, sep = fixed("::")),
                Replicate = word(cn, 4, sep = fixed("::")),
                pert_well = word(cn, 5, sep = fixed("::"))) %>% 
  dplyr::left_join(dplyr::mutate(inst_meta, pert_dose = round(pert_dose, 5))) %>% 
  dplyr::distinct(Label, SampleID, pert_dose, pert_dose_unit, Replicate, CompoundPlate, screen) %>%
  dplyr::rename(Dose = pert_dose, DoseUnit = pert_dose_unit) %>%
  dplyr::group_by(CompoundPlate, SampleID, Dose, DoseUnit) %>% 
  dplyr::arrange(SampleID, Replicate) %>% 
  dplyr::mutate(Replicate = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(dplyr::distinct(compound_list, CompoundName, SampleID))

Keep <- Condition.annotations %>% 
  dplyr::distinct(SampleID, screen, CompoundPlate, CompoundName) %>% 
  dplyr::group_by(CompoundName) %>% 
  dplyr::top_n(1,readr::parse_number(screen)) %>%  # This is just a heuristic
  dplyr::top_n(1,readr::parse_number(CompoundPlate)) %>%
  dplyr::ungroup() 

Condition.annotations <- Condition.annotations %>% 
  dplyr::semi_join(Keep) 
  
rownames(M) <- 0:(dim(M)[1]-1)
M <- M[as.character(Condition.annotations$Label),]

DRC %>% 
  dplyr::left_join(analyte_meta) %>%  
  dplyr::rename(ModelID = depmap_id,
                EC50 = Inflection,
                LowerAsymptote = Lower_Limit,
                UpperAsymptote = Upper_Limit) %>% 
  dplyr::semi_join(dplyr::distinct(Condition.annotations, CompoundPlate, SampleID)) %>% 
  dplyr::group_by(ModelID, SampleID, cellset) %>%  
  dplyr::top_n(1, frac_var_explained) %>% 
  dplyr::group_by(ModelID, SampleID) %>%  
  dplyr::top_n(1, cellset) %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct(ModelID, SampleID, CompoundPlate, EC50, LowerAsymptote, UpperAsymptote, Slope) %>% 
  write_csv("data/24Q2/PRISMOncologyReferenceResponseCurves.csv")


DRC %>% 
  dplyr::semi_join(dplyr::distinct(Condition.annotations, CompoundPlate, SampleID)) %>% 
  dplyr::left_join(analyte_meta) %>%  
  dplyr::filter(is.finite(AUC)) %>% 
  dplyr::group_by(depmap_id, SampleID, cellset) %>%  
  dplyr::top_n(1, frac_var_explained) %>% 
  dplyr::group_by(depmap_id, SampleID) %>%  
  dplyr::top_n(1, cellset) %>% 
  dplyr::ungroup() %>% 
  reshape2::acast(depmap_id ~ SampleID,
                  value.var = "AUC") %>%
  write.csv("data/24Q2/PRISMOncologyReferenceAUCMatrix.csv")


DRC %>% 
  dplyr::semi_join(dplyr::distinct(Condition.annotations, CompoundPlate, SampleID)) %>% 
  dplyr::left_join(analyte_meta) %>%  
  dplyr::filter(is.finite(log2.IC50)) %>% 
  dplyr::group_by(depmap_id, SampleID, cellset) %>%  
  dplyr::top_n(1, frac_var_explained) %>% 
  dplyr::group_by(depmap_id, SampleID) %>%  
  dplyr::top_n(1, cellset) %>% 
  dplyr::ungroup() %>% 
  reshape2::acast(depmap_id ~ SampleID,
                  value.var = "log2.IC50") %>%
  write.csv("data/24Q2/PRISMOncologyReferenceLog2IC50Matrix.csv")


Condition.annotations %>% 
  dplyr::select(-screen, -CompoundName) %>% 
  dplyr::distinct() %>% 
  write_csv("data/24Q2/PRISMOncologyReferenceViabilityConditions.csv")

t(2^M) %>% 
  write.csv("data/24Q2/PRISMOncologyReferenceViabilityMatrix.csv")

MM <- M %>% 
  reshape2::melt() %>%
  dplyr::rename(Label = Var1,
                ModelID = Var2,
                LFC = value) %>% 
  dplyr::mutate(Label = as.numeric(Label)) %>% 
  dplyr::left_join(Condition.annotations) %>% 
  dplyr::group_by(ModelID, SampleID, Dose, DoseUnit, CompoundName) %>% 
  dplyr::summarise(LFC = median(LFC, na.rm = T)) %>% 
  dplyr::ungroup()


Collapsed.condition.annotations <- MM %>%
  dplyr::distinct(CompoundName, SampleID, Dose, DoseUnit) %>% 
  dplyr::arrange(SampleID, Dose) %>% 
  dplyr::mutate(Label = paste0(CompoundName," (", SampleID, ") @", Dose, " ", DoseUnit), 
                .before = CompoundName) %>%
  dplyr::select(-CompoundName)

MM <- MM %>% 
  dplyr::left_join(Collapsed.condition.annotations) %>% 
  reshape2::acast(Label ~ ModelID, value.var = "LFC")
  

Collapsed.condition.annotations %>% 
  write_csv("data/24Q2/PRISMOncologyReferenceLog2ViabilityCollapsedConditions.csv")

t(MM) %>% 
  write.csv("data/24Q2/PRISMOncologyReferenceLog2ViabilityCollapsedMatrix.csv")



# One column for each pool

Confounder_pool <- analyte_meta %>% 
  dplyr::distinct(depmap_id, pool_id) %>%
  dplyr::mutate(dummy = 1) %>% 
  dplyr::filter(!is.na(depmap_id)) %>% 
  reshape2::acast(depmap_id ~ pool_id, value.var = "dummy", fill = 0) 

Confounder_cellset <- analyte_meta %>% 
  dplyr::distinct(depmap_id, cellset) %>%
  dplyr::mutate(dummy = 1) %>% 
  dplyr::filter(!is.na(depmap_id)) %>% 
  reshape2::acast(depmap_id ~ cellset, value.var = "dummy", fill = 0) 

Confounder_QC <- QC %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::filter(!is.na(depmap_id), PASS) %>%
  dplyr::distinct(screen, cellset, depmap_id,
                  NC.median, NC.mad, PC.median, PC.mad, SSMD) %>% 
  tidyr::pivot_longer(c(NC.median, NC.mad, PC.median, PC.mad, SSMD),
                      names_to = "QC.col", values_to = "QC.val") %>% 
  reshape2::acast(depmap_id ~ screen  + QC.col, value.var = "QC.val", 
                  fill = 0, fun.aggregate = function(x) median(x , na.rm = T)) 

cl <-rownames(Confounder_QC)
Confounder <- cbind(Confounder_QC[cl, ], Confounder_pool[cl, ], Confounder_cellset[cl,])
Confounder <- Confounder[,!duplicated(t(Confounder))]

write.csv(Confounder, "data/24Q2/PRISMOncologyReferenceConfounderMatrix.csv")

