#Wrapper for SpATS
f_spats <- function(data, response = response, random = random, fixed = fixed, 
                    genotype.as.random = genotype.as.random, genotype = genotype) {
  SpATS(response = response, random = as.formula(random), 
        fixed = as.formula(fixed),
        spatial = ~PSANOVA(RangeBL, RowBL, nseg = c(20,20), nest.div=c(2,2)), 
        genotype = genotype, genotype.as.random = genotype.as.random, data = data,
        control = list(maxit = 100, tolerance = 1e-03, monitoring = 0))
}

get_h2 <- function(obj){
  ifelse(length(unique(obj$data$Rep))> 1, SpATS::getHeritability(obj), NA)
}

#helper function
make_design_matrix <- function(geno, names) {
  Nrow <- length(geno)
  Ncol <- length(names)
  col <- match(geno, names)
  frame <- data.frame(i = c(1:Nrow), j = col, v = rep(1,Nrow))
  frame <- subset(frame, is.na(col) == FALSE)
  L <- as.list(frame)
  X <- spam::spam(L, nrow = Nrow, ncol = Ncol)
  return(X)
}

#helper function
construct_genotype_prediction_matrix <- function(object, newdata) {
  Z_geno = make_design_matrix(newdata[,object$model$geno$genotype], object$terms$geno$geno_names)
  if(object$model$geno$as.random)
    Z_geno <- Z_geno[,object$terms$geno$ndx]
  else
    Z_geno <- Z_geno[, object$terms$geno$ndx[2:length(object$terms$geno$ndx)]]
  Z_geno	
}

#extract BLUE from SpATS
get_BLUE_spats <- function(object) {
  fitted <- object$fitted
  intercept <- object$coeff['Intercept']
  gen_mod_mat <- construct_genotype_prediction_matrix(object, object$data)
  gen_coeff <- object$coeff[1:ncol(gen_mod_mat)]
  geno_pred <- as.vector(gen_mod_mat %*% gen_coeff)
  BLUE <- as.data.frame(intercept + gen_coeff) %>% 
    rownames_to_column() %>% 
    as_tibble() %>% 
    dplyr::rename(Gen_Name = 1, BLUE = 2)
}

#extract BLUE from SpATS
spat_corr_spats <- function(object) {
  fitted <- object$fitted
  intercept <- object$coeff['Intercept']
  gen_mod_mat <- construct_genotype_prediction_matrix(object, object$data)
  gen_coeff <- object$coeff[1:ncol(gen_mod_mat)]
  geno_pred <- as.vector(gen_mod_mat %*% gen_coeff)
  residual <- object$residuals
  plot_corr <- as.data.frame(intercept + geno_pred + residual) %>% rename(spat_corr = 1) %>% 
    bind_cols(object$data, .) %>% as_tibble() %>% dplyr::select(-value)
}

#extract spatial trend
get_spatial <- function(object){
  fitted <- object$fitted
  intercept <- object$coeff['Intercept']
  gen_mod_mat <- construct_genotype_prediction_matrix(object, object$data)
  gen_coeff <- object$coeff[1:ncol(gen_mod_mat)]
  geno_pred <- as.vector(gen_mod_mat %*% gen_coeff)
  spatial <- fitted-geno_pred-intercept
  spatial <- cbind(object$data, spatial) %>% as_tibble()
  return(spatial)
}

#calculate heritability using BLUEs
get_h2_years <- function(data, fixed = fixed, random = random){
  #check if random is factor and convert if required
  if(!is.factor(data$Gen_Name)){
    print("Gen_Name converted to factor")
    data$Gen_Name <- as.factor(data$Gen_Name)
  }
  mod <- asreml(fixed = as.formula(paste("BLUE ~", fixed)),
                random = as.formula(paste("~", random)),
                data = data,
                trace = FALSE)
  GenVar <- summary(mod)$varcomp["Gen_Name", "component"]
  ErrVar <- summary(mod)$varcomp["units!R","component"]
  H2 <- round(GenVar/(GenVar + ErrVar/2), 3)
}

#get spatially corrected Plot values
# > either by adding the residual to the BLUE, if both Lots measured
# > OR by subtracting the design effect, if only 1 Lot measured
spat_corr_asreml <- function(data, fixed, random, residual){
  if(!is.factor(data$Gen_Name)){
    data$Gen_Name <- as.factor(data$Gen_Name)
  }
  if(length(unique(data$Lot)) > 1){
    #get a spatially corrected value per plot,
    #by adding the residual to the BLUE
    as <- asreml(fixed = as.formula(fixed),
                 random = as.formula(paste("~", random)),
                 residual = as.formula(paste("~", residual)),
                 data = data,
                 trace = FALSE)
    pv1 <- predict(as, "Gen_Name")
    residuals <- as$residuals %>% as.data.frame()
    out <- full_join(data, pv1$pvals[c("Gen_Name", "predicted.value")], by = "Gen_Name") %>% 
      mutate(Gen_Name = as.factor(Gen_Name)) %>% 
      bind_cols(., residuals) %>%
      mutate(spat_corr = predicted.value + e) %>% 
      dplyr::select(-predicted.value, -e)
  } else {
    #get a spatially corrected value per plot,
    #by subtracting the design effect from the observed value
    message("No reps. Removing Lot and Genotype effect from model")
    #adjusted model, only spatial effects
    as <- asreml(fixed = as.formula("value ~ Row + Range"),
                 random = as.formula(paste("~", random)),
                 residual = as.formula("~ar1(Range):ar1(Row)"),
                 data = data,
                 trace = FALSE)
    fitted <- fitted(as) %>% as_tibble() %>% dplyr::rename(fitted="value")
    nranges <- length(unique(data$Range))
    nrows <- length(unique(data$Row))
    coeffs <- as.vector(as$coefficients[[1]])[1:(nranges+nrows)]
    Range <- tibble("Range" = c(1:nranges), "RangeEffect" = coeffs[1:nranges])
    Row <- tibble("Row" = c(1:nrows), "RowEffect" = coeffs[(nranges+1):(nranges+nrows)])
    design <- tibble("Range" = rep(1:nranges, each = nrows), "Row" = rep(1:nrows, nranges))
    Effects <- full_join(design, Range) %>% full_join(., Row) %>% 
      mutate(DesignEffect = RangeEffect + RowEffect) %>% 
      mutate_if(is.integer, as.factor)
    out <- bind_cols(data, fitted) %>% full_join(., Effects) %>% 
      mutate(spat_corr = value-DesignEffect) %>% 
      dplyr::select(-fitted, -contains("Effect"))
  }
}

#get Year-wise BLUE
get_BLUE_asreml <- function(data, fixed, random, residual){
  if(!is.factor(data$Gen_Name)){
    data$Gen_Name <- as.factor(data$Gen_Name)
  }
  if(length(unique(data$Lot)) > 1){
    as <- asreml(fixed = as.formula(fixed),
                 random = as.formula(paste("~", random)),
                 residual = as.formula(paste("~", residual)),
                 data = data,
                 trace = FALSE)
    pv <- predict(as, "Gen_Name")
    out <- pv$pvals[c("Gen_Name", "predicted.value")] %>% data.frame() %>% as_tibble() %>% 
      dplyr::rename(BLUE = predicted.value)
  } else {
    message("No reps. Removing Lot and Genotype effect from model")
    #adjusted model, only spatial effects
    as <- asreml(fixed = as.formula("value ~ Row + Range"),
                 random = as.formula(paste("~", random)),
                 residual = as.formula("~ar1(Range):ar1(Row)"),
                 data = data,
                 trace = FALSE)
    fitted <- fitted(as) %>% as_tibble() %>% dplyr::rename(fitted="value")
    nranges <- length(unique(data$Range))
    nrows <- length(unique(data$Row))
    coeffs <- as.vector(as$coefficients[[1]])[1:(nranges+nrows)]
    Range <- tibble("Range" = c(1:nranges), "RangeEffect" = coeffs[1:nranges])
    Row <- tibble("Row" = c(1:nrows), "RowEffect" = coeffs[(nranges+1):(nranges+nrows)])
    design <- tibble("Range" = rep(1:nranges, each = nrows), "Row" = rep(1:nrows, nranges))
    Effects <- full_join(design, Range) %>% full_join(., Row) %>% 
      mutate(DesignEffect = RangeEffect + RowEffect) %>% 
      mutate_if(is.integer, as.factor)
    out <- bind_cols(data, fitted) %>% full_join(., Effects) %>% 
      mutate(spat.corr = value-DesignEffect) %>% 
      #average, so that each check gets one "pseudo-BLUE"
      group_by(Gen_Name) %>% 
      dplyr::summarize(spat.corr.mean = mean(spat.corr, na.rm = TRUE)) %>% 
      dplyr::select(Gen_Name, spat.corr.mean) %>% dplyr::rename(BLUE = spat.corr.mean) %>% 
      unique()
  }
}

#calculate repeatability for traits measured on > 1 Lots
get_w2_asreml <- function(data, fixed, random, residual, cullis = TRUE){
  if(!is.factor(data$Gen_Name)){
    data$Gen_Name <- as.factor(data$Gen_Name)
  }
  if(length(unique(data$Lot)) > 1){
    as <- asreml(fixed = as.formula(fixed),
                 random = as.formula(paste("~", random)),
                 residual = as.formula(paste("~", residual)),
                 data = data,
                 trace = FALSE)
    if(cullis) {
      R <- unlist(strsplit(random, " "))[1]
      GenVar <- summary(as)$varcomp[R,"component"]
      #Predict genotypic mean values to estimate avsed
      pv.G <- predict(as, classify = "Gen_Name")
      avsed.G<-pv.G$avsed
      #The Cullis version of heritability based on avsed^2 and V.G
      w2 <- 1-(avsed.G^2/(2*GenVar))
    } else {
      R <- unlist(strsplit(random, " "))[1]
      GenVar <- summary(as)$varcomp[R,"component"]
      err_srcs <- grep(pattern = "R$", rownames(summary(as)$varcomp), value = TRUE)
      ErrVar_comp <- NULL
      for(i in 1:length(err_srcs)){
        ErrVar_comp[i] <- summary(as)$varcomp[err_srcs[i],"component"]
      }
      ErrVar <- mean(ErrVar_comp)
      w2 <- GenVar/(GenVar + ErrVar/length(ErrVar_comp))
    }
  } else {
    message("No reps. Repeatability not estimated.")
    w2 <- NA
  }
}

#calculate heritability, using spatially corrected plot values
get_h2_asreml1 <- function(data, fixed, random, residual, cullis = TRUE){
  if(!is.factor(data$Gen_Name)){
    data$Gen_Name <- as.factor(data$Gen_Name)
  }
  as <- asreml(fixed = as.formula(fixed),
               random = as.formula(paste("~", random)),
               residual = as.formula(paste("~", residual)),
               data = data,
               trace = FALSE)
  if(cullis){
    GenVar <- summary(as)$varcomp["Gen_Name","component"]
    #Predict genotypic mean values to estimate avsed
    pv.G <- predict(as, classify = "Gen_Name")
    pv.G$pvals
    avsed.G<-pv.G$avsed
    #The Cullis version of heritability based on avsed^2 and V.G
    h2<- 1-(avsed.G^2/(2*GenVar))
  } else {
    #number of Years
    rep <- strsplit(fixed, "~") %>% lapply(., "[[", 2) %>% unlist() %>% gsub(" ", "", .)
    n1 <- as.character(as.data.frame(data)[,rep]) %>% unique() %>% length()
    #number of Year:Replicate combinations
    n2 <- data %>% group_by(Exp, Lot) %>% nest() %>% nrow
    GenVar <- summary(as)$varcomp["Gen_Name","component"]
    ErrVar <- summary(as)$varcomp["units!R","component"]
    GenExpVar <- summary(as)$varcomp["Gen_Name:Exp","component"]
    h2 <- GenVar/(GenVar +  GenExpVar/n1 + ErrVar/n2)
  }
}

#calculate heritability, using year-wise BLUEs
get_h2_asreml2 <- function(data, fixed, random, residual, cullis = TRUE){
  if(!is.factor(data$Gen_Name)){
    data$Gen_Name <- as.factor(data$Gen_Name)
  }
  as <- asreml(fixed = as.formula(fixed),
               random = as.formula(random),
               residual = as.formula(residual),
               data = data,
               trace = TRUE)
  if(cullis){
    GenVar <- summary(as)$varcomp["Gen_Name","component"]
    #Predict genotypic mean values to estimate avsed
    pv.G <- predict(as, classify = "Gen_Name")
    avsed.G<-pv.G$avsed
    #The Cullis version of heritability based on avsed^2 and V.G
    h2<- 1-(avsed.G^2/(2*GenVar))
  } else {
    #number of Years
    rep <- strsplit(fixed, "~") %>% lapply(., "[[", 2) %>% unlist() %>% gsub(" ", "", .)
    n <- as.character(as.data.frame(data)[,rep]) %>% unique() %>% length()
    GenVar <- summary(as)$varcomp["Gen_Name","component"]
    ErrVar <- summary(as)$varcomp["units!R","component"]
    h2 <- GenVar/(GenVar + ErrVar/n)
  }
}

#get BLUEs
get_BLUE_years_asreml <- function(data, fixed, random, residual){
  #check if random is factor and convert if required
  if(!is.factor(data$Gen_Name)){
    data$Gen_Name <- as.factor(data$Gen_Name)
  }
  as <- asreml(fixed = as.formula(fixed),
               random = as.formula(paste("~", random)),
               residual = as.formula(paste("~", residual)),
               data = data,
               trace = FALSE)
  pv <- predict(as, "Gen_Name")
  out <- pv$pvals[c("Gen_Name", "predicted.value")] %>% data.frame() %>% as_tibble() %>% 
    dplyr::rename(BLUE = predicted.value)
}

#helper function
remove_reps <- function(data, value, missing){
  remove <- ifelse(sum(is.na(data %>% pull(value)), na.rm = TRUE) > missing, TRUE, FALSE)
}

