# Helper functions ----

extract_covars_from_nested <- function(tbl, from, vars)
{
  dat_vars <- list()
  for(i in vars){
    dat_vars[[i]] <- do.call("c", lapply(tbl[[from]], "[[", i))
  }
  if(length(vars)>1){
    vars <- do.call("cbind", dat_vars)
    out <- cbind(tbl, vars)
    out <- as_tibble(out)
  } else {
    out <- as_tibble(dat_vars)
    out <- bind_cols(tbl, out)
  }
  return(out)
}

f_scale <- function(x)
{
  (x - min(x, na.rm = TRUE))/diff(range(x, na.rm = TRUE))
}

findOutlier <- function(tbl, var, cutoff = 3) 
{
  ## Calculate the sd
  vardata <- pull(tbl[,var])
  sd <- sd(vardata, na.rm = TRUE)
  mean <- mean(vardata, na.rm = TRUE)
  ## Identify the cells with value greater than cutoff * sd (column wise)
  outliers_idx <- which(vardata > mean + cutoff * sd | vardata < mean - cutoff * sd)
  tbl[outliers_idx, var] <- NA
  return(tbl)
}

chauvenet <- function(tbl, var)
{
  vec <- tbl %>% pull(var)
  mean <- mean(vec)
  sd <- sd(vec)
  n <- length(vec)
  crit <- 1.0/(2*n)
  d = abs(vec-mean)/sd
  prob = erfc(d)
  idx <- which(prob < crit)
  tbl_out <- tbl[-idx,]
  return(tbl_out)
}

get_lm_pars <- function(obj, pred)
{
  slope <- unname(coef(obj)[pred])
}

get_RMSE <- function(obj)
{
  RSS <- c(crossprod(obj$residuals))
  MSE <- RSS / length(obj$residuals)
  RMSE <- sqrt(MSE)
}

#calculate correlation and p-values
do_cor_test <- function(data, x, y, use = "pairwise.complete.obs", 
                        method = "pearson", return = "estimate") 
{
  cor <- cor.test(data %>% pull(x), data %>% pull(y), 
                  use = use, method = method) %>% 
    broom::tidy(.) %>% pull(return)
  return(cor)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#extract BLUE from SpATS
get_spat_corr_spats <- function(
  object, # a fitted SpATS object
  response, #the modelled trait
  element_return = "minimal") #either "minimal" (only corrected values per plot) or "full" (including design) 
{
  fitted <- object$fitted
  intercept <- object$coeff['Intercept']
  gen_mod_mat <- construct_genotype_prediction_matrix(object, object$data)
  gen_coeff <- object$coeff[1:ncol(gen_mod_mat)]
  geno_pred <- as.vector(gen_mod_mat %*% gen_coeff)
  residual <- object$residuals
  if(element_return == "full"){
    plot_corr <- as.data.frame(intercept + geno_pred + residual) %>% rename(spat_corr = 1) %>% 
      bind_cols(object$data, .) %>% as_tibble() %>% dplyr::select(-one_of(response), -weights)
  } else if(element_return == "minimal"){
    plot_corr <- as.data.frame(intercept + geno_pred + residual) %>% rename(spat_corr = 1) %>% 
      as_tibble()
  }
}