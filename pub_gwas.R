#' pUb GWAS analysis
#'
#' @param pub_cohort Data frame. Cohort data and covariates
#' @param model_dir String. File path for results to be saved to
#' @param imputed_pca Data frame. PCA analysis results
#' @param cwow Data frame. Hippocampus staining information
#'
#' @return String. Returns `model_dir` invisibly
run_pub_gwas <- function(pub_cohort, model_dir, imputed_pca, cwow) {

  model_path <- function(path) fs::path(model_dir, path)

  pheno_data <- pub_cohort %>%
    select(FID = fid, IID = iid, PUB = pub_log)

  pheno_file <- model_path("pheno.txt")

  write_delim(pheno_data, file = pheno_file, delim = " ")

  model_covars <- pub_cohort %>%
    left_join(cwow, by = "npid") %>%
    select(FID = fid, IID = iid, AGE = aad_imputed, HIPPO = hippo) %>%
    inner_join(imputed_pca, by = c("FID", "IID")) %>%
    mutate(HIPPO = as.numeric(HIPPO) - 1)

  covar_file <- model_path("covars.txt")

  write_delim(model_covars, file = covar_file, delim = " ")

  plink_file <- here::here("data/imputed/merged_qc/merged_qc2")

  model_outfile <- model_path("model")

  plink <- here::here("plink")

  maf_run <- glue::glue(
    "{plink} --bfile {plink_file} --freq --out {model_outfile} --keep {pheno_file} --hide-covar",
  )

  system(maf_run)


  model_run <- glue::glue(
    "{plink} --bfile {plink_file} --adjust --linear sex --covar {covar_file} --ci 0.95 --out {model_outfile} --pheno {pheno_file} --hide-covar"
  )

  system(model_run)

  invisible(model_dir)

}
