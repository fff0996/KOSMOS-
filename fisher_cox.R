library(dplyr)
library(survival)

## ========== STEP 1: 데이터 준비 ==========

## Binary alteration (0 vs 1/2/3)
alter_binary <- alter_df_full
alter_binary[alter_binary > 0] <- 1

## Sample matching
matched_samples <- intersect(rownames(alter_binary), anno$SubjectNo)

alter_matched <- alter_binary[matched_samples, ]
anno_matched <- anno %>%
  filter(SubjectNo %in% matched_samples) %>%
  arrange(match(SubjectNo, matched_samples))

cat("Matched samples:", nrow(alter_matched), "\n")
cat("Genes:", ncol(alter_matched), "\n")

## ========== STEP 2: Fisher's Exact Test ==========

fisher_results <- lapply(colnames(alter_matched), function(gene) {
  
  test_data <- data.frame(
    res.non = anno_matched$res.non,
    altered = alter_matched[, gene]
  ) %>%
    filter(!is.na(res.non), !is.na(altered))
  
  # 2x2 table
  tab <- table(test_data$res.non, test_data$altered)
  
  # Table이 2x2가 아니면 skip
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(NULL)
  }
  
  # Fisher test
  ft <- fisher.test(tab)
  
  data.frame(
    gene = gene,
    non_res_alt = tab["non_res", "1"],
    non_res_total = sum(tab["non_res", ]),
    res_alt = tab["res", "1"],
    res_total = sum(tab["res", ]),
    freq_non_res = tab["non_res", "1"] / sum(tab["non_res", ]),
    freq_res = tab["res", "1"] / sum(tab["res", ]),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
})

fisher_res <- bind_rows(fisher_results) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    diff = freq_res - freq_non_res
  ) %>%
  arrange(p_value)

cat("\n=== Fisher's Exact Test Results ===\n")
cat("Total genes tested:", nrow(fisher_res), "\n")
cat("Significant (p < 0.05):", sum(fisher_res$p_value < 0.05, na.rm=TRUE), "\n")
cat("FDR significant (p_adj < 0.05):", sum(fisher_res$p_adj < 0.05, na.rm=TRUE), "\n")

## Top results
cat("\nTop 10 genes:\n")
print(head(fisher_res, 10))

write.csv(fisher_res, "fisher_test_response_fixed.csv", row.names = FALSE)

## ========== STEP 3: Cox Regression ==========

cox_results <- lapply(colnames(alter_matched), function(gene) {
  
  cox_data <- data.frame(
    pfs = anno_matched$pfs,
    PFS_event = anno_matched$PFS_event,
    altered = alter_matched[, gene]
  ) %>%
    filter(!is.na(pfs), !is.na(PFS_event), !is.na(altered))
  
  # Minimum 3 altered samples
  if (sum(cox_data$altered) < 3) {
    return(NULL)
  }
  
  # Cox regression
  tryCatch({
    fit <- coxph(Surv(pfs, PFS_event) ~ altered, data = cox_data)
    
    data.frame(
      gene = gene,
      n_altered = sum(cox_data$altered),
      n_wt = sum(cox_data$altered == 0),
      hazard_ratio = exp(coef(fit)),
      HR_lower = exp(confint(fit)[1]),
      HR_upper = exp(confint(fit)[2]),
      p_value = summary(fit)$coefficients[5],
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    return(NULL)
  })
})

cox_res <- bind_rows(cox_results) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

cat("\n=== Cox Regression Results ===\n")
cat("Total genes tested:", nrow(cox_res), "\n")
cat("Significant (p < 0.05):", sum(cox_res$p_value < 0.05, na.rm=TRUE), "\n")
cat("FDR significant (p_adj < 0.05):", sum(cox_res$p_adj < 0.05, na.rm=TRUE), "\n")

## Top results
cat("\nTop 10 prognostic genes:\n")
print(head(cox_res, 10))

write.csv(cox_res, "cox_regression_pfs_fixed.csv", row.names = FALSE)

## ========== STEP 4: Summary ==========

cat("\n=== Summary ===\n")
cat("Fisher test:", nrow(fisher_res), "genes\n")
cat("  - p < 0.05:", sum(fisher_res$p_value < 0.05), "\n")
cat("  - p_adj < 0.05:", sum(fisher_res$p_adj < 0.05), "\n")

cat("\nCox regression:", nrow(cox_res), "genes\n")
cat("  - p < 0.05:", sum(cox_res$p_value < 0.05), "\n")
cat("  - p_adj < 0.05:", sum(cox_res$p_adj < 0.05), "\n")

## Both significant
both_sig <- intersect(
  fisher_res %>% filter(p_value < 0.05) %>% pull(gene),
  cox_res %>% filter(p_value < 0.05) %>% pull(gene)
)

cat("\nGenes significant in BOTH tests:", length(both_sig), "\n")
print(both_sig)
