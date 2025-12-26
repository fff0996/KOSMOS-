library(dplyr)
library(survival)

## ========== STEP 1: 데이터 준비 ==========

gene_list <- unique(dat_somatic$Hugo_Symbol)  # 558 genes

dat_somatic <- dat_somatic %>%
  mutate(Sample_clean = Tumor_Sample_Barcode %>%
           as.character() %>%
           str_remove("_snvindel") %>%
           str_remove("_CancerScan") %>%
           str_remove("_RNA") %>%
           gsub("_", "-", .))

mut_matrix <- dat_somatic %>%
  select(Sample_clean, Hugo_Symbol) %>%
  distinct() %>%
  mutate(mutated = 1) %>%
  pivot_wider(names_from = Hugo_Symbol, 
              values_from = mutated, 
              values_fill = 0) %>%
  column_to_rownames("Sample_clean")


all_samples <- anno$SUBJ_ID

# Missing 샘플들 (mutation 없는 25명)
missing_samples <- setdiff(all_samples, rownames(mut_matrix))

cat("Samples with mutations:", nrow(mut_matrix), "\n")
cat("Samples without mutations:", length(missing_samples), "\n")

# Missing 샘플들을 0으로 채워서 추가
if (length(missing_samples) > 0) {
  missing_matrix <- matrix(0, 
                          nrow = length(missing_samples), 
                          ncol = ncol(mut_matrix))
  rownames(missing_matrix) <- missing_samples
  colnames(missing_matrix) <- colnames(mut_matrix)
  
  mut_matrix <- rbind(mut_matrix, missing_matrix)
}

mut_matrix_final <- mut_matrix

# Sample matching
matched_samples <- intersect(rownames(mut_matrix_final), anno$SUBJ_ID)

mut_matched <- mut_matrix_final[matched_samples, ]
anno_matched <- anno %>%
  filter(SUBJ_ID %in% matched_samples) %>%
  arrange(match(SUBJ_ID, matched_samples))



## ========== STEP 2: Fisher's Exact Test ==========

fisher_results <- lapply(colnames(mut_matched), function(gene) {
  
  test_data <- data.frame(
    response = anno_matched$Responder.vs...Non.responder,
    mutated = mut_matched[, gene]
  ) %>%
    filter(!is.na(response), !is.na(mutated))
  
  tab <- table(test_data$response, test_data$mutated)
  
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(NULL)
  }
  
  # Fisher test (p-value용)
  ft <- fisher.test(tab)
  
  # Odds ratio: 0.5 correction
  tab_corrected <- tab + 0.5
  or_corrected <- (tab_corrected["1", "1"] * tab_corrected["0", "0"]) / 
                  (tab_corrected["1", "0"] * tab_corrected["0", "1"])
  
  data.frame(
    gene = gene,
    non_res_mut = tab["0", "1"],
    non_res_total = sum(tab["0", ]),
    res_mut = tab["1", "1"],
    res_total = sum(tab["1", ]),
    freq_non_res = tab["0", "1"] / sum(tab["0", ]),
    freq_res = tab["1", "1"] / sum(tab["1", ]),
    odds_ratio = or_corrected,
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

cox_results <- lapply(colnames(mut_matrix), function(gene) {
  
  cox_data <- data.frame(
    pfs = anno$pfs,
    PFS_event = anno_$PFS_event,
    altered = mut_matrix[, gene]
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



##cox sig heatmap 암종별 
library(pheatmap)
library(dplyr)

## ========== Cox 유의한 gene의 암종별 alteration frequency ==========

## 1) Cox에서 p < 0.05인 gene
prog_genes <- cox_res %>%
  filter(p_value < 0.05) %>%
  pull(gene)

cat("Prognostic genes (Cox p < 0.05):", length(prog_genes), "\n")

## 2) 이 gene들로 alteration matrix
alter_cox <- alter_binary[matched_samples, prog_genes]

## 3) Cancer type 리스트
cancer_types <- sort(unique(anno_matched$dx))
cat("Cancer types:", length(cancer_types), "\n")

## 4) Frequency matrix 생성
freq_matrix_cox <- matrix(0, nrow = length(prog_genes), ncol = length(cancer_types))
rownames(freq_matrix_cox) <- prog_genes
colnames(freq_matrix_cox) <- cancer_types

for (gene in prog_genes) {
  for (dx in cancer_types) {
    samples_dx <- anno_matched %>% filter(dx == !!dx) %>% pull(SubjectNo)
    if (length(samples_dx) > 0 && gene %in% colnames(alter_cox)) {
      freq_matrix_cox[gene, dx] <- sum(alter_cox[samples_dx, gene]) / length(samples_dx)
    }
  }
}

## 5) 확인
cat("\nFrequency matrix:", dim(freq_matrix_cox), "\n")
cat("Max frequency:", max(freq_matrix_cox), "\n")
cat("Min frequency:", min(freq_matrix_cox), "\n")

## Genes with max freq = 0 확인
gene_max_freq <- apply(freq_matrix_cox, 1, max)
cat("Genes with max freq = 0:", sum(gene_max_freq == 0), "\n")

if (sum(gene_max_freq == 0) > 0) {
  zero_genes <- names(gene_max_freq[gene_max_freq == 0])
  cat("\nZero frequency genes:\n")
  print(zero_genes)
  
  # 왜 0인지 확인
  for (gene in head(zero_genes, 3)) {
    cat("\n", gene, ":\n")
    cat("  In alter_binary:", sum(alter_binary[, gene]), "\n")
    cat("  In alter_cox:", ifelse(gene %in% colnames(alter_cox), 
                                   sum(alter_cox[, gene]), "NOT IN MATRIX"), "\n")
  }
}

## 6) Cancer type 이름 정리
colnames(freq_matrix_cox) <- gsub("cancer", "", colnames(freq_matrix_cox))
colnames(freq_matrix_cox) <- gsub("\\(.*\\)", "", colnames(freq_matrix_cox))
colnames(freq_matrix_cox) <- trimws(colnames(freq_matrix_cox))

## 7) Max frequency 기준 legend
max_freq_cox <- ceiling(max(freq_matrix_cox) * 10) / 10

## 8) Heatmap
pdf("cancer_type_heatmap_prognostic_genes_fixed.pdf", width = 12, height = 10)
pheatmap(
  freq_matrix_cox,
  color = colorRampPalette(c("white", "#FFB6A3", "#E64B35", "#C0201C"))(100),
  breaks = seq(0, max_freq_cox, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 9,
  border_color = "gray80",
  main = "Cancer type-specific alteration frequency of prognostic genes",
  legend_breaks = seq(0, max_freq_cox, by = 0.2),
  angle_col = 45
)
dev.off()

write.csv(freq_matrix_cox, "cancer_type_freq_prognostic_genes.csv", row.names = TRUE)

cat("\nFiles saved!\n")

library(ComplexHeatmap)
library(dplyr)
library(grid)

## 1) Pathways
pathways <- list(
  `DNA Damage Response (DDR)` = c("TP53", "BRCA2", "ATM", "ATR", "BRCA1", "PALB2", "MDM2", "BAP1"),
  `Cell Cycle` = c("CDK12", "RB1", "CDKN2A", "CCND1", "CCNE1"),
  `Epigenetic Modifier` = c("ARID1A", "ARID1B", "CREBBP", "KDM5A", "ARID2", "PBRM1", "IDH2", "IDH1"),
  `RAS/MAPK` = c("ERBB2", "FGFR3", "ALK", "STK11", "NF1", "PTEN", "BRAF", "PIK3CA", "FGFR2", "ERBB3", "KRAS"),
  `TGF-β` = c("TGFBR2", "MYC", "FBXW7", "SMAD4"),
  `WNT` = c("APC", "EP300", "CTNNB1")
)
pathway_genes <- unique(unlist(pathways))

## 2) Samples
matched_samples <- intersect(rownames(alter_df_full_updated), anno$SubjectNo)
anno_matched <- anno %>%
  filter(SubjectNo %in% matched_samples) %>%
  arrange(res.non)
sample_order <- anno_matched$SubjectNo

## 3) Alteration matrix
alter_pathway <- alter_df_full_updated[sample_order, pathway_genes]

## 4) Gene frequency
gene_freq <- colSums(alter_pathway > 0) / nrow(alter_pathway)
gene_order <- names(sort(gene_freq, decreasing = TRUE))

## 5) Oncoprint matrix
onco_mat <- matrix("", nrow = length(pathway_genes), ncol = length(sample_order))
rownames(onco_mat) <- pathway_genes
colnames(onco_mat) <- sample_order

for (i in 1:length(sample_order)) {
  for (j in 1:length(pathway_genes)) {
    val <- alter_pathway[i, j]
    if (val == 1) {
      onco_mat[j, i] <- "MUT"
    } else if (val == 2) {
      onco_mat[j, i] <- "CNV"
    } else if (val == 3) {
      onco_mat[j, i] <- "MUT_CNV"
    }
  }
}

## 6) Reorder
onco_mat <- onco_mat[gene_order, ]
gene_freq <- gene_freq[gene_order]

## 7) Row split
row_split <- rep(NA, length(gene_order))
names(row_split) <- gene_order
for (pw in names(pathways)) {
  for (gene in pathways[[pw]]) {
    if (gene %in% gene_order) {
      row_split[gene] <- pw
    }
  }
}
row_split <- factor(row_split, levels = names(pathways))

## 8) Top annotation - 깔끔하게
gap_pos <- sum(anno_matched$res.non == "non_res")

# Top barplot for number of altered genes
altered_counts <- colSums(onco_mat != "")

top_ha <- HeatmapAnnotation(
  `No. of altered\ngenes` = anno_barplot(
    altered_counts,
    gp = gpar(fill = "#4DBBD5"),
    height = unit(2, "cm"),
    border = FALSE,
    axis_param = list(
      side = "left",
      labels_rot = 0
    )
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  gap = unit(1, "mm")
)

## 9) Plot
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = "#EEEEEE", col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#3C5488", col = NA))
  },
  CNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#F39B7F", col = NA))
  },
  MUT_CNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#8B4789", col = NA))
  }
)

pdf("oncoprint_pathway_clean.pdf", width = 14, height = 10)
oncoPrint(
  onco_mat,
  alter_fun = alter_fun,
  col = c("MUT" = "#3C5488", "CNV" = "#F39B7F", "MUT_CNV" = "#8B4789"),
  top_annotation = top_ha,
  right_annotation = rowAnnotation(
    ` ` = anno_barplot(
      gene_freq * 100, 
      gp = gpar(fill = "#3C5488"),
      border = FALSE,
      width = unit(3, "cm"),
      axis_param = list(
        at = c(0, 50, 100),
        labels = c("0", "50", "100"),
        side = "top"
      )
    ),
    annotation_name_gp = gpar(fontsize = 0)
  ),
  row_split = row_split,
  column_split = factor(
    c(rep("Non-Responder", gap_pos), 
      rep("Responder", length(sample_order) - gap_pos)),
    levels = c("Non-Responder", "Responder")
  ),
  show_column_names = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9),
  pct_side = "right",
  pct_gp = gpar(fontsize = 8),
  column_title = "Genomic Alterations (n = 121)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10),
  border = TRUE
)
dev.off()
