# ======================================================
# 300_marker_selection_chisq_fullMD.R — χ² top-300 on FULL MD (no inner fold)
# ======================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(readr); library(purrr)
})

# ---------- 사용자 파라미터 ----------
FNAME       <- "10Xv2"
PROC_DIR    <- file.path("C:/test/scRNAseq_ensemble/Data/processed", FNAME)
OUT_DIR     <- file.path(PROC_DIR, "markers"); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TOP_N       <- 300     # 클래스별 선별 개수
MIN_POS_CT  <- 5       # class+ & gene+ 최소 셀 수(너무 희소한 유전자 제외)
USE_YATES   <- TRUE    # Yates 연속성 보정 사용

# ---------- 로드 ----------
seu    <- readRDS(file.path(PROC_DIR, paste0(FNAME, "_seurat_sct.rds")))
splits <- readRDS(file.path(PROC_DIR, paste0(FNAME, "_splits.rds")))
stopifnot("cell_type" %in% colnames(seu@meta.data))

meta <- seu@meta.data %>% mutate(cell_id = rownames(.))
train_ids <- splits$train_ids
train_idx <- match(train_ids, colnames(seu))
stopifnot(!any(is.na(train_idx)))

# counts(원시) → presence/absence
bin <- GetAssayData(seu, slot = "counts") > 0  # genes x cells (logical sparse)
bin_tr <- bin[, train_idx, drop = FALSE]
y_tr   <- factor(meta$cell_type[train_idx])
classes <- levels(y_tr)

# ---------- 유틸: 2x2 χ²(p), logOR, stat ----------
chisq_p_logor <- function(a, b, c, d, yates = TRUE) {
  # a=class+ & gene+, b=class- & gene+, c=class+ & gene-, d=class- & gene-
  n <- a + b + c + d
  # Haldane-Anscombe 0.5 보정으로 logOR 안정화
  logOR <- log(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))
  num <- abs(a * d - b * c)
  if (yates) num <- pmax(0, num - n * 0.5)     # 연속성 보정
  den <- (a + b) * (c + d) * (a + c) * (b + d)
  stat <- (n * (num^2)) / (den + 1e-12)
  p <- pchisq(stat, df = 1, lower.tail = FALSE)
  list(p = p, logOR = logOR, stat = stat)
}

# ---------- MD 전체에서 클래스별 top-300 ----------
marker_list <- vector("list", length(classes))
names(marker_list) <- classes

for (cls in classes) {
  y <- as.integer(y_tr == cls)           # 1 vs rest
  pos_idx <- which(y == 1)
  neg_idx <- which(y == 0)
  
  # 2x2 카운트 (벡터화)
  a <- Matrix::rowSums(bin_tr[, pos_idx, drop = FALSE])  # class+ & gene+
  b <- Matrix::rowSums(bin_tr[, neg_idx, drop = FALSE])  # class- & gene+
  c <- sum(y == 1) - a                                   # class+ & gene-
  d <- sum(y == 0) - b                                   # class- & gene-
  
  a <- as.numeric(a); b <- as.numeric(b); c <- as.numeric(c); d <- as.numeric(d)
  
  # 너무 희소한 유전자 제거
  keep <- (a >= MIN_POS_CT)
  genes_kept <- rownames(bin_tr)[keep]
  
  st <- chisq_p_logor(a[keep], b[keep], c[keep], d[keep], yates = USE_YATES)
  
  df <- tibble(
    gene = genes_kept,
    a = a[keep], b = b[keep], c = c[keep], d = d[keep],
    stat = st$stat, p = st$p, logOR = st$logOR,
    pos_rate = a[keep] / (a[keep] + c[keep]),
    neg_rate = b[keep] / (b[keep] + d[keep])
  ) %>%
    arrange(desc(stat), p, desc(abs(logOR))) %>%   # χ² 큰 순(= p 작은 순과 동일)
    slice_head(n = TOP_N) %>%
    mutate(class = cls)
  
  marker_list[[cls]] <- df
}

markers_fullMD <- bind_rows(marker_list)

# ---------- 저장 ----------
saveRDS(markers_fullMD, file.path(OUT_DIR, sprintf("%s_markers_fullMD.rds", FNAME)))
write_csv(markers_fullMD, file.path(OUT_DIR, sprintf("%s_markers_fullMD.csv", FNAME)))

# (옵션) 클래스별/유니온 패널도 함께 저장
panel_union <- tibble(gene = unique(markers_fullMD$gene))
write_csv(panel_union, file.path(OUT_DIR, sprintf("%s_marker_union_fullMD.csv", FNAME)))

message("✓ Full-MD marker selection done. Genes: ", nrow(panel_union),
        "  →  ", file.path(OUT_DIR, sprintf("%s_markers_fullMD.*", FNAME)))
