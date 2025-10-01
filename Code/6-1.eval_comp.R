# ======================================================
# 6-1.eval_comp.R — Compare Stacking vs Base(LR/SVM): UMAP & per-class F1
# ======================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(Seurat); library(tibble)
})

# ---------- 사용자 파라미터 ----------
FNAME    <- "10Xv2"
PROC_DIR <- file.path("C:/test/scRNAseq_ensemble/Data/processed", FNAME)
PRED_DIR <- file.path(PROC_DIR, "predictions")
L0_DIR   <- file.path(PROC_DIR, "level0")
OUT_DIR  <- file.path(PROC_DIR, "eval_comp"); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------- 입력 파일 경로 ----------
seu_file    <- file.path(PROC_DIR, paste0(FNAME, "_seurat_sct.rds"))
split_file  <- file.path(PROC_DIR, paste0(FNAME, "_splits.rds"))

# stacking (from step 5)
stack_pred_csv <- file.path(PRED_DIR, sprintf("%s_test_predictions.csv", FNAME))
stack_prob_rds <- file.path(PRED_DIR, sprintf("%s_test_meta_probs_df.rds", FNAME))
stack_prob_csv <- file.path(PRED_DIR, sprintf("%s_test_meta_probs.csv", FNAME))

# base refit probs (from step 5)
lr_prob_rds  <- file.path(L0_DIR, "test_pred_lr_marker_refit.rds")
svm_prob_rds <- file.path(L0_DIR, "test_pred_svm_marker_refit.rds")

# ---------- 로드 ----------
seu    <- readRDS(seu_file)
splits <- readRDS(split_file)
meta   <- seu@meta.data
stopifnot("cell_type" %in% colnames(meta))

test_ids <- splits$test_ids
classes  <- levels(factor(meta$cell_type))

# ---------- 유틸: 확률행렬 → 하드라벨 ----------
hard_label_from_probs <- function(M, classes, row_ids = NULL) {
  stopifnot(!is.null(colnames(M)))
  use_classes <- intersect(classes, colnames(M))
  if (length(use_classes) == 0) stop("클래스 교집합이 없습니다.")
  M <- as.matrix(M[, use_classes, drop = FALSE])
  # 행 정렬
  if (!is.null(row_ids)) {
    common <- intersect(row_ids, rownames(M))
    M <- M[common, , drop = FALSE]
  }
  stopifnot(!is.null(rownames(M)))
  idx <- max.col(M, ties.method = "first")
  tibble(cell_id = rownames(M),
         pred    = use_classes[idx])
}

# ---------- 1) stacking pred_df 만들기 ----------
if (file.exists(stack_pred_csv)) {
  stack_pred <- readr::read_csv(stack_pred_csv, show_col_types = FALSE) %>%
    dplyr::select(cell_id, pred)
} else {
  # 확률에서 즉석 생성
  if (file.exists(stack_prob_rds)) {
    probs <- readRDS(stack_prob_rds)           # data.frame
    if (is.null(rownames(probs))) {
      # rownames가 없으면 CSV에서 복구
      if (file.exists(stack_prob_csv)) {
        probs_raw <- readr::read_csv(stack_prob_csv, show_col_types = FALSE)
        rn <- dplyr::coalesce(probs_raw$cell_id, probs_raw[[1]])
        probs <- as.data.frame(probs_raw[, setdiff(names(probs_raw), "cell_id"), drop = FALSE])
        rownames(probs) <- rn
      } else stop("Stacking 확률에서 rownames를 찾을 수 없습니다.")
    }
  } else if (file.exists(stack_prob_csv)) {
    probs_raw <- readr::read_csv(stack_prob_csv, show_col_types = FALSE)
    rn <- dplyr::coalesce(probs_raw$cell_id, probs_raw[[1]])
    probs <- as.data.frame(probs_raw[, setdiff(names(probs_raw), "cell_id"), drop = FALSE])
    rownames(probs) <- rn
  } else stop("Stacking 예측/확률 파일을 찾을 수 없습니다.")
  
  probs <- probs[intersect(test_ids, rownames(probs)), , drop = FALSE]
  stack_pred <- hard_label_from_probs(probs, classes)
  readr::write_csv(stack_pred, stack_pred_csv)
}

# ---------- 2) base(LR/SVM) 하드라벨 ----------
stopifnot(file.exists(lr_prob_rds), file.exists(svm_prob_rds))
lr_probs  <- readRDS(lr_prob_rds)    # N_test x K (matrix)
svm_probs <- readRDS(svm_prob_rds)   # N_test x K (matrix)

lr_probs  <- lr_probs [intersect(test_ids, rownames(lr_probs)),  , drop = FALSE]
svm_probs <- svm_probs[intersect(test_ids, rownames(svm_probs)), , drop = FALSE]

lr_pred  <- hard_label_from_probs(lr_probs,  classes)
svm_pred <- hard_label_from_probs(svm_probs, classes)

# ---------- 3) truth ----------
truth_df <- tibble(
  cell_id = test_ids,
  truth   = as.character(meta[test_ids, "cell_type"])
)

# ---------- 4) per-class F1 계산 함수 (옵션 A) ----------
per_class_f1_df <- function(df, classes) {
  df <- tibble(truth = factor(df$truth, levels = classes),
               pred  = factor(df$pred,  levels = classes))
  mat <- as.data.frame.matrix(table(df$truth, df$pred))  # zero-filled by factor levels
  out <- lapply(classes, function(cls){
    TP <- mat[cls, cls]
    FN <- sum(mat[cls, ]) - TP
    FP <- sum(mat[, cls]) - TP
    precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
    recall    <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
    f1        <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
    tibble(class = cls, TP=TP, FP=FP, FN=FN,
           precision=precision, recall=recall, f1=f1)
  })
  dplyr::bind_rows(out)
}

# ---------- 5) 각 모델의 per-class F1 ----------
stack_f1 <- truth_df %>%
  dplyr::inner_join(stack_pred, by = "cell_id") %>%
  per_class_f1_df(classes) %>%
  dplyr::mutate(model = "Stacking")

lr_f1 <- truth_df %>%
  dplyr::inner_join(lr_pred, by = "cell_id") %>%
  per_class_f1_df(classes) %>%
  dplyr::mutate(model = "Base-LR")

svm_f1 <- truth_df %>%
  dplyr::inner_join(svm_pred, by = "cell_id") %>%
  per_class_f1_df(classes) %>%
  dplyr::mutate(model = "Base-SVM")

f1_long <- dplyr::bind_rows(stack_f1, lr_f1, svm_f1) %>%
  dplyr::mutate(class = factor(class, levels = classes),
                model = factor(model, levels = c("Base-LR", "Base-SVM", "Stacking")))

readr::write_csv(f1_long, file.path(OUT_DIR, sprintf("%s_per_class_f1_compare.csv", FNAME)))

# ---------- 6) per-class F1 비교 그래프 ----------
p_f1_comp <- ggplot(f1_long, aes(x = class, y = f1, fill = model)) +
  geom_col(position = position_dodge(width = 0.75)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = sprintf("%s — Per-class F1: Base(LR/SVM) vs Stacking", FNAME),
       x = NULL, y = "F1")
ggsave(file.path(OUT_DIR, sprintf("%s_per_class_f1_compare.png", FNAME)),
       p_f1_comp, width = 10, height = 5, dpi = 300)

# ---------- 7) UMAP: Truth / Stacking / LR / SVM (각각 별도 PNG) ----------
# Seurat meta에 예측 라벨 주입
seu@meta.data$pred_label_stack <- NA_character_
seu@meta.data$pred_label_lr    <- NA_character_
seu@meta.data$pred_label_svm   <- NA_character_

seu@meta.data[stack_pred$cell_id, "pred_label_stack"] <- stack_pred$pred
seu@meta.data[lr_pred$cell_id,    "pred_label_lr"]    <- lr_pred$pred
seu@meta.data[svm_pred$cell_id,   "pred_label_svm"]   <- svm_pred$pred

# factor level 통일 (색상 일관성)
seu@meta.data$cell_type        <- factor(seu@meta.data$cell_type,        levels = classes)
seu@meta.data$pred_label_stack <- factor(seu@meta.data$pred_label_stack, levels = classes)
seu@meta.data$pred_label_lr    <- factor(seu@meta.data$pred_label_lr,    levels = classes)
seu@meta.data$pred_label_svm   <- factor(seu@meta.data$pred_label_svm,   levels = classes)

# Truth
p_umap_truth <- DimPlot(seu, reduction = "umap", group.by = "cell_type",
                        cells = test_ids, label = FALSE) +
  ggtitle("UMAP — Truth (Test only)")
ggsave(file.path(OUT_DIR, sprintf("%s_umap_truth.png", FNAME)),
       p_umap_truth, width = 8, height = 6, dpi = 300)

# Stacking
p_umap_stack <- DimPlot(seu, reduction = "umap", group.by = "pred_label_stack",
                        cells = test_ids, label = FALSE) +
  ggtitle("UMAP — Stacking (Test only)")
ggsave(file.path(OUT_DIR, sprintf("%s_umap_stack.png", FNAME)),
       p_umap_stack, width = 8, height = 6, dpi = 300)

# Base LR
p_umap_lr <- DimPlot(seu, reduction = "umap", group.by = "pred_label_lr",
                     cells = test_ids, label = FALSE) +
  ggtitle("UMAP — Base LR (Test only)")
ggsave(file.path(OUT_DIR, sprintf("%s_umap_base_lr.png", FNAME)),
       p_umap_lr, width = 8, height = 6, dpi = 300)

# Base SVM
p_umap_svm <- DimPlot(seu, reduction = "umap", group.by = "pred_label_svm",
                      cells = test_ids, label = FALSE) +
  ggtitle("UMAP — Base SVM (Test only)")
ggsave(file.path(OUT_DIR, sprintf("%s_umap_base_svm.png", FNAME)),
       p_umap_svm, width = 8, height = 6, dpi = 300)

message("✓ Saved per-class F1 comparison and UMAPs → ", OUT_DIR)
