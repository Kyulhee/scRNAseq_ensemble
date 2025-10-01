# ======================================================
# 80_eval_visualize_test.R — Test 평가/시각화 (CTISL 재현 평가)
# ======================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(Seurat); library(patchwork); library(tibble)
})

# ---------- 사용자 파라미터 ----------
FNAME    <- "10Xv2"
PROC_DIR <- file.path("C:/test/scRNAseq_ensemble/Data/processed", FNAME)
PRED_DIR <- file.path(PROC_DIR, "predictions")
OUT_DIR  <- file.path(PROC_DIR, "eval"); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------- 입력 경로 ----------
seu_file   <- file.path(PROC_DIR, paste0(FNAME, "_seurat_sct.rds"))
split_file <- file.path(PROC_DIR, paste0(FNAME, "_splits.rds"))
pred_csv   <- file.path(PRED_DIR, sprintf("%s_test_predictions.csv", FNAME))         # cell_id,pred,conf,second,margin
prob_rds   <- file.path(PRED_DIR, sprintf("%s_test_meta_probs_df.rds", FNAME))       # data.frame rownames=cell_id
prob_csv   <- file.path(PRED_DIR, sprintf("%s_test_meta_probs.csv", FNAME))          # CSV with cell_id col

# ---------- 로드 ----------
seu    <- readRDS(seu_file)
splits <- readRDS(split_file)

# 메타데이터(정답 라벨)
meta <- seu@meta.data
stopifnot("cell_type" %in% colnames(meta))
test_ids <- splits$test_ids

# ---------- pred_df 생성 (우선순위: pred_csv → prob_rds → prob_csv) ----------
if (file.exists(pred_csv)) {
  pred_df <- readr::read_csv(pred_csv, show_col_types = FALSE) %>%
    dplyr::select(cell_id, pred)
} else {
  # 확률 행렬 읽기
  if (file.exists(prob_rds)) {
    probs <- readRDS(prob_rds)              # data.frame (rownames=cell_id, cols=classes)
  } else if (file.exists(prob_csv)) {
    probs_raw <- readr::read_csv(prob_csv, show_col_types = FALSE)
    rn <- dplyr::coalesce(probs_raw$cell_id, probs_raw[[1]])
    probs <- as.data.frame(probs_raw[, setdiff(names(probs_raw), "cell_id"), drop = FALSE])
    rownames(probs) <- rn
  } else {
    stop("test_predictions.csv, test_meta_probs_df.rds, test_meta_probs.csv 중 아무 것도 찾지 못했습니다.")
  }
  
  # test_ids 순서로 subset
  common <- intersect(test_ids, rownames(probs))
  if (length(common) == 0) stop("확률 테이블에서 test_ids를 찾지 못했습니다.")
  probs  <- probs[common, , drop = FALSE]
  
  # 클래스 목록
  classes <- levels(factor(meta$cell_type))
  pred_classes <- colnames(probs)
  if (is.null(pred_classes)) stop("확률 테이블의 열 이름(클래스)이 없습니다.")
  use_classes <- intersect(classes, pred_classes)
  if (length(use_classes) == 0) stop("클래스 교집합이 없습니다.")
  
  # 하드 라벨 산출
  M <- as.matrix(probs[, use_classes, drop = FALSE])
  top_idx <- max.col(M, ties.method = "first")
  pred_df <- tibble(cell_id = rownames(M),
                    pred    = use_classes[top_idx])
  
  # 재사용 편의를 위해 저장
  readr::write_csv(pred_df, pred_csv)
}

# ---------- truth_df 구성 ----------
truth_df <- tibble(
  cell_id = test_ids,
  truth   = as.character(meta[test_ids, "cell_type"])
)

eval_df <- truth_df %>%
  inner_join(pred_df, by = "cell_id") %>%
  mutate(correct = truth == pred)

# ---------- 1) 교차표 (counts / row-normalized) ----------
classes <- levels(factor(meta$cell_type))
eval_df <- eval_df %>%
  mutate(truth = factor(truth, levels = classes),
         pred  = factor(pred,  levels = classes))

tab_counts <- as.data.frame.matrix(table(eval_df$truth, eval_df$pred))
# row-normalized (정답 기준 비율)
tab_rowprop <- prop.table(as.matrix(tab_counts), margin = 1) %>% as.data.frame.matrix()

# 저장 (※ rownames_to_tibble -> rownames_to_column 이 정답)
tab_counts_df  <- tab_counts  %>% rownames_to_column("truth")
tab_rowprop_df <- tab_rowprop %>% rownames_to_column("truth")
write_csv(tab_counts_df,  file.path(OUT_DIR, sprintf("%s_contingency_counts.csv", FNAME)))
write_csv(tab_rowprop_df, file.path(OUT_DIR, sprintf("%s_contingency_rowprop.csv", FNAME)))

# ---------- 2) 정량 지표 (accuracy, macro/micro F1) ----------
accuracy <- mean(eval_df$correct)

# per-class metrics
per_class <- eval_df %>%
  count(truth, pred, name = "n") %>%
  complete(truth = classes, pred = classes, fill = list(n = 0)) %>%
  group_by(truth) %>%
  summarise(
    TP = sum(n[truth == pred]),
    FN = sum(n[truth != pred]),
    .groups = "drop"
  ) %>%
  # FP는 pred 기준으로 다시 집계
  left_join(
    eval_df %>%
      count(pred, truth, name = "n2") %>%
      complete(pred = classes, truth = classes, fill = list(n2 = 0)) %>%
      group_by(pred) %>% summarise(FP = sum(n2[pred != truth]), .groups = "drop"),
    by = c("truth" = "pred")
  ) %>%
  mutate(
    precision = ifelse((TP + FP) > 0, TP / (TP + FP), 0),
    recall    = ifelse((TP + FN) > 0, TP / (TP + FN), 0),
    f1        = ifelse((precision + recall) > 0, 2 * precision * recall / (precision + recall), 0)
  ) %>%
  rename(class = truth)

macro_f1 <- mean(per_class$f1)
# micro-F1 (단일 라벨 다중분류에서는 정확도와 동일)
sum_TP <- sum(per_class$TP); sum_FP <- sum(per_class$FP); sum_FN <- sum(per_class$FN)
micro_precision <- ifelse(sum_TP + sum_FP > 0, sum_TP / (sum_TP + sum_FP), 0)
micro_recall    <- ifelse(sum_TP + sum_FN > 0, sum_TP / (sum_TP + sum_FN), 0)
micro_f1        <- ifelse(micro_precision + micro_recall > 0,
                          2 * micro_precision * micro_recall / (micro_precision + micro_recall), 0)

metrics_summary <- tibble(
  metric = c("accuracy", "macro_f1", "micro_f1"),
  value  = c(accuracy, macro_f1, micro_f1)
)
write_csv(metrics_summary, file.path(OUT_DIR, sprintf("%s_metrics_summary.csv", FNAME)))
write_csv(per_class,      file.path(OUT_DIR, sprintf("%s_metrics_per_class.csv", FNAME)))

# ---------- 3) Confusion heatmap (정답 기준 정규화) ----------
heat_df <- tab_rowprop %>%
  rownames_to_column("truth") %>%
  pivot_longer(-truth, names_to = "pred", values_to = "prop") %>%
  mutate(truth = factor(truth, levels = classes),
         pred  = factor(pred,  levels = classes))

p_heat <- ggplot(heat_df, aes(pred, truth, fill = prop)) +
  geom_tile(color = "white") +
  scale_fill_gradient(limits = c(0, 1), name = "Row prop") +
  labs(x = "Predicted", y = "Truth",
       title = sprintf("%s — Confusion (row-normalized)", FNAME)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, sprintf("%s_confusion_heatmap.png", FNAME)), p_heat, width = 7.5, height = 6, dpi = 300)

# ---------- 4) Per-class F1 barplot ----------
p_f1 <- per_class %>%
  mutate(class = factor(class, levels = classes)) %>%
  ggplot(aes(x = class, y = f1)) +
  geom_col() +
  geom_hline(yintercept = macro_f1, linetype = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = sprintf("%s — Per-class F1 (macro=%.3f, acc=%.3f)", FNAME, macro_f1, accuracy),
       x = NULL, y = "F1")
ggsave(file.path(OUT_DIR, sprintf("%s_per_class_f1_bar.png", FNAME)), p_f1, width = 7.5, height = 4.5, dpi = 300)

# ---------- 5) UMAP: Truth vs Pred on Test cells ----------
# 메타데이터에 예측 라벨 넣기(테스트 셀만)
pred_map <- pred_df %>% select(cell_id, pred)
seu@meta.data$pred_label <- NA_character_
seu@meta.data[ pred_map$cell_id, "pred_label"] <- pred_map$pred

# Test 셀만 시각화
p_umap_truth <- DimPlot(seu, reduction = "umap", group.by = "cell_type",
                        cells = test_ids, label = FALSE) + ggtitle("UMAP — Truth (Test only)")
p_umap_pred  <- DimPlot(seu, reduction = "umap", group.by = "pred_label",
                        cells = test_ids, label = FALSE) + ggtitle("UMAP — Pred (Test only)")

#p_umap <- p_umap_truth + p_umap_pred + plot_layout(ncol = 2)
ggsave(file.path(OUT_DIR, sprintf("%s_umap_pred.png", FNAME)),
       p_umap_pred, width = 8, height = 6, dpi = 300)
ggsave(file.path(OUT_DIR, sprintf("%s_umap_truth.png", FNAME)),
       p_umap_truth, width = 8, height = 6, dpi = 300)


message("✓ Saved: contingency CSVs, metrics, heatmap, F1 bars, UMAP — → ", OUT_DIR)
