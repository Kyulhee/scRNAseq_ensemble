# ======================================================
# level0_marker_models.R — 1-vs-rest (marker=Full-MD 고정) OOF: LR & SVM
# ======================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(readr)
  library(glmnet); library(e1071); library(purrr); library(tidyr)
})

# ---------- 사용자 파라미터 ----------
FNAME      <- "10Xv2"
PROC_DIR   <- file.path("C:/test/scRNAseq_ensemble/Data/processed", FNAME)
M_DIR      <- file.path(PROC_DIR, "markers")
OUT_DIR    <- file.path(PROC_DIR, "level0"); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
MODEL_DIR  <- file.path(OUT_DIR, "models_marker"); dir.create(MODEL_DIR, showWarnings = FALSE)

USE_WEIGHTS <- FALSE   # 클래스 불균형 가중 (원하면 TRUE)
RIDGE_ALPHA <- 0       # glmnet ridge (alpha=0)
set.seed(42)

# ---------- 데이터 로드 ----------
seu    <- readRDS(file.path(PROC_DIR, paste0(FNAME, "_seurat_sct.rds")))
splits <- readRDS(file.path(PROC_DIR, paste0(FNAME, "_splits.rds")))
stopifnot("cell_type" %in% colnames(seu@meta.data))

meta    <- seu@meta.data %>% mutate(cell_id = rownames(.))
classes <- levels(factor(meta$cell_type))

# 표현 행렬 (SCTransform: slot="data"), gene x cell (dgCMatrix)
X_all <- GetAssayData(seu, slot = "data")
stopifnot(all(colnames(X_all) == rownames(meta)))

# ---------- Full-MD 마커 읽기 (모든 fold에 공통 사용) ----------
mk_full_path <- file.path(M_DIR, sprintf("%s_markers_fullMD.rds", FNAME))
stopifnot(file.exists(mk_full_path))
mk_full <- readRDS(mk_full_path)
mk_list <- split(mk_full$gene, mk_full$class)

# ---------- MD(Train) 서브셋 준비 ----------
train_ids <- splits$train_ids
test_ids  <- splits$test_ids  # 여기서는 사용하지 않음(OOF만 생성)
train_idx <- match(train_ids, colnames(X_all))
stopifnot(!any(is.na(train_idx)))

X_md <- X_all[, train_idx, drop = FALSE]
y_md <- factor(meta$cell_type[train_idx], levels = classes)

# ---------- OOF 저장 객체 (MD 크기) ----------
N <- ncol(X_md); T <- length(classes)
oof_lr  <- matrix(NA_real_, nrow = N, ncol = T, dimnames = list(train_ids, classes))
oof_svm <- matrix(NA_real_, nrow = N, ncol = T, dimnames = list(train_ids, classes))

# ---------- 유틸 ----------
mk_weights <- function(y_bin) {
  if (!USE_WEIGHTS) return(NULL)
  pos <- sum(y_bin == 1); neg <- sum(y_bin == 0)
  if (pos == 0 || neg == 0) return(NULL)
  ifelse(y_bin == 1, 0.5/pos, 0.5/neg)
}

# ---------- MD 내부 CV folds ----------
vfold <- splits$vfold$splits  # 이미 train에서 만든 3-fold

for (i in seq_along(vfold)) {
  sp <- vfold[[i]]
  tr_ids <- rsample::analysis(sp)$cell_id    # MD 내부의 fold-train ids
  va_ids <- rsample::assessment(sp)$cell_id  # MD 내부의 fold-valid ids
  
  # MD 인덱스 기준으로 매핑
  tr_idx <- match(tr_ids, colnames(X_md))
  va_idx <- match(va_ids, colnames(X_md))
  message(sprintf("[OOF] Fold %d: train=%d, valid=%d", i, length(tr_idx), length(va_idx)))
  
  for (cls in classes) {
    genes_cls <- intersect(mk_list[[cls]], rownames(X_md))
    if (length(genes_cls) == 0) {
      warning(sprintf("Fold %d, class %s: no markers in assay. Skipped.", i, cls)); next
    }
    
    Xtr <- t(X_md[genes_cls, tr_idx, drop = FALSE])  # cell × gene
    Xva <- t(X_md[genes_cls, va_idx, drop = FALSE])
    ytr <- as.integer(y_md[tr_idx] == cls)
    
    w <- mk_weights(ytr)
    
    # ---- LR (glmnet, ridge, 5-fold 내부 CV) ----
    cvfit <- cv.glmnet(x = Xtr, y = ytr, family = "binomial",
                       alpha = RIDGE_ALPHA, standardize = TRUE,
                       weights = w, nfolds = 5, type.measure = "deviance")
    prob_lr <- as.numeric(predict(cvfit, newx = Xva, type = "response", s = "lambda.min"))
    
    # ---- SVM (RBF, prob=TRUE) ----
    ytr_fac <- factor(ytr, levels = c(0,1))
    svmfit <- e1071::svm(x = as.matrix(Xtr), y = ytr_fac,
                         kernel = "radial", probability = TRUE, scale = TRUE)
    pr_va <- predict(svmfit, as.matrix(Xva), probability = TRUE)
    pr_mat <- attr(pr_va, "probabilities")
    prob_svm <- if (!is.null(pr_mat) && "1" %in% colnames(pr_mat)) pr_mat[, "1"] else
      if (!is.null(pr_mat)) pr_mat[, ncol(pr_mat)] else
        as.numeric(pr_va)
    
    # ---- OOF 저장 (MD 행렬에 채움) ----
    oof_lr [va_ids, cls]  <- prob_lr
    oof_svm[va_ids, cls]  <- prob_svm
    
    # (선택 저장: 디버깅/재현용)
    saveRDS(list(class = cls, fold = i, genes = genes_cls, cvfit = cvfit),
            file = file.path(MODEL_DIR, sprintf("lr_oof_fold%02d_%s.rds", i, cls)))
    saveRDS(list(class = cls, fold = i, genes = genes_cls, svmfit = svmfit),
            file = file.path(MODEL_DIR, sprintf("svm_oof_fold%02d_%s.rds", i, cls)))
  }
}

# ---------- 산출 ----------
saveRDS(oof_lr,  file = file.path(OUT_DIR, "oof_pred_lr_marker.rds"))
saveRDS(oof_svm, file = file.path(OUT_DIR, "oof_pred_svm_marker.rds"))

message(sprintf("OOF saved: LR non-NA=%d / SVM non-NA=%d  (N_MD=%d, T=%d)",
                sum(!is.na(oof_lr)), sum(!is.na(oof_svm)), nrow(oof_lr), ncol(oof_lr)))
message("→ ", OUT_DIR)
