# ======================================================
# 70_test_predict_with_meta.R — Use refit base + trained meta to predict Test
# ======================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(readr)
  library(glmnet); library(purrr); library(tidyr)
})

# ---------- 사용자 파라미터 ----------
FNAME      <- "10Xv2"
PROC_DIR   <- file.path("C:/test/scRNAseq_ensemble/Data/processed", FNAME)
L0_DIR     <- file.path(PROC_DIR, "level0")
L1_DIR     <- file.path(PROC_DIR, "level1")
MODEL_DIR  <- file.path(L0_DIR, "models_marker_refit")
OUT_DIR    <- file.path(PROC_DIR, "predictions"); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

USE_LOGIT_INPUT <- FALSE  # meta가 logit으로 학습되었다면 TRUE로

# ---------- 로드 ----------
seu    <- readRDS(file.path(PROC_DIR, paste0(FNAME, "_seurat_sct.rds")))
splits <- readRDS(file.path(PROC_DIR, paste0(FNAME, "_splits.rds")))
meta   <- seu@meta.data %>% mutate(cell_id = rownames(.))
classes <- levels(factor(meta$cell_type))

# 메타 모델 & 스펙
cvfit_meta <- readRDS(file.path(L1_DIR, "meta_glmnet_multinomial.rds"))
feat_spec  <- readRDS(file.path(L1_DIR, "meta_feature_spec.rds"))
feat_names <- feat_spec$feature_names
stopifnot(length(feat_names) > 0)

# Test IDs
test_ids <- splits$test_ids
stopifnot(length(test_ids) > 0)

# 표현 행렬
X_all <- GetAssayData(seu, slot = "data")  # genes x cells
te_idx <- match(test_ids, colnames(X_all))
stopifnot(!any(is.na(te_idx)))
X_te <- X_all[, te_idx, drop = FALSE]

# ---------- Refit base 모델 불러와 Test 확률 생성 ----------
# (현재는 marker_LR / marker_SVM 두 블록만 가정; 임베딩 추가 시 동일 패턴으로 확장)
pred_lr  <- matrix(NA_real_, nrow = length(test_ids), ncol = length(classes),
                   dimnames = list(test_ids, classes))
pred_svm <- matrix(NA_real_, nrow = length(test_ids), ncol = length(classes),
                   dimnames = list(test_ids, classes))

for (cls in classes) {
  # LR
  lr_obj_path <- file.path(MODEL_DIR, sprintf("lr_full_%s.rds", cls))
  if (file.exists(lr_obj_path)) {
    lr_obj <- readRDS(lr_obj_path)
    genes  <- lr_obj$genes
    genes  <- intersect(genes, rownames(X_te))
    if (length(genes)) {
      Xte <- t(X_te[genes, , drop = FALSE])
      pred_lr[, cls] <- as.numeric(predict(lr_obj$cvfit, newx = Xte, type = "response", s = "lambda.min"))
    } else {
      warning(sprintf("LR test: no genes found for class %s", cls))
      pred_lr[, cls] <- 0.0
    }
  } else {
    warning(sprintf("LR model missing for class %s", cls))
    pred_lr[, cls] <- 0.0
  }
  
  # SVM
  svm_obj_path <- file.path(MODEL_DIR, sprintf("svm_full_%s.rds", cls))
  if (file.exists(svm_obj_path)) {
    svm_obj <- readRDS(svm_obj_path)
    genes   <- svm_obj$genes
    genes   <- intersect(genes, rownames(X_te))
    if (length(genes)) {
      Xte <- t(X_te[genes, , drop = FALSE])
      pr  <- predict(svm_obj$svmfit, as.matrix(Xte), probability = TRUE)
      prp <- attr(pr, "probabilities")
      pred_svm[, cls] <- if (!is.null(prp) && "1" %in% colnames(prp)) prp[, "1"] else
        if (!is.null(prp)) prp[, ncol(prp)] else
          as.numeric(pr)
    } else {
      warning(sprintf("SVM test: no genes found for class %s", cls))
      pred_svm[, cls] <- 0.0
    }
  } else {
    warning(sprintf("SVM model missing for class %s", cls))
    pred_svm[, cls] <- 0.0
  }
}

# 저장(선택)
saveRDS(pred_lr,  file = file.path(L0_DIR, "test_pred_lr_marker_refit.rds"))
saveRDS(pred_svm, file = file.path(L0_DIR, "test_pred_svm_marker_refit.rds"))

# ---------- 메타 입력행렬 구성 (블록명::클래스명) ----------
# 메타 학습 때 썼던 네이밍 규칙과 동일하게 컬럼명 생성
colnames(pred_lr)  <- paste("marker_LR",  colnames(pred_lr),  sep = "::")
colnames(pred_svm) <- paste("marker_SVM", colnames(pred_svm), sep = "::")

X_meta_test <- cbind(pred_lr, pred_svm)

# (옵션) 임베딩 브랜치가 있으면 같은 네이밍으로 cbind 후 정렬
# e.g., pred_embed_lr <- readRDS(...); colnames <- paste("embed_LR", classes, sep="::"); X_meta_test <- cbind(X_meta_test, pred_embed_lr)

# **특징 순서 맞추기**
missing_cols <- setdiff(feat_names, colnames(X_meta_test))
if (length(missing_cols)) {
  # 없으면 0으로 채움(혹은 훈련시 중앙값 사용)
  fill <- matrix(0, nrow = nrow(X_meta_test), ncol = length(missing_cols),
                 dimnames = list(rownames(X_meta_test), missing_cols))
  X_meta_test <- cbind(X_meta_test, fill)
}
# 훈련 시 사용한 순서로 정렬
X_meta_test <- X_meta_test[, feat_names, drop = FALSE]

# (옵션) 메타가 logit으로 학습되었으면 동일 변환
if (isTRUE(feat_spec$use_logit) || USE_LOGIT_INPUT) {
  eps <- 1e-6
  X_meta_test <- pmin(pmax(X_meta_test, eps), 1 - eps)
  X_meta_test <- log(X_meta_test / (1 - X_meta_test))
}

# ---------- 메타 예측 (robust reshape) ----------
pr_raw <- predict(
  cvfit_meta,
  newx = as.matrix(X_meta_test),
  s    = "lambda.min",
  type = "response"
)

# 메타 학습때 쓰던 클래스 순서
classes <- if (!is.null(feat_spec$classes)) feat_spec$classes else levels(factor(meta$cell_type))
n_test  <- length(test_ids)
n_cls   <- length(classes)

if (is.array(pr_raw) && length(dim(pr_raw)) == 3L) {
  # pr_raw: [N, K, |s|], 여기선 |s|=1
  dn2 <- dimnames(pr_raw)[[2]]
  pred_classes <- if (!is.null(dn2)) dn2 else classes
  ord <- match(classes, pred_classes)
  
  if (any(is.na(ord))) {
    stop("Class name mismatch. Missing in predict(): ",
         paste(classes[is.na(ord)], collapse=", "))
  }
  
  # 3D → 2D (λ=1 축 제거). drop=TRUE 로 꼭 2D로 만듦
  probs_mat <- pr_raw[, ord, 1, drop = TRUE]  # now N x K matrix
  
} else if (is.matrix(pr_raw)) {
  # 이미 2D인 경우: [N x K] 또는 [K x N]
  if (nrow(pr_raw) == n_cls && ncol(pr_raw) == n_test) {
    probs_mat <- t(pr_raw)
  } else {
    probs_mat <- pr_raw
  }
  
  # 열 이름 보정
  if (is.null(colnames(probs_mat)) || length(colnames(probs_mat)) != n_cls) {
    colnames(probs_mat) <- classes
  }
  
} else {
  # 벡터로 납작해진 드문 케이스 (길이 N*K)
  v <- as.numeric(pr_raw)
  if (length(v) != n_test * n_cls) {
    stop("Unexpected shape from glmnet predict(): length=", length(v),
         " expected=", n_test * n_cls)
  }
  probs_mat <- matrix(v, nrow = n_test, ncol = n_cls, byrow = FALSE)
  colnames(probs_mat) <- classes
}
# --- after you have probs_mat ready ---

# 0) dimnames 안전하게
if (is.null(rownames(probs_mat))) rownames(probs_mat) <- test_ids
if (is.null(colnames(probs_mat))) colnames(probs_mat) <- classes

# 1) 하드 라벨 및 부가 지표
top_idx    <- max.col(probs_mat, ties.method = "first")
pred_label <- colnames(probs_mat)[top_idx]
pred_conf  <- probs_mat[cbind(seq_len(nrow(probs_mat)), top_idx)]

# 2nd best + margin(확신도 격차)
ord2          <- t(apply(probs_mat, 1, order, decreasing = TRUE))
second_idx    <- ord2[, 2]
second_label  <- colnames(probs_mat)[second_idx]
second_conf   <- probs_mat[cbind(seq_len(nrow(probs_mat)), second_idx)]
margin        <- pred_conf - second_conf

pred_tbl <- tibble::tibble(
  cell_id = rownames(probs_mat),
  pred    = pred_label,
  conf    = pred_conf,
  second  = second_label,
  margin  = margin
)

# 2) 저장
readr::write_csv(pred_tbl,
                 file.path(OUT_DIR, sprintf("%s_test_predictions.csv", FNAME))
)

# (이미 만든 probs_df 저장도 유지)
readr::write_csv(tibble::as_tibble(as.data.frame(probs_mat), rownames = "cell_id"),
                 file.path(OUT_DIR, sprintf("%s_test_meta_probs.csv", FNAME))
)
saveRDS(as.data.frame(probs_mat),
        file.path(OUT_DIR, sprintf("%s_test_meta_probs_df.rds", FNAME))
)

message("✓ Saved test_predictions.csv and prob matrices → ", OUT_DIR)
