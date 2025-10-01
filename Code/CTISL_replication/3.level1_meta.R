# ======================================================
# 50_level1_meta.R — Train the meta (stacking) model on OOF predictions
# ======================================================
suppressPackageStartupMessages({
  library(glmnet); library(Matrix); library(dplyr); library(readr)
})

# ---------- 사용자 파라미터 ----------
FNAME     <- "10Xv2"
PROC_DIR  <- file.path("C:/test/scRNAseq_ensemble/Data/processed", FNAME)
L0_DIR    <- file.path(PROC_DIR, "level0")
OUT_DIR   <- file.path(PROC_DIR, "level1"); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ALPHA_META <- 0      # ridge multinomial
NFOILDS    <- 5
USE_LOGIT  <- FALSE  # 확률 대신 logit 사용하려면 TRUE

# ---------- 로드 (라벨/ID) ----------
seu    <- readRDS(file.path(PROC_DIR, paste0(FNAME, "_seurat_sct.rds")))
splits <- readRDS(file.path(PROC_DIR, paste0(FNAME, "_splits.rds")))
meta   <- seu@meta.data %>% mutate(cell_id = rownames(.))
classes <- levels(factor(meta$cell_type))
train_ids <- splits$train_ids

# ---------- OOF 블록 로딩 함수 ----------
read_oof <- function(path, block_name, classes, train_ids) {
  if (!file.exists(path)) return(NULL)
  M <- readRDS(path)                   # N x T, rownames = ids, colnames = classes
  # 정렬/부분집합
  M <- M[train_ids, classes, drop = FALSE]
  # NA 가드: 열 중앙값으로 대체 후, 행별로 0–1 클램프
  if (anyNA(M)) {
    for (j in seq_len(ncol(M))) {
      cj <- M[, j]; med <- median(cj[!is.na(cj)])
      cj[is.na(cj)] <- ifelse(is.finite(med), med, 0)
      M[, j] <- cj
    }
  }
  # (옵션) logit 변환
  if (USE_LOGIT) {
    eps <- 1e-6; M <- pmin(pmax(M, eps), 1 - eps)
    M <- log(M / (1 - M))
  }
  colnames(M) <- paste(block_name, classes, sep = "::")
  M
}

# ---------- OOF 블록 읽기 ----------
blk_list <- list(
  read_oof(file.path(L0_DIR, "oof_pred_lr_marker.rds"),  "marker_LR",  classes, train_ids),
  read_oof(file.path(L0_DIR, "oof_pred_svm_marker.rds"), "marker_SVM", classes, train_ids)
)

# (옵션) 임베딩 브랜치가 있으면 자동 결합
emb_lr  <- read_oof(file.path(L0_DIR, "oof_pred_lr_embed.rds"),  "embed_LR",  classes, train_ids)
emb_svm <- read_oof(file.path(L0_DIR, "oof_pred_svm_embed.rds"), "embed_SVM", classes, train_ids)
if (!is.null(emb_lr))  blk_list <- c(blk_list, list(emb_lr))
if (!is.null(emb_svm)) blk_list <- c(blk_list, list(emb_svm))

# 유효 블록 체크
blk_list <- blk_list[!vapply(blk_list, is.null, logical(1))]
stopifnot(length(blk_list) >= 1)

# ---------- 메타 입력행렬/타깃 구성 ----------
X <- do.call(cbind, blk_list)               # (N_MD) x (B*T)
y <- factor(meta$cell_type[match(train_ids, rownames(meta))], levels = classes)

# 저장(재현/디버깅용)
saveRDS(X, file = file.path(OUT_DIR, "X_oof_meta.rds"))

# ---------- 메타 학습 (multinomial ridge) ----------
cvfit_meta <- cv.glmnet(
  x = as.matrix(X),
  y = y,
  family = "multinomial",
  alpha = ALPHA_META,
  nfolds = NFOILDS,
  type.measure = "deviance",
  standardize = TRUE
)

saveRDS(cvfit_meta, file = file.path(OUT_DIR, "meta_glmnet_multinomial.rds"))

# ---------- 메타데이터(특징 순서/블록) 저장 ----------
feat_spec <- list(
  FNAME        = FNAME,
  classes      = classes,
  train_ids    = train_ids,
  use_logit    = USE_LOGIT,
  feature_names = colnames(X),
  blocks        = names(blk_list) # just to hint how many were combined
)
saveRDS(feat_spec, file = file.path(OUT_DIR, "meta_feature_spec.rds"))

# ---------- 메시지 ----------
message("✓ Meta model trained on OOF. Features: ", ncol(X), "  Samples: ", nrow(X))
message("λ.min = ", tryCatch(cvfit_meta$lambda.min, error = function(e) NA))
message("→ Saved to: ", OUT_DIR)
