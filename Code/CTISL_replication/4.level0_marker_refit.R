# ======================================================
# 60_level0_marker_refit.R — Refit base (LR/SVM) on FULL MD markers
# ======================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(readr)
  library(glmnet); library(e1071)
})

# ---------- 사용자 파라미터 ----------
FNAME      <- "10Xv2"
PROC_DIR   <- file.path("C:/test/scRNAseq_ensemble/Data/processed", FNAME)
M_DIR      <- file.path(PROC_DIR, "markers")
OUT_DIR    <- file.path(PROC_DIR, "level0"); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
MODEL_DIR  <- file.path(OUT_DIR, "models_marker_refit"); dir.create(MODEL_DIR, showWarnings = FALSE)

USE_WEIGHTS <- FALSE      # 클래스 불균형 가중
RIDGE_ALPHA <- 0          # glmnet ridge
NFOILDS     <- 5
set.seed(42)

# ---------- 로드 ----------
seu    <- readRDS(file.path(PROC_DIR, paste0(FNAME, "_seurat_sct.rds")))
splits <- readRDS(file.path(PROC_DIR, paste0(FNAME, "_splits.rds")))
stopifnot("cell_type" %in% colnames(seu@meta.data))

meta    <- seu@meta.data %>% mutate(cell_id = rownames(.))
classes <- levels(factor(meta$cell_type))
X_all   <- GetAssayData(seu, slot = "data")  # genes x cells
stopifnot(all(colnames(X_all) == rownames(meta)))

# MD(train) 서브셋
train_ids <- splits$train_ids
tr_idx    <- match(train_ids, colnames(X_all))
stopifnot(!any(is.na(tr_idx)))

X_md <- X_all[, tr_idx, drop = FALSE]
y_md <- factor(meta$cell_type[tr_idx], levels = classes)

# Full-MD 마커
mk_full_path <- file.path(M_DIR, sprintf("%s_markers_fullMD.rds", FNAME))
stopifnot(file.exists(mk_full_path))
mk_full <- readRDS(mk_full_path)
mk_list <- split(mk_full$gene, mk_full$class)

mk_weights <- function(y){
  if (!USE_WEIGHTS) return(NULL)
  pos <- sum(y == 1); neg <- sum(y == 0)
  if (pos == 0 || neg == 0) return(NULL)
  ifelse(y == 1, 0.5/pos, 0.5/neg)
}

# ---------- 클래스별 refit ----------
for (cls in classes) {
  genes_cls <- intersect(mk_list[[cls]], rownames(X_md))
  if (length(genes_cls) == 0) {
    warning(sprintf("[REFIT] class %s: no markers in assay. Skipped.", cls)); next
  }
  
  Xtr <- t(X_md[genes_cls, , drop = FALSE])   # cells x genes
  ytr <- as.integer(y_md == cls)
  w   <- mk_weights(ytr)
  
  # ---- LR (ridge, CV) ----
  cvfit <- cv.glmnet(x = Xtr, y = ytr, family = "binomial",
                     alpha = RIDGE_ALPHA, standardize = TRUE,
                     weights = w, nfolds = NFOILDS, type.measure = "deviance")
  saveRDS(list(class = cls, genes = genes_cls, cvfit = cvfit),
          file = file.path(MODEL_DIR, sprintf("lr_full_%s.rds", cls)))
  
  # ---- SVM (RBF, prob=TRUE) ----
  ytr_fac <- factor(ytr, levels = c(0,1))
  svmfit <- e1071::svm(x = as.matrix(Xtr), y = ytr_fac,
                       kernel = "radial", probability = TRUE, scale = TRUE)
  saveRDS(list(class = cls, genes = genes_cls, svmfit = svmfit),
          file = file.path(MODEL_DIR, sprintf("svm_full_%s.rds", cls)))
  
  message(sprintf("[REFIT] class=%s  | genes=%d", cls, length(genes_cls)))
}

message("✓ Base refit models saved → ", MODEL_DIR)
