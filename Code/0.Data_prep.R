# ======================================================
# 01_preprocess_flex.R   —  Count MTX + tsv 세트 전처리
# ======================================================

# 패키지 로드 -------------------------------------------------------------
suppressPackageStartupMessages({
  library(Matrix);  library(Seurat);  library(dplyr);  library(readr)
  library(rsample); library(stringr); library(ggplot2); library(patchwork)
})

# ---------- 사용자 파라미터 ----------
FNAME      <- "10Xv2"          # 바꿔가며 사용
path <- "C:/test/scRNAseq_ensemble/Data/"
INPUT_DIR  <- paste0(path,"raw/", FNAME)
OUTPUT_DIR <- paste0(path,"processed/", FNAME)


SPECIES    <- "human"          # ^MT- vs ^mt-
SEED       <- 42
MIN_GENES  <- 200
#MAX_GENES <- 6000
#MAX_MT_PCT <- 15
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------- 파일 경로 ----------
mtx_file  <- file.path(INPUT_DIR, paste0(FNAME, "_counts.mtx"))
if (!file.exists(mtx_file))          # counts가 없으면 expr 사용
  mtx_file <- file.path(INPUT_DIR, paste0(FNAME, "_expr.mtx"))
genes_tsv <- file.path(INPUT_DIR, paste0(FNAME, "_genes.tsv"))
bars_tsv  <- file.path(INPUT_DIR, paste0(FNAME, "_barcodes.tsv"))
meta_csv  <- file.path(INPUT_DIR, paste0(FNAME, "_meta.csv"))

stopifnot(file.exists(mtx_file), file.exists(genes_tsv),
          file.exists(bars_tsv), file.exists(meta_csv))

# ---------- 1) 데이터 읽기 ----------
counts <- readMM(mtx_file)
genes  <- read_tsv(genes_tsv, col_names = FALSE, show_col_types = FALSE)$X1
cells  <- read_tsv(bars_tsv,  col_names = FALSE, show_col_types = FALSE)$X1
meta   <- read_csv(meta_csv,  show_col_types = FALSE)

MIN_CELLS <- 3   # CTISL 기본

keep_genes <- Matrix::rowSums(counts > 0) >= MIN_CELLS
counts     <- counts[keep_genes, , drop = FALSE]
length(keep_genes)
genes <- genes[keep_genes]

message(sprintf("After gene filter (≥%d cells): %d genes remain",
                MIN_CELLS, nrow(counts)))

rownames(counts) <- genes
colnames(counts) <- cells
if (!"cell_id" %in% names(meta)) meta$cell_id <- meta[[1]]
meta <- meta %>% distinct(cell_id, .keep_all = TRUE)
meta <- meta %>% filter(cell_id %in% cells)
rownames(meta) <- meta$cell_id
meta <- meta[cells, , drop = FALSE]   # order to counts

message(sprintf("Loaded matrix: %d genes × %d cells", nrow(counts), ncol(counts)))

# ---------- 2) Seurat 객체 & QC ----------
seu <- CreateSeuratObject(counts = counts, meta.data = meta, project = FNAME)

# 미토 퍼센트
pat <- if (tolower(SPECIES) == "human") "^MT-" else "^mt-"
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = pat)

# 기본 QC violin 저장
png(file.path(OUTPUT_DIR, paste0(FNAME, "_qc_violin.png")), width=900, height=350)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0.1) + NoLegend()
dev.off()

# 필터링
seu <- subset(seu,
              subset = nFeature_RNA >= MIN_GENES)
                # nFeature_RNA <= MAX_GENES &
                # percent.mt   <= MAX_MT_PCT)
message(sprintf("After QC: %d genes × %d cells", nrow(seu), ncol(seu)))

# ---------- 3) 정규화 & 차원축소 ----------
seu <- SCTransform(seu, vars.to.regress = "percent.mt", verbose = FALSE)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)

png(file.path(OUTPUT_DIR, paste0(FNAME, "_umap.png")), width=700, height=500)
DimPlot(seu, reduction = "umap", group.by = if ("cell_type" %in% colnames(seu@meta.data)) "cell_type" else "orig.ident",
        label = TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()

# ---------- 4) Train/Test split & 3-fold 내부 CV ----------
set.seed(SEED)
meta_split <- seu@meta.data %>% mutate(cell_id = rownames(.))

meta_split

split  <- initial_split(meta_split, prop = 0.8,
                        strata = if ("cell_type" %in% names(meta_split)) "cell_type" else NULL)
train_md <- training(split);    test_md <- testing(split)
vfold    <- vfold_cv(train_md, v = 3,
                     strata = if ("cell_type" %in% names(train_md)) "cell_type" else NULL)

saveRDS(list(train_ids = train_md$cell_id,
             test_ids  = test_md$cell_id,
             vfold     = vfold),
        file = file.path(OUTPUT_DIR, paste0(FNAME, "_splits.rds")))
saveRDS(seu, file = file.path(OUTPUT_DIR, paste0(FNAME, "_seurat_sct.rds")))
saveRDS(GetAssayData(seu, slot = "data"),
        file = file.path(OUTPUT_DIR, paste0(FNAME, "_expr_norm.rds")))

write_csv(tibble(feature = VariableFeatures(seu)),
          file.path(OUTPUT_DIR, paste0(FNAME, "_variable_features.csv")))

table(test_md$cell_type)

message("✓ 전처리 완료 → 출력 폴더: ", OUTPUT_DIR)
