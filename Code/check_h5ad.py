import scanpy as sc, pandas as pd, scipy.sparse as sp
from scipy.io import mmwrite
from pathlib import Path

# ---------- 사용자 경로 ----------
DATA_DIR = Path(r"C:/test/scRNAseq_ensemble/CTISL/data/10Xv2")  # 폴더 이름 = FNAME
RAW_ROOT = Path(r"C:/test/scRNAseq_ensemble/Data/raw")
# ----------------------------------

FNAME = DATA_DIR.name
h5ad_file  = DATA_DIR / f"{FNAME}.h5ad"
label_file = DATA_DIR / f"{FNAME}label.txt"
out_dir    = RAW_ROOT / FNAME
out_dir.mkdir(parents=True, exist_ok=True)

def read_labels(label_path):
    with open(label_path, encoding="utf-8") as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]
    if lines and lines[0].lower() in {"x", "cell_type", "label"}:
        lines = lines[1:]
    return pd.Series(lines, name="cell_type")

adata = sc.read_h5ad(h5ad_file)
adata.var_names_make_unique()
adata.obs_names.name = adata.obs_names.name or "cell_id"
adata.var_names.name = adata.var_names.name or "gene"

labels = read_labels(label_file)
assert len(labels) == adata.n_obs, "라벨 개수 ≠ 셀 개수!"
adata.obs["cell_type"] = labels.values

mat  = adata.layers["counts"] if "counts" in adata.layers else (
       adata.raw.X if adata.raw is not None else adata.X)
kind = "counts" if "counts" in adata.layers or adata.raw is not None else "expr"
if not sp.issparse(mat):
    mat = sp.csr_matrix(mat)

mmwrite(out_dir / f"{FNAME}_{kind}.mtx", mat.T)
pd.Series(adata.var_names).to_csv(out_dir / f"{FNAME}_genes.tsv",
                                  sep="\t", index=False, header=False)
pd.Series(adata.obs_names).to_csv(out_dir / f"{FNAME}_barcodes.tsv",
                                  sep="\t", index=False, header=False)
adata.obs.reset_index().to_csv(out_dir / f"{FNAME}_meta.csv", index=False)

print("✓ saved to", out_dir)
