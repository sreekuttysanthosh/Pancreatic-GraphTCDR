"""
Checks to run BEFORE building GraphTCDR. Reads saved files only; changes nothing except writing
data/processed/integrated_with_auc.csv.

Run:  python pre_graph_checks.py --ae_holdout_mse 0.014964
"""
import argparse, os
import numpy as np, pandas as pd
from scipy.stats import spearmanr

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=r"D:\Pancreatic_GraphTCDR")
ap.add_argument("--ae_holdout_mse", type=float, default=0.014964)   # from your checkpoint (epoch 388)
ap.add_argument("--seed", type=int, default=42)
a = ap.parse_args()
RAW = os.path.join(a.root, "data", "raw"); PROC = os.path.join(a.root, "data", "processed")
line = "=" * 70

# ---------------------------------------------------------------- A. autoencoder vs simple baselines
print(line, "\nA. AUTOENCODER vs BASELINES (held-out cell lines)\n" + line)
try:
    import torch
    from sklearn.model_selection import train_test_split
    X = torch.cat([torch.as_tensor(torch.load(os.path.join(RAW, f), map_location="cpu")).float()
                   for f in ("GeneExpressionFeature402.pt", "CNVFeature402.pt", "miRNAFeature402.pt")], 1).numpy()
    tr, ho = train_test_split(np.arange(X.shape[0]), test_size=40, random_state=a.seed, shuffle=True)  # same call as your notebook
    Xtr, Xho = X[tr], X[ho]
    mu = Xtr.mean(0)
    mean_mse = float(np.mean((Xho - mu) ** 2))
    k = min(512, len(tr) - 1)
    _, _, Vt = np.linalg.svd(Xtr - mu, full_matrices=False)
    V = Vt[:k].T
    pca_mse = float(np.mean((((Xho - mu) @ V) @ V.T + mu - Xho) ** 2))
    print(f"held-out MSE  predict-train-mean : {mean_mse:.6f}")
    print(f"held-out MSE  PCA-{k:<3d}            : {pca_mse:.6f}")
    print(f"held-out MSE  autoencoder (yours): {a.ae_holdout_mse:.6f}")
    if a.ae_holdout_mse >= mean_mse * 0.95:
        print("-> WARNING: autoencoder is barely better than the training mean. Embeddings may carry little information.")
    elif a.ae_holdout_mse > pca_mse:
        print("-> Autoencoder beats the mean but not PCA. PCA is a legitimate (simpler) cell representation.")
    else:
        print("-> Autoencoder beats both baselines on held-out cell lines.")
except Exception as e:
    print("Skipped part A:", repr(e))

# ---------------------------------------------------------------- B. graph diagnostics
print("\n" + line, "\nB. GRAPH DIAGNOSTICS\n" + line)
sp = os.path.join(PROC, "splits")
train_df = pd.read_csv(os.path.join(sp, "cold_drug_train.csv"))
test_df = pd.read_csv(os.path.join(sp, "cold_drug_test.csv"))
train_drugs, test_drugs = set(train_df["drug_index"]), set(test_df["drug_index"])
print("split check: train drugs", len(train_drugs), "| test drugs", len(test_drugs), "| overlap", len(train_drugs & test_drugs))

ca = np.load(os.path.join(PROC, "pancreatic28_cell_adjacency_methylation_threshold075.npy"))
n = ca.shape[0]
print(f"\ncell-cell graph: {n} nodes | density {ca.sum() / (n * (n - 1)):.3f} | min/median/max degree "
      f"{int(ca.sum(1).min())}/{int(np.median(ca.sum(1)))}/{int(ca.sum(1).max())}")
if ca.sum() / (n * (n - 1)) > 0.9:
    print("-> Cell-cell graph is almost complete: it carries almost no structure at threshold 0.75.")

da = np.load(os.path.join(PROC, "drug_adjacency_fingerprint_threshold075.npy"))
deg = da.sum(1)
print(f"\ndrug-drug graph: {da.shape[0]} nodes | undirected edges {int(np.triu(da, 1).sum())} | "
      f"isolated nodes {int((deg == 0).sum())} ({100 * (deg == 0).mean():.1f}%) | median degree {np.median(deg):.0f}")
tr_idx = np.array(sorted(train_drugs)); te_idx = np.array(sorted(test_drugs))
te_with_train_nb = int((da[np.ix_(te_idx, tr_idx)].sum(1) > 0).sum())
print(f"test drugs with >=1 TRAINING-drug neighbour: {te_with_train_nb} of {len(te_idx)} "
      f"({100 * te_with_train_nb / len(te_idx):.1f}%)")
if te_with_train_nb / len(te_idx) < 0.5:
    print("-> Most test drugs have no similarity link into the training graph. At test time they will be isolated;"
          "\n   consider Tanimoto k-nearest-neighbour edges instead of Pearson >= 0.75 on Morgan bits.")

edges = pd.read_csv(os.path.join(PROC, "train_response_cell_drug_edges_lower50.csv"))
eset = set(zip(edges["omics_index"], edges["drug_index"]))
train_df["has_edge"] = [(c, d) in eset for c, d in zip(train_df["omics_index"], train_df["drug_index"])]
rho = spearmanr(train_df["has_edge"].astype(int), train_df["log_ic50"])[0]
print(f"\nresponse edges: {len(edges)} | share of TRAIN pairs that have an edge: {train_df.has_edge.mean():.3f}")
print(f"Spearman(has_edge, log_ic50) on training pairs: {rho:.3f}")
print("-> This is how much a training pair's OWN edge reveals about its label. Test drugs have no such edges,\n"
      "   so a model that reads this edge sees different inputs at train and test time (self-leakage).")

# ---------------------------------------------------------------- C. AUC at pair level
print("\n" + line, "\nC. AUC TARGET ON THE EXISTING PAIRS\n" + line)
tp = os.path.join(PROC, "targets_raw_v1.csv")
if os.path.exists(tp):
    t = pd.read_csv(tp, usecols=["ccle_name", "name", "auc"])
    auc_pair = t.groupby(["ccle_name", "name"], as_index=False)["auc"].median().rename(columns={"auc": "auc_pair"})
    integ = pd.read_csv(os.path.join(PROC, "integrated_modeling_table.csv"))
    m = integ.merge(auc_pair, on=["ccle_name", "name"], how="left")
    print("integrated pairs:", len(integ), "| with AUC:", int(m.auc_pair.notna().sum()))
    print("Spearman(log_ic50, AUC) at pair level:", round(spearmanr(m["log_ic50"], m["auc_pair"], nan_policy="omit")[0], 3))
    m.to_csv(os.path.join(PROC, "integrated_with_auc.csv"), index=False)
    print("saved integrated_with_auc.csv")
    for nm, d in (("train", train_df), ("test", test_df)):
        mm = d.merge(auc_pair, on=["ccle_name", "name"], how="left") if {"ccle_name", "name"} <= set(d.columns) else None
        if mm is not None:
            print(f"AUC {nm}: median {mm.auc_pair.median():.3f} | share >= 1: {(mm.auc_pair >= 1).mean():.3f}")
else:
    print("targets_raw_v1.csv not found; skipping.")
print("\nDone.")