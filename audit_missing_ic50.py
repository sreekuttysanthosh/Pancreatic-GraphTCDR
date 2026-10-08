import os
import numpy as np
import pandas as pd

ROOT = r"D:\Pancreatic_GraphTCDR"

RAW = os.path.join(
    ROOT,
    "data",
    "raw",
    "secondary-screen-dose-response-curve-parameters.csv"
)

MAP = os.path.join(
    ROOT,
    "data",
    "processed",
    "pancreatic28_cell_mapping.csv"
)

print("Loading PRISM response data...")

df = pd.read_csv(
    RAW,
    usecols=["ccle_name", "ic50", "auc", "lower_limit"]
)

print(f"Loaded records: {len(df):,}")

for column in ["ic50", "auc", "lower_limit"]:
    df[column] = pd.to_numeric(df[column], errors="coerce")

df["ic50_present"] = (
    np.isfinite(df["ic50"]) &
    (df["ic50"] > 0)
)

mapping = pd.read_csv(MAP)

pancreatic_names = set(
    mapping["cellLineName"].astype(str)
)

pancreatic_df = df[
    df["ccle_name"].isin(pancreatic_names)
].copy()


def report(title, data):

    print("\n" + "=" * 70)
    print(title)
    print(f"Records: {len(data):,}")
    print("=" * 70)

    rows = []

    groups = [
        ("IC50 present", data["ic50_present"]),
        ("IC50 missing/non-finite", ~data["ic50_present"])
    ]

    for label, mask in groups:

        group = data[mask]

        rows.append({
            "group": label,
            "N": len(group),
            "pct_of_records": round(
                100 * len(group) / len(data), 2
            ),
            "median_AUC": round(
                group["auc"].median(), 3
            ),
            "mean_AUC": round(
                group["auc"].mean(), 3
            ),
            "pct_AUC>=1": round(
                100 * (group["auc"] >= 1).mean(), 2
            ),
            "pct_lower_limit>0.5": round(
                100 * (group["lower_limit"] > 0.5).mean(), 2
            ),
            "median_lower_limit": round(
                group["lower_limit"].median(), 3
            ),
            "pct_AUC_missing": round(
                100 * group["auc"].isna().mean(), 2
            )
        })

    result = pd.DataFrame(rows)

    print(result.to_string(index=False))


report(
    "ALL PRISM CELL LINES",
    df
)

report(
    "28 PANCREATIC CELL LINES",
    pancreatic_df
)

print("\n" + "=" * 70)
print("SANITY CHECKS")
print("=" * 70)

print(f"Pancreatic records: {len(pancreatic_df):,}")

print(
    f"Pancreatic cell lines found: "
    f"{pancreatic_df['ccle_name'].nunique()}"
)

print(
    f"IC50 present: "
    f"{pancreatic_df['ic50_present'].sum():,}"
)

print(
    f"IC50 missing/non-finite: "
    f"{(~pancreatic_df['ic50_present']).sum():,}"
)

print(
    f"AUC missing/non-finite: "
    f"{(~np.isfinite(pancreatic_df['auc'])).sum():,}"
)

print("\nDONE")
