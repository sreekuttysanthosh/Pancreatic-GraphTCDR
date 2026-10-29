# Do Drug-Response Models Generalize to Unseen Drugs?
### A long-tail and similarity audit in pancreatic cancer cell lines (PRISM)

> **Research prototype. Preclinical cell-line data only. Not for clinical use.**

A GraphTCDR-inspired study (Zhang et al., *Neural Networks*, 2026). It is an **evaluation study**, not a reproduction of the paper and not a new architecture.

## Question
How reliably do drug-response models predict responses to drugs they have never seen, and how much does performance fall when chemically similar drugs are excluded from training?

## Key findings
Same model (RF + autoencoder embeddings + Morgan fingerprints), same data, three ways of testing (3 split seeds each):

| Evaluation | R² | Spearman |
|---|---|---|
| Random pair split (97–98% of test drugs also in training) | 0.584 ± 0.011 | 0.733 |
| Random drug-grouped split | 0.153 ± 0.072 | 0.380 |
| Bemis–Murcko scaffold split | 0.066 ± 0.009 | 0.332 |
| Cell-line-mean baseline (grouped/scaffold) | about 0.01 | |

1. **The evaluation design changes the apparent performance.** Apparent performance falls sharply once test drugs are unseen (0.584 to 0.153), and again when scaffold analogues are removed (to 0.066). Notebooks 13 and 14.
2. **Performance depends on similarity to the training drugs.** For test drugs with maximum Tanimoto similarity below about 0.4 to any training drug, R² is near zero in both grouped splits. At matched similarity the two splits give similar R² in the lower bins. Exploratory, 3 seeds (14.12).
3. **Most of what is predicted is drug potency.** Between 62% and 79% of label variance lies between drugs (14.13–14.14). A fingerprint-only RF reaches R² 0.116 (random drug) and 0.037 (scaffold); adding a cell-line offset gives 0.128 and 0.050, against 0.153 and 0.066 for RF + AE (14.17).
4. **Drug-specific cell-line signal is small.** Within-drug Spearman is 0.28 (random drug) and 0.25 (scaffold) for RF, against 0.21 for a cell-line-mean baseline. The mean gain is 0.07 and 0.04; cluster-bootstrap intervals (resampling drugs) exclude zero overall, but not for scaffold-split drugs with similarity below 0.3 (14.15).
5. **Long tail.** Extreme IC50 values (|log10 IC50| > 3) are rare but have poor curve fits (Notebook 06). They were kept and clipped, not removed. Quality and tail weighting gave small gains (R² +0.01 to +0.02) and shifted error between response bins rather than removing it (Notebook 11).
6. **Same pattern on the frozen cold-drug split** (Notebook 11): median per-drug Spearman 0.33 (RF) vs 0.24 (cell-line mean).

## Data
- PRISM secondary repurposing screen; 28 pancreatic cell lines matched to GraphTCDR multi-omics data.
- 1,228 drugs, 18,037 drug–cell-line pairs (median of repeated measurements); only pairs with a finite fitted IC50.
- Cell features: 512-d autoencoder embedding of expression, CNV and miRNA (PCA 40-d as alternative). Methylation is used only for the cell-similarity graph.
- Drug features: RDKit Morgan fingerprints, 2,048 bits.
- Target `y_clip`: log10(IC50) clipped to the training 1st/99th percentiles [-3.2443, 1.2946] (`data/processed/target_v2/target_policy.json`). `y_raw` (unclipped) is reported separately.

## Evaluation designs
| Design | What is held out |
|---|---|
| Frozen cold-drug split (Notebooks 07, 11) | 246 drugs (3,802 pairs); validation grouped by drug (99 drugs) |
| Random pair split (Notebook 13) | Random pairs; drugs are shared between train and test |
| Random drug-grouped split, 3 seeds (Notebooks 13–14) | About 21% of drugs |
| Bemis–Murcko scaffold split, 3 seeds (Notebook 14) | Whole scaffold groups (983 scaffolds, 42 acyclic drugs as own groups) |

## Frozen-split model comparison (`y_clip`, test, 3-seed mean)
| Model | R² | MAE | Spearman |
|---|---|---|---|
| RF + AE (LDS + quality weighted) | 0.204 | 0.529 | 0.413 |
| RF + AE (tail + quality weighted) | 0.197 | 0.525 | 0.418 |
| RF + AE (quality weighted) | 0.192 | 0.521 | 0.429 |
| RF + AE | 0.182 | 0.521 | 0.420 |
| RF + PCA | 0.168 | 0.522 | 0.418 |
| CNN + AE | 0.125 | 0.553 | 0.343 |
| MLP + AE | 0.114 | 0.557 | 0.329 |
| XGBoost + AE | 0.104 | 0.563 | 0.290 |
| Graph A + residual + kNN drug graph | 0.058 | 0.574 | 0.259 |
| Floor: cell-line mean | 0.005 | 0.598 | 0.123 |

RF + AE bootstrap (246 test drugs, 1,000 resamples): R² 0.18, 95% CI [0.09, 0.28]. RF beats the graph model, MLP and XGBoost on paired differences; RF vs CNN and weighted vs unweighted RF are not distinguishable at 95%. In this small cold-drug setting the GraphTCDR-inspired model did not outperform RF; this is not a replication of the paper's results.

Full tables: `results/v2/FINAL_headline_table_yclip.csv`, `FINAL_all_scores.csv`, `bootstrap_ci.csv`. New analyses: `scaffold_split_comparison.csv`, `similarity_binned_results.csv`, `gain_cluster_bootstrap.csv`, `fingerprint_baselines.csv`.

## Integrity checks (Notebook 14, Sections 14.10–14.11)
- No missing targets and no drug–cell-line pair overlap across partitions.
- Drug and scaffold disjointness asserted for every split seed.
- The autoencoder was trained on unlabeled omics only (reconstruction loss; 402 cell lines, 362 train / 40 holdout; checkpoint chosen on holdout reconstruction error). 23 of the 28 pancreatic lines were in its training set. No response labels were used.

## Differences from the GraphTCDR paper
| Aspect | Paper | This project |
|---|---|---|
| Cell lines | 402 (pan-cancer) | 28 pancreatic (pan-cancer training as an extra experiment) |
| Drug fingerprint | PubChem 881-bit | RDKit Morgan 2,048-bit |
| Outliers | Removed (IQR) | Kept; training target clipped |
| Split | Random 8:1:1, 5-fold CV | Frozen cold-drug, random-drug and scaffold splits |
| Target | Normalised IC50 | log10(IC50), clipped (primary) and raw |

Scores are not comparable with the paper's. The paper lists the long-tailed PRISM response distribution as an open issue for future work; this project examines that tail and unseen-drug generalization.

## Limitations
- Results apply to unseen drugs on known cell lines, not unseen cell lines or patients.
- Only 28 cell lines; scaffold and similarity analyses cover RF + AE only (CNN and graph models were not rerun), so no model ranking under these splits is claimed.
- Three split seeds; test sets differ in size and composition; similarity bins were chosen by inspection and are exploratory.
- Similarity may be confounded with drug class or mechanism; no mechanism-of-action check was done.
- Pairs without a finite IC50 (21,206 of 40,119 records) were excluded; they had higher AUC, so results do not apply to non-responding pairs.
- Extreme values are treated as likely fit artefacts on the basis of association with poor curve fits, which is not proof.
- Clipping bounds were fixed from the original training split.
- Several paper baselines (DeepCDR, DeepTTA, GraphDRP, BANDRP, DRN-CDR) were not re-implemented. No external validation. Preclinical only.

## Repository structure
- `notebooks/`: 01–04 data preparation; 05/05b baselines and embeddings; 06 tail diagnosis; 07 target policy; 08 autoencoder and graph preparation; 09 old graph notebook (superseded); 10 GraphTCDR v2; 11 final evaluation; 12 pan-cancer; 13 split comparison; 14 scaffold, similarity and audit analyses.
- `data/processed/target_v2/`: splits and target policy.
- `results/v2/`: predictions, tables, figures.
- Raw PRISM data and model checkpoints are not included; consult the data providers' terms.

## Reproduce
Install `requirements.txt`, run notebooks in order. Use Notebook 11 for frozen-split numbers, 13 for the random-pair comparison and 14 for scaffold, similarity and baseline analyses.