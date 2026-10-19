\# Quality-Aware Assessment of Long-Tailed Drug-Response Prediction for Unseen Drugs in Pancreatic Cancer Cell Lines



A GraphTCDR-inspired study using PRISM. Preclinical cell-line prediction only; no clinical claims.



\## Question

For drugs never seen in training, how well can drug fingerprints and multi-omics predict response

in pancreatic cancer cell lines? Do a heterogeneous graph model and quality-aware / tail-aware

training improve on simpler baselines?



\## Data

\- PRISM secondary repurposing screen, 28 pancreatic cell lines matched to GraphTCDR omics

\- 1,228 drugs, 18,037 cell-drug pairs (median of repeated measurements)

\- Cell features: autoencoder (512-d, from expression, CNV, miRNA) or PCA (40-d)

\- Drug features: RDKit Morgan fingerprint, 2,048 bits

\- Methylation is used only to build the cell-similarity graph



\## Evaluation design

| Split | Pairs | Drugs |

|---|---|---|

| Train | 12,929 | 883 |

| Validation | 1,306 | 99 |

| Test (frozen) | 3,802 | 246 |



Cold-drug split: zero drug overlap (asserted in code). Validation is drug-grouped. The test set

was never used for training, weighting, early stopping or tuning.



\## Target policy

\- Raw log10(IC50) has implausible extremes (about -110 to +288). 59 test pairs (1.5%) have |y| > 3 and mostly poor curve fits.

\- \*\*Primary target `y\_clip`\*\*: clipped to training 1st/99th percentiles \[-3.2443, 1.2946]

\- \*\*Secondary target `y\_raw`\*\*: unclipped, reported separately

\- Bounds computed from training data only (`data/processed/target\_v2/target\_policy.json`)



\## Results (y\_clip, test, mean over 3 seeds)

| Model | R² | MAE | Spearman |

|---|---|---|---|

| RF + AE (tail+quality weighted) | 0.197 | 0.525 | 0.418 |

| RF + AE (quality weighted) | 0.192 | 0.521 | 0.429 |

| RF + AE | 0.182 | 0.521 | 0.420 |

| RF + PCA | 0.168 | 0.522 | 0.418 |

| CNN + AE (quality weighted) | 0.145 | 0.556 | 0.340 |

| CNN + AE | 0.125 | 0.553 | 0.343 |

| MLP + AE | 0.114 | 0.557 | 0.329 |

| XGBoost + AE | 0.104 | 0.563 | 0.290 |

| Graph A + residual + kNN drug graph | 0.058 | 0.574 | 0.259 |

| Graph A + residual | 0.051 | 0.580 | 0.234 |

| Graph B (response edges, half data) | 0.035 (SD 0.059) | 0.597 | 0.192 |

| Floor: cell-line mean | 0.005 | 0.598 | 0.123 |

| Floor: global mean | -0.006 | 0.604 | n/a |



Full tables: `results/v2/FINAL\_headline\_table\_yclip.csv`, `FINAL\_headline\_table\_yraw.csv`, `FINAL\_all\_scores.csv`.



!\[Leaderboard](results/v2/fig\_leaderboard.png)



\## Findings

1\. Every model beats the constant-prediction floor on unseen drugs. Random Forest is best and most stable.

2\. Reweighting (quality, tail, or both) gives a small, consistent gain in RF (+0.01 to +0.015 R²) and in the CNN. Quality-only and tail-only effects cannot be separated.

3\. The GraphTCDR-inspired model did not beat simpler baselines here. Adding the paper's residual connection, 16 heads, and a Tanimoto kNN drug graph moved it only slightly. Response edges did not help on average.

4\. Long tail: all models shrink predictions toward the centre (overpredict the sensitive tail, underpredict the resistant tail). Weighting shifts error between bins rather than removing it. Models still enrich the most sensitive 10% of responses about 3-fold over chance (RF precision about 0.30 vs 0.10 random; graph about 0.20).

5\. On `y\_raw`, all models score R² near 0 (best 0.026) because 59 extreme labels dominate the error.



!\[Bias by bin](results/v2/fig\_bias\_by\_bin.png)



\## Differences from the GraphTCDR paper

| Aspect | Paper | This project |

|---|---|---|

| Cell lines | 402 | 28 pancreatic |

| Drug fingerprint | PubChem 881-bit | RDKit Morgan 2,048-bit |

| Outliers | Removed (IQR) before splitting | Kept in test; training target clipped |

| Split | Random 8:1:1, 5-fold CV | Fixed cold-drug split |

| Target | Normalised IC50 | log10(IC50), clipped (primary) and raw |



Scores are not comparable with the paper's. This is a GraphTCDR-inspired adaptation, not a reproduction.



\## Limitations

\- Only 28 cell lines; seeds capture model randomness only (fixed drug split, no split-uncertainty estimate)

\- Extreme values treated as likely fit artefacts based on association with poor curve fit, not proven

\- Pairs without a finite IC50 (21,206 records) were excluded, likely biasing toward drugs that reached 50% inhibition

\- Autoencoder trained on 402 lines that may include the 28 (unsupervised, no response labels)

\- Several paper baselines (DeepCDR, DeepTTA, GraphDRP, BANDRP, DRN-CDR) were not re-implemented

\- No external validation; preclinical only



\## Repository

\- `notebooks/`: 01-04 data prep; 05 baselines (v2 section); 06 tail diagnosis; 07 target policy (v2 section); 08 autoencoder and graph prep; 09 old graph notebook (superseded); 10 graph v2; 11 final evaluation

\- `data/processed/target\_v2/`: splits and target policy

\- `results/v2/`: predictions, tables, figures

\- `requirements.txt`: package versions



\## Reproduce

Run notebooks in order; use `11\_Final\_Evaluation.ipynb` for all reported numbers.

