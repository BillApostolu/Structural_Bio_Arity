# Output file reference

Every file below is written under `Data/` (created by the scripts, ignored by git). Placeholders in file names: `{r}` radius in Å, `{k}` sequence gap, `{W}` smoothing window, `{weight}` = `weighted` or `unweighted`.

Contents:

1. [Pair list and chain mapping](#1-pair-list-and-chain-mapping) (`data_acquisition.py`, `phase1.py`)
2. [Phase 2 signatures](#2-phase-2-signatures) (`phase2.py`)
3. [Phase 3 validation](#3-phase-3-validation) (`phase3.py`)
4. [Phase 4 summary statistics](#4-phase-4-summary-statistics) (`phase4.py`)
5. [Phase 4 confusion matrices](#5-phase-4-confusion-matrices) (`phase4.py`)

---

## 1. Pair list and chain mapping

### `target_proteins.txt` (tab-separated, written by `data_acquisition.py`)

| Column | Description |
|---|---|
| `UniProt_ID` | UniProt accession of the human protein; also the AF2 model ID (`AF-{UniProt_ID}`) |
| `PDB_ID` | Experimental structure used for this protein (first PDB cross-reference in the UniProt entry) |

### `mapping_log.csv` (written by `phase1.py`)

One row per protein that passed Phase 1.

| Column | Description |
|---|---|
| `UniProt_ID` | UniProt accession |
| `PDB_ID` | PDB entry (lower case) |
| `Matched_Chain` | Chain of the PDB entry that aligned best to the AF2 sequence and was used for all comparisons |

---

## 2. Phase 2 signatures

### `phase2_results/Signatures_r_{r}_k_{k}_{weight}_W{W}.csv`

One row per protein and model (AF2 and PDB), computed on the profile smoothed with window W (W = 1 is the raw, unsmoothed profile). NaN positions (residues unresolved in the PDB structure) are ignored.

| Column | Description |
|---|---|
| `Protein_ID` | UniProt accession |
| `Model_Type` | `AF2` or `PDB` |
| `Window_Size` | W |
| `p_10`, `p_25`, `p_50`, `p_75`, `p_90` | Arity signature: 10th–90th percentiles of the arity profile (`numpy.nanpercentile`, linear interpolation) |
| `Entropy_Hp` | Persistence entropy: Shannon entropy of the normalised persistences divided by log(number of merges); 0–1 |
| `Mean_Persistence` | Mean persistence of the merged components, in arity units (how far the arity rises between a component's birth and its merge) |
| `Fraction_CP` | Number of merges with positive persistence divided by the number of residues |
| `Max_CC` | Maximum number of separate connected components reached during the filtration |
| `PLM` | Persistent Local Maxima: number of persistences greater than T<sub>ν</sub> = 2.0 |

---

## 3. Phase 3 validation

### `phase3_results/Phase3_Validation_r_{r}_k_{k}_unweighted_W{W}.csv`

One row per protein. Phase 3 uses unweighted profiles only, restricted to the residues resolved in both AF2 and the PDB.

| Column | Description |
|---|---|
| `Protein_ID` | UniProt accession |
| `Window_Size` | W |
| `EMD_Score` | 1D Earth Mover's (Wasserstein-1) distance between the AF2 and PDB arity value distributions |
| `DTW_Error` | pLDDT-weighted Dynamic Time Warping cost between the two profiles, divided by the number of residues compared |
| `Mean_pLDDT` | Mean AF2 pLDDT over the residues compared |
| `Cross_Domain_PAE` | Mean AF2 PAE (Å) over all residue pairs more than 50 positions apart in sequence; NaN if the PAE file could not be downloaded |
| `Hallucination_Flag` | `True` when `DTW_Error` ≤ 1.0, `Cross_Domain_PAE` > 15 Å and `Mean_pLDDT` ≥ 70: good local agreement with the experiment, but AF2 reports high uncertainty for the placement of distant regions |
| `AF2_Fold_Angle`, `PDB_Fold_Angle` | Hinge angle θ in degrees, computed for flagged proteins only; NaN otherwise, or if no hinge with prominence ≥ 2.0 is found, or if fewer than 50 residues are compared |
| `Angle_Deviation` | \|`AF2_Fold_Angle` − `PDB_Fold_Angle`\| in degrees |

---

## 4. Phase 4 summary statistics

### `phase4_results/Global_Grid_Search_Stats.csv`, `Sweep_W.csv`, `Sweep_k.csv`, `Sweep_r.csv`

All four files have the same columns, one row per configuration, based on unweighted profiles.

| File | Rows |
|---|---|
| `Global_Grid_Search_Stats.csv` | r ∈ {10, 15} × k ∈ {0, 1, 2, 5} × W ∈ {1, 5, 11, 15, 25, 31} |
| `Sweep_W.csv` | r = 10, k = 0, W ∈ {1, 5, 11, 15, 25, 31} |
| `Sweep_k.csv` | r = 10, W = 1, k ∈ {0, 1, 2, 5} |
| `Sweep_r.csv` | k = 0, W = 1, r ∈ {10, 15} |

### 4.1 Configuration

| Column | Description |
|---|---|
| `r` | Neighbourhood radius in Å |
| `k` | Sequence gap: neighbours with \|i − j\| ≤ k are ignored. `k = 0` excludes only the residue itself (standard arity of Cazals & Sarti); `k > 0` removes local contacts such as those inside helices |
| `W` | Smoothing window; `W = 1` is the raw profile, larger W averages over longer stretches (≈ secondary structure at 5, domains at 25–31) |
| `Total_Proteins` | Proteins present in both the Phase 2 and Phase 3 outputs for this configuration |

### 4.2 Signature distance

| Column | Description |
|---|---|
| `D_sig_2D_Mean` | Mean over proteins of the Euclidean distance between the AF2 and PDB signatures in (p<sub>25</sub>, p<sub>75</sub>). Lower means closer packing distributions |
| `D_sig_3D_Mean` | Same in (p<sub>25</sub>, p<sub>50</sub>, p<sub>75</sub>) |

Absolute values depend on r and k, because the arity scale does (see *Known limitations* in the README).

### 4.3 Error by reference packing mode

The four reference modes are defined once, by Ward clustering of the PDB signatures (p<sub>25</sub>, p<sub>75</sub>) at r = 10, k = 0, W = 1, and ordered by packing: Mode 0 (lowest, labelled IDPs in the report), Mode 1 (molten globules), Mode 2 (dynamic globulars), Mode 3 (highest, rigid globulars). The same assignment is used for every configuration. For each `c` = 0, 1, 2, 3:

| Column | Description |
|---|---|
| `DTW_Mode{c}_Med` | Median `DTW_Error` of the proteins in mode c |
| `EMD_Mode{c}_Med` | Median `EMD_Score` of the proteins in mode c |
| `Mode{c}_Size` | Number of proteins of this configuration in mode c |

These are medians only; no statistical test between modes is computed.

### 4.4 Paired AF2 vs PDB comparison

For each metric `m` in `p_10`, `p_25`, `p_50`, `p_75`, `p_90`, `Entropy_Hp`, `Mean_Persistence`, `Fraction_CP`, `Max_CC`, `PLM` (defined in section 2), using the proteins that have both an AF2 and a PDB value:

| Column | Description |
|---|---|
| `{m}_AF2_Med`, `{m}_PDB_Med` | Median of the metric for AF2 and for the PDB |
| `{m}_pval` | Two-sided Wilcoxon signed-rank p-value for AF2 vs PDB (`zero_method="zsplit"`), not corrected for multiple testing |
| `{m}_Cliffs` | Cliff's delta = P(AF2 > PDB) − P(AF2 < PDB), from −1 to 1. Negative means AF2 values tend to be lower than the PDB values. Common rule of thumb for \|δ\|: < 0.147 negligible, < 0.33 small, < 0.474 medium, otherwise large |

All four columns are NaN when 10 or fewer proteins have both values.

### 4.5 Flagged inter-domain orientation

| Column | Description |
|---|---|
| `Total_Hallucinations` | Number of proteins with `Hallucination_Flag = True` in Phase 3 |
| `Mean_Angle_Dev` | Mean `Angle_Deviation` (degrees) of the flagged proteins for which a hinge was found in both AF2 and the PDB; NaN if none |
| `Max_Angle_Dev` | Maximum `Angle_Deviation` (degrees) among the same proteins |

---

## 5. Phase 4 confusion matrices

### `phase4_results/Confusion_Matrix_r_10_k_{k}_W_{W}.csv`

A 4 × 4 table of protein counts for one configuration (k ∈ {0, 1, 2, 5}, all six W). Unlike the fixed reference modes in section 4.3, the PDB and AF2 signatures (p<sub>25</sub>, p<sub>75</sub>) are clustered independently for this configuration (Ward, 4 clusters, each ordered by packing).

| Part | Description |
|---|---|
| Rows `True PDB 0 (IDP)` … `True PDB 3 (Rigid)` | Mode of the protein from clustering the PDB signatures |
| Columns `Pred AF2 0` … `Pred AF2 3` | Mode of the protein from clustering the AF2 signatures |

The diagonal counts proteins placed in the same mode by both; the agreement printed by `phase4.py` is the diagonal sum divided by the total.

The k-ablation results (`phase4.py`, "Ablation sweep") are printed to the console and not saved to a file.