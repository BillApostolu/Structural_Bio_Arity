# Structural_Bio_Arity

**Multiscale, arity-based comparison of AlphaFold 2 predictions with experimental PDB structures**

This repository contains the analysis pipeline for *Project 10: AlphaFold understanding it through usage*. It builds on the Cα arity packing signatures of Cazals & Sarti (2024) and extends them with:

- a **sequence-gap parameter k**, which ignores contacts between residues that are close in sequence, so that long-range (tertiary) packing can be separated from local secondary structure;
- **sliding-window smoothing** (window W) of the per-residue arity profile at several structural scales;
- **1D persistence (union–find filtration)** metrics computed on the arity profile;
- **per-protein profile errors** between AlphaFold 2 (AF2) and the PDB: Earth Mover's Distance and pLDDT-weighted Dynamic Time Warping;
- a **flag for uncertain inter-domain orientation** based on AF2's Predicted Aligned Error (PAE), checked against the hinge angle in the experimental structure.

The analysis in the report covers 8,050 human proteins that have both an AF2 model and an experimental structure. The full write-up, figures and discussion are in [`Project10_Report.pdf`](Project10_Report.pdf).

---

## Repository layout

```
.
├── scripts/
│   ├── data_acquisition.py   # build the UniProt/PDB pair list, download AF2 models and PDB entries
│   ├── phase1.py             # parse, align, compute arity profiles, download PAE
│   ├── phase2.py             # smoothing, arity signatures, persistence metrics
│   ├── map_generation.py     # 2D/3D global arity maps
│   ├── phase3.py             # EMD, DTW, inter-domain flag, hinge angles
│   └── phase4.py             # clustering, statistics, parameter sweeps, ablation
├── csv_columns.md            # every column of every output CSV
├── Project10_Report.pdf      # project report
├── requirements.txt
└── LICENSE
```

All data and results are written to `Data/`, which is created by the scripts and ignored by git.

## Installation

Python 3.10 or newer.

```bash
git clone https://github.com/BillApostolu/Structural_Bio_Arity.git
cd Structural_Bio_Arity
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Running the pipeline

The scripts use paths relative to the current directory, so run them **from the repository root**, in this order:

| Step | Command | What it does | Main outputs (under `Data/`) |
|---|---|---|---|
| 1a | `python scripts/data_acquisition.py` | Builds the UniProt/PDB pair list and downloads AF2 models and PDB entries | `target_proteins.txt`, `raw/AlphaFold2/`, `raw/PDB/` |
| 1b | `python scripts/phase1.py` | Parsing, sequence alignment, arity profiles, PAE download | `processed_1d_signals_new/`, `mapping_log.csv` |
| 2 | `python scripts/phase2.py` | Smoothing, arity signatures, persistence metrics | `phase2_results/` |
| 2b | `python scripts/map_generation.py` | 2D and 3D global arity maps | `phase2_results/maps/` |
| 3 | `python scripts/phase3.py` | EMD, DTW, inter domain flag, hinge angles | `phase3_results/` |
| 4 | `python scripts/phase4.py` | Reference modes, statistics, sweeps, ablation, confusion matrices | `phase4_results/` and a console summary |

Notes:

- `data_acquisition.py` and `phase1.py` need internet access (UniProt, AlphaFold DB, RCSB PDB). `data_acquisition.py` can be stopped and restarted: pairs already listed in `Data/target_proteins.txt` are skipped. The number of pairs to collect is set by `TARGET_PAIRS_NEEDED` (default 20,000); fewer pass the alignment filter in Phase 1.
- **Disk usage.** Phase 1 writes 192 small `.npy` files per protein (6 radii × 8 gaps × weighted/unweighted × AF2/PDB) and Phase 2 a further 1,152 (× 6 windows). For about 8,000 proteins this is more than 10 million files. Phase 4 only uses r ∈ {10, 15} and k ∈ {0, 1, 2, 5}; to save space and time you can trim `R_VALUES` and `K_VALUES` in `phase1.py` to those values.
- Re-running `phase1.py` rebuilds `Data/mapping_log.csv` from scratch.
- Every column of every output CSV is described in [`csv_columns.md`](csv_columns.md).

---

## Method

### Phase 1 — Data acquisition, parsing and synchronisation

**Target selection (`data_acquisition.py`).** UniProt is queried for reviewed human entries (Swiss-Prot, taxonomy 9606) with at least one PDB cross-reference. For each entry the first PDB ID in the cross-reference list is used. The AF2 model is downloaded through the AlphaFold DB API and the PDB entry in legacy `.pdb` format from RCSB. A pair is kept only if both downloads succeed, so entries that are distributed only as mmCIF (very large structures) drop out.

**Parsing (`phase1.py`).** Cα coordinates are read from the first model of each file, keeping standard residues and the modified residues MSE, CSO, SEC and PYL. For AF2 the B-factor column (pLDDT) is kept as well. The AF2 model has a single chain; for the PDB entry every chain is parsed.

**Sequence synchronisation.** The AF2 sequence is aligned to each PDB chain with Smith–Waterman local alignment (BLOSUM62, gap open −10, gap extend −1) and the best-scoring chain is kept. The pair is rejected if the aligned residues cover less than 40% of that PDB chain. PDB Cα coordinates are then placed into AF2 residue numbering; AF2 positions with no aligned, resolved PDB residue are set to NaN.

**Arity.** For residue i with Cα–Cα distances d<sub>ij</sub>, the unweighted arity at radius r and sequence gap k is the number of neighbours

```math
a_i(r,k) = \left|\left\{\, j \;:\; d_{ij} \le r,\ |i-j| > k \,\right\}\right|
```

and the weighted arity applies a Gaussian distance-decay kernel to the same neighbours:

```math
a^{w}_i(r,k) = \sum_{j \,:\, d_{ij} \le r,\ |i-j| > k} \exp\!\left(-\left(\frac{d_{ij}}{r}\right)^{2}\right)
```

k = 0 excludes only the residue itself and gives the standard arity of Cazals & Sarti; k > 0 also discards neighbours within k positions in sequence. Arities are computed for r ∈ {5, 8, 10, 12, 15, 20} Å and k ∈ {0, 1, 2, 3, 4, 5, 7, 10}, for the AF2 model and for the synchronised PDB chain. In the PDB chain only resolved residues count as neighbours, and unresolved positions stay NaN. Each profile is saved as a 1D array in AF2 numbering, one file per (r, k), weighting and model.

The AF2 PAE matrix is downloaded from the AlphaFold DB and stored together with the Cα coordinates and pLDDT for Phase 3. The matched PDB chain of every protein is recorded in `Data/mapping_log.csv`.

### Phase 2 — Smoothing, arity signatures and persistence (`phase2.py`, `map_generation.py`)

Every Phase 1 profile (all (r, k), weighted and unweighted, AF2 and PDB) is processed at window sizes W ∈ {1, 5, 11, 15, 25, 31}.

**Smoothing.** A centred moving average of width W; W = 1 is the raw profile. The mean is taken over the resolved residues inside the window, and positions that are unresolved in the PDB structure stay NaN. The window sizes roughly follow structural scales: W = 5 for secondary structure (helices, strands), W = 11–15 for super-secondary motifs and W = 25–31 for domains.

**Arity signature.** The 10th, 25th, 50th, 75th and 90th percentiles of the (smoothed) profile, computed with `numpy.nanpercentile` (linear interpolation, NaNs ignored). Cazals & Sarti define the signature as the inverse CDF, i.e. the smallest arity whose cumulative frequency reaches the percentile; `np.nanpercentile(..., method="inverted_cdf")` reproduces that definition exactly.

**Persistence (arity filtration).** Residues are added in increasing order of arity (a sublevel-set filtration of the 1D profile) and union–find tracks stretches of consecutive residues. When two stretches meet, the younger one (the one born at the higher arity) dies; its persistence is the arity at the merge minus its birth arity. NaN positions are removed first. The metrics are:

- **Persistence entropy H<sub>p</sub>:** Shannon entropy of the normalised persistences, divided by log m, where m is the number of merges with positive persistence (range 0–1);
- **Mean persistence;**
- **Fraction of positive critical points:** m divided by the number of residues;
- **Maximum number of connected components** reached during the filtration;
- **Persistent Local Maxima (PLM):** number of persistences greater than T<sub>ν</sub> = 2.0.

Results are written to one CSV per (r, k), weighting and W, and the smoothed profiles are saved for Phase 3.

**Global arity maps (`map_generation.py`).** For every signatures file, the 2D map plots p<sub>25</sub> against p<sub>75</sub> as density histograms for the PDB and AF2 and as a scatter overlay; the 3D map adds p<sub>50</sub> and colours points by local density. Densely packed proteins lie towards the upper right of the maps and sparsely packed or disordered proteins towards the lower left.

### Phase 3 — Profile alignment and inter-domain orientation (`phase3.py`)

Phase 3 uses the **unweighted** smoothed profiles for every (r, k) and W. For each protein the AF2 and PDB profiles are restricted to the positions resolved in both; unresolved stretches are removed from both profiles before comparison.

**EMD.** The 1D Earth Mover's (Wasserstein-1) distance between the AF2 and PDB arity value distributions (`scipy.stats.wasserstein_distance`):

```math
\mathrm{EMD} = \int_{0}^{\infty} \left| F_{\mathrm{PDB}}(x) - F_{\mathrm{AF2}}(x) \right| \, dx
```

**pLDDT-weighted DTW.** Dynamic Time Warping between the two profiles with local cost

```math
C(i,j) = \frac{\mathrm{pLDDT}_{\mathrm{AF2},\,j}}{100}\,\left(a_{\mathrm{PDB},\,i} - a_{\mathrm{AF2},\,j}\right)^{2}
```

so that deviations where AF2 is confident cost more and low-confidence regions are down-weighted. The accumulated cost along the optimal warping path is divided by the number of residues compared.

**Inter-domain orientation flag.** The long-range PAE is the mean AF2 PAE over all residue pairs more than 50 positions apart in sequence. A protein is flagged (column `Hallucination_Flag`) when its arity profile agrees well with the experiment and AF2 is confident locally, but AF2 reports high error for the relative placement of distant regions:

- DTW error ≤ 1.0, **and**
- mean long-range PAE > 15 Å, **and**
- mean pLDDT over the compared residues ≥ 70.

**Hinge angle.** For flagged proteins, a hinge is located separately in the AF2 and PDB profiles as the most prominent minimum of the smoothed profile (`scipy.signal.find_peaks` on the inverted profile, prominence ≥ 2.0, the same value as the PLM threshold in Phase 2). Proteins with fewer than 50 compared residues or without such a minimum get no angle. With hinge Cα H and vectors **v**<sub>A</sub>, **v**<sub>B</sub> from H to the centroids of the Cα atoms before and after it,

```math
\theta = \arccos\left(\frac{\vec{v}_A \cdot \vec{v}_B}{\lVert\vec{v}_A\rVert\,\lVert\vec{v}_B\rVert}\right)
```

The output contains θ for AF2 and for the PDB and the absolute deviation between them.

### Phase 4 — Statistics, sweeps and ablation (`phase4.py`)

Phase 4 merges the Phase 2 signatures and Phase 3 results (unweighted) for each (r, k, W) and computes per protein the signature distances

```math
D^{2D}_{\mathrm{sig}} = \sqrt{\left(p^{\mathrm{PDB}}_{25} - p^{\mathrm{AF2}}_{25}\right)^{2} + \left(p^{\mathrm{PDB}}_{75} - p^{\mathrm{AF2}}_{75}\right)^{2}}
```

and D<sub>sig</sub><sup>3D</sup>, which adds the p<sub>50</sub> term.

- **Reference packing modes.** Ward agglomerative clustering into 4 clusters of the PDB signatures (p<sub>25</sub>, p<sub>75</sub>) at r = 10, k = 0, W = 1. Clusters are ordered by mean packing, from Mode 0 (lowest) to Mode 3 (highest); the report labels them IDPs, molten globules, dynamic globulars and rigid globulars. These assignments are held fixed for every configuration.
- **Per-configuration statistics.** For each (r, k, W): mean D<sub>sig</sub>; median DTW and EMD per mode; for ten metrics (the five percentiles and the five persistence metrics) a two-sided Wilcoxon signed-rank test of AF2 against the PDB, with Cliff's delta as effect size; and the number of flagged proteins with their mean and maximum hinge-angle deviation.
- **Sweeps and grid search.** One-at-a-time sweeps over W (r = 10, k = 0), k ∈ {0, 1, 2, 5} (r = 10, W = 1) and r ∈ {10, 15} (k = 0, W = 1), then a full grid over r ∈ {10, 15}, k ∈ {0, 1, 2, 5} and W ∈ {1, 5, 11, 15, 25, 31}. The configuration with the lowest mean D<sub>sig</sub><sup>2D</sup> is printed.
- **Ablation of k.** At r = 10 and W ∈ {1, 5, 15, 25}, k = 0 is compared with k ∈ {1, 2, 5} by paired Wilcoxon tests and Cliff's delta on DTW, EMD and both D<sub>sig</sub> (printed to the console).
- **Confusion matrices.** At r = 10, for each k ∈ {0, 1, 2, 5} and each W, the PDB and AF2 signatures are clustered independently (Ward, 4 clusters, ordered by packing) and cross-tabulated; agreement is the fraction of proteins on the diagonal.

---

## Known limitations

- **AF2 and PDB profiles cover different residues in Phase 2.** AF2 profiles span the full UniProt sequence, whereas PDB profiles contain only the aligned, resolved residues (the 40% coverage filter is relative to the PDB chain, so short constructs of long proteins pass). The Phase 2 signatures and persistence metrics, and therefore D<sub>sig</sub> and the paired tests in Phase 4, compare the full AF2 chain with the resolved part of the experimental chain. Count-based metrics (`Max_CC`, `PLM`) grow with chain length. AF2 arities also count neighbours that are missing from the PDB chain, and in the persistence step removing NaN positions joins the resolved segments on either side of a gap. Phase 3 (EMD, DTW) restricts both profiles to shared positions, but its AF2 arities are still computed with all AF2 neighbours.
- **Arity scale depends on r and k.** Larger k and smaller r give smaller arity values and therefore smaller absolute DTW, EMD and D<sub>sig</sub>. Comparisons of these errors across configurations (sweeps, ablation, choosing the configuration with the lowest D<sub>sig</sub>) and the fixed DTW threshold of 1.0 are affected by this. Normalising by the arity scale of each configuration would separate scale from agreement.
- **Independent hinges.** The hinge is located separately in the AF2 and PDB profiles, so the two angles can be measured around different residues.
- **PDB entry selection.** The first PDB ID listed in each UniProt entry is used, rather than the structure with the best coverage or resolution.
- **Multiple testing.** p-values are reported without correction across metrics and configurations. With thousands of proteins almost every test is significant, so Cliff's delta is the more informative number.
- **Terminology.** The output column is named `Hallucination_Flag`, but the flag marks proteins whose AF2 model itself reports a high PAE between distant regions, i.e. cases where AF2 signals its own uncertainty about inter-domain placement.

## Data sources and references

Structures and PAE files are downloaded at run time and are not redistributed in this repository; they remain subject to the terms of their providers.

- Cazals, F., & Sarti, E. (2024). AlphaFold predictions on whole genomes at a glance: a coherent view on packing properties, pLDDT values, and disordered regions. *bioRxiv*. https://doi.org/10.1101/2024.11.16.623929
- Jumper, J., Evans, R., Pritzel, A., et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature*, 596, 583–589. https://doi.org/10.1038/s41586-021-03819-2
- Varadi, M., Bertoni, D., Magana, P., et al. (2024). AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. *Nucleic Acids Research*, 52(D1), D368–D375. https://doi.org/10.1093/nar/gkad1011
- Worldwide Protein Data Bank (wwPDB) — https://www.wwpdb.org; structures downloaded from RCSB PDB — https://www.rcsb.org
- UniProt — https://www.uniprot.org

## License

The code is released under the MIT License; see [`LICENSE`](LICENSE).