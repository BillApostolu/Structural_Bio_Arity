Breakdown of every column generated in the final CSVs

### 1. Hyperparameters & Metadata

* **r:** The spatial radius in Angstroms used to define the hydration shell or packing neighborhood (e.g., 10 or 15).
* **k:** The sequence gap parameter used to isolate long-range contacts. `k=0` is the standard Cazals baseline; `k>0` filters out local secondary structure.
* **W:** The sliding window size used as a low-pass filter to smooth the 1D signal (e.g., W=1 for atoms, W=25 for domains).
* **Total_Proteins:** The number of valid, paired AF2 and PDB proteins successfully evaluated for this specific configuration.

### 2. Macroscopic Geometric Distance

* **D_sig_2D_Mean / D_sig_3D_Mean:** The overall average Euclidean distance between AlphaFold's arity percentiles and the PDB's arity percentiles across the whole dataset. A lower number means AF2 globally captured the experimental packing density more accurately.

### 3. Structural Class Benchmarking (IDP vs. Globular)

*These columns test whether AlphaFold's predictive error is systematically worse for IDPs compared to Rigid Globular domains using the Mann-Whitney U test.*

* **DTW_IDP_Med / DTW_Glob_Med:** The median Dynamic Time Warping 1D alignment error for the IDP cluster vs. the Globular cluster.
* **EMD_IDP_Med / EMD_Glob_Med:** The median Earth Mover's Distance global distribution error for the IDP cluster vs. the Globular cluster.
* **D_sig_2D_IDP_Med / D_sig_2D_Glob_Med:** The median 2D signature distance for the IDP cluster vs. the Globular cluster. (Same logic applies to the **D_sig_3D** columns).
* **[Metric]_MW_pval:** The p-value from the Mann-Whitney U test proving if the difference between IDPs and Globulars is statistically significant.
* **[Metric]_Cliffs:** Cliff's Delta effect size. A positive number indicates the error is systematically higher for IDPs than for Globular proteins.

### 4. Paired Topological Wiring (AF2 vs. PDB)

*These columns pair AF2 predictions directly against their experimental PDB counterparts to evaluate internal topology using the Wilcoxon signed rank test.*

* **Entropy_Hp_AF2_Med / PDB_Med:** The median Persistence Entropy of the topological diagrams.
* **Mean_Persistence_AF2_Med / PDB_Med:** The average "lifespan" of topological features (how long independent structural islands survive before merging).
* **Fraction_CP_AF2_Med / PDB_Med:** Fraction of positive Critical Points (the frequency of domain-merging events).
* **Max_CC_AF2_Med / PDB_Med:** The Maximum number of Connected Components (how heavily fragmented the protein's folding topology gets at its peak).
* **PLM_AF2_Med / PDB_Med:** Persistent Local Maxima. Count of major structural peaks that survive the noise threshold.
* **p_10_... through p_90_AF2_Med / PDB_Med:** The raw arity percentiles. (As you noticed, AF2 systematically under-packs, so the AF2 medians will generally be lower here).
* **[Metric]_pval:** The Wilcoxon p-value proving if AlphaFold systematically deviates from experimental reality for that specific metric.
* **[Metric]_Cliffs:** Cliff's Delta effect size. A negative number means AlphaFold systematically predicts *lower* values than the PDB.

### 5. Quaternary Hallucination Validation

* **Total_Hallucinations:** The absolute count of proteins flagged as having a perfect local fold but an uncertain, unphysical global domain orientation (based on your DTW and PAE thresholds).
* **Mean_Angle_Dev:** For the hallucinated proteins only, this is the average deviation (in degrees) of the macroscopic hinge angle between AlphaFold's prediction and the experimental PDB reality.

---