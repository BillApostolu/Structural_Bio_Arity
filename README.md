# The Project Pipeline: Multiscale Topological, Sequential & Quaternary Assessment of AlphaFold 2 Predictions

## Phase 1: Data Acquisition, Geometric Parsing & Synchronization

**The goal is to translate 3D Cartesian coordinates into raw, perfectly aligned, and heavily parameterized topological data.**

* **Data Acquisition & Target Selection:** Curate a dataset of roughly 7,500 proteins that have experimentally resolved structures in the Protein Data Bank (PDB) alongside predictions from AlphaFold 2 (AF2). Concurrently fetch all auxiliary files required for downstream phases, including the AlphaFold PAE JSONs.
* **Coordinate Parsing:** Parse the coordinate files to isolate the $C_{\alpha}$ atom coordinates.
* **Mathematical Sequence Synchronization:** Strictly map the experimental PDB sequence to the full AlphaFold sequence using an alignment algorithm, padding missing unresolved PDB indices with null values to maintain perfect spatial-sequential alignment.
* **Multiscale & Long-Range Arity Calculation:** Calculate both standard local arity and long-range arity for each $C_{\alpha}$ atom. Formally, let the arity of a $C_{\alpha}$ carbon be the number of neighboring $C_{\alpha}$ atoms within a distance range $r$. Assuming the polypeptide chain has $n$ amino acids, define $L_{a}=\{a_{1},...,a_{n}\}$ as the arities of the $n$ $C_{\alpha}$ carbons. Define long-range arity by isolating residues separated by a sequence distance of $|i-j| > k$ that fall within a spatial radius of $d \le r$. Apply distance-decay kernel weights to these connections to capture nuanced packing density alongside unweighted binary thresholds.
* **Data Structuring:** Instead of a monolithic matrix, save the sequence of arity values as independent, synchronized 1D spatial time-series arrays. Each array is uniquely parameterized across multiple pairs of sequence separation ($k$) and spatial radius ($r$), generating separate files for both weighted and unweighted calculations.

---

## Phase 2: Topological Signatures, Mapping & Local Filtration

**The goal is to extract geometric baselines, expanding the distribution metrics and introducing local topological smoothing to distinguish ordered from disordered regions.**

* **Expanded Distribution Generation:** Compute the normalized arity histogram and use an expanded Cumulative Distribution Function (CDF) suite. Let $L=[A_{1},...,A_{m}]$ be the unique arities sorted by increasing value. The number of $C_{\alpha}$ with a given arity defines the histogram $H_{C_{\alpha}}$. We define the $C_{\alpha}$ Arity Signature for a suite of specific percentiles $P_{K}=\{10th, 25th, 50th, 75th, 90th\}$ as the smallest arity such that the CDF is larger than $p_{k}$:

$$arity\_C_{\alpha}(p_{k}) = CDF_{C_{\alpha}}^{-1}(p_{k}) \equiv \arg\min_{j} (CDF_{C_{\alpha}}(A_{j}) \ge p_{k})$$


* **Global Arity Mapping (2D & 3D):** Plot the 25th and 75th percentiles to create a 2D Global Arity Map. Extend this to a 3D representation by including the 50th percentile. At this global level, proteins physically cluster based on their macroscopic packing density: rigid globular folds inherently congregate in the high-arity top-right region, while Intrinsically Disordered Proteins (IDPs) isolate themselves in the low-arity bottom-left region.
* **Sliding Window Locality (1D Signal Smoothing):** Apply a moving average filter across the synchronized 1D sequence arrays. The filter begins sliding dynamically based on the sequence gap $k$. We utilize varied window sizes ($W \in \{5, 11, 15, 25, 31\}$) that directly map to physical protein hierarchy:
* $W=5$: Captures local secondary structure ($\alpha$-helices and $\beta$-strands).
* $W=11 \text{ to } 15$: Targets super-secondary motifs.
* $W=25 \text{ to } 31$: Acts as a heavy low-pass filter ideal for isolating macro-domains.


Mathematically, the smoothed arity at residue $i$ is $S_i = \frac{1}{W} \sum_{j} a_j$. At this local level, deep topological troughs isolate Intrinsically Disordered Regions (IDRs), while sustained high-arity plateaus map to rigid domains.
* **Arity Filtration:** Sort the $C_{\alpha}$ atoms by increasing arity and run the Union-Find algorithm on this sorted sequence to track connected components. Using a filtration $\mathcal{G}_{u}$ encoding a sequence of nested subgraphs, let $N_{cc}(u)$ represent the number of connected components at each value of $u$. Extract topological quantifications such as Persistence Entropy, Mean Persistence, Fraction of Positive Critical points, Maximum number of connected components, and Persistent Local Maxima.

* **Global Map Evaluation (Signature Distance):** To systematically evaluate the vast hyperparameter space $(r, k, W)$ without relying on visual inspection, calculate the Euclidean Signature Distance ($D_{sig}$) between the AF2 and PDB distributions. For a 2D map, this is $D_{sig} = \sqrt{(p_{25, PDB} - p_{25, AF2})^2 + (p_{75, PDB} - p_{75, AF2})^2}$. The parameter pair yielding the lowest mean $D_{sig}$ across the dataset identifies where AF2 most accurately captures global packing density.

---

## Phase 3: Advanced Structural Alignment & Quaternary Validation

**The goal is to physically validate spatial arrangements using combined metrics and efficient mathematical formulations.**

* **Global Packing Similarity (1D EMD):** Use the Earth Mover’s Distance (EMD) to compare the full CDFs between AF2 and PDB models. The distance is calculated as the exact "work" required to match the distributions:

$$EMD(PDB, AF2) = \int_{0}^{\infty} |CDF_{PDB}(x) - CDF_{AF2}(x)| \, dx$$


* **Multivariate Spatial Motif Alignment (DTW):** Deploy Dynamic Time Warping on the heavily smoothed macrotopology arrays from Phase 2. To handle unresolved PDB coordinates without distorting the alignment, missing sequence blocks are uniformly excised from both the experimental and predicted arrays prior to DTW, evaluating only the shared resolved subsets. The cost matrix incorporates a penalty multiplier $\gamma=2.0$ driven by AlphaFold's pLDDT, mathematically forcing the algorithm to forgive structural deviations in intrinsically disordered regions where AF2's confidence is low.

$$Cost(i, j) = \left(1 + \gamma \frac{100 - pLDDT_{AF2, j}}{100}\right) (Arity_{PDB, i} - Arity_{AF2, j})^2$$


* **Hallucination Flagging:** Use AlphaFold's cross-domain Predicted Aligned Error (PAE) to flag models with low DTW error (perfect local folding) but high cross-domain PAE (uncertain quaternary orientation). The pipeline triggers a hallucination flag under the following conditions:

$$(\text{DTW Error} < \text{Threshold}) \land (\text{Cross-Domain PAE} > \text{Threshold})$$


* **Quaternary Physical Validation:** Validate flagged hallucinations by calculating the rigid-body angle between the centers of mass of the distinct domains. Instead of using arbitrary sequence cutoffs to find the hinge, locate the inter domain linkers using **Topological Prominence** same as Persistence (via signal processing on the inverted smoothed arity array). By enforcing a prominence threshold of $2.0$, this continuous 1D minimum perfectly mirrors the discrete Persistent Local Maxima ($T_{\nu} = 2.0$) from Phase 2, isolating the physical hinge coordinate $H$ to compute the macroscopic fold angle $\theta$. Then flagged hallucinations are validated by calculating the rigid body angle between the centers of mass of the distinct domains in the AF2 model, benchmarking it directly against the ground truth angle derived from the experimental PDB coordinates in Phase 4. The macroscopic fold angle $\theta$ is computed via the dot product of domain vectors originating from a shared hinge point $H$:

$$\theta = \arccos \left( \frac{\vec{v}_A \cdot \vec{v}_B}{|\vec{v}_A| |\vec{v}_B|} \right)$$


---

## Phase 4: Statistical Benchmarking & Structural Ablation

**The goal is to statistically validate AlphaFold's topological deviations from experimental reality and prove the efficacy of the pipeline's novel parameters.**

- **Geometric K-Means Clustering:** Deploy unsupervised K-Means clustering ($K=2$) strictly on the geometric Arity Signatures ($p_{25}, p_{75}$) of the dataset. This mathematically partitions the protein universe into Intrinsically Disordered Proteins (IDPs) and Rigid Globular domains based purely on packing density.
- **Statistical Benchmarking (Distance Metrics):** Apply independent non-parametric tests (Mann-Whitney U test) to the 1D EMD and DTW error distributions to compare alignment accuracy across the distinct structural classifications. This statistically determines if AlphaFold's predictive error is systematically higher in intrinsically disordered regions compared to rigid domains.
- **Statistical Benchmarking (Paired Topological Metrics):** Group the Filtration metrics (Persistence Entropy, Mean Persistence, Persistent Local Maxima) into paired sets (AF2 vs. PDB). Apply non-parametric paired tests (Wilcoxon signed-rank), pairing every $p$-value with an Effect Size metric (Cliff’s Delta) to systematically quantify where AlphaFold's internal topological wiring deviates from experimental ground truth.
- **Quaternary Hallucination Analysis:** Isolate the subset of proteins flagged as hallucinations in Phase 3. Statistically analyze the distribution of their macroscopic rigid-body angle deviations ($\theta$) to quantify the severity of AlphaFold's quaternary misalignments.
- **Methodological Ablation Study:** Systematically test whether the addition of multiscale long-range arity yields statistically superior structural classification compared to standard local arity. Using a paired Wilcoxon signed-rank test on the DTW alignment errors, compare predictions processed with standard arity ($k=0$) against those processed with long-range isolating gaps ($k>0$) to quantify the magnitude of information gain.