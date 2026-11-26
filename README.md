# The Proteomics-informed Renal tumor Integrative Subtyping Model (PRISM)

We developed the Proteomics-informed Renal tumor Integrative Subtyping Model (PRISM), which uses subtype-specific enrichment scores of each omics to classify non-clear cell renal tumor (nccRCC) NMF subtypes, and evaluated individual omics layers and their combinations using leave-one-out threshold fitting.

## File Structure Overview

`calculate_enrichment_scores.R`: Calculate the normalized enrichment scores of each omics for each sample using the Mann–Whitney–Wilcoxon Gene-Set Test (MWW–GST).

`calculate_probability.py`: Calculate the probability of each sample belonging to each nccRCC NMF subtype based on the enrichment scores of different omics.

`experiment.ipynb`: Experiment to select the optimal combination of omics layers for classifying nccRCC NMF subtypes.

## 📖 Citation

If you find this work useful, please consider citing:

```
@article{}
```

---

> **⚠️ Reproducibility Note:** The code in this repository may undergo minor modifications during analysis. For methodological details, please refer to the *Methods* section of the original publication.
