# AMP-structural-bias-audit
Datasets and codes for auditing structural bias on AMP ML classifiers

# Data availability
- validation_datasets include the FASTA sequences of external validation datasets GRAMPA and non-GRAMPA.
- training_datasets include the FASTA sequences of training datasets AMPlify (balanced), PepNet, AMP scanner v2, and ampPEPpy.
- fold_predictions include the predicted structural classes or folds (1, 5, 6, 7) associated with GRAMPA and non-GRAMPA datasets.
- activity_predictions include the predicted activity (AMP or non-AMP) with class probability associated with GRAMPA and non-GRAMPA datasets from the 16 classifiers.
- ProstT5-embeddings are the calculated ProsT5 embeddings for GRAMPA, non-GRAMPA, AMPlify (balanced), PepNet, AMP scanner v2, and ampPEPpy datasets.

# Code availability
- R script to project GRAMPA and non-GRAMPA datasets onto a ternary plot representation.
- Python script to produce UMAP projections of ProstT5-derived 3Di structural embeddings for the training and validation datasets.
