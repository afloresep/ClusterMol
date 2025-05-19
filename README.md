# ClusterMol

ClusterMol is a Python tool designed to cluster molecules based on structural similarity, using UMAP dimensionality reduction and centroid-based cluster analysis. It reads molecular data (e.g., from DrugBank) with SMILES strings, performs clustering, computes centroids, and outputs a processed Excel file with cluster assignments and distances to centroids.

---

## Features

- Takes an Excel file with molecular data and SMILES as input.
- Converts SMILES to molecular fingerprints.
- Reduces dimensionality with UMAP.
- Clusters molecules (e.g., k-means).
- Calculates cluster centroids with optional iteration for stability.
- Outputs an Excel file with:
  - Cluster ID for each molecule.
  - Distance to cluster centroid.

---

## Requirements

- Python 3.8+
- RDKit
- pandas
- scikit-learn
- umap-learn
- openpyxl

You can install dependencies via pip:

```bash
pip install -r requirements.txt
```
