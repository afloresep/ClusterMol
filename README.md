## ClusterMol

ClusterMol is a Python command-line tool for clustering molecules based on their structural features. It supports both fingerprint-based clustering (Butina algorithm) and property-based clustering (K-Means on MQN fingerprints), with UMAP for dimensionality reduction and visualization.

---

### Features

* **Input:** Reads molecular data from CSV/Excel files with SMILES strings.
* **Fingerprint generation:** Computes molecular fingerprints (e.g., RDKit) or MQN descriptors.
* **Dimensionality reduction:** Applies UMAP to project high-dimensional fingerprints into lower dimensions for visualization or clustering.
* **Clustering algorithms:**

  * **Butina:** Gold-standard, threshold-based clustering using Tanimoto similarity.
  * **K-Means:** Centroid-based clustering on MQN descriptors, with user-defined cluster count.
* **Centroid calculation:** Optionally iterates to refine centroids for stability.
* **Output:** Exports an Excel file with:

  * Original molecular data.
  * Assigned cluster ID for each molecule.
  * Distance of each molecule to its cluster centroid.

---

### Requirements

* **Python:** 3.8 or higher
* **Libraries:**

  * RDKit
  * pandas
  * scikit-learn
  * umap-learn
  * openpyxl

Install dependencies with:

  ```bash
  pip install -r requirements.txt
  ```

---

### Installation

1. Clone the repository:

 ```bash 
  git clone [https://github.com/your-org/ClusterMol.git](https://github.com/your-org/ClusterMol.git)
  cd ClusterMol
 ``` 

2. Install requirements:

   ```bash
    pip install -r requirements.txt
    ```

---

### Usage

Run the main script from the `scripts/` directory:

```bash
python scripts/main.py [OPTIONS]
```

#### Required arguments

* `--dataframe PATH`
  Path to the input file (.csv or .xlsx) containing molecular data. Must include a SMILES column.

* `--cluster {Butina,KMeans}`
  Choose clustering algorithm.

#### Optional arguments

* `--smilescol NAME`
  Column name for SMILES strings (default: `SMILES`).
* `--threshold FLOAT`
  (Butina only) Tanimoto similarity cutoff (0.1–1.0). Default: 0.5.
* `--num-clusters INT`
  (K-Means only) Number of clusters to generate.
* `--output DIR`
  Directory for the output file. Default: current working directory.
* `--umap-components INT`
  Number of UMAP dimensions (default: 2).
* `--centroid-iterations INT`
  Number of refinement iterations for centroids (default: 1).

---

### Examples

1. **K-Means** with 400 clusters, SMILES in column `smiles`:

  ```bash
  python main.py --dataframe ../data/DB_corrected.csv --cluster KMeans --num-clusters 400 --smilescol smiles --output ../results/
   ```

2. **Butina** with similarity threshold 0.7, SMILES in column `canonical_smiles`:
  ```bash
    python main.py  --dataframe ../data/DB_corrected.csv  --cluster Butina  --threshold 0.7  --smilescol canonical_smiles --output ../results/
  ```


### Visualization

* A helper script generates TMAP visualizations for cluster inspection.
<!--TODO: * UMAP plots can also be produced by passing `--umap-plot` to the main script. -->


### License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
