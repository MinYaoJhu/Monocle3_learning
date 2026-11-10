# Monocle 3: study notes <img src="assets/logo.svg" align="right" width="130" alt="Monocle 3 tutorial notes badge"/>

**My worked tutorials for single-cell clustering, trajectory and pseudotime analysis with Monocle 3**

[![Monocle 3](https://img.shields.io/badge/Monocle%203-cole--trapnell--lab.github.io-0F3B3A)](https://cole-trapnell-lab.github.io/monocle3/)
![R](https://img.shields.io/badge/R-4.5-276DC3?logo=r&logoColor=white)
![Seurat](https://img.shields.io/badge/Seurat-v5-blue)

In September 2025 I worked through the [Monocle 3](https://cole-trapnell-lab.github.io/monocle3/) documentation from the Trapnell Lab, using the *C. elegans* example datasets, and then the SeuratWrappers vignette that connects Seurat v5 objects to Monocle 3. Each notebook (`.Rmd`) follows a published tutorial with my notes and edits; the rendered `.md` files can be read directly on GitHub.

## 📚 Notebooks

| # | Notebook | Data | Topics |
|---|---|---|---|
| 0 | [Getting started](0_Getting_started_with_Monocle_3.md) | *C. elegans* L2 larvae (Cao et al. 2017) | Installing Monocle 3, the `cell_data_set` class, building a CDS object |
| 1 | [Clustering and classifying cells](1_Clustering_and_classifying_your_cells.md) · [part 2](1_Clustering_and_classifying_your_cells-2.md) | *C. elegans* L2 larvae | Pre-processing, UMAP, batch correction, clustering, marker genes, cell-type annotation, selecting cell subsets; part 2: automated annotation with Garnett (writes [`marker_file.txt`](marker_file.txt)) |
| 2 | [Constructing single-cell trajectories](2_Constructing_single-cell_trajectories.md) · [part 2](2_Constructing_single-cell_trajectories-2.md) | *C. elegans* embryo (Packer et al. 2019) | Learning the principal graph, ordering cells in pseudotime, choosing the root programmatically, 3D trajectories |
| 3 | [Differential expression analysis](3_Differential_expression_analysis.md) · [part 2](3_Differential_expression_analysis-2.md) | *C. elegans* embryo and L2 larvae | Regression models, controlling for batch, model comparison, graph-autocorrelation across clusters; part 2: genes that change along pseudotime and gene modules in the AFD neuron lineage |
| 00 | [Monocle 3 with Seurat](00_Calculating-Trajectories-with-Monocle-3-and-Seurat.md) | HCA bone marrow 40k (SeuratData) | Converting a Seurat v5 object to a CDS with SeuratWrappers, trajectory and pseudotime on a human dataset |

## ▶️ Running the notebooks

Run the notebooks in order (0 → 3); each part-2 notebook continues from its part 1.

- **Notebook 1:** part 1 selects a subset of cells interactively with `choose_cells()` and saves `1_cds.rds` and `1_cds_subset.rds`, which part 2 loads. These files are not stored in the repository; re-run part 1 to create them.
- **Notebooks 2 and 3:** the tutorials choose the trajectory root interactively with `order_cells(cds)`. To make the results exactly reproducible, the ordered objects were saved with `save_monocle_objects()` and are included here: [`2_cds_after_ordering_monocle_save/`](2_cds_after_ordering_monocle_save), [`3_cds_after_ordering_monocle_save/`](3_cds_after_ordering_monocle_save) and [`3_AFD_lineage_cds_after_ordering_monocle_save/`](3_AFD_lineage_cds_after_ordering_monocle_save). The part-2 notebooks load them with `load_monocle_objects()`. Notebook 2 part 2 also shows how to choose the root programmatically (`get_earliest_principal_node()`), which avoids the interactive step.
- **Notebook 00:** install the dataset package from SeuratData:
  ```r
  install.packages("https://seurat.nygenome.org/src/contrib/hcabm40k.SeuratData_3.0.0.tar.gz",
                   repos = NULL, type = "source")   # 84 MB
  ```

## 🧠 Skills practiced

- **Single-cell workflow in Monocle 3:** pre-processing, dimensionality reduction, batch correction, clustering and marker detection
- **Cell-type annotation:** marker-based and automated (Garnett) classification
- **Trajectory inference:** principal graphs, pseudotime, root selection, 3D trajectories
- **Differential expression:** generalized linear models, graph-autocorrelation tests and co-regulated gene modules
- **Interoperability:** moving between Seurat v5 and Monocle 3 with SeuratWrappers

## 📂 Data

- ***C. elegans* datasets** (notebooks 0–3) are downloaded automatically by the notebooks from the Trapnell Lab website (`cao_l2_*` and `packer_embryo_*` `.rds` files).
- **HCA bone marrow 40k** (notebook 00) is installed from [SeuratData](https://seurat.nygenome.org/src/contrib/) as shown above.

## License and credits

- Notebooks 0–3 follow the [Monocle 3 documentation](https://cole-trapnell-lab.github.io/monocle3/) by the Trapnell Lab (Monocle 3, MIT License).
- Notebook 00 is adapted from the SeuratWrappers vignette "Calculating Trajectories with Monocle 3 and Seurat" by the Satija Lab and remains under **GPL-3.0**, the [seurat-wrappers](https://github.com/satijalab/seurat-wrappers) license.
- My own notes and code are released under the MIT License. See [`NOTICE.md`](NOTICE.md) for details.

## AI assistance

The README, repository clean-up, path fixes and the badge were prepared with the assistance of Claude Opus 5.5 (Anthropic), used through Claude Code. The study notes were written by Min-Yao Jhu. All changes were reviewed and approved by the author.

## Contact

**Min-Yao Jhu**, Senior Research Scientist, Donald Danforth Plant Science Center · [minyaojhu.com](https://www.minyaojhu.com) · [GitHub](https://github.com/MinYaoJhu)
