# Master Thesis Code: Prediction and Visualization of Molecular Substructures using Mass Spectrometry

This repository contains the code, benchmarks, visualizations, and experiments for my master’s thesis in Media Informatics at Hochschule Düsseldorf:

> **Prediction and Visualization of Molecular Substructures using Mass Spectrometry**  
> Kevin R. Zielke, 2025  
> Supervised by Prof. Dr. Florian Huber, Advisors David Sokalski, M.Sc.  
> DOI: [10.13140/RG.2.2.20175.34725](https://doi.org/10.13140/RG.2.2.20175.34725)

The work focuses on improving the interpretation of molecular similarity in mass spectrometry by introducing **Stacked Similarity Maps (SSM)** and related metrics for group-based similarity assessment in chemical space exploration.  

This repository primarily serves as:
- A transparent record of the experimental methodology
- A reference implementation for the concepts introduced in the thesis
- Supplementary material for readers of the thesis

---

## Reproducibility & Usage
The notebooks have some dependencies to code I implemented in the backend of the [ms_chemical_space_explorer](https://github.com/florian-huber/ms_chemical_space_explorer) web application (which contains a cleaner and more modular implementation than this exploratory research code). Unfortunately this project is currently under development and private. Therefore you can't run the notebooks yourself. However this might change in the future, so here is how you would go about it:

### Setup
Download `compounds_ms2structures.csv` and `biostructures_combined.csv` from [Zenodo](https://doi.org/10.5281/zenodo.15650045) and place them under `data/datasets/`.

To create and activate the conda environment:

```bash
conda env create -f chemspace.yml
conda activate chemspace
```

This environment includes the required Python version, RDKit, and the scientific Python stack used in the experiments. The notebooks create some files on disk and are intended to run sequentially from `00_biostructures_dataset_selection.ipynb` to `06_example_vis.ipynb`. 
