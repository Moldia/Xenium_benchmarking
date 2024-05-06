# Xenium_benchmarking

This is a public repository to reproduce the analysis presented in the study available here:

Marco Salas et al. **Optimizing Xenium In Situ data utility by quality assessment and best practice analysis workflows** , 2024.

**Abstract of the study**

The Xenium In Situ platform is a new spatial transcriptomics product commercialized by 10X Genomics capable of mapping hundreds of genes in situ at a subcellular resolution. Given the multitude of commercially available spatial transcriptomics technologies, recommendations in choice of platform and analysis guidelines are increasingly important. Herein, we explore 25 Xenium datasets generated from multiple tissues and species comparing scalability, resolution, data quality, capacities and limitations with eight other spatially resolved transcriptomics technologies and commercial platforms. In addition, we benchmark the performance of multiple open source computational tools, when applied to Xenium datasets, in tasks including preprocessing, cell segmentation, selection of spatially variable features and domain identification. This study serves as the first independent analysis of the performance of Xenium, and provides best practices and recommendations for analysis of such datasets.

***

## Datasets
The Xenium datasets used in though this study combined datasets provided by 10X Genomics, datasets published elsewhere and datasets generated specifically for this project. The original files can be found in:
- **10X Genomics datasets**: Most of the datasets used in this study were publically available in https://www.10xgenomics.com/datasets [03.05.2024]
- **Published datasets**: the spinal chord datasets used were originally published in Kukanja, Langseth et al. (Cell, 2024) (  https://doi.org/10.1016/j.cell.2024.02.030)
- **Freshly generated datasets**: four mouse brain sections, labelled as "hm" through the study, were profiled for this analysis. Their original data can be downloaded from: https://doi.org/10.5281/zenodo.10566172
 
To facilitate the reproducibility of our analysis, we also provided already pre-formated AnnData files for each dataset used in the study, available at [ADD ZENODO REPO UPON PUBLICAION]

An **example dataset** (human spinal chord) used as an input for the end-to-end pipeline developed can be downloaded from  https://doi.org/10.5281/zenodo.11120922

***

## Folder structure


***

## Cloning and adding
Please first clone the environment by running in the terminal:

```git clone https://github.com/Moldia/Xenium_benchmarking.git```

Navigate to the folder:

```cd Xenium_benchmarking```

create a conda environment using the provided .yml file by: 

```conda env create --name xb --file=xenium_benchmarking.yml```

Acivate the conda environment by: 

```conda activate xb```

And install using pip:

```pip install -e . ```

