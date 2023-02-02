# Xenium_benchmarking

This is a public repository for the following study available here:

Marco Salas et al. **Evaluating spatial gene expression heterogeneity in Xenium In Situ data** (bioRxiv), 2023.

**Abstract**

Xenium is a new spatial transcriptomics product commercialized by 10X Genomics capable of generating maps of hundreds of transcripts in situ at a subcellular resolution. Herein, we explore  5 Xenium datasets of the mouse brain and  2 of human breast cancer by comparing scalability, resolution, data quality, capacities and limitations, with other spatially resolved transcriptomics (SRT) technologies. In addition, we benchmarked the performance of multiple open source computational tools when applied to Xenium datasets in tasks including cell segmentation, segmentation-free analysis, selection of spatially variable genes and domain identification, among others. To our knowledge, this study serves as the first independent analysis of the performance of Xenium.

***

## Datasets
The Xenium datasets used in this analysis were provided by 10X Genomics. 
The three mouse brain coronal sections (“ms brain multisection”) are publically available datasets and can be downloaded at from the 10X website (https://www.10xgenomics.com/xenium-preview-data , https://www.10xgenomics.com/resources/datasets). 
The mouse brain full coronal and half coronal sections (named as “ms brain coronal” and “ms brain ROI” in Figure 1B), as well as the human breast sections are available upon request to 10X Genomics.

***

## Folder structure


***

## Cloning and adding
In a clean conda environment with pip installed, run in the terminal:

```git clone https://github.com/Moldia/Xenium_benchmarking.git```

Navigate to the folder:

```cd Xenium_benchmarking```

And install using pip:

```pip install -e . ```

