
# 2 Spatial annotation

***

This folder contains files and scripts used to perform spatial annotation of
different spatially-resolved datasets, including:
- Xenium replicate 1
- Xenium replicate 2
- Xenium replicate 3
- Vizgen
- Starmap
- Visium on FF tissue
- Visium with FISH
- Visium on FFPE tissue
- smFISH
- osmFISH
- exseq
- Hybriss
- baristaseq
- merFISH

All datasets were annotated using the spatial map of cell, colored by cell class. Annotation was done using polygons drawn using bézier curves and saved as scalable vector graphics (SVG). The folder `images` thus contain the cell maps, while the folder `SVG` contain such annotation polygons for each dataset, respectivelly.

Delimitation and naming of spatial regions of the brain was accomplished by manual comparisson to coronal P56 Allen brain reference at different sectioning depths (68, 75, 76, 77, 78), which are available in the folder `allen_brain_reference`.

Naming conventions from the Allen brain used for annotation herein are listed in the file `spatial_annotation.csv`, at different hierachical levels.

Finally, assigment of spatial annotation to cells and spots for each dataset was done matching relative spatial coordiantes from the original datasets to their respective SVG polygons (see code in `annotate_coordenates.Rmd`). Minor adjustments for this assignemnt were stored in the file `spatial_adjustement.csv`, for each technology respectivelly. 


