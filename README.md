##Single-cell expression profiling and spatial mapping into tissue of origin
You will find in this repository the R source code used to perform the analysis presented in the manuscript. The different function and the global work flow will be presented using simple (made up) example datasets. The data used in the manuscript is also available but requires more computing power.

###Example datasets
The example dataset is located in the folder "example_dataset", you can obtain them by downloading this project as a .zip file or by cloning the git repository on your machine. The example dataset is divided in 3 files :
 - `example_data_RNA_seq.csv` : contains some fake RNA-seq counts values for 10 cells (rows) and 10 genes (columns), when applying the method to your data before this step your data should have been cleaned up with only the interesting genes overlapping the atlas selected, the cells selected following some sequencing quality controls, and potentially normalized.
 - `example_data_atlas.csv` : contains some fake binary expression data values for 1000 voxels (rows) and 10 genes (columns. Note: the genes are the same as the ones and in the same order as in the RNA-seq, this has to be the case)
 - `example_3D_coordinates_atlas` : contains some spatial coordinates (a sphere) for the 1000 voxels in the atlas. NOTE: there is no atlas cell ID, this means that the coordinates in this file need to be in the same order as the row in the `example_data_atlas.csv`

###Work flow in R
The following section will describe the entire work flow, assuming you have cloned the repository on your pc, simply open a terminal, go to the nbt_spatial_backmapping directory and start `R`

#####Loading the data and the needed functions
```R
#loading data
rna_seq <- read.table("example_dataset/example_data_RNA_seq.csv",header=TRUE,sep="\t")
atlas <- read.table("example_dataset/example_data_atlas.csv",header=TRUE,sep="\t")
coordinates <- read.table("example_dataset/example_3D_coordinates_atlas",header=TRUE,sep=",")

#Loading the analysis functions
source("spatial_mapping.R")
```

#####Computing the specificity scores on the RNA-seq data
```R
specificity_matrix <- specificity_scores(rna_seq)
```

#####Mapping for one cell
```R
#Take the first sequenced cell informations
specificity_score_cell_1 <- specificity_matrix[1,]

#You need a binary expression vector for this cell, the count threshold is for you to decide given your RNA-seq data
example_count_threshold <- 10
binary_expression_cell_1 <- sapply(rna_seq[1,],function(x){ifelse(x>example_count_threshold,1,0)})

#Optionally you can compute the specificity scores on the atlas side too to penalize the mismatches between the RNA-seq and the atlas in both ways
number_of_points_in_altas <- 1000
atlas_specificity_score <- 1/apply(atlas,2,function(x){length(x[x>0])/number_of_points_in_altas})
#If you decide not to, simply set it as a vector of 0s
atlas_specificity_score <- as.numeric(vector(length=number_of_points_in_altas))

#Launch the mapping against the atlas for cell 1
mapping_result_cell_1 <- spatial_map_scoring(specificity_score_cell_1,binary_expression_cell_1,bin.ratio2)
```


