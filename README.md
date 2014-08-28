##Single-cell expression profiling and spatial mapping into tissue of origin
You will find in this repository the R source code used to perform the analysis presented in the manuscript. The different function and the global work flow will be presented using simple (made up) example datasets. The data used in the manuscript is also available but requires more computing power, some `bash` scripts to launch the computations in parallel on a LSF cluster system are also made available.

###Example datasets
The example dataset is located in the folder "example_dataset", you can obtain them by downloading this project as a .zip file or by cloning the git repository on your machine. The example dataset is divided in 3 files :
 - example_data_RNA_seq.csv : contains some fake RNA-seq counts values for 10 cells (rows) and 10 genes (columns), when applying the method to your data before this step your data should have been cleaned up with only the interesting genes overlapping the atlas selected, the cells selected following some sequencing quality controls, and potentially normalized.
 - example_data_atlas.csv : contains some fake binary expression data values for 1000 voxels (rows) and 10 genes (columns. Note: the genes are the same as the ones and in the same order as in the RNA-seq, this has to be the case)
 - example_3D_coordinates_atlas : contains some spatial coordinates (a sphere) for the 1000 voxels in the atlas. NOTE: there is no atlas cell ID, this means that the coordinates in this file need to be in the same order as the row in the `example_data_atlas.csv`

###Work flow in R
The following section will describe the entire work flow, assuming you have cloned the repository on your pc, simply open a terminal, go to the nbt_spatial_backmapping directory and start `R`

#####Loading the data
```R
rna_seq = read.table("example_dataset/example_data_RNA_seq.csv",header=TRUE,sep="\t")
atlas = read.table("example_dataset/example_data_atlas.csv",header=TRUE,sep="\t")
coordinates = read.table("example_dataset/example_data_RNA_seq.csv",header=TRUE,sep="\t")
```


