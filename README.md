##Single-cell expression profiling and spatial mapping into tissue of origin
You will find in this repository the R source code used to perform the analysis presented in the manuscript. The different function and the global work flow will be presented using simple (made up) example datasets. The data used in the manuscript is also available but requires more computing power.

The workflow presented below doesn't present the functions used in details, for that, please see the annotated code in the file `spatial_mapping.R`

###Example datasets
The example dataset is located in the folder `example_dataset`, you can obtain them by downloading this project as a .zip file or by cloning the git repository on your machine. The example dataset is divided in 3 files :
 - `example_data_RNA_seq.csv` : contains some fake RNA-seq counts values for 10 cells (rows) and 10 genes (columns), when applying the method to your data before this step your data should have been cleaned up with only the interesting genes overlapping the atlas selected, the cells selected following some sequencing quality controls, and potentially normalized.
 - `example_data_atlas.csv` : contains some fake binary expression data values for 1000 voxels (rows) and 10 genes (columns. Note: the genes in this file need to be in the same order as the ones  in the RNA-seq file.)
 - `example_3D_coordinates_atlas` : contains some spatial coordinates (a sphere) for the 1000 voxels in the atlas. NOTE: there is no atlas cell ID, this means that the rows in this file need to be in the same order as the row in the `example_data_atlas.csv`

###Work flow in R
The following section will describe the entire work flow, assuming you have cloned the repository on your pc, simply open a terminal, go to the nbt_spatial_backmapping directory and start `R`

#####Loading the data and the needed functions
```R
#loading data
rna_seq <- read.table("example_dataset/example_data_RNA_seq.csv",header=TRUE,sep="\t")
atlas <- read.table("example_dataset/example_data_atlas.csv",header=TRUE,sep="\t")
coordinates <- read.table("example_dataset/example_3D_coordinates_atlas",header=FALSE,sep=",")

#Loading the analysis functions
source("spatial_mapping.R")
```

#####Computing the specificity scores on the RNA-seq data
```R
specificity_matrix <- specificity_scores(rna_seq)
```

#####Mapping the 10 cells
```R
#Optionally you can compute the specificity scores on the atlas side too to penalize the mismatches between the RNA-seq and the atlas in both ways
number_of_points_in_altas <- 1000
atlas_specificity_score <- 1/apply(atlas,2,function(x){length(x[x>0])/number_of_points_in_altas})
#If you decide not to, simply set it as a vector of 0s
atlas_specificity_score <- as.numeric(vector(length=number_of_points_in_altas))

#Iterate over the sequenced cells
example_results_scores = sapply(seq_along(rna_seq[,1]),function(cell_num) {

	print(paste("Mapping cell",cell_num))
	cell = rna_seq[cell_num,]	

	specificity_score <- specificity_matrix[cell_num,]

	#You need a binary expression vector for this cell, the count threshold is for you to decide given your RNA-seq data
	example_count_threshold <- 10
	binary_expression_cell <- sapply(rna_seq[1,],function(x){ifelse(x>example_count_threshold,1,0)})

	#Launch the mapping against the atlas for the cell
	mapping_result_cell <- spatial_map_scoring(specificity_score,binary_expression_cell,atlas_specificity_score,atlas)

	#Save the result in a file
	write.table(file=paste("mapping_results_example/mapping_result_cell_",cell_num,sep=""),mapping_result_cell)
	
	mapping_result_cell
})

```
NOTE: Everything here is done sequencially but of course the scores for each cell can be computed separately in parallel, which would drastically decrease the computational time.

####Simulating data to find the confidence thresholds
Now you have the object `example_results_scores` which contains the score for each sequenced cell (columns) against every reference voxel in the atlas (rows). You still need to find the thresholds above which you can assume that a voxel is a match for a cell. To do this, you will generate a null distribution from your data by simulating random sequenced cells. First create a directory to store the generated simulated data, for instance `/example_simluated_data`.

Then to generate simulated cells from you RNA-seq data simply run in R to generate 100 cells for each of the 10 real cells
```R
generate_simulated_data(specificity_matrix,100,"example_simulated_data/")
```

This command will create 10 datasets containing 100 simulated cells each. Each dataset has two file:
 - n.data : is the table of specificity scores
 - n_bin.data : is the table of binary expression infered from the specificity scores

You can then map the simulated datasets. For example to map all the simulated cells from all the generated datasets (which corresponds to 100 simulated cells generated from a random sampling of the specificity score of all the sequenced cells)
```R
#Iterate over the 10 datasets
cell_counter = 1
mapping_simulated_list = list()
for(dataset in 1:10) {
	randomCells = read.table(paste("example_simulated_data/",dataset,"_bin.data",sep=""))
	randomRatios = read.table(paste("example_simulated_data/",dataset,".data",sep=""))

	#map each cell
	for(cell_i in 1:100) {
		print(cell_i)
		mapping_simulated_list[[cell_counter]] <- spatial_map_scoring(randomRatios[cell_i,],randomCells[cell_i,],atlas_specificity_score,atlas)
		cell_counter = cell_counter +1
	}
}

#convert to a dataframe for easier use
mapping_simulated = data.frame(mapping_simulated_list)
colnames(mapping_simulated)=paste("cell_",1:length(mapping_simulated_list),sep="")
```

#####Choosing the threshold
From the previous step you will able to plot the proportion of simulated cells scoring higher than different thresholds, which will allow you you choose sensible threshold(s) for your dataset. For instance you can see the proportion of your simulated above the thresholds (-2,-1,0,1,2) having at least (1,6,11) match in the atlas with the follwing command
```R
#creating the different thresholds
threshold_score <- seq(-2,2,1)
threshold_number <- seq(1,11,5)


#Computing the number of simulated cells meeting the requirements for every combination of threshold_score and threshold_number
number_cells_above = matrix(nrow=length(threshold_score),ncol=length(threshold_number))
for(j in 1:length(threshold_number)){
	print(j)
	number_cells_above[,j] = unlist(lapply(seq_along(threshold_score),function(i){
		sup = apply(mapping_simulated,2,function(x){
			length(x[x>threshold_score[i]])
		})
		length(sup[sup>threshold_number[j]])
	}))
}
number_cells_above = data.frame(number_cells_above)
rownames(number_cells_above) = as.character(threshold_score)
colnames(number_cells_above) = as.character(threshold_number)

#Compute the proportions
proportion_cells_thresholds = number_cells_above/(dim(mapping_simulated)[2])
proportion_cells_thresholds
```

The following functions assume you chose 3 thresholds (composed of a score and a minimum number of mapped voxels in the atlas), one for high confidence mapping, one for medium confidence and the last one for simple hypothesis

Once the threshold are chosen you can output a summary of the results for the real cells by using the following command
```R
#defining the thresholds as you chose them
#high confidence minimum score
h_thres = 2
#high confidence minimum numbers of mapped voxels
h_thres_num = 11

#medium confidence minimum score
m_thres = 1
#medium confidence minimum numbers of mapped voxels
m_thres_num = 6

#hypothesis confidence minimum score
l_thres = 0
#hypothesis confidence minimum numbers of mapped voxels
l_thres_num = 6


#Summarizing the results
results <- summary_results(example_results_scores,h_thres,m_thres,l_thres,h_thres_num,m_thres_num,l_thres_num)
```


####Visualizing the results in 3D
With the coordinates of the points in the atlas in 3D in CSV separated by commas as shown in the example file `example_3D_coordinates_atlas` you can directly see the results in [bioWeb3D](http://www.ebi.ac.uk/~jbpettit/map_viewer). Just start by importing the CSV file as a dataset file. Then you can use the follwing command to create the mapping file in R :
```R
scaled_results <- scale_res(example_results_scores,h_thres,m_thres,l_thres,h_thres_num,m_thres_num,l_thres_num,rownames(rna_seq))
write.table(file="example_scaled.csv",sep=",",scaled_results,row.names=FALSE,quote=FALSE)
```

Then you simply add this newly created file as "cluster data".

screenshot

![screenshot](http://www.ebi.ac.uk/~jbpettit/map_viewer/screen.png)


The full workflow on the example dataset provided is listed below:
```R
rna_seq <- read.table("example_dataset/example_data_RNA_seq.csv",header=TRUE,sep="\t")
atlas <- read.table("example_dataset/example_data_atlas.csv",header=TRUE,sep="\t")
coordinates <- read.table("example_dataset/example_3D_coordinates_atlas",header=FALSE,sep=",")

#Loading the analysis functions
source("spatial_mapping.R")
specificity_matrix <- specificity_scores(rna_seq)


number_of_points_in_altas <- 1000
atlas_specificity_score <- 1/apply(atlas,2,function(x){length(x[x>0])/number_of_points_in_altas})
#If you decide not to, simply set it as a vector of 0s
atlas_specificity_score <- as.numeric(vector(length=number_of_points_in_altas))

#Iterate over the sequenced cells
example_results_scores = sapply(seq_along(rna_seq[,1]),function(cell_num) {

	print(paste("Mapping cell",cell_num))
	cell = rna_seq[cell_num,]	

	specificity_score <- specificity_matrix[cell_num,]

	#You need a binary expression vector for this cell, the count threshold is for you to decide given your RNA-seq data
	example_count_threshold <- 10
	binary_expression_cell <- sapply(rna_seq[1,],function(x){ifelse(x>example_count_threshold,1,0)})

	#Launch the mapping against the atlas for the cell
	mapping_result_cell <- spatial_map_scoring(specificity_score,binary_expression_cell,atlas_specificity_score,atlas)
	
	mapping_result_cell
})

generate_simulated_data(specificity_matrix,100,"example_simulated_data/")


cell_counter = 1
mapping_simulated_list = list()
for(dataset in 1:10) {
	randomCells = read.table(paste("example_simulated_data/",dataset,"_bin.data",sep=""))
	randomRatios = read.table(paste("example_simulated_data/",dataset,".data",sep=""))

	#map each cell
	for(cell_i in 1:100) {
		print(cell_i)
		mapping_simulated_list[[cell_counter]] <- spatial_map_scoring(randomRatios[cell_i,],randomCells[cell_i,],atlas_specificity_score,atlas)
		cell_counter = cell_counter +1
	}
}

#convert to a dataframe for easier use
mapping_simulated = data.frame(mapping_simulated_list)
colnames(mapping_simulated)=paste("cell_",1:length(mapping_simulated_list),sep="")

#creating the different thresholds
threshold_score <- seq(-2,2,1)
threshold_number <- seq(1,11,5)

number_cells_above = matrix(nrow=length(threshold_score),ncol=length(threshold_number))
for(j in 1:length(threshold_number)){
	print(j)
	number_cells_above[,j] = unlist(lapply(seq_along(threshold_score),function(i){
		sup = apply(mapping_simulated,2,function(x){
			length(x[x>threshold_score[i]])
		})
		length(sup[sup>threshold_number[j]])
	}))
}
number_cells_above = data.frame(number_cells_above)
rownames(number_cells_above) = as.character(threshold_score)
colnames(number_cells_above) = as.character(threshold_number)

#Compute the proportions
proportion_cells_thresholds = number_cells_above/(dim(mapping_simulated)[2])
proportion_cells_thresholds


#defining the thresholds as you chose them
#high confidence minimum score
h_thres = 2
#high confidence minimum numbers of mapped voxels
h_thres_num = 11

#medium confidence minimum score
m_thres = 1
#medium confidence minimum numbers of mapped voxels
m_thres_num = 6

#hypothesis confidence minimum score
l_thres = 0
#hypothesis confidence minimum numbers of mapped voxels
l_thres_num = 6


#Summarizing the results
results <- summary_results(example_results_scores,h_thres,m_thres,l_thres,h_thres_num,m_thres_num,l_thres_num)

scaled_results <- scale_res(example_results_scores,h_thres,m_thres,l_thres,h_thres_num,m_thres_num,l_thres_num,rownames(rna_seq))
write.table(file="example_scaled.csv",sep=",",scaled_results,row.names=FALSE,quote=FALSE)
```








