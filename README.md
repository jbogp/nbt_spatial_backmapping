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
mapping_result_cell_1 <- spatial_map_scoring(specificity_score_cell_1,binary_expression_cell_1,atlas_specificity_score)
```

Of course this is just for the first sequenced cell, to get the mapping results for the other cells, you simply have to iterate over this. Note that each mapping computation is independent from the others, meaning that you can easily parallelize these computations depending on the system you are using (multi-core plugins in R or different jobs if you have access to a cluster). If your expression atlas is quite large, such parallelization will become necessary as the scoring process for each cell will take time. 

In this case you will probably save each `mapping_result_cell_n` on your system. We provide a function to get such results together in an easy to use object. 

#Getting results together
Assuming I have saved the mapping scores for each of my 10 sequenced cells in the folder `mapping_results_example` with the file names:
 - `mapping_result_cell_1`
 - `mapping_result_cell_2`
 - ...
 - `mapping_result_cell_10`

in R I can simply get all the results together by doing
```R
example_results_scores = together("mapping_results_example/mapping_result_cell_",10)
```

####Simulating data to find the confidence thresholds
Now you have the object `example_results_scores` which contains the score for each sequenced cell against every reference voxel in the atlas. You still need to find the thresholds above which you can assume that a voxel is a match for a cell. To do this, you will generate a null distribution from your data by simulating random sequenced cells. First create a directory to store the generated simulated data, for instance `/example_simluated_data`.

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
for(dataset in 1:10) {
	randomCells = read.table(paste("example_simulated_data/",dataset,"_bin.data",sep=""))
	randomRatios = read.table(paste("example_simulated_data/",dataset,".data",sep=""))

	mapping_simulated_list = list()

	#map each cell
	for(cell_i in 1:100) {
		mapping_simulated_list[[cell_counter]] <- spatial_map_scoring(randomRatios[cell_i,],randomCells[cell_i,],atlas_specificity_score)
		cell_counter = cell_counter +1
	}
}
```

#####Choosing the threshold
From the previous step you will able to plot the proportion of simulated cells scoring higher than different thresholds, which will allow you you choose sensible threshold(s) for your dataset. The following functions assume you chose 3 thresholds (composed of a score and a minimum number of mapped voxels in the atlas), one for high confidence mapping, one for medium confidence and the last one for simple hypothesis

Once the threshold are chose you can output a summary of the results for the real cells by using the following command
```R
#defining the thresholds as you chose them
#high confidence minimum score
h_tres = 2
#high confidence minimum numbers of mapped voxels
h_tres_num = 30

#medium confidence minimum score
m_tres = 2
#medium confidence minimum numbers of mapped voxels
m_tres_num = 30

#hypothesis confidence minimum score
l_tres = 2
#hypothesis confidence minimum numbers of mapped voxels
l_tres_num = 30


#Summarizing the results
results <- summary_results(example_results_scores,h_tres,m_tres,l_tres,h_tres_num,m_tres_num,l_tres_num)
```


####Visualizing the results in 3D
With the coordinates of the points in the atlas in 3D in CSV separated by commas as shown in the example file `example_3D_coordinates_atlas` you can directly see the results in (bioWeb3D)[http://www.ebi.ac.uk/~jbpettit/map_viewer]. Just start by importing the CSV file as a dataset file. Then you can use the follwing command to create the mapping file in R :
```R
scaled_results <- scale_res(example_results_scores,h_tres,m_tres,l_tres,h_tres_num,m_tres_num,l_tres_num)
write.table(file="example_scaled.csv",sep=",",scaled_results,row.names=FALSE)
```

Then you simply add this newly created file as "cluster data"



