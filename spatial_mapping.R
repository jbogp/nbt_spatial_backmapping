################################################################################################
#Computes the specificity scores for the RNA-seq dataset
#
#Parameters:
# - rna_seq_counts : Matrix containing the count number of read for every sequenced cells (rows) and every considered gene (columns)
#
################################################################################################
specificity_scores <- function(rna_seq_counts) {
	#Divide each row of the data by the mean of each column
	t(apply(rna_seq_counts,1,function(x){x/colMeans(rna_seq_counts)}))
}

################################################################################################
#Gets results together for simualted data (same as together with a subset structure)
#Use this function to put together the results from each simulated cell in an easy to use data structure
#
#Parameters: 
# - folder : the folder to look into e.g "results/"
# - num_set : number of sets of simulated cells
# - num : number of files to look for
# - num_units_per_cell : number of rows for each file (number of units in the reference dataset)
################################################################################################
simul_together <- function(folder, num_set, num, num_units_per_cell) {
	result <- matrix(ncol=(num*num_set), nrow=num_units_per_cell)
	j=1
	for(k in 1:num_set) {
		for(i in 1:num) {
			print(i)
			if(file.exists(paste(folder,"set_",k,"_cell_",i,sep=""))){
					result[,j] = as.numeric(read.table(paste(folder,"set_",k,"_cell_",i,sep=""),header=TRUE)[,1])
					j=j+1
			}
		}
	}
	result
}


################################################################################################
#Generates the mapping summary object from res
#
#Parameters:
# - res : matrix containing the sequenced cell by columns, each row is the mapping score for every unit in the reference dataset
# - h_thres: high confidence threshold
# - m_thres: medium confidence threshold
# - l_thres: hypothesis confidence threshold
# - h_thres_res: high confidence minumum number of units
# - m_thres_res: medium confidence minumum number of units
# - l_thres_res: hypothesis confidence minumum number of units
################################################################################################
summary_results <- function(res,h_thres,m_thres,l_thres,h_thres_num,m_thres_num,l_thres_num){


	results = apply(res,2,function(x){
		if(length(x[x>=h_thres])>=h_thres_num){
			c("high",length(x[x>=h_thres]))
		}
		else if(length(x[x>=m_thres])>=m_thres_num){
			c("medium",length(x[x>=m_thres]))
		}
		else if(length(x[x>=l_thres])>=l_thres_num){
			c("hypothesis",length(x[x>=l_thres]))
		}
		else{
			c("not_mapped",0)
		}
	})

	results = t(data.frame(results))

	results

}


################################################################################################
#Scale the results to create a file for 3D visualization in bioweb3D
#
#Parameters:
# - res : matrix containing the sequenced cell by columns, each row is the mapping score for every unit in the reference dataset
# - h_thres: high confidence threshold
# - m_thres: medium confidence threshold
# - l_thres: hypothesis confidence threshold
# - h_thres_res: high confidence minumum number of units
# - m_thres_res: medium confidence minumum number of units
# - l_thres_res: hypothesis confidence minumum number of units
# - names_vec: vector with the names of the cells
################################################################################################
scale_res <- function(res,h_tres,m_tres,l_tres,h_tres_num,m_tres_num,l_tres_num,names_vec) {
	### Iterate on every sequenced cell
	sca = apply(res,2,function(x) {
		vec = unlist(x)

		h_above = which(vec>=h_thres)
		h_below = which(vec<h_thres)

		m_above = which(vec>=m_thres)
		m_below = which(vec<m_thres)

		l_above = which(vec>=l_thres)
		l_below = which(vec<l_thres)

		#### Depending on the confidence we have in the mapping, we apply a different scaling, to have readable results in every case
		if(length(h_above)>= h_thres_num){
			if(max(vec[vec>=h_thres]) == min(vec[vec>=h_thres])){
				vec[h_above] = 1
			}
			else{
				vec[h_above] = ceiling((1-(vec[vec>=h_thres]-min(vec[vec>=h_thres]))/(max(vec[vec>=h_thres])-min(vec[vec>=h_thres])))*4)+1
			}
			vec[h_below]= 100
		}
		else if(length(m_above)>= m_thres_num){
			if(max(vec[vec>=m_thres]) == min(vec[vec>=m_thres])){
				vec[m_above] = 1
			}
			else{
				vec[m_above] = ceiling((1-(vec[vec>=m_thres]-min(vec[vec>=m_thres]))/(max(vec[vec>=m_thres])-min(vec[vec>=m_thres])))*4)+1
			}
			vec[m_below]= 100
		}
		else if(length(l_above)>= l_thres_num){
			if(max(vec[vec>=l_thres]) == min(vec[vec>=l_thres])){
				vec[l_above] = 1
			}
			else{
				vec[l_above] = ceiling((1-(vec[vec>=l_thres]-min(vec[vec>=l_thres]))/(max(vec[vec>=l_thres])-min(vec[vec>=l_thres])))*4)+1
			}
			vec[l_below]= 100
		}
		else{
			vec[] = 100
		}
		
		vec

	
	})
 	colnames(sca) = names_vec

	sca


}

################################################################################################
#Generate simulated sequenced cell from the specificity score of the real sequenced cells
#
#Parameters:
# - specificity_matrix : contains the specificity score for each cell (rows) and each gene (columns)
# - simulated_cells_per_cell : Number of simulated cells to create for each sequenced cell
# - output_folder : folder (has to exist) in whih to create the resulting files
#
################################################################################################
generate_simulated_data <- function(specificity_matrix,simulated_cells_per_cell,output_folder) {

	#itrating over the sequenced cells
	for(h  in 1:length(rownames(specificity_matrix))) {

		randomCells <- matrix(ncol=100, nrow=length(colnames(specificity_matrix)))

		for(i in 1:simulated_cells_per_cell) {
			tempVec = NULL
			tempVec = as.numeric(sample(specificity_matrix[h,]))
			randomCells[,i] = tempVec
		}
		randomCells = data.frame(randomCells)
		rownames(randomCells) = colnames(specificity_matrix)
		randomCells = t(randomCells)

		#generating the binary expression too
		randomCells.bin = randomCells
		randomCells.bin[randomCells.bin>0] = 1

		#outputing the files
		write.table(randomCells, file=paste(output_folder,h,".data",sep=""), quote=F)
		write.table(randomCells.bin, file=paste(output_folder,h,"_bin.data",sep=""), quote=F)
	}
}



################################################################################################
#Main scoring function. Computes the scores of one sequenced cell's specificity scores against every unit in the reference atlas
#
#Parameters:
# - scores: vector of size M (number of genes considered) containing the specificity scores for one sequenced cell 
# - vec: binarized vector of size M of the expression from the sequenced cell for the considered genes
# - scores_atlas : vector of size M containing the prortion of units in the reference Atlas expressing each gene (or 0s if you watn 1 way penalization)
# - bin : matrix of binary expression reference data (atlas) (columns are genes, rows are volume or surface binary units)
#
################################################################################################
spatial_map_scoring <- function(scores,vec,scores_atlas,bin) {
	ret = NULL

	vec = as.numeric(vec)
	scores = as.numeric(scores)
	#Add the sequenced cell binary expression to each row of the reference atlas
	sumMat = apply(bin,1,function(x){x+vec})
	#Apply to each row of sumMat
	scoreMat = sapply(1:nrow(t(sumMat)),function(j) {
		unlist(lapply(seq_along(sumMat[,j]),function(i) {
			#mismatch, we will penalise			
			if(sumMat[i,j]==1){
				#this means the gene is not expressed in the rna-seq
				if(scores[i] == 0){
					-scores_atlas[i]/(1+scores_atlas[i])
				}
				else {
					-scores[i]/(1+scores[i])
				}
			}
			else if(sumMat[i,j] ==2){
				scores[i]/(1+scores[i])
			}
			else {
				0
			}
		}))

	})
	

	colSums(scoreMat)

}
