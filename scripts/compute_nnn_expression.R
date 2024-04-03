# Extract expression values and next nearest neighbours - Resolve data ----

library(FNN)
library(data.table)

# read file
transcript_coord <- read.table("~/your_mc.Resolve_mols_coords.csv" , sep= ",", header=TRUE)

# order via gene names
transcript_coord_sorted <- transcript_coord[order(transcript_coord$molecule),]
transcript_coord_sorted_matrix <- data.matrix(atranscript_coord_sorted)

# count number of called molecules per gene and extract expression values
# return lines when gene name changes, add value of zero and length of data matrix
# substract value from previous value in data.frame which returns the total gene expression count per transcript
lines <- which(c(FALSE, tail(transcript_coord_sorted$molecule,-1) != head(transcript_coord_sorted$molecule,-1)))
dim <- as.matrix(dim(atranscript_coord_sorted))
lines <- append(lines,(dim[1,1]+1))#add length of matrix +1 to bottom row to count also last transcript
lines <- append(1,lines)#add one as top row to count the first transcript
exp <- as.data.table(lines)
exp[, expression := lines - shift(lines, fill = first(lines))]

# extract gene names and combine with expression values
gene_names <- unique(transcript_coord_sorted$transcript_coord_sorted.gene)
dim_expression <- as.matrix(dim(exp))
new <- cbind(gene_names, exp[2:dim_expression[1],2])

# create vector with median next nearest distance per transcript
# compute nnn with k=1 among the same transcripts and return median to empty vector
nnn_dis=c()
for (i in seq(from=1, to=length(lines), by=1)) {
genes <- transcript_coord_sorted_matrix[lines[i]:(lines[i+1]-1),1:2] # lines given return when value for gene name changes, thus a value of one has to be subtracted
distance <- knn.dist(genes, k=1, algorithm=c("kd_tree"))
nnn_dis <- c(nnn_dis,median(distance))
}

# combine with gene name and expression values and save file
combine <- cbind(new,nnn_dis)
write.table(combine,"~/your_mc_exp_nnn.txt")


# Extract expression values and next nearest neighbours - Vizgen data ----
# for vizgen data the follow adaptations are required
library(FNN)
library(data.table)

# read file
transcript_coord <- read.table("~/your_viz_detected_transcripts.csv" , sep= ",", header=TRUE)

# order via gene names
transcript_coord_sorted <- transcript_coord[order(transcript_coord$gene),]
transcript_coord_sorted <- data.frame(transcript_coord_sorted$global_x,transcript_coord_sorted$global_y,transcript_coord_sorted$gene)
transcript_coord_sorted_matrix <- data.matrix(transcript_coord_sorted)

# count number of called molecules per gene and extract expression values
# return lines when gene name changes, add value of zero and length of data matrix
# substract value from previous value in data.frame which returns the total gene expression count per transcript
lines <- which(c(FALSE, tail(transcript_coord_sorted$transcript_coord_sorted.gene,-1) != head(transcript_coord_sorted$transcript_coord_sorted.gene,-1)))
dim <- as.matrix(dim(transcript_coord_sorted))
lines <- append(lines,(dim[1,1]+1))#add length of matrix +1 to bottom row to count also last transcript
lines <- append(1,lines)#add one as top row to count the first transcript
exp <- as.data.table(lines)
exp[, expression := lines - shift(lines, fill = first(lines))]

# extract gene names and combine with expression values
gene_names <- unique(transcript_coord_sorted$transcript_coord_sorted.gene)
dim_expression <- as.matrix(dim(exp))
new <- cbind(gene_names, exp[2:dim_expression[1],2])

# create vector with median next nearest distance per transcript
# compute nnn with k=1 among the same transcripts and return median to empty vector
nnn_dis=c()
for (i in seq(from=1, to=length(lines), by=1)) {
genes <- transcript_coord_sorted_matrix[lines[i]:(lines[i+1]-1),1:2] # lines given return when value for gene name changes, thus a value of one has to be subtracted
distance <- knn.dist(genes, k=1, algorithm=c("kd_tree"))
nnn_dis <- c(nnn_dis,median(distance))
}


#combine gene names, expression and nnn distance
#write file
combine <- cbind(new,nnn_dis)
write.table(combine,"~/_merscope_exp_nnn.txt")


# Extract expression values and next nearest neighbours - Xenium data ---- 
library(FNN)
library(data.table)

# read file
transcript_coord  <- data.table::fread("~/your_xenium_transcripts.csv.gz", sep= ",")

# order via gene names
transcript_coord_sorted <- transcript_coord[order(transcript_coord$feature_name),]
transcript_coord_sorted <- data.frame(transcript_coord_sorted$feature_name,transcript_coord_sorted$x_location,transcript_coord_sorted$y_location)
transcript_coord_sorted_matrix <- data.matrix(transcript_coord_sorted)

# count number of called molecules per gene and extract expression values
# return lines when gene name changes, add value of zero and length of data matrix
# substract value from previous value in data.frame which returns the total gene expression count per transcript
lines <- which(c(FALSE, tail(transcript_coord_sorted$transcript_coord_sorted.feature_name,-1) != head(transcript_coord_sorted$transcript_coord_sorted.feature_name,-1)))
dim <- as.matrix(dim(transcript_coord_sorted))
lines <- append(lines,(dim[1,1]+1))#add length of matrix +1 to bottom row to count also last transcript
lines <- append(1,lines)#add one as top row to count the first transcript
exp <- as.data.table(lines)
exp[, expression := lines - shift(lines, fill = first(lines))]



# extract gene names and combine with expression values
gene_names <- unique(transcript_coord_sorted$transcript_coord_sorted.feature_name)
dim_expression <- as.matrix(dim(exp))
new <- cbind(gene_names, exp[2:dim_expression[1],2])


# create vector with median next nearest distance per transcript
# compute nnn with k=1 among the same transcripts and return median to empty vector
nnn_dis=c()
for (i in seq(from=1, to=length(lines), by=1)) {
genes <- transcript_coord_sorted_matrix[lines[i]:(lines[i+1]-1),2:3] # lines given return when value for gene name changes, thus a value of one has to be subtracted
distance <- knn.dist(genes, k=1, algorithm=c("kd_tree"))
nnn_dis <- c(nnn_dis,median(distance))
}

#combine gene names, expression and nnn distance
#write file
combine <- cbind(new,nnn_dis)
write.table(combine,"~/_xenium_exp_nnn.txt")
