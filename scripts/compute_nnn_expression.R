
# extract expression values and next nearest neighbours from resolve data

library(FNN)
library(data.table)

# read file
a <- read.table("~/your_mc.Resolve_mols_coords.csv" , sep= ",", header=TRUE)

# order via gene names
as <- a[order(a$molecule),]
ams <- data.matrix(as)

# count number of called molecules per gene and extract expression values
# return lines when gene name changes, add value of zero and length of data matrix
# substract value from previous value in data.frame which returns the total gene expression count per transcript
jump <- which(c(FALSE, tail(as$molecule,-1) != head(as$molecule,-1)))
xy <- as.matrix(dim(as))
jump <- append(jump,xy[1,1])
jump <- append(0,jump)
jump3 <- as.data.table(jump)
jump3[, expression := jump - shift(jump, fill = first(jump))]

# extract gene names and combine with expression values
genes <- unique(as$molecule)
yz <- as.matrix(dim(jump3))
new <- cbind(genes, jump3[2:yz[1],2])

# create vector with median next nearest distance per transcript
nnn_dis=c()
for (i in seq(from=1, to=length(jump), by=1)) {
gene2 <- ams[jump[i]:jump[i+1],1:2]
dist2 <- knn.dist(gene2, k=1, algorithm=c("kd_tree"))
nnn_dis <- c(nnn_dis,median(dist2))
}

# combine with gene name and expression values and save file
comp <- cbind(new,nnn_dis)
write.table(comp,"~/your_mc_exp_nnn.txt")


# for vizgen data the follow adaptations are required

library(FNN)
library(data.table)

# read file
a <- read.table("~/your_viz_detected_transcripts.csv" , sep= ",", header=TRUE)

# order via gene names
as <- a[order(a$gene),]
as <- data.frame(as$global_x,as$global_y,as$gene)

ams <- data.matrix(as)

# count number of called molecules per gene and extract expression values
jump <- which(c(FALSE, tail(as$as.gene,-1) != head(as$as.gene,-1)))
xy <- as.matrix(dim(as))
jump <- append(jump,xy[1,1])
jump <- append(0,jump)
jump3 <- as.data.table(jump)
jump3[, expression := jump - shift(jump, fill = first(jump))]

# extract gene names and combine with expression values
genes <- unique(as$as.gene)
yz <- as.matrix(dim(jump3))
new <- cbind(genes, jump3[2:yz[1],2])

# create vector with median next nearest distance per transcript
nnn_dis=c()
for (i in seq(from=1, to=length(jump), by=1)) {
gene2 <- ams[jump[i]:jump[i+1],1:2]
dist2 <- knn.dist(gene2, k=1, algorithm=c("kd_tree"))
nnn_dis <- c(nnn_dis,median(dist2))
}

# combine with gene name and expression values and save file
comp <- cbind(new,nnn_dis)
write.table(comp,"~/your_viz_exp_nnn.txt")


# for Xenium 
library(FNN)
library(data.table)

# read file
a <- data.table::fread("~/your_xenium_transcripts.csv.gz", sep= ",")

# order via gene names
as <- a[order(a$feature_name),]
ass <- data.frame(as$feature_name,as$x_location,as$y_location,as$qv)
ams <- data.matrix(ass)

# count number of called molecules per gene and extract expression values
jump <- which(c(FALSE, tail(ass$as.feature_name,-1) != head(ass$as.feature_name,-1)))
xy <- as.matrix(dim(ass))
jump <- append(0,jump)
jump <- append(jump,xy[1,1])
jump3 <- as.data.table(jump)
jump3[, expression := jump - shift(jump, fill = first(jump))]
genes <- unique(ass$as.feature_name)
yz <- as.matrix(dim(jump3))
new <- cbind(genes, jump3[2:yz[1],2])

nnn_dis=c()

for (i in seq(from=1, to=length(jump), by=1)) {
gene2 <- ams[jump[i]:jump[i+1],2:3]
dist2 <- knn.dist(gene2, k=1, algorithm=c("kd_tree"))
nnn_dis <- c(nnn_dis,median(dist2))
}

comp <- cbind(new,nnn_dis)
write.table(comp,"~/your_xenium_exp_nnn.txt")

