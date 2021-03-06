################################################################################
# Running PATHIFIER (Drier et al., 2013)
# Author: Miguel Angel Garcia-Campos - Github: https://github.com/AngelCampos
################################################################################

# Install/load package ##############
source("http://bioconductor.org/biocLite.R")
if (!require("pathifier")){biocLite("pathifier", ask=FALSE); library(pathifier)}

# Load expression data for PATHIFIER
exp.matrix <- read.delim(file =file.choose(), as.is = T, row.names = 1)

# Load Genesets annotation 
gene_sets <- as.matrix(read.delim(file = file.choose(), header = F, sep = "\t",
                                  as.is = T))


#  Generate a list that contains genes in genesets
gs <- list()
for (i in 1:nrow(gene_sets)){
  a <- as.vector(gene_sets[i,3:ncol(gene_sets)])
  a <- na.omit(a)
  a <- a[a != ""]
  a <- matrix(a, ncol = 1)
  gs[[length(gs)+1]] <- a
  rm(a,i)
}

# Generate a list that contains the names of the genesets used
pathwaynames <- as.list(gene_sets[,1])

# Generate a list that contains the previos two lists: genesets and their names
PATHWAYS <- list(); PATHWAYS$gs <- gs; PATHWAYS$pathwaynames <- pathwaynames

# Prepare data and parameters ##################################################
# Extract information from binary phenotypes. 1 = Normal, 0 = Tumor
normals <- as.vector(as.logical(exp.matrix[1,]))
exp.matrix <- as.matrix(exp.matrix[-1, ])

# Calculate MIN_STD
N.exp.matrix <- exp.matrix[,as.logical(normals)]
rsd <- apply(N.exp.matrix, 1, sd)
min_std <- quantile(rsd, 0.25)

# Calculate MIN_EXP 
min_exp <- quantile(as.vector(exp.matrix), 0.1) # Percentile 10 of data

# Filter low value genes. At least 10% of samples with values over min_exp
# Set expression levels < MIN_EXP to MIN_EXP
over <- apply(exp.matrix, 1, function(x) x > min_exp)
G.over <- apply(over, 2, mean)
G.over <- names(G.over)[G.over > 0.1]
exp.matrix <- exp.matrix[G.over,]
exp.matrix[exp.matrix < min_exp] <- min_exp

# Set maximum 5000 genes with more variance
V <- names(sort(apply(exp.matrix, 1, var), decreasing = T))[1:5000]
V <- V[!is.na(V)]
exp.matrix <- exp.matrix[V,]
genes <- rownames(exp.matrix) # Checking genes
allgenes <- as.vector(rownames(exp.matrix))

# Generate a list that contains previous data: gene expression, normal status,
# and name of genes
DATASET <- list(); DATASET$allgenes <- allgenes; DATASET$normals <- normals
DATASET$data <- exp.matrix

# Run Pathifier
PDS <- quantify_pathways_deregulation(DATASET$data, 
                                      DATASET$allgenes,
                                      PATHWAYS$gs,
                                      PATHWAYS$pathwaynames,
                                      DATASET$normals, 
                                      maximize_stability = T,
                                      attempts = 10,
                                      logfile="logfile.txt",
                                      min_std = min_std,
                                      min_exp = min_exp)
# Remove unnecesary data
rm(gene_sets, exp.matrix, allgenes, DATASET, PATHWAYS, rsd, V, over, G.over, 
   N.exp.matrix, gs, genes, min_exp, min_std, pathwaynames)

# Save image into working directory
save.image("PDS.RData")

###############################################################################################################
###############################################################################
## Heatmap for Pathifier results in R
### Author: Angel Garcia-Campos https://github.com/AngelCampos
#### Base by wonderful: Sebastian Raschka https://github.com/rasbt
###############################################################################

###############################################################################
### Installing and/or loading required packages
###############################################################################

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

###############################################################################
## Load Pathifier results and turn into a matrix
###############################################################################

load ("PDS.RData")
PDSmatrix <- t(mapply(FUN = c, PDS$scores))

###############################################################################
## Creating Custom Palette
###############################################################################

# creates a own color palette passing from blue, green yellow to dark red
my_palette <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(n = 1000))

###############################################################################
## Clustering Methods
###############################################################################

# If you want to change the default clustering method (complete linkage method
# with Euclidean distance measure), this can be done as follows: For non-square
# matrix, we can define the distance and cluster based on our matrix data by

row.distance = dist(PDSmatrix, method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(PDSmatrix), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")

# Arguments for the dist() function are: euclidean (default), maximum, canberra,
# binary, minkowski, manhattan
# And arguments for hclust(): complete (default), single, average, mcquitty,
# median, centroid, ward.D2

###############################################################################
## Assign Column labels (Optional)
###############################################################################

colLabels <- as.character(normals)
colLabels[colLabels == "TRUE"] <- "#377EB8"
colLabels[colLabels == "FALSE"] <- "#E41A1C"

###############################################################################
## Plotting the Heatmap!! (where all colorful things happen...)
###############################################################################

png("heatmap.png", # Name of png file
    width = 6 * 500,      # Easier scaling 6*500 = 3000 pixels
    height = 6 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # font size

heatmap.2(PDSmatrix,
          main = "PDS-Heatmap",   # heat map title
          density.info = "none",  # turns off density plot inside color legend
          trace = "none",         # turns off trace lines inside the heat map
          margins = c(10,21),     # widens margins around plot
          col = my_palette,       # use on color palette defined earlier
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize = 0.8,          # size of color key
          #Additional Options
          ## Color labeled columns
          ColSideColors = colLabels
)
## Legend for ColumnSide color labeling
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Normals", "Tumors"), # category labels
       col = c("dodgerblue", "firebrick1"),  # color key
       lty= 1,          # line style
       lwd = 5, unit    # line width
)
dev.off() # close the PNG device

###############################################################################################################
################################################################################
# Plotting Principal Cuves for Pathifier Results
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################
# Installing and/or loading required packages
if (!require("plotly")) {
  install.packages("plotly", dependencies = TRUE)
  library(plotly)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

################################################################################
# Load Pathifier results and turn into a matrix
load ("PDS.RData")
PDSmatrix <- t(mapply(FUN = c, PDS$scores))

################################################################################
# Plotting Principal Curves
curve.data.PDS <- function(pw){
  x <- PDS$curves[[pw]][which(normals),] # Healthy samples
  y <- PDS$curves[[pw]][which(!normals),] # Tumoral samples
  v <- PDS$z[[pw]][which(normals),]
  z <- PDS$z[[pw]][which(!normals),]
  p <- plot_ly(x = x[,1], y = x[,2],z = x[,3], type = "scatter3d", 
               name = "Healthy", marker = list(size = 4))
  add_trace(p, x = y[,1], y = y[,2], z = y[,3], type = "scatter3d", 
            name = "Tumors",marker = list(size = 4)) %>%
    add_trace(p, x = v[,1], y = v[,2], z = v[,3], type = "scatter3d", 
              name = "Data-Healthy", marker = list(size = 2)) %>%
    add_trace(p, x = z[,1], y = z[,2], z = z[,3], type = "scatter3d", 
              name = "Data-Tumoral", marker = list(size = 2)) %>%
    layout(title = paste("Principal curve of", rownames(PDSmatrix)[pw]))
}

curve.PDS <- function(pw){
  x <- PDS$curves[[pw]][which(normals),] # Healthy samples
  y <- PDS$curves[[pw]][which(!normals),] # Tumoral samples
  p <- plot_ly(x = x[,1], y = x[,2],z = x[,3], type = "scatter3d", 
               name = "Healthy", marker = list(size = 4))
  add_trace(p, x = y[,1], y = y[,2], z = y[,3], type = "scatter3d", 
            name = "Tumors",marker = list(size = 4)) %>%
    layout(title = paste("Principal curve of", rownames(PDSmatrix)[pw]))
}

################################################################################
# Examples:
# Plotting the samples on their principal curves of pathway # 1
curve.data.PDS(1)
# Plotting the data of the samples and their projected points to the principal
# curve of pathway # 2
curve.PDS(2)




