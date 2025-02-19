#Working copies of the original datasets.
xCIS_0h <- CIS_0h
xCIS_12h <- CIS_12h
xCIS_24h <- CIS_24h
xCIS_48h <- CIS_48h

#Normalize the gene expression counts without log transformation
xCIS_0h$mat <- MERINGUE::normalizeCounts(xCIS_0h$gexp, log=FALSE)
xCIS_12h$mat <- MERINGUE::normalizeCounts(xCIS_12h$gexp, log=FALSE)
xCIS_24h$mat <- MERINGUE::normalizeCounts(xCIS_24h$gexp, log=FALSE)
xCIS_48h$mat <- MERINGUE::normalizeCounts(xCIS_48h$gexp, log=FALSE)

#Identifies genes present in all four time points
genes.shared <- Reduce(intersect, list(rownames(xCIS_0h$mat), 
                                       rownames(xCIS_12h$mat), 
                                       rownames(xCIS_24h$mat), 
                                       rownames(xCIS_48h$mat)))

#Creates combined gene expression matrix for the shared genes across all time points.
CIS_gexp <- data.frame(cbind(
  as.matrix(xCIS_0h$mat[genes.shared,]),
  as.matrix(xCIS_12h$mat[genes.shared,]),
  as.matrix(xCIS_24h$mat[genes.shared,]),
  as.matrix(xCIS_48h$mat[genes.shared,])
))

#Creates a metadata dataframe indicating the time points (0h, 12h, 24h, 48h) 
#for each sample.
meta <- data.frame(time = c(
  rep(0, ncol(xCIS_0h$mat)),
  rep(12, ncol(xCIS_12h$mat)), ## unit of days
  rep(24, ncol(xCIS_24h$mat)),
  rep(48, ncol(xCIS_48h$mat))
))  

#Assigns row names corresponding to the column names of CIS_gexp
rownames(meta) <- colnames(CIS_gexp)
table(meta)

#Initializes an annotation vector (annot) where all samples are labeled as 'other'.
annot <- rep('other', nrow(meta))
#Assigns row names from meta to annot
names(annot) <- rownames(meta)
cortex <- readRDS('/Users/deev/Documents/R/Kidney App/AKI-CIS-irl-ctrl_Cortex_spots.rds')
names(cortex) <- gsub('-', '.',names(cortex))
annot[rownames(meta) %in% names(cortex)] <- 'cortex'
interface <- readRDS('/Users/deev/Documents/R/Kidney App/AKI-CIS-irl-ctrl_Interface_spots.rds')
names(interface) <- gsub('-', '.',names(interface))
annot[rownames(meta) %in% names(interface)] <- 'interface'
medulla <- readRDS('/Users/deev/Documents/R/Kidney App/AKI-CIS-irl-ctrl_Medulla_spots.rds')
names(medulla) <- gsub('-', '.',names(medulla))
annot[rownames(meta) %in% names(medulla)] <- 'medulla'
table(annot)
meta$annot <- annot
head(meta)
colnames(meta) <- c('time', 'compartment')

#Retrieves the range of X and Y coordinates for the 0h time point.
#Sets a translation factor (t = 2500) to separate time points visually.
range(xCIS_0h$pos[,1])
range(xCIS_0h$pos[,2])
t <- 2500

#Adds an increasing offset to the X-coordinates of later time points.
#Ensures that different time points are displayed side-by-side rather than 
#overlapping in spatial visualizations
xCIS_12h$pos[,1] <- xCIS_12h$pos[,1] + t
range(xCIS_12h$pos[,1])
range(xCIS_12h$pos[,2])

xCIS_24h$pos[,1] <- xCIS_24h$pos[,1] + 2*t
range(xCIS_24h$pos[,1])
range(xCIS_24h$pos[,2])

xCIS_48h$pos[,1] <- xCIS_48h$pos[,1] + 3*t

CIS_gexp.sub <- t(CIS_gexp[medulla_up_genes,rownames(meta)])
dim(CIS_gexp.sub)
dim(meta)

#Combines all time points’ spatial positions into a single dataframe
pos.all <- rbind(xCIS_0h$pos, xCIS_12h$pos, xCIS_24h$pos, xCIS_48h$pos)
#Replaces hyphens in row names
rownames(pos.all) <- gsub('-', '.', rownames(pos.all))
head(pos.all)
plot(pos.all)
## factors (time and compartment)
head(meta)

## make polygon hack
#creates a circle around a given (x, y) position in the pos matrix
generate_circle <- function(pos, radius = 10, n_points = 10) {
  #Converts pos into a matrix
  pos <- as.matrix(pos)
  x = pos[1]
  y = pos[2]
  theta <- seq(0, 2 * pi, length.out = n_points + 1)  # Create angles for the circle
  #Generates x and y coordinates for n_points around the circle
  return(data.frame(
    X = x + radius * cos(theta),
    Y = y + radius * sin(theta)
  ))
}

#Applies generate_circle() to each row (cell) in pos.all, 
#creating a list of circular polygons.
pos.poly <- lapply(rownames(pos.all), function(i) { generate_circle(pos.all[i,]) })
names(pos.poly) <- rownames(pos.all)
plot(pos.poly[[1]])

#Iterates over each cell and structures the data into:
#genes: Gene expression data for the cell.
#poly: The polygon coordinates defining the cell shape.
#factors: Metadata, including time point and kidney compartment.
#xy: The spatial coordinates of the cell.
cells <- names(pos.poly)
final <- lapply(cells, function(cell) {
  out <- list(
    genes = as.list(CIS_gexp.sub[cell,]),
    poly = as.matrix(pos.poly[[cell]]),
    factors = list(time = as.character(meta[cell,1]), compartment = as.character(meta[cell,2])),
    xy = as.numeric(pos.all[cell,])
  )
  return(out)
})
names(final) <- cells

exportJSON <- toJSON(final, auto_unbox = TRUE)
write(exportJSON, "/Users/deev/Documents/R/Kidney App/CI_cells_Medulla_Up.json")

#Prepares a clusters.json file, containing:
#rows: Gene names
#cols: Cell names
#matrix: A scaled gene expression matrix, normalizing each gene 
# to its maximum value across samples.
cols = rownames(CIS_gexp.sub)
rows = colnames(CIS_gexp.sub)
matrix = lapply(1:ncol(CIS_gexp.sub), function(i) { as.numeric(CIS_gexp.sub[,i]/max(CIS_gexp.sub[,i])) }) ## need to be scaled?
final2 <- list(rows=rows, cols=cols, matrix=matrix)
exportJSON <- toJSON(final2)
write(exportJSON, "/Users/deev/Documents/R/Kidney App/CI_clusters_Medulla_Up.json")

###
base_url <- "http://localhost:8000/"

vc <- VitessceConfig$new(schema_version = "1.0.16", name = "Murine Kidney Cold Ischemia Medulla Up")
dataset <- vc$add_dataset("CI")
dataset$add_file(
  url = paste0(base_url, "CI_clusters_Medulla_Up.json"), 
  file_type = "clusters.json"
)
dataset$add_file(
  url = paste0(base_url, "CI_cells_Medulla_Up.json"),
  file_type = "cells.json"
)


desc <- vc$add_view(dataset, Component$DESCRIPTION)
desc <- desc$set_props(description = "Murine Kidney Cold Ischemia 2 with Cortex Up Genes")

spatial <- vc$add_view(dataset, Component$SPATIAL)
spatial_layers <- vc$add_view(dataset, Component$LAYER_CONTROLLER)
gene_list <- vc$add_view(dataset, Component$FEATURE_LIST)

vc$layout(vconcat(
  hconcat(desc, spatial_layers, gene_list),
  hconcat(vconcat(spatial))
))

vc$widget(theme = "light", width = "100%")

vc$export(with_config = TRUE, out_dir = "/Users/deev/Documents/R/Kidney App/vitessce-app")

