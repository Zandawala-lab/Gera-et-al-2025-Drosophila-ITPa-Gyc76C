# --- Figure 13 Supplement 3 (S3)- Single-cell transcriptome analyses -----------
# ------------------------------------------------------------------------------
# This is the code for:
# Figure 13 Supplement 3: Panels A-D
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Files to download before running this script:
# ------------------------------------------------------------------------------
# loom file from https://scope.aertslab.org/#/Davie_et_al_Cell_2018/Davie_et_al_Cell_2018%2FGoodwin_Fly_AdultVNC_elife54074.loom/gene

# Put file in the 'input' directory: ITP scRNA analysis/input/

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load packages
# ------------------------------------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

PATH_input = "./input/"
PATH_output = "./output/"

# Check if folders exist, if not, make them
if(!dir.exists(file.path(PATH_input))){
  dir.create(file.path(PATH_input))
}
if(!dir.exists(file.path(PATH_output))){
  dir.create(file.path(PATH_output))
}

# ------------------------------------------------------------------------------
# Options for saving plots
# ------------------------------------------------------------------------------

# When running script for the first time, set to TRUE:
write_plots = TRUE           # TRUE - save/replicate figure plots
                             # FALSE - plots not saved outside of R
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# Import the list of genes of interest
load(paste0(PATH_input,"neuropeptides.rda"))               # neuropeptidelist
load(paste0(PATH_input,"neuropeptide_receptors.rda"))      # NPRlist
load(paste0(PATH_input,"monoamine_receptors.rda"))         # MARlist
load(paste0(PATH_input,"neurotransmitter_receptors.rda"))  # NTRlist

# ------------------------------------------------------------------------------
# Start from here if you want to start with the original data
# Otherwise, skip to line 125
# ------------------------------------------------------------------------------

# Load VNC scRNA data
VNC_loom <- Connect(filename = paste0(PATH_input,"Goodwin_Fly_AdultVNC_elife54074.loom"), mode = 'r')
VNC_loom
VNC_loom[["/matrix"]]
#download loom file from:
#https://scope.aertslab.org/#/Davie_et_al_Cell_2018/Davie_et_al_Cell_2018%2FGoodwin_Fly_AdultVNC_elife54074.loom/gene

# Gene * cell matrix 
VNC_mat <- VNC_loom[["/matrix"]][,]
VNC_mat <- Matrix::Matrix(VNC_mat, sparse=T)
VNC_mat <- Matrix::t(VNC_mat)

# Extract the cell id and gene id and set them as the column and row names of the matrix
VNC_cellid <- VNC_loom[["/col_attrs/CellID"]][]
VNC_geneid <- VNC_loom[["/row_attrs/Gene"]][]

colnames(VNC_mat) <- VNC_cellid
rownames(VNC_mat) <- VNC_geneid

save(VNC_mat, file = paste0(PATH_input,"VNC_matrix.rda"))

# Create seurat object
VNC_seurat <- CreateSeuratObject(counts = VNC_mat, project = "adultVNC")
View(VNC_seurat@meta.data)

# Check data quality
plots <- lapply(c("nFeature_RNA", "nCount_RNA"), function(x) {
  VlnPlot(VNC_seurat, features = x) + theme_minimal()
})
wrap_plots(plots, ncol = 2)

FeatureScatter(VNC_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')

# Normalize data
VNC_seurat <- NormalizeData(VNC_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Scale data
all.genes <- rownames(VNC_seurat)
VNC_seurat <- ScaleData(VNC_seurat, features = all.genes)

save(VNC_seurat, file = paste0(PATH_input,"VNC_seurat.rda"))
load(paste0(PATH_input,"VNC_seurat.rda"))

# Extract IAG neurons
IAG <- subset(x=VNC_seurat, subset = AstA > 0 & CCAP > 0 & ITP > 1 & Phm > 0 & amon > 0)
Cells(IAG)
IAG[["new_cluster"]] <- "IAG"
save(IAG, file = paste0(PATH_input,"IAG.rda"))

# Extract non_IAG neurons
nonIAG <- subset(x=VNC_seurat, subset = AstA == 0 & CCAP == 0 & ITP > 1 & Phm > 0 & amon > 0)
Cells(nonIAG)
nonIAG[["new_cluster"]] <- "non_IAG"
save(nonIAG, file = paste0(PATH_input,"nonIAG.rda"))

# Merge both ITP cells
VNCITPcells <- merge(nonIAG, IAG)
save(VNCITPcells, file = paste0(PATH_input,"VNCITPcells.rda"))

# ------------------------------------------------------------------------------
# Start from here if you want to analyze ITP cells or only replicate the figures
# ------------------------------------------------------------------------------
load(paste0(PATH_input,"VNCITPcells.rda"))

# Figure 13 S3A - Marker genes for ITP cells -----------------------------------
p <- DotPlot(object = VNCITPcells, features = c("amon","svr","Pal2", "Phm","Cadps",
                                                "ITP", "CCAP", "AstA","Gpb5","Ddc"), 
        dot.min = 0.05, scale = FALSE, group.by = "new_cluster") + 
        RotatedAxis() + 
        scale_size(range = c(1,8)) +
        labs(x = NULL, y = NULL)

if(write_plots){
ggsave((paste0(PATH_output,"Figure 13S3A.pdf")), 
       plot = p, width = 16, height = 10, units = "cm")
}

# Figure 13 S3B - Monoamine receptors in ITP cells -----------------------------
NPexpression <- FetchData(object = VNCITPcells, layer = "count", vars = c(MARlist))
NPmeans <- colMeans(NPexpression)
NPmeans <- names(NPmeans[NPmeans > 0.05])
q <- DotPlot(object = VNCITPcells, features = c(NPmeans), dot.min = 0.05, 
        scale = FALSE, group.by = "new_cluster") + 
        RotatedAxis() + 
        scale_size(range = c(1,8)) +
        labs(x = NULL, y = NULL)
if(write_plots){
ggsave((paste0(PATH_output,"Figure 13S3B.pdf")), 
       plot = q, width = 18, height = 10, units = "cm")
}

# Figure 13 S3C - Neuropeptide receptors in ITP cells --------------------------
NPexpression <- FetchData(object = VNCITPcells, layer = "count", vars = c(NPRlist))
NPmeans <- colMeans(NPexpression)
NPmeans <- NPmeans[NPmeans > 0.05]
NPmeansorted <- as.matrix(names(sort(NPmeans,decreasing = TRUE)))
r <- DotPlot(object = VNCITPcells, features = c(NPmeansorted), dot.min = 0.05, 
        scale = FALSE, group.by = "new_cluster") + 
        RotatedAxis() + 
        scale_size(range = c(1,8)) +
        labs(x = NULL, y = NULL)

if(write_plots){
ggsave((paste0(PATH_output,"Figure 13S3C.pdf")), 
         plot = r, width = 20, height = 10, units = "cm")
}

# Figure 13 S3D -Neurotransmitter receptors in ITP cells -----------------------
NPexpression <- FetchData(object = VNCITPcells, layer = "count", vars = c(NTRlist))
NPmeans <- colMeans(NPexpression)
NPmeans <- names(NPmeans[NPmeans > 0.05])
s <- DotPlot(object = VNCITPcells, features = c(NPmeans), dot.min = 0.05, 
             scale = FALSE, group.by = "new_cluster") + 
    RotatedAxis() + 
    scale_size(range = c(1,8)) + 
    labs(x = NULL, y = NULL)

if(write_plots){
ggsave((paste0(PATH_output,"Figure 13S3D.pdf")), 
       plot = s, width = 25, height = 10, units = "cm")
}
