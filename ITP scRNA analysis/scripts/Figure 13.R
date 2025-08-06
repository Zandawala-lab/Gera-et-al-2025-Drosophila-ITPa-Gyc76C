# --- Figure 13 - Single-cell transcriptome analyses ---------------------------
# ------------------------------------------------------------------------------
# This is the code for:
# Figure 13: Panels G-J
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Files to download before running this script:
# ------------------------------------------------------------------------------
# loom file from https://scope.aertslab.org/#/Davie_et_al_Cell_2018/*/welcome

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
# Otherwise, skip to line 131
# ------------------------------------------------------------------------------
# Load brain scRNA data
Brain_loom <- Connect(filename = paste0(PATH_input,"Aerts_Fly_AdultBrain_Filtered_57k.loom"), mode = 'r')
Brain_loom
Brain_loom[["/matrix"]]
# Download loom file from https://scope.aertslab.org/#/Davie_et_al_Cell_2018/*/welcome

# Gene * cell matrix 
Brain_mat <- Brain_loom[["/matrix"]][,]
Brain_mat <- Matrix::Matrix(Brain_mat, sparse=T)
Brain_mat <- Matrix::t(Brain_mat)

# Extract the cell id and gene id and set them as the column and row names of the matrix
Brain_cellid <- Brain_loom[["/col_attrs/CellID"]][]
Brain_geneid <- Brain_loom[["/row_attrs/Gene"]][]

colnames(Brain_mat) <- Brain_cellid
rownames(Brain_mat) <- Brain_geneid

save(Brain_mat, file = paste0(PATH_input,"Brain_matrix.rda"))

# Create seurat object
Brain_seurat <- CreateSeuratObject(counts = Brain_mat, project = "adultbrain")
View(Brain_seurat@meta.data)

# Check data quality
plots <- lapply(c("nFeature_RNA", "nCount_RNA"), function(x) {
  VlnPlot(Brain_seurat, features = x) + theme_minimal()
})
wrap_plots(plots, ncol = 2)

FeatureScatter(Brain_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')

# Normalize data
Brain_seurat <- NormalizeData(Brain_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Scale data
all.genes <- rownames(Brain_seurat)
Brain_seurat <- ScaleData(Brain_seurat, features = all.genes)

save(Brain_seurat, file = paste0(PATH_input,"Brain_seurat.rda"))
load(paste0(PATH_input,"Brain_seurat.rda"))

# Extract L_NSC_DH31
L_NSC_DH31 <- subset(x=Brain_seurat, subset = ITP > 2 & Dh31 > 4 & amon > 0 & Phm > 0)
Cells(L_NSC_DH31)
L_NSC_DH31[["new_cluster"]] <- "L_NSC_DH31"
save(L_NSC_DH31, file = paste0(PATH_input,"L_NSC_DH31.rda"))

# Extract LN_ITP
LN_ITP <- subset(x=Brain_seurat, subset = ITP > 1 & NPF > 1 & cry > 0 & Phm > 0)
Cells(LN_ITP)
LN_ITP[["new_cluster"]] <- "LN_ITP"
save(LN_ITP, file = paste0(PATH_input,"LN_ITP.rda"))

# Extract L_NSC_ITP
L_NSC_ITP <- subset(x=Brain_seurat, subset = Tk > 1 & sNPF > 1 & ITP > 1 & ImpL2 > 1 & Crz == 0)
Cells(L_NSC_ITP)
L_NSC_ITP[["new_cluster"]] <- "L_NSC_ITP"
save(L_NSC_ITP, file = paste0(PATH_input,"L_NSC_ITP.rda"))

# Merge all ITP cell types
brainITPcells <- merge(L_NSC_ITP, y = c(L_NSC_DH31, LN_ITP), 
                       add.cell.ids = c("L_NSC_ITP", "L_NSC_DH31", "LN_ITP"), project = "ITP_brain")
                       save(brainITPcells, file = paste0(PATH_input,"brainITPcells.rda"))

# ------------------------------------------------------------------------------
# Start from here if you want to analyze ITP cells or only replicate the figures
# ------------------------------------------------------------------------------
load(paste0(PATH_input,"brainITPcells.rda"))
                      
# Figure 13G - marker genes for ITP cells --------------------------------------
p <- DotPlot(object = brainITPcells, features = c("amon","svr","Pal2", "Phm","Cadps",
                                             "ImpL2","cry","Trpm","ITP","NPF", 
                                             "Dh31","sNPF","Tk","Lk","Gpb5"), 
        dot.min = 0.05, scale = FALSE, group.by = "new_cluster") + 
        RotatedAxis() + 
        scale_size(range = c(1,8)) +
        labs(x = NULL, y = NULL)
                      
if(write_plots){
ggsave((paste0(PATH_output,"Figure 13G.pdf")), 
       plot = p, width = 20, height = 10, units = "cm")
}


# Figure 13H - Monoamine receptors in ITP cells  -------------------------------
NPexpression <- FetchData(object = brainITPcells, layer = "count", vars = c(MARlist))
NPmeans <- colMeans(NPexpression)
NPmeans <- names(NPmeans[NPmeans > 0.05])
q <- DotPlot(object = brainITPcells, features = c(NPmeans), dot.min = 0.05, 
        scale = FALSE, group.by = "new_cluster") + RotatedAxis() + 
        scale_size(range = c(1,8)) +
        labs(x = NULL, y = NULL)

if(write_plots){
ggsave((paste0(PATH_output,"Figure 13H.pdf")), 
       plot = q, width = 20, height = 10, units = "cm")
}

# Figure 13I - Neuropeptide receptors in ITP cells -----------------------------
NPexpression <- FetchData(object = brainITPcells, layer = "count", vars = c(NPRlist))
NPmeans <- colMeans(NPexpression)
NPmeans <- NPmeans[NPmeans > 0.05]
NPmeansorted <- as.matrix(names(sort(NPmeans,decreasing = TRUE)))
r <- DotPlot(object = brainITPcells, features = c(NPmeansorted), dot.min = 0.05, 
        scale = FALSE, group.by = "new_cluster") + RotatedAxis() + 
        scale_size(range = c(1,8)) +
        labs(x = NULL, y = NULL)

if(write_plots){
ggsave((paste0(PATH_output,"Figure 13I.pdf")), 
         plot = r, width = 30, height = 10, units = "cm")
} 


# Figure 13J - Neurotransmitter receptors in ITP cells -------------------------
NPexpression <- FetchData(object = brainITPcells, layer = "count", vars = c(NTRlist))
NPmeans <- colMeans(NPexpression)
NPmeans <- names(NPmeans[NPmeans > 0.05])
s <- DotPlot(object = brainITPcells, features = c(NPmeans), dot.min = 0.05, 
             scale = FALSE, group.by = "new_cluster") + 
    RotatedAxis() + scale_size(range = c(1,8)) + 
    labs(x = NULL, y = NULL)

if(write_plots){
ggsave((paste0(PATH_output,"Figure 13J.pdf")), 
       plot = s, width = 30, height = 10, units = "cm")
}
