suppressPackageStartupMessages({
    library(Seurat)
    library(celldex)
    library(ShinyCell)
    library(dplyr)
    library(SingleR)
})

setwd("/data/nsiegel/projects/amp_data/phase2")

# Load the main object
load(file=file.path("objects", "06.rda"))

# Update metadata
obj@meta.data$Stim <- obj@meta.data$Type
obj@meta.data$Tissue_Type <- obj@meta.data$broad.group
obj@meta.data$Fine_Cell_Type <- obj@meta.data$fine.type
obj@meta.data$Broad_Cell_Type <- obj@meta.data$broad.type
obj@meta.data$Donor <- gsub("AMPSLEkid_cells", "Donor",  obj@meta.data$sample)

columns_to_keep <-  c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'Donor', 'Stim', 'Tissue_Type', 'Fine_Cell_Type', 'Broad_Cell_Type')
obj@meta.data <- obj@meta.data[, columns_to_keep, drop = FALSE]

# Define a function for PCA, UMAP, and clustering
process_and_save <- function(seu, name) {
    print(paste("Processing", name))
    seu <- RunPCA(seu, npcs = 30, verbose = TRUE)
    seu <- RunUMAP(seu, reduction = "pca", dims = 1:20)
    seu <- FindNeighbors(seu, reduction = "pca", dims = 1:20)
    seu <- FindClusters(seu, resolution = 0.5)
    save(seu, file=file.path("objects", paste0(name, ".rda")))
    print(paste("Saved", name, "to file"))
    return(seu)
}

# Subset and process each cell type
cell_types <- c("B Cell", "DN", "GLOM", "INTL", "LOH", "Myeloid Cell", "NK Cell", "Plasma Cell", "PT", "T Cell", "UND")
for (cell_type in cell_types) {
    subset_seu <- subset(obj, subset = Broad_Cell_Type == cell_type)
    subset_seu <- process_and_save(subset_seu, gsub(" ", "", cell_type))
}

rm(obj)

# # Load annotation references
# monaco.ref <- celldex::MonacoImmuneData()
# blueprint.ref <- celldex::BlueprintEncodeData()

# # Define a function for SingleR annotation
# annotate_and_save <- function(seu, name, reference, label_type, suffix) {
#     print(paste("Annotating", name, "with", label_type))
#     load(file=file.path("objects", paste0(name, ".rda")))
#     sce <- LayerData(seu)
#     annotation <- SingleR(test = sce, assay.type.test = 1, ref = reference, labels = reference[[label_type]])
#     seu@meta.data[[paste0(name, ".", suffix)]] <- annotation$pruned.labels
#     save(seu, file=file.path("objects", paste0(name, suffix, ".rda")))
#     print(paste("Saved", name, suffix, "annotations"))
# }

# # Annotate and save each cell type
# annotations <- list(
#     Blueprint = "label.main",
#     Monaco = "label.fine"
# )

# for (cell_type in cell_types) {
#     seu_name <- gsub(" ", "", cell_type)
#     for (ann_type in names(annotations)) {
#         annotate_and_save(seu = get(seu_name), name = seu_name, reference = get(paste0(tolower(ann_type), ".ref")), label_type = annotations[[ann_type]], suffix = paste0(tolower(ann_type), "05"))
#     }
# }

# # Write Shiny files
# shiny_dir <- "kidney_amp2_sle_scRNAseq/"
# shiny_prefixes <- c("sc1", "sc2", "sc3", "sc4", "sc5", "sc6", "sc7", "sc8", "sc9", "sc10", "sc11")
# cell_types_labels <- c("B Cell", "DN", "GLOM", "INTL", "LOH", "Myeloid", "NK", "Plasma Cell", "PT", "T Cell", "UND")

# for (i in seq_along(cell_types)) {
#     seu <- get(gsub(" ", "", cell_types[i]))
#     scConf <- createConfig(seu)
#     makeShinyFiles(seu, scConf,
#                    gene.mapping = TRUE,
#                    shiny.prefix = shiny_prefixes[i],
#                    shiny.dir = shiny_dir,
#                    gex.assay = "SCT",
#                    default.multigene = present_genes)
#     print(paste("Shiny files for", cell_types[i], "created"))
# }

# # Create Shiny codes
# makeShinyCodesMulti(
#     shiny.title = "AMP Phase II SLE", shiny.footnotes = NULL,
#     shiny.prefix = shiny_prefixes,
#     shiny.headers = cell_types_labels, 
#     shiny.dir = shiny_dir
# )
# print("Shiny codes created")