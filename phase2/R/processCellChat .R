library(Seurat)
library(CellChat)
library(future)

processCellChat <- function(obj, 
                            save_path = "objects/00_cellchat.rds", 
                            num_workers = 8, 
                            search = "Secreted Signaling") {
  # Reformat Seurat object
  data.input <- GetAssayData(obj)
  meta <- data.frame(labels = Idents(obj), row.names = names(Idents(obj)))
  
  # Create CellChat object
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels", assay = "SCT")
  
  # Set the ligand-receptor interaction database
  CellChatDB <- CellChatDB.human
  
  # Subset the database based on the search flag
  CellChatDB.use <- subsetDB(CellChatDB, search = search, key = "annotation")
  cellchat@DB <- CellChatDB.use
  
  # Subset database
  cellchat <- CellChat::subsetData(cellchat)
  
  # Parallelize analysis
  future::plan("multisession", workers = num_workers)
  
  # Identify genes and interactions
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  
  # Calculate the average gene expression per cell group
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  
  # Filter min cells for interaction network
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # Extract communication dataframe
  df.net <- subsetCommunication(cellchat)
  
  # Compute pathway probability
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Calculate cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  
  # Add embeddings
  cellchat <- addReduction(object = cellchat, seu.obj = obj)
  
  # Save the CellChat object
  tryCatch({
    print(paste("Saving CellChat object to:", save_path))
    saveRDS(cellchat, file = save_path)
  }, error = function(e) {
    print(paste("Error saving file:", e$message))
  })
  
  return(cellchat)
}

# Example usage:
# seurat_obj <- readRDS("path_to_your_seurat_object.rds")
# cellchat_result <- processCellChat(seurat_obj, search = "Secreted Signaling")
