# run_cellchat.R
library(CellChat)  # Load the CellChat library

# Load the cellchat object; adjust the path as necessary
cellchat <- readRDS("objects/01_cellchat.rds")

# Run the function with the loaded cellchat object
runCellChatApp(cellchat)