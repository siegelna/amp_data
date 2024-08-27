cleanCellChatNames <- function(cellchat) {
  # Ensure the input is a CellChat object with a net$count matrix
  if (!("CellChat" %in% class(cellchat)) || is.null(cellchat@net$count)) {
    stop("The provided object is not a valid CellChat object or does not contain a 'net$count' matrix.")
  }
  
  # Get column and row names
  col_names <- colnames(cellchat@net$count)
  row_names <- rownames(cellchat@net$count)
  
  # Remove numbers and periods from column names
  clean_col_names <- gsub("[0-9]+\\.", "", col_names)
  # Remove numbers and periods from row names
  clean_row_names <- gsub("[0-9]+\\.", "", row_names)
  
  # Trim any leading or trailing whitespace
  clean_col_names <- trimws(clean_col_names)
  clean_row_names <- trimws(clean_row_names)
  
  # Update the column and row names of the matrix
  colnames(cellchat@net$count) <- clean_col_names
  rownames(cellchat@net$count) <- clean_row_names
  
  # Return the modified CellChat object
  return(cellchat)
}

# Example usage:
# cleaned_cellchat <- cleanCellChatNames(cellchat)
