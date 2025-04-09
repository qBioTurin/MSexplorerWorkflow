createFolder <- function(folder_path) {
  # Controlla se la cartella esiste
  if (!dir.exists(folder_path)) {
    # Crea la cartella
    dir.create(folder_path, recursive = TRUE)
    message("Folder created: ", folder_path)
  } else {
    message("The folder already exists: ", folder_path)
  }
}
