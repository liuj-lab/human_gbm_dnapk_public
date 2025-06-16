save_to_metadata_file <- function(metadata_file, text) {
  file_conn <- file(metadata_file, "a")
  writeLines(text, file_conn)
  close(file_conn)
}

create_dir <- function(path) {
  if (!file.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}
