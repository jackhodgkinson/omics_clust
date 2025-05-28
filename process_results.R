# Combine all result files
result_files <- list.files("results", pattern = "\\.csv$", full.names = TRUE)

# Function to extract output_name from filename
extract_output_name <- function(filepath) {
  fname <- basename(filepath)
  sub("_\\d+\\.csv$", "", fname)
}

for (file in result_files) {
  output_name <- extract_output_name(file)
  
  # Read the file and add output_name column
  df <- read.csv(file)
  df$output_name <- output_name
  
  # Save to new CSV file as <output_name>.csv
  output_csv <- paste0(output_name, ".csv")
  write.csv(df, output_csv, row.names = FALSE)
  cat("Saved:", output_csv, "\n")
  
  # If "task_id" column exists, find all unique task_ids in this file
  if ("task_id" %in% colnames(df)) {
    task_ids <- unique(df$task_id)
    
    for (task_id in task_ids) {
      for (ext in c("out", "err")) {
        log_file <- file.path("logs", paste0(output_name, "_", task_id, ".", ext))
        if (file.exists(log_file)) {
          file.remove(log_file)
          cat("Deleted:", log_file, "\n")
        } else {
          cat("Not found:", log_file, "\n")
        }
      }
    }
  } else {
    cat("No task_id column found in", file, "\n")
  }
}

