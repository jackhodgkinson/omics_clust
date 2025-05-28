# Combine all results files 
# Get task_id from command line argument 
args <- commandArgs(trailingOnly = TRUE)
task_id <- args[1]

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
  
  # If task_id is in the data, delete matching .out/.err files
  if ("task_id" %in% colnames(df) && task_id %in% df$task_id) {
    for (ext in c("out", "err")) {
      log_file <- file.path("logs", paste0(output_name, "_", task_id, ".", ext))
      if (file.exists(log_file)) {
        file.remove(log_file)
        cat("Deleted:", log_file, "\n")
      } else {
        cat("Not found:", log_file, "\n")
      }
    }
  } else {
    cat("Task ID", task_id, "not found in", file, "\n")
  }
}

