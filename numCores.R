# numCores.R
numCores <- function(fallback = 1) {
  # Try multiple known environment variables used by different schedulers
  env_vars <- c(
    "SLURM_CPUS_PER_TASK",     # Slurm
    "PBS_NUM_PPN",             # PBS
    "NSLOTS",                  # SGE
    "OMP_NUM_THREADS",         # OpenMP or general hint
    "NUM_THREADS",             # Generic cloud var
    "NUMBER_OF_PROCESSORS"     # Windows-specific sometimes
  )
  
  # Try each environment variable
  for (var in env_vars) {
    val <- Sys.getenv(var, unset = NA)
    if (!is.na(val) && nzchar(val)) {
      cores <- suppressWarnings(as.integer(val))
      if (!is.na(cores) && cores > 0) return(cores)
    }
  }
  
  # Fall back to local detection
  cores <- parallel::detectCores(logical = FALSE)
  if (is.na(cores) || cores < 1) return(fallback)
  return(cores)
}