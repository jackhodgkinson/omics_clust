# MDS Threshold Results Analysis 
# Load packages
library(tidyverse)

log_con <- file("mds_threshold.txt", open = "wt")

# Start logging both output and messages
sink(log_con)                     # Standard output
sink(log_con, type = "message")

# Code
mds_results <- read.csv("mds_results_combined.csv")

# Reorder data 
mds_results <- mds_results %>%
  arrange(Parameter.ID, Columns, True.Number.of.Groups,RunTime)

# Remove duplicates
mds_results <- mds_results[!duplicated(mds_results), ]

# Count NA
sum(is.na(mds_results$RunTime))
print(paste("Percentage of Failed Runs:", round(sum(is.na(mds_results$RunTime))/nrow(mds_results)*100, 2),"%"))

# Remove NA 
mds_results <- na.omit(mds_results)

# Mclust success 
mclust_success <- mds_results %>%
  filter(str_detect(Method, "mclust")) %>%
  mutate(Method =  str_replace(Method, "mclustBootstrapLRT: ", "")) 

mclust_summary <- mclust_success %>%
  mutate(Match = True.Number.of.Groups == Determined.Number.of.Groups) %>%
  group_by(Method) %>%
  summarise(Total = n(),
    Matches = sum(Match, na.rm = TRUE),
    MatchRate = round(Matches / Total * 100, 2))

print(mclust_summary)

# OTRIMLE Success 
otrimle_success <- mds_results %>%
  filter(str_detect(Method, "OTRIMLE")) %>%
  select(-Method)

otrimle_summary <- otrimle_success %>%
  mutate(Match = True.Number.of.Groups == Determined.Number.of.Groups) %>%
  group_by(Parameter.Label, Columns) %>%
  summarise(Total = n(),
            Matches = sum(Match, na.rm = TRUE),
            MatchRate = round(Matches / Total * 100, 2))

print(otrimle_summary)

sink(type = "message")
sink()
close(log_con) 