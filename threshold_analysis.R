# MDS Threshold Results Analysis 
# Load packages
library(tidyverse)

# Code
mds_results <- read.csv("mds_results.csv")

# Change data 
mds_results <- mds_results %>%
  mutate(TrueGroupLabel = ifelse(True.Number.of.Groups == 1, "1", ">1"),
         DeterminedGroupLabel = Determined.Number.of.Groups) %>%
  select(-True.Number.of.Groups, -Determined.Number.of.Groups)

# Mclust success 
mclust_success <- mds_results %>%
  filter(str_detect(Method, "mclust")) %>%
  mutate(Method =  str_replace(Method, "mclustBootstrapLRT: ", "")) 


mclust_summary <- mclust_success %>%
  mutate(Match = TrueGroupLabel == DeterminedGroupLabel) %>%
  group_by(Method) %>%
  summarise(Total = n(),
    Matches = sum(Match, na.rm = TRUE),
    MatchRate = round(Matches / Total * 100, 2))

print(mclust_summary)

mclust_success %>%
  filter(Method %in% c("EEE","EEI")) %>%
  group_by(Parameter.Label, Columns, Method) %>%
  mutate(Match = TrueGroupLabel == DeterminedGroupLabel) %>%
  summarise(Total = n(),
            Matches = sum(Match, na.rm = TRUE),
            MatchRate = round(Matches / Total * 100, 2))

# OTRIMLE Success 
otrimle_success <- mds_results %>%
  filter(str_detect(Method, "OTRIMLE")) %>%
  select(-Method)

otrimle_summary <- otrimle_success %>%
  mutate(Match = TrueGroupLabel == DeterminedGroupLabel) %>%
  group_by(Parameter.Label, Columns) %>%
  summarise(Total = n(),
            Matches = sum(Match, na.rm = TRUE),
            MatchRate = round(Matches / Total * 100, 2))

print(otrimle_summary)