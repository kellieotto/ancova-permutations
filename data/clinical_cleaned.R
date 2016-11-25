# Kellie Ottoboni
# November 14, 2016
# Script to privatize the RB dataset

require(dplyr)
require(reshape2)



ga <- read.csv("GA1102.csv", header = TRUE, stringsAsFactors = FALSE)

# Remove missing
ga <- ga %>% filter(tr != ".")

# Change treatment values to drug A vs drug B
ga <- ga %>% mutate(tr = factor(tr, labels = c("A", "B")))

# Rename subject ids
ga <- ga %>% mutate(SUBJID = as.numeric(factor(SUBJID)))

# Rename site ids
ga <- ga %>% mutate(SITEID = as.numeric(factor(SITEID)))

# Take means
ga <- ga %>% group_by(SUBJID, SITEID, VISITNUM, tr, country) %>%
  summarise_each(funs(mean)) %>% 
  ungroup()

write.csv(ga, file = "clinical_cleaned.csv", row.names = FALSE)
