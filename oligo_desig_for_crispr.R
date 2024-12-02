# OLIGO DESIGN FOR CRISPR
# AUTOMATED TABLE

setwd("./Projects/MGB2_screen/followup_screen/")

# Load spreadsheet
csv = read.csv("Oligo design for CRISPR.csv")

csv = csv[,c(2,3,8,11)]

# Forward and reverse gRNA sequences are in different columns
# Need to consolidate them into one column for ordering

# Install and load tidyr if not already installed
# install.packages("tidyr")
library(tidyr)
library(dplyr)

# Use pivot_longer to combine the specified columns
csv_long <- csv %>%
  pivot_longer(cols = c(Sense.Add.5..overhang, Anti.sense.Add.5..overhang),
               names_to = "Type",
               values_to = "Sequence") %>%
  group_by(Gene, Guide) %>%
  arrange(Gene, Guide, desc(Type)) %>%  # Arrange to alternate Type
  ungroup() %>%
  mutate(Guide = case_when(
    Type == "Sense.Add.5..overhang" ~ paste0(Guide, "_F"),
    Type == "Anti.sense.Add.5..overhang" ~ paste0(Guide, "_R")
  ))


# Display the result
print(csv_long)
# Check
head(arrange(csv, Gene))


write.csv(csv_long, "./oligo_order_table.csv", 
          quote = FALSE, row.names = FALSE)
