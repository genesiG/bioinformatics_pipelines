## Read data table with each replicate (well) measurement as a column
viabilityNTC = read.csv("C:/Users/thephillipslab/Desktop/mock_values.csv")

## Convert data to long format using the tidyr package
if(!require("tidyr")){
  install.packages("tidyr")
  library(tidyr)
}
viability_longer = pivot_longer(viabilityNTC,
                                cols = c("Rep1", "Rep2", "Rep3", "Rep4"),
                                names_to = NULL) 

## Turn your conditions into factors
viability_longer$Pre.treatment <- factor(viability_longer$Pre.treatment)
viability_longer$Treatment <- factor(viability_longer$Treatment)
viability_longer$Genotype <- factor(viability_longer$Genotype)

## Make sure your control is the first level
viability_longer$Pre.treatment <- relevel(viability_longer$Pre.treatment, "DMSO")
viability_longer$Treatment <- relevel(viability_longer$Treatment, "DMSO")
viability_longer$Genotype <- relevel(viability_longer$Genotype, "NTC")

## Group data by each condition then create summary table to make bar plots

# Use the piping operator (%>%) from the dplyr package
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}

viability_summary <- viability_longer %>%
  group_by(Genotype) %>%
  group_by(Treatment, .add = T) %>%
  group_by(Pre.treatment, .add = T) %>%
  summarise(
    sd = sd(value, na.rm = TRUE),
    value = mean(value)
  )

viability_summary

## Make your control the first level
viability_summary$Pre.treatment <- relevel(factor(viability_summary$Pre.treatment), "DMSO")
viability_summary$Treatment <- relevel(factor(viability_summary$Treatment), "DMSO")
viability_summary$Genotype <- relevel(factor(viability_summary$Genotype), "NTC")

## Plot data
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}

library(RColorBrewer)
cols <- brewer.pal(9,"BuPu")[5:9]

# Create ggplot object with the data table in the long format
prettyPlot <- ggplot(viability_longer, 
                     aes(x = Pre.treatment, 
                         # Turn relative (normalized) values into percentages
                         y = value*100, 
                         # Color data by genotype
                         fill = Genotype)) + 
  
  # Add/modify the theme
  theme_light() +
  theme(strip.text = element_text(color = "black")) +
  
  # Use the "BuPu" color palette (or any other palette) from RColorBrewer
  scale_fill_brewer(palette = "BuPu") +
  
  # Label axes
  ylab("% Relative cell viability at day 5 post-treatment (mean +/- SD)") +
  xlab("Pre-treatment (18 days)") +
  
  # Add geom column (or bar) using the SUMMARY TABLE
  geom_col(data = viability_summary, 
           # Make the border of the columns black
           color = "black",
           # Dodge the columns' position so it looks prettier
           position=position_dodge(0.9)) + 
  
  # Add jitter points if you want
   # geom_jitter(position = position_jitterdodge(jitter.width = 0.5, 
   #                                             dodge.width = 0.9), 
   #             color = "black") + 
  
  # Add error bars 
  geom_errorbar(aes(ymin = (value-sd)*100, 
                    ymax = (value+sd)*100),
                # Adjust width of error bars
                width = 0.2, 
                # For grouped columns, dodge the error bars by the same amount as for the columns
                position=position_dodge(0.9), 
                # use the SUMMARY table of mean and sd
                data = viability_summary
  ) +
  # Each sample was pre-treated by 18 days then each pre-treated sample was treated by another 5 days with both concentrations
  # We set our x-axis as pre-treatment, now we facet data by treatment
  facet_wrap(vars(Treatment), ncol = 1) 

# Export as pdf so it looks BEAUTIFUL
pdf("C:/Users/thephillipslab/Desktop/prettyPlot.pdf", 
    width = 6, 
    height = 7)
prettyPlot
dev.off()

