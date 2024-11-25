# Clear environment # 

rm(list = ls())

# Setting up... #

library(ggplot2)
library(dplyr)
library(reshape2)
library(readxl)
library(broom)
library(lme4)
library(insight)
library(lmerTest)
library(stringr)


###################
# 1. Loading data #
###################

# Open data and platemaps #

pb <- read.delim("rnai_analysis_v1.txt")

plates <- read_excel("platemap_rnai_analysis_v1.xlsx")

# Annotate plate by adding platemap values through merging, IGNORE WARNING #

plates$MEASUREMENT.SET.ID <- as.character(plates$MEASUREMENT.SET.ID)
pb <- inner_join(pb, plates, by = c("MEASUREMENT.SET.ID", "Well.Name"))
head(pb)

# Remove missing values # 

pb_clean <- pb[which(pb$Cell..Average.Intensity_Average..Paramecium. != ""), ]
pb_clean <- pb_clean[which(pb_clean$MEASUREMENT.SET.ID != ""), ]
pb_clean <- pb_clean[which(pb_clean$Condition != "NA"), ]
pb_clean <- pb_clean[which(pb_clean$Rep != "NA"), ]

# Convert characters to numeric # 

pb_clean$Cell..Average.Intensity_Average..Paramecium. <- as.numeric(pb_clean$Cell..Average.Intensity_Average..Paramecium.)
pb_clean$Cell..Area_Average..Paramecium. <- as.numeric(pb_clean$Cell..Area_Average..Paramecium.)
pb_clean$Cell..Shape.Factor_Average..Paramecium. <- as.numeric(pb_clean$Cell..Shape.Factor_Average..Paramecium.)
pb_clean$Cell..Length_Average..Paramecium. <- as.numeric(pb_clean$Cell..Length_Average..Paramecium.)
pb_clean$Cell..Breadth_Average..Paramecium. <- as.numeric(pb_clean$Cell..Breadth_Average..Paramecium.)
pb_clean$Cell..Minimum.Intensity_Average..Paramecium. <- as.numeric(pb_clean$Cell..Minimum.Intensity_Average..Paramecium.)
pb_clean$Cell..Maximum.Intensity_Average..Paramecium. <- as.numeric(pb_clean$Cell..Maximum.Intensity_Average..Paramecium.)


############################
# 2. Subsetting and gating #
############################

# Day 6 cycloheximide data only # 

pb_clean_6 <- pb_clean[which(pb_clean$Day =="6"),] 
pb_clean_6_cycloheximide <- pb_clean_6[which(pb_clean_6$Treatment =="cycloheximide"),] 

# Gate dataset: live Pb cells between 0.6 - 0.9 # 

pb_gated <- pb_clean_6_cycloheximide[which(pb_clean_6_cycloheximide$Cell..Shape.Factor_Average..Paramecium. > 0.6), ]
pb_gated <- pb_gated[which(pb_gated$Cell..Shape.Factor_Average..Paramecium. < 0.90), ]
pb_gated <- pb_gated[which(pb_gated$Cell..Area_Average..Paramecium. < 10000), ]
#pb_gated <- pb_gated[which(pb_gated$Cell..Average.Intensity_Average..Paramecium. > 50), ]
#pb_gated <- pb_gated[which(pb_gated$Cell..Average.Intensity_Average..Paramecium. < 150), ]

h1 <- ggplot(pb_gated, aes(x = Cell..Average.Intensity_Average..Paramecium.)) +
  geom_histogram()


######################
# 3. Processing data #
######################

# Average fluorescence per well #

pb_f <- pb_clean_6_cycloheximide %>% group_by(Condition, Rep, Well.Name) %>% 
  summarize_at(vars(Cell..Average.Intensity_Average..Paramecium.), list(Mean = mean, sd = sd))

# Group raw data (for fluorescence) by condition, then calculate average #

pb_fl_average <- pb_gated %>% 
  group_by(Condition) %>% 
  summarize_at(vars(Cell..Average.Intensity_Average..Paramecium.), list(av_fluo = mean, sd_fluo = sd))
head(pb_fl_average)


############
# 4. Stats #
############

pb_fl_average <- pb_gated %>% 
  group_by(Condition, Rep, Well.Name) %>% 
  summarize_at(vars(Cell..Average.Intensity_Average..Paramecium.), list(av_fluo = mean, sd_fluo = sd))
head(pb_fl_average)

pbf <- as.data.frame(pb_fl_average)

# Split Well.Name into row and column - then convert these to factors
pbf$well_row <- as.factor(str_sub(pb_fl_average$Well.Name, 1, 1))
pbf$well_col <- as.factor(str_sub(pb_fl_average$Well.Name, 2, 3))

# function for extracting % variance explained by mixed model terms

percentvariance <- function(model){
  fix <- get_variance(model)$var.fixed
  ran <- get_variance(model)$var.random
  res <- get_variance(model)$var.residual
  total <- fix + ran + res
  result <- list(fixed = (fix/total)*100, random = (ran/total)*100, resid = (res/total)*100)
  return(result)
}

# Fluorescence data

# add column for experimental block
pbf$expblock <- rep(1, times=(dim(pbf)[1]))
pbf$expblock[pbf$Rep>6] <- 2
pbf$expblock <- as.factor(pbf$expblock)

# Re-level dataset so that comparisons are relative to scramble_150 (-ive control)

pbf$Condition <- as.factor(pbf$Condition)
pbf$Condition <- relevel(pbf$Condition, ref="scramble_150")

pbf$Condition <- factor(pbf$Condition, levels = c("scramble_150", "u2af", "33", "34", "35", "53", "54", 
                                                  "67", "60", "22", "61", "63", "43",
                                                  "24", "25", "26", "62", "28", "29", "30", "31", "32",
                                                  "37", "66", "44", "50", "55", "56"))

# Compare models with different structure 

glmm1a <- lmer(av_fluo ~ Condition + (1|well_row) + (1|well_col) + (1|expblock), data=pbf) 
glmm1b <- lmer(av_fluo ~ Condition + (1|Well.Name) + (1|expblock), data=pbf) 
glmm1c <- lmer(av_fluo ~ Condition + (1|expblock), data=pbf) 

AIC(glmm1a, glmm1b, glmm1c)
# glmm1a has best fit

percentvariance(glmm1a)
summary(glmm1a)

# Plot estimates from the model #

library(sjPlot)

sjPlot::plot_model(glmm1a) # in default plot red are negative compared to control, blue are positive