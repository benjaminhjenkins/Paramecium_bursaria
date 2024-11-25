# Clear environment # 

rm(list = ls())

# Setting up... #

library(ggplot2)
library(dplyr)
library(reshape2)
library(readxl)
line <- "black"


# Open data and platemaps #

pb <- read.delim("glycan_feeding_highlight_v1.txt")

platemap <- read_excel("platemap_glycan_feeding_highlight_v1.xlsx")
head(platemap)

ylab = expression(paste("cell number ", "(ml"^"-1", ")"))

# Replace plate ID with day #

pb$MEASUREMENT.SET.ID <- ifelse(pb$MEASUREMENT.SET.ID == "809", "0u1",
                                ifelse(pb$MEASUREMENT.SET.ID == "810", "0u2",
                                       ifelse(pb$MEASUREMENT.SET.ID == "812", "0h1", 
                                              ifelse(pb$MEASUREMENT.SET.ID == "813", "0h2",
                                                     ifelse(pb$MEASUREMENT.SET.ID == "824", "2u1",
                                                            ifelse(pb$MEASUREMENT.SET.ID == "825", "2u2",
                                                                  ifelse(pb$MEASUREMENT.SET.ID == "827", "2h1", 
                                                                         ifelse(pb$MEASUREMENT.SET.ID == "828", "2h2", NA ))))))))

# Annotate well.name by adding env value through merging, IGNORE WARNING #

names(pb)[names(pb) == 'vars'] <- 'Well.Name'

pb_annotated <- inner_join(pb, platemap, by = "Well.Name")
head(pb_annotated)

# Split dataset #

pb_annotated_1 <- pb_annotated[which(pb_annotated$MEASUREMENT.SET.ID %in% c("0u1", "0h1", "2u1", "2h1")), ]
pb_annotated_2 <- pb_annotated[which(pb_annotated$MEASUREMENT.SET.ID %in% c("0u2", "0h2", "2u2", "2h2")), ]

# Replace Env with name #

pb_annotated_1$Env <- ifelse(pb_annotated_1$Env == "1", "- 1",
                           ifelse(pb_annotated_1$Env == "2", "Glc 0.1",
                                  ifelse(pb_annotated_1$Env == "3", "Glc 1",
                                         ifelse(pb_annotated_1$Env == "4", "GlcN 0.1", 
                                                ifelse(pb_annotated_1$Env == "5", "GlcN 1",
                                                       ifelse(pb_annotated_1$Env == "6", "GlcNAc 0.1", 
                                                              ifelse(pb_annotated_1$Env == "7", "GlcNAc 1", 
                                                                     ifelse(pb_annotated_1$Env == "8", "GlcNAc-6P 0.1", 
                                                                            ifelse(pb_annotated_1$Env == "9", "GlcNAc-6P 1",
                                                                                   ifelse(pb_annotated_1$Env == "10", "GlcN-6P 0.1", NA))))))))))

pb_annotated_2$Env <- ifelse(pb_annotated_2$Env == "1", "- 2",
                             ifelse(pb_annotated_2$Env == "2", "GlcN-6P 1",
                                    ifelse(pb_annotated_2$Env == "3", "CH2 0.1",
                                           ifelse(pb_annotated_2$Env == "4", "CH2 1", 
                                                  ifelse(pb_annotated_2$Env == "5", "CH3 0.1",
                                                         ifelse(pb_annotated_2$Env == "6", "CH3 1", 
                                                                ifelse(pb_annotated_2$Env == "7", "CH4 0.1", 
                                                                       ifelse(pb_annotated_2$Env == "8", "CH4 1", 
                                                                              ifelse(pb_annotated_2$Env == "9", "CH5 0.1",
                                                                                     ifelse(pb_annotated_2$Env == "10", "CH5 1", NA))))))))))

# Merge #

pb_annotated_merged <- rbind(pb_annotated_1, pb_annotated_2)

# Remove missing values # 

pb_clean <- pb_annotated_merged[which(pb_annotated_merged$Cell..Average.Intensity_Average..Paramecium. != ""), ]
pb_clean <- pb_annotated_merged[which(pb_annotated_merged$MEASUREMENT.SET.ID != ""), ]

# Convert characters to numeric # 

pb_clean$Cell..Average.Intensity_Average..Paramecium. <- as.numeric(pb_clean$Cell..Average.Intensity_Average..Paramecium.)
pb_clean$Cell..Area_Average..Paramecium. <- as.numeric(pb_clean$Cell..Area_Average..Paramecium.)
pb_clean$Cell..Shape.Factor_Average..Paramecium. <- as.numeric(pb_clean$Cell..Shape.Factor_Average..Paramecium.)
pb_clean$Cell..Length_Average..Paramecium. <- as.numeric(pb_clean$Cell..Length_Average..Paramecium.)
pb_clean$Cell..Breadth_Average..Paramecium. <- as.numeric(pb_clean$Cell..Breadth_Average..Paramecium.)
pb_clean$Cell..Minimum.Intensity_Average..Paramecium. <- as.numeric(pb_clean$Cell..Minimum.Intensity_Average..Paramecium.)
pb_clean$Cell..Maximum.Intensity_Average..Paramecium. <- as.numeric(pb_clean$Cell..Maximum.Intensity_Average..Paramecium.)

# Gate dataset #    # Average Pb cell between 0.5 - 0.8 # 

pb_gated <- pb_clean[which(pb_clean$Cell..Shape.Factor_Average..Paramecium. > 0.5), ]
pb_gated <- pb_gated[which(pb_gated$Cell..Shape.Factor_Average..Paramecium. < 0.8), ]
pb_gated <- pb_gated[which(pb_gated$Cell..Average.Intensity_Average..Paramecium. > 50), ]
#pb_gated <- pb_gated[which(pb_gated$Cell..Average.Intensity_Average..Paramecium. < 300), ]


##########################
# Change in fluorescence # 
##########################

# Average fluorescence per well #

pb_f <- pb_gated %>% group_by(MEASUREMENT.SET.ID, Env, Replicate) %>% 
  summarize(Mean = mean(Cell..Average.Intensity_Average..Paramecium.))

# Change #

pb_f_change_a <- dcast(pb_f, Env + Replicate ~ MEASUREMENT.SET.ID, value.var = 'Mean')
pb_f_change_a[is.na(pb_f_change_a)] <- 0
pb_f_change_a$diff <- pb_f_change_a$"2h1" - pb_f_change_a$"0h1"
pb_f_change_a$pct_change <- (pb_f_change_a$"2h1" - pb_f_change_a$"0h1")/pb_f_change_a$"0h1" *100
pb_f_change_a <- na.omit(pb_f_change_a)
# Add 100 to give positive values for stats #
pb_f_change_a$pct_change_stats <- pb_f_change_a$pct_change + 100

pb_f_change_b <- dcast(pb_f, Env + Replicate ~ MEASUREMENT.SET.ID, value.var = 'Mean')
pb_f_change_b[is.na(pb_f_change_b)] <- 0
pb_f_change_b$diff <- pb_f_change_b$"2h2" - pb_f_change_b$"0h2"
pb_f_change_b$pct_change <- (pb_f_change_b$"2h2" - pb_f_change_b$"0h2")/pb_f_change_b$"0h2" *100
pb_f_change_b <- na.omit(pb_f_change_b)
# Add 100 to give positive values for stats #
pb_f_change_b$pct_change_stats <- pb_f_change_b$pct_change + 100

# Merge #

pb_f_change <- rbind(pb_f_change_a, pb_f_change_b)
head(pb_f_change)

# Aesthetics #

pb_f_change$Env <- as.factor(pb_f_change$Env)
pb_f_change$pct_change <- as.numeric(pb_f_change$pct_change)


################
# Plot results #
################

p1 <- ggplot(pb_f_change, aes(x = Env, y = pct_change)) +
  geom_boxplot(colour = line, outlier.shape = NA, size =1) + 
  geom_point(aes(x = Env, y = pct_change, colour = Env),
             position = position_jitter(width =.2), size = 2.5) +
  scale_y_continuous(name = "change in fluorescence (%)",
                     breaks = seq(-60, 30, 10),
                     limits = c(-70, 30)) +
  scale_x_discrete(limits = c("- 1", "CH5 1", "CH4 1", "CH3 1", "CH2 1", "GlcNAc 1", "GlcNAc-6P 1", "GlcN 1", "GlcN-6P 1", "Glc 1"),
                   name = "condition") +
  theme_bw() +
  theme(legend.position = "", 
        axis.text = element_text(size = 20, colour = line),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 50, hjust = 1))


#########
# Stats #
#########

# Linear model #

pb_f_change$Env <- factor(pb_f_change$Env, levels = c("- 1", "CH5 1", "CH4 1", "CH3 1", "CH2 1", "GlcNAc 1", "GlcNAc-6P 1", "GlcN 1", "GlcN-6P 1", "Glc 1"))
pb_f_change$pct_change_log <- log(pb_f_change$pct_change_stats)
a = pb_f_change

s1 <- summary(lm(pct_change_log ~ Env, data = a))