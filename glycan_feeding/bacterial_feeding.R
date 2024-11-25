# Clear environment # 

rm(list = ls())

# Setting up... #

library(ggplot2)
library(dplyr)
library(reshape2)
library(readxl)
line <- "black"


# Open data and platemaps #

pb <- read.delim("bacterial_feeding_v1.txt")

plates <- read_excel("platemap_bacterial_feeding_v1.xlsx")

# Annotate plate by adding platemap values through merging, IGNORE WARNING #

plates$MEASUREMENT.SET.ID <- as.character(plates$MEASUREMENT.SET.ID)
pb <- inner_join(pb, plates, by = c("MEASUREMENT.SET.ID", "Well.Name"))
head(pb)

# Remove missing values # 

pb_clean <- pb[which(pb$Cell..Average.Intensity_Average..Paramecium. != ""), ]
pb_clean <- pb[which(pb$MEASUREMENT.SET.ID != ""), ]
pb_clean <- pb[which(pb$Condition != "NA"), ]

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
pb_gated <- pb_gated[which(pb_gated$Cell..Average.Intensity_Average..Paramecium. < 300), ]


###################
# Processing data #
###################

# Count per day # 

pb_count <- pb_gated %>% group_by(Well.Name, MEASUREMENT.SET.ID, Condition, Gram, Treatment, Day, Rep) %>% tally()
head(pb_count)

# Group count data by condition and day, then calculate average #

pb_average <- pb_count %>% 
  group_by(Condition, Gram, Treatment, Day) %>% 
  summarize_at(vars(n), list(av_count = mean, sd_count = sd))
head(pb_average)

# Group raw data (for fluorescence) by condition and day, then calculate average #

pb_fl_average <- pb_gated %>% 
  group_by(Condition, Gram, Treatment, Day) %>% 
  summarize_at(vars(Cell..Average.Intensity_Average..Paramecium.), list(av_fluo = mean, sd_fluo = sd))
head(pb_average)

# Merge average values for count and fluorescence data #

pb_total <- merge(pb_average, pb_fl_average, by = c("Condition", "Treatment", "Day"))
head(pb_total)


### Average fluorescence per well ####

pb_f <- pb_gated %>% group_by(Condition, Gram, Treatment, Day, Rep, Well.Name) %>% 
  summarize_at(vars(Cell..Average.Intensity_Average..Paramecium.), list(Mean = mean))



########################
# Change compared to 0 #
########################

# Change in cell number #

pb_1 <- pb_count %>% group_by(Condition, Gram, Treatment, Rep, Day, n) %>% tally()
pb_2 <- pb_1[which(pb_1$Day == "0"), ]

pb_change <- merge(pb_2, pb_1, by=c("Treatment", "Gram", "Condition", "Rep"))
pb_change$pct_change <- ((pb_change$"n.y" - pb_change$"n.x")/pb_change$"n.x") * 100
# Add lowest negative value for stats #
pb_change$pct_change_stats <- pb_change$pct_change + 100
head(pb_change)

# Aesthetics #

pb_change$Condition.y <- as.factor(pb_change$Day.y)
pb_change$pct_change <- as.numeric(pb_change$pct_change)


# Change in fluorescence #

pb_f1 <- pb_gated %>% group_by(Condition, Gram, Treatment, Rep, Day) %>% 
  summarize(Mean = mean(Cell..Average.Intensity_Average..Paramecium.))
pb_f2 <- pb_f1[which(pb_f1$Day == "0"), ]

pb_f_change <- merge(pb_f2, pb_f1, by=c("Treatment", "Gram", "Condition", "Rep"))
pb_f_change$pct_change <- ((pb_f_change$"Mean.y" - pb_f_change$"Mean.x")/pb_f_change$"Mean.x") * 100
# Add lowest negative value for stats #
pb_f_change$pct_change_stats <- pb_f_change$pct_change + 100
head(pb_f_change)

# Aesthetics #

pb_f_change$Condition.y <- as.factor(pb_f_change$Day.y)
pb_f_change$pct_change <- as.numeric(pb_f_change$pct_change)


#########
# Plots #
#########

pb_change_1 <- pb_change[which(pb_change$Day.y == "1"), ]
pb_f_change_1 <- pb_f_change[which(pb_f_change$Day.y == "1"), ]

pb_change_2 <- pb_change[which(pb_change$Day.y == "2"), ]
pb_f_change_2 <- pb_f_change[which(pb_f_change$Day.y == "2"), ]


p1 <- ggplot(subset(pb_f_change_1, Treatment %in% c("heatkilled")),
             aes(x = Gram, y = pct_change)) +
  geom_boxplot(colour = line, outlier.shape = NA, size =1) + 
  geom_point(aes(x = Gram, y = pct_change, colour = Condition),
             position = position_jitter(width =.2), size = 2.5) +
  scale_y_continuous(name = "change in per host fluorescence (%)",
                     breaks = seq(-20, 20, 10),
                     limits = c(-25, 25)) +
  scale_x_discrete(limits = c("-", "negative", "positive"),
                   name = "condition") +
  theme_bw() +
  theme(text = element_text(size=25, colour = line),
        legend.position = "", 
        plot.title = element_text(size = 18, hjust = 0.5), 
        axis.text = element_text(size = 18, colour = line),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(angle = 70, hjust = 1))


p2 <- ggplot(subset(pb_change_1, Treatment %in% c("heatkilled")),
             aes(x = Condition, y = pct_change)) +
  geom_boxplot(colour = line, outlier.shape = NA, size =1) + 
  geom_point(aes(x = Condition, y = pct_change, colour = Condition),
             position = position_jitter(width =.2), size = 2.5) +
  scale_y_continuous(name = "change in host cell number (%)",
                     breaks = seq(-100, 100, 50),
                    limits = c(-100, 100)) +
  scale_x_discrete(limits = c("-", "E. coli", "Klebsiella", "Citrobacter", "Serratia", 
                              "Acinetobacter", "Pseudomonas", "Enterococcus", "Lactobacillus",
                              "Staphylococcus", "Bacillus", "Bifidobacterium", "Leucobacter"),
                   name = "condition") +
  theme_bw() +
  theme(text = element_text(size=25, colour = line),
        legend.position = "", 
        plot.title = element_text(size = 18, hjust = 0.5), 
        axis.text = element_text(size = 18, colour = line),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(angle = 70, hjust = 1))


p3 <- ggplot(subset(pb_f_change_1, Treatment %in% c("heatkilled")),
              aes(x = Condition, y = pct_change)) +
  geom_boxplot(colour = line, outlier.shape = NA, size =1) + 
  geom_point(aes(x = Condition, y = pct_change, colour = Condition),
             position = position_jitter(width =.2), size = 2.5) +
  scale_y_continuous(name = "change in per host fluorescence (%)",
                     breaks = seq(-20, 20, 10),
                    limits = c(-20, 20)) +
  scale_x_discrete(limits = c("-", "E. coli", "Klebsiella", "Citrobacter", "Serratia", 
                              "Acinetobacter", "Pseudomonas", "Enterococcus", "Lactobacillus",
                              "Staphylococcus", "Bacillus", "Bifidobacterium", "Leucobacter"),
                   name = "condition") +
  theme_bw() +
  theme(text = element_text(size=25, colour = line),
        legend.position = "", 
        plot.title = element_text(size = 18, hjust = 0.5), 
        axis.text = element_text(size = 18, colour = line),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(angle = 70, hjust = 1))


p4 <- ggplot(subset(pb_change_1, Treatment %in% c("heatkilled")),
              aes(x = Condition, y = pct_change)) +
  geom_boxplot(colour = line, outlier.shape = NA, size =1) + 
  geom_point(aes(x = Condition, y = pct_change, colour = Condition),
             position = position_jitter(width =.2), size = 2.5) +
  scale_y_continuous(name = "change host cell number (%)",
                     breaks = seq(-100, 100, 50),
                    limits = c(-100, 100)) +
  scale_x_discrete(limits = c("-", "E. coli", "E. coli + GlcNAc", "Klebsiella", "Klebsiella + GlcNAc", 
                              "Citrobacter", "Citrobacter + GlcNAc", "Serratia", "Serratia + GlcNAc",
                              "Acinetobacter", "Acinetobacter + GlcNAc", "Pseudomonas", "Pseudomonas + GlcNAc",
                              "Enterococcus", "Enterococcus + GlcNAc", "Staphylococcus", "Staphylococcus + GlcNAc",
                              "Bacillus", "Bacillus + GlcNAc"),
                   name = "condition") +
  theme_bw() +
  theme(text = element_text(size=25, colour = line),
        legend.position = "", 
        plot.title = element_text(size = 18, hjust = 0.5), 
        axis.text = element_text(size = 18, colour = line),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(angle = 70, hjust = 1))


p5 <- ggplot(subset(pb_f_change_1, Treatment %in% c("heatkilled")),
              aes(x = Condition, y = pct_change)) +
  geom_boxplot(colour = line, outlier.shape = NA, size =1) + 
  geom_point(aes(x = Condition, y = pct_change, colour = Condition),
             position = position_jitter(width =.2), size = 2.5) +
  scale_y_continuous(name = "change in per host fluorescence (%)",
                     breaks = seq(-20, 20, 10),
                    limits = c(-20, 20)) +
  scale_x_discrete(limits = c("-", "E. coli", "E. coli + GlcNAc", "Klebsiella", "Klebsiella + GlcNAc", 
                              "Citrobacter", "Citrobacter + GlcNAc", "Serratia", "Serratia + GlcNAc",
                              "Acinetobacter", "Acinetobacter + GlcNAc", "Pseudomonas", "Pseudomonas + GlcNAc",
                              "Enterococcus", "Enterococcus + GlcNAc", "Staphylococcus", "Staphylococcus + GlcNAc",
                              "Bacillus", "Bacillus + GlcNAc"),
                   name = "condition") +
  theme_bw() +
  theme(text = element_text(size=25, colour = line),
        legend.position = "", 
        plot.title = element_text(size = 18, hjust = 0.5), 
        axis.text = element_text(size = 18, colour = line),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(angle = 70, hjust = 1))


#########
# Stats #
#########

# Log Linear model #

pb_f_change_1h <- pb_f_change_1[which(pb_f_change_1$Treatment == "heatkilled"), ]
pb_f_change_1h$pct_change_log <- log(pb_f_change_1h$pct_change_stats)

pb_f_change_1h$Gram <- factor(pb_f_change_1h$Gram, levels = c("negative", "positive", "-"))
pb_f_change_1h$pct_change_log <- log(pb_f_change_1h$pct_change_stats)

a = pb_f_change_1h

s1 <- summary(lm(pct_change_log ~ Gram, data = a))


# Generalized linear models for supplementary data #

pb_change_1h <- pb_change_1[which(pb_change_1$Treatment == "heatkilled"), ]
pb_change_1h$pct_change_log <- log(pb_change_1h$pct_change_stats)

pb_change_1h$Condition <- factor(pb_change_1h$Condition, levels = c("Bacillus", "Staphylococcus", "Enterococcus", "Pseudomonas", "Acinetobacter", "Serratia", "Citrobacter", "Klebsiella", "E. coli", "-", "E. coli + GlcNAc", "Klebsiella + GlcNAc", 
                                                                    "Citrobacter + GlcNAc", "Serratia + GlcNAc", 
                                                                    "Acinetobacter + GlcNAc", "Pseudomonas + GlcNAc",
                                                                    "Enterococcus + GlcNAc", "Lactobacillus", 
                                                                    "Staphylococcus + GlcNAc", "Bacillus + GlcNAc", "Bifidobacterium", 
                                                                    "Leucobacter"))

pb_f_change_1h$Condition <- factor(pb_f_change_1h$Condition, levels = c("Bacillus", "Staphylococcus", "Enterococcus", "Pseudomonas", "Acinetobacter", "Serratia", "Citrobacter", "Klebsiella", "E. coli", "-", "E. coli + GlcNAc", "Klebsiella + GlcNAc", 
                                                                        "Citrobacter + GlcNAc", "Serratia + GlcNAc", 
                                                                        "Acinetobacter + GlcNAc", "Pseudomonas + GlcNAc",
                                                                        "Enterococcus + GlcNAc", "Lactobacillus", 
                                                                        "Staphylococcus + GlcNAc", "Bacillus + GlcNAc", "Bifidobacterium", 
                                                                        "Leucobacter"))


q1 <- glm(pct_change_log ~ Condition, data = subset(pb_change_1h, Condition %in% c("-", "E. coli", "E. coli + GlcNAc", "Klebsiella", "Klebsiella + GlcNAc", 
                                                                                                              "Citrobacter", "Citrobacter + GlcNAc", "Serratia", "Serratia + GlcNAc", 
                                                                                                              "Acinetobacter", "Acinetobacter + GlcNAc", "Pseudomonas", "Pseudomonas + GlcNAc",
                                                                                                              "Enterococcus", "Enterococcus + GlcNAc", "Lactobacillus", "Staphylococcus", 
                                                                                                              "Staphylococcus + GlcNAc", "Bacillus", "Bacillus + GlcNAc", "Bifidobacterium", 
                                                                                                              "Leucobacter")))

summary(q1)


q2 <- glm(pct_change_log ~ Condition, data = subset(pb_f_change_1h, Condition %in% c("-", "E. coli", "E. coli + GlcNAc", "Klebsiella", "Klebsiella + GlcNAc", 
                                                                                                                "Citrobacter", "Citrobacter + GlcNAc", "Serratia", "Serratia + GlcNAc", 
                                                                                                                "Acinetobacter", "Acinetobacter + GlcNAc", "Pseudomonas", "Pseudomonas + GlcNAc",
                                                                                                                "Enterococcus", "Enterococcus + GlcNAc", "Lactobacillus", "Staphylococcus", 
                                                                                                                "Staphylococcus + GlcNAc", "Bacillus", "Bacillus + GlcNAc", "Bifidobacterium", 
                                                                                                                "Leucobacter")))

summary(q2)
