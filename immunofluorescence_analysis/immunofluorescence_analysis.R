# Clear environment # 

rm(list = ls())

library(ggplot2)
library(dplyr)
library(readxl)
line <- "black"

# Update path to correct directory # 

scan <- read_excel("scanR_all_merged_output_v1.xlsx")

# Format data # 

scan %>% group_by(Condition2, Channel, CellID, Population)
scanA <- scan[which(scan$Channel == "Algae"), ]
scanC <- scan[which(scan$Channel == "CLP"), ]


# ALGAL FLUORESCENCE #

p1 <- ggplot(subset(scanA, Population %in% c("P", "PL", "L")),
             aes(x = Condition2, y = Mean)) +
  geom_boxplot(colour = line, outlier.shape = NA, size =1) + 
  geom_point(aes(x = Condition2, y = Mean, colour = Population),
             position = position_jitter(width =.2), size = 2.5) +
  scale_y_continuous(name = "mean algal autofluorescence",
                     breaks = seq(0, 30000, 5000),
                     limits = c(5000, 30000)) +
  scale_x_discrete(limits = c("Untreated", "Cyclo_63"),
                   name = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 20, colour = line),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 50, hjust = 1))


p1.1 <- ggplot(subset(scanA, Population %in% c("P")),
             aes(x = Condition2, y = Mean)) +
  geom_boxplot(colour = line, outlier.shape = NA, size =1) + 
  geom_point(aes(x = Condition2, y = Mean, colour = Population),
             position = position_jitter(width =.2), size = 2.5) +
  scale_y_continuous(name = "mean algal autofluorescence",
                     breaks = seq(0, 30000, 5000),
                     limits = c(5000, 30000)) +
  scale_x_discrete(limits = c("Untreated", "Cyclo_63"),
                   name = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 20, colour = line),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 50, hjust = 1))


p1.2 <- ggplot(subset(scanA, Population %in% c("PL", "L")),
               aes(x = Condition2, y = Mean)) +
  geom_boxplot(colour = line, outlier.shape = NA, size =1) + 
  geom_point(aes(x = Condition2, y = Mean, colour = Population),
             position = position_jitter(width =.2), size = 2.5) +
  scale_y_continuous(name = "mean algal autofluorescence",
                     breaks = seq(0, 30000, 5000),
                     limits = c(5000, 30000)) +
  scale_x_discrete(limits = c("Untreated", "Cyclo_63"),
                   name = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 20, colour = line),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 50, hjust = 1))


# CLP FLUORESCENCE #

p2 <- ggplot(subset(scanC, Population %in% c("P", "PL", "L")),
             aes(x = Condition2, y = Max)) +
  geom_boxplot(colour = line, outlier.shape = NA, size =1) + 
  geom_point(aes(x = Condition2, y = Max, colour = Population),
             position = position_jitter(width =.2), size = 2.5) +
  scale_y_continuous(name = "max CLP immunofluorescence",
                     breaks = seq(0, 2000000, 500000),
                     limits = c(0, 2250000)) +
  scale_x_discrete(limits = c("Untreated", "Cyclo_63"),
                   name = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 20, colour = line),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 50, hjust = 1))



p2.1 <- ggplot(subset(scanC, Population %in% c("P")),
               aes(x = Condition2, y = Max)) +
  geom_boxplot(colour = line, outlier.shape = NA, size =1) + 
  geom_point(aes(x = Condition2, y = Max, colour = Population),
             position = position_jitter(width =.2), size = 2.5) +
  scale_y_continuous(name = "max CLP immunofluorescence",
                     breaks = seq(0, 2000000, 500000),
                     limits = c(0, 2250000)) +
  scale_x_discrete(limits = c("Untreated", "Cyclo_63"),
                   name = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 20, colour = line),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 50, hjust = 1))


p2.2 <- ggplot(subset(scanC, Population %in% c("PL", "L")),
               aes(x = Condition2, y = Max)) +
  geom_boxplot(colour = line, outlier.shape = NA, size =1) + 
  geom_point(aes(x = Condition2, y = Max, colour = Population),
             position = position_jitter(width =.2), size = 2.5) +
  scale_y_continuous(name = "max CLP immunofluorescence",
                     breaks = seq(0, 2000000, 500000),
                     limits = c(0, 2250000)) +
  scale_x_discrete(limits = c("Untreated", "Cyclo_63"),
                   name = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 20, colour = line),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 50, hjust = 1))



#########
# Stats #
#########

# Prepare datasets # 

scanA_PPLL <- scanA[which(scanA$Population == "P" |
                            scanA$Population == "PL" |
                            scanA$Population == "L"), ]

scanA_P <- scanA[which(scanA$Population == "P"), ]

scanA_PLL <- scanA[which(scanA$Population == "PL" |
                           scanA$Population == "L"), ]


scanC_PPLL <- scanC[which(scanC$Population == "P" |
                            scanC$Population == "PL" |
                            scanC$Population == "L"), ]

scanC_P <- scanC[which(scanC$Population == "P"), ]

scanC_PLL <- scanC[which(scanC$Population == "PL" |
                           scanC$Population == "L"), ]


# Set order, so that 'untreated' is the first comparison #

scanA$Condition2 <- factor(scanA$Condition2, levels = c("Untreated", "Cyclo_63"))
scanC$Condition2 <- factor(scanC$Condition2, levels = c("Untreated", "Cyclo_63"))

a = scanA_PPLL
b = scanA_P
c = scanA_PLL

d = scanC_PPLL
e = scanC_P
f = scanC_PLL


# Stats summaries #

s1 <- summary(lm(Mean ~ Condition2, data = a))
s2 <- summary(lm(Mean ~ Condition2, data = b))
s3 <- summary(lm(Mean ~ Condition2, data = c))

s4 <- summary(lm(Max ~ Condition2, data = d))
s5 <- summary(lm(Max ~ Condition2, data = e))
s6 <- summary(lm(Max ~ Condition2, data = f))