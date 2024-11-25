# Setting up... #

rm(list = ls())

library(ggplot2)
library(dplyr)
library(readxl)
line <- "black"

# Update path to correct directory #

pb <- read_excel("chitinase_assay_data_v1.xlsx")

pb_grouped <- group_by(pb, Name, Substrate)
head(pb_grouped)

pb_grouped$Name = factor(pb_grouped$Name, levels = c("CLP","CLP(d)", "HEXO"))
pb_grouped$Substrate = factor(pb_grouped$Substrate, levels = c("1mer", "2mer", "3mer"))

p1 <- ggplot(pb_grouped, aes(x = Name, y = Av, fill = Name)) +
  geom_bar(stat = "identity", position = position_dodge(), colour = line, size = 1) +
  geom_errorbar(aes(ymin = Av - SE, ymax = Av + SE, width = .2),
                position = position_dodge(.9), width = 0.15, size = 1, colour = line) +
  scale_y_continuous(name = "% activity", 
                     breaks = seq(0, 10, 2),
                     limits = c(0, 10)) +
  scale_x_discrete(name = "") +
  facet_grid(. ~ Substrate) +
  theme_classic()

p2 <- ggplot(pb_grouped, aes(x = Name, y = Activity)) +
  geom_boxplot(colour = line, outlier.shape = NA, size =1) + 
  geom_point(aes(x = Name, y = Activity, colour = Name),
             position = position_jitter(width =.2), size = 2.5) +
  scale_y_continuous(name = "% activity", 
                     breaks = seq(0, 10, 2),
                     limits = c(-1, 10)) +
  scale_x_discrete(name = "") +
  facet_grid(. ~ Substrate) +
  theme_classic()


pb_1 <- pb_grouped[which(pb_grouped$Substrate == "1mer"), ]
pb_1$pctactivity <- pb_1$Activity / 100
pb_1$Name <- factor(pb_1$Name, levels = c("CLP", "CLP(d)", "HEXO"))

q1 <- glm(pctactivity ~ Name, family = quasibinomial(link = "logit"), data = pb_1)
summary(q1)


pb_2 <- pb_grouped[which(pb_grouped$Substrate == "2mer"), ]
pb_2$pctactivity <- pb_2$Activity / 100
pb_2$Name <- factor(pb_2$Name, levels = c("CLP", "CLP(d)", "HEXO"))

q2 <- glm(pctactivity ~ Name, family = quasibinomial(link = "logit"), data = pb_2)
summary(q2)


pb_3 <- pb_grouped[which(pb_grouped$Substrate == "3mer"), ]
pb_3$pctactivity <- pb_3$Activity / 100
pb_3$Name <- factor(pb_3$Name, levels = c("CLP", "CLP(d)", "HEXO"))

q3 <- glm(pctactivity ~ Name, family = quasibinomial(link = "logit"), data = pb_3)
summary(q3)