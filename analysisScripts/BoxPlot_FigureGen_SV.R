library(ggplot2)
library(dplyr)
#Clean up everything
rm(list = ls())

setwd('~/RNAILabs/users/ej276743/spec_loop/')
data <- read.csv("Induced_Unbiased_renameligands.csv")

names(data)[names(data) == 'Temperature..K.'] <- 'Temperature'

data$Group <- paste(data$Structure,data$Temperature)

Shape <- c("2KHY 300K" = 17, "2KHY 350K" = 19, "4LCK 300K" = 16, "4LCK 350K" = 15)

average_scores <- data %>%
  group_by(Name) %>%
  summarise(average_score = mean(S))

data$Name <- factor(data$Name, levels = average_scores$Name[order(average_scores$average_score, decreasing = TRUE)])

#png(filename = "Biased_docking_24.png")

# Create the plot
ggplot(data, aes(x = S, y = Name, group = Name)) +
  geom_boxplot(fill = "azure", outlier.shape = NA, width=.5)+
  geom_point(aes(shape = Group,
                 color = Group),
             position = position_jitterdodge(jitter.width = .5, dodge.width = 1),
             size = 2.5, alpha=.7) +
  scale_shape_manual(values = Shape) +
  scale_color_manual(values = c("#D55E00", "#009E73", "#E2DA6C", "#56DDFD")) +
  labs(title = "Docking Score of PKZ18 Analogues when Docked to 2KHY and 4LCK",
       y = "PKZ18 Analogue", x = "Binding Affinity / S Score (kcal / mol)",
       shape = "Structure and Temperature") +
 # theme_minimal()+
  scale_x_continuous(breaks = seq(-5, -8, by = -.5))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, face = "bold", size = 12), 
        axis.text.x = element_text(face = "bold", size = 12), 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),  # # Remove major grid lines
        panel.grid.major.y = element_line(linewidth = 0.25),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="none") +
  xlim(-5.5, -8.5)

t=  theme(axis.text.x = element_text(size = 15, angle = 90, face="bold",vjust=0.25, hjust=1),
          axis.text.y = element_text(size = 16,face="bold"), #, vjust=-0.5),
          legend.text = element_text(size = 16, face = "bold"),
          #   legend.title = element_blank(),
          plot.title=element_text(hjust=0.5,face="bold", size=14),
          axis.title=element_text(size=15, face="bold"),
          #axis.text=element_text(size=12,face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="none")
#dev.off()
