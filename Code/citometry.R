library(tidyverse)
library(viridis)

fcm.data <- read.csv("/Users/Artemis/Documents/GitHub/Remember-R/Data/abundancia_eventos_citometria_pch23.csv", header = TRUE, sep = ";", dec = ",", skip = 0, row.names = 1)
fcm.data <- fcm.data[,-c(4,10,13,16,19,22,25,28)]
fcm <- as.data.frame(t(fcm.data))
fcm <- fcm[,-6]
fcm <- add_column(fcm, 
                  Time = c(rep("t0", 3), rep("t1", 9), rep("t3", 8)),
                  Condition = c(rep("In.Cond", 3), rep("Control", 3), rep("Cu", 2), rep("Fe", 2), rep("Cu:Fe", 2), rep("Control", 2), rep("Cu", 2), rep("Fe", 2), rep("Cu:Fe", 2)),
                  .before = "PPE", 
                  .name_repair = "minimal")
fcm$Time <- factor(fcm$Time, levels = unique(fcm$Time))
fcm$Condition <- factor(fcm$Condition, levels = unique(fcm$Condition))
long <- gather(fcm[,c(1:3,5)], "Event", "Count", 3:4)
long$Event <- factor(long$Event, levels = unique(long$Event))

t0 <- ggplot(long, aes(x = Condition, 
                       y = Count, 
                       fill = Event)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("#16DB93", "#F94144")) + 
  hrbrthemes::theme_ipsum(base_size = 23, axis_col = "#000000", axis_title_size = 30, plot_title_size = 40, axis_text_size = 20) + 
  theme() + 
  facet_wrap(~Time)
t0
ggsave(plot = t0, 
       filename = "t0.tiff",
       path = "Products/", 
       width = 21, 
       height = 13, 
       dpi = 600)

t1 <- long %>% 
  filter(Time == "t1") %>% 
  ggplot(aes(x = Condition, 
             y = Count, 
             fill = Event)) + 
  geom_violin() +
  geom_boxplot() + 
  scale_fill_viridis(discrete = TRUE, option = "plasma", direction = -1) + 
  hrbrthemes::theme_ipsum(base_size = 20, axis_col = "#000000", axis_title_size = 30, plot_title_size = 40, axis_text_size = 20) + 
  theme()

ggsave(plot = t1, 
       filename = "t1.tiff",
       path = "Products/", 
       width = 21, 
       height = 13, 
       dpi = 600)

t3 <- long %>% 
  filter(Time == "t3") %>% 
  ggplot(aes(x = Condition, 
             y = Count, 
             fill = Event)) + 
  geom_boxplot() + 
  scale_fill_viridis(discrete = TRUE, option = "plasma", direction = -1) + 
  hrbrthemes::theme_ipsum(base_size = 20, axis_col = "#000000", axis_title_size = 30, plot_title_size = 40, axis_text_size = 20) + 
  theme()

ggsave(plot = t3, 
       filename = "t3.tiff",
       path = "Products/", 
       width = 21, 
       height = 13, 
       dpi = 600)

chao1 %>% 
  rename("Chao1" = "S.chao1") %>% 
  ggplot(aes(x = Site, y = Chao1, fill = Site)) + 
  geom_violin(width = 1.2) + 
  geom_boxplot(width = 0.2, color = "grey", alpha = 0.2) + 
  scale_fill_viridis(discrete = TRUE) + 
  geom_jitter(color = "black", size = 0.5, alpha = 0.5) +
  theme_ipsum() + 
  theme(legend.position = "none") + 
  ggtitle("Chao1 richness")



long.fcm <- gather(rltv.fcm, "Event", "Count", c(2,3,4,5,6,8))
fcmplot <- ggplot(long.fcm, aes(x = Sample, 
                                y = Count, 
                                fill = Event)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Condition", 
       y = "Number of events", 
       title = "Number of events per bottles (conditions)") + 
  scale_fill_viridis(discrete = TRUE, begin = 0, end = 1, direction = 1, option = "inferno") + 
  hrbrthemes::theme_ipsum(base_size = 25, axis_col = "#000000", axis_title_size = 35, plot_title_size = 40) + 
  theme(legend.position = "bottom", legend.direction = "horizontal")
legend <- cowplot::get_legend(fcmplot)
grapho <- fcmplot + theme(legend.position = "none")
ggsave(plot = grapho, 
       filename = "number_events_conditions.tiff",
       path = "Products/", 
       width = 31, 
       height = 19, 
       dpi = 600)
ggsave(ggpubr::as_ggplot(legend), 
       filename = "legend.tiff", 
       path = "Products/", 
       width = 21,
       height = 10,
       dpi = 900)
