# Title     : to plot the pearson correlation
# Objective :
# Created by: yiquan
# Created on: 1/26/21
library(ggplot2)

site_count_df <- read.csv('results/site_count.csv',header=TRUE)
nep_df <- read.csv('results/one_mut.csv',header=TRUE)

site_count_plot <- ggplot(site_count_df, aes(1:nrow(site_count_df))) +
  geom_line(aes(y = input1), color = "darkred") +
  geom_line(aes(y = input2), color="red") +
  geom_line(aes(y = output1), color = "steelblue") +
  geom_line(aes(y = output2), color="lightblue") +
  labs(x = "position", y = "sample frequency")+
  theme_linedraw() +
  theme(plot.title=element_blank(), legend.title=element_blank())
ggsave("graph/site_mut_freq.png",site_count_plot, width = 16, height = 8, units = "cm", dpi = 300)
input_plot <- ggplot(nep_df, aes(input1, input2))+
  geom_point()+
  labs(x = "input1", y = "input2")+
  xlim(1, 150) +ylim(1, 150)+ theme_linedraw()+
  theme(plot.title=element_blank(),
        legend.title=element_blank()
  )

output_plot <- ggplot(nep_df, aes(output1, output2))+ geom_point() + labs(x = "output1", y = "output2")+ xlim(1, 150) +ylim(1, 150)+
  theme_linedraw() +
  theme(plot.title=element_blank(),
        legend.title=element_blank()
        )
ggsave("graph/Correlation_input.png",input_plot, width = 16, height = 14, units = "cm", dpi = 300)
ggsave("graph/Correlation_output.png",output_plot, width = 16, height = 14, units = "cm", dpi = 300)