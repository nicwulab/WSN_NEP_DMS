# Title     : to plot the pearson correlation
# Objective :
# Created by: yiquan
# Created on: 1/26/21
library(ggplot2)


nep_df <- read.csv('results/nep_mut.csv',header=TRUE)

input_plot <- ggplot(nep_df, aes(input1, input2))+ geom_point() + labs(x = "input1", y = "input2")+ xlim(1, 300) +ylim(1, 300)+
  theme_linedraw() +
  theme(plot.title=element_blank(),
        legend.title=element_blank()
        )

cor.test(nep_df$input1, nep_df$input2, method = "pearson", conf.level = 0.95)
cor.test(nep_df$output1, nep_df$output2, method = "pearson", conf.level = 0.95)
output_plot <- ggplot(nep_df, aes(output1, output2))+ geom_point() + labs(x = "input1", y = "input2")+ xlim(1, 300) +ylim(1, 300)+
  theme_linedraw() +
  theme(plot.title=element_blank(),
        legend.title=element_blank()
        )
ggsave("graph/Correlation_input.png",input_plot,
    width = 16, height = 14, units = "cm", dpi = 300)
ggsave("graph/Correlation_output.png",output_plot,
    width = 16, height = 14, units = "cm", dpi = 300)