#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(stringr)
require(cowplot)

plot_enrich_heatmap <- function(norm_enrich_table, WTresibox){
  textsize <- 7
  p <-  ggplot() +
          geom_tile(data=norm_enrich_table,aes(x=resi,y=aa,fill=average_RF)) +
          labs("Dose (mg)") +
          scale_fill_gradientn(colours=c("blue", "white", "red"),
                limits=c(-5.5,1.5),
                values=rescale(c(-5.5, 0, 1.5)),
                breaks=c(-5,-4,-3,-2,-1,0,1),
                labels=c('-5','-4','-3','-2','-1','0','1'),
                guide="colorbar",
                na.value="grey") +
          theme_cowplot(12) +
          theme(axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title=element_text(size=textsize,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1)) +
          guides(fill = guide_colorbar(title.theme=element_text(size=textsize,face="bold",colour='black',hjust=0.5),
                                       label.theme=element_text(size=textsize,face="bold",colour='black'),
                                       frame.colour="black",
                                       frame.linewidth = 1,
                                       ticks = TRUE,
                                       ticks.colour = "black",
                                       barwidth = 0.5, barheight = 6, title="Relative\naffinity")) +
          geom_point(data=WTresibox, aes(x=x, y=y), color='gray', size=0.5) +
          xlab("Position") +
          ylab("Amino acid")
  ggsave('graph/norm_affinity_heatmap.png',p,width=7, height=2.2, dpi=1200)
  }

aa_level <- rev(c('E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W'))
WTresibox  <- read_tsv('script/WT_heatmap.tsv')
Analyzed_onemut_df  <- read_csv('results/Analyzed_onemut.csv')

arranged_onemut_df <- Analyzed_onemut_df %>%
                       filter(NEP_Mutation!='WT') %>%
                       mutate(Pos=str_sub(NEP_Mutation,2,-2)) %>%
                       mutate(resi=str_sub(NEP_Mutation,1,-2)) %>%
                       mutate(aa=str_sub(NEP_Mutation,-1,-1)) %>%
                       filter(aa %in% aa_level) %>%
                       mutate(aa=factor(aa,levels=aa_level)) %>%
                       mutate(Pos=factor(Pos,levels=as.character(seq(2,113)))) %>%
                       arrange(Pos) %>%
                       mutate(resi=factor(resi,levels=unique(resi))) %>%
                       mutate(Pos=as.numeric(as.character(Pos))) %>%
                       select(NEP_Mutation, resi, Pos, aa, average_RF)
plot_enrich_heatmap(arranged_onemut_df, WTresibox)
