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
  p1 <-  ggplot() +
          geom_tile(data=norm_enrich_table,aes(x=resi,y=aa,fill=output1_RF)) +
          labs("Dose (mg)") +
          scale_fill_gradientn(colours=c("blue","blue", "white", "red"),
                limits=c(-2,1.5),
                values=rescale(c(-2,-0.3, 0, 1.5)),
                breaks=c(-2,-1,0,1),
                labels=c('-2','-1','0','1'),
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
                                       barwidth = 0.5, barheight = 6, title="Relative\n Fitness")) +
          geom_point(data=WTresibox, aes(x=x, y=y), color='black', size=0.5) +
          xlab("Position") +
          ylab("Amino acid")
  ggsave('graph/NEP1_RF_heatmap.png',p1,width=12, height=2.7, dpi=1200)

  p2 <-  ggplot() +
          geom_tile(data=norm_enrich_table,aes(x=resi,y=aa,fill=output2_RF)) +
          labs("Dose (mg)") +
          scale_fill_gradientn(colours=c("blue","blue", "white", "red"),
                limits=c(-2,1.5),
                values=rescale(c(-2,-0.3, 0, 1.5)),
                breaks=c(-2,-1,0,1),
                labels=c('-2','-1','0','1'),
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
                                       barwidth = 0.5, barheight = 6, title="Relative\n Fitness")) +
          geom_point(data=WTresibox, aes(x=x, y=y), color='black', size=0.5) +
          xlab("Position") +
          ylab("Amino acid")
  ggsave('graph/NEP2_RF_heatmap.png',p2,width=12, height=2.7, dpi=1200)


  p_min <-  ggplot() +
          geom_tile(data=norm_enrich_table,aes(x=resi,y=aa,fill=min_RF)) +
          labs("Dose (mg)") +
          scale_fill_gradientn(colours=c("blue","blue", "white", "red"),
                limits=c(-2,1.5),
                values=rescale(c(-2,-0.3, 0, 1.5)),
                breaks=c(-2,-1,0,1),
                labels=c('-2','-1','0','1'),
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
                                       barwidth = 0.5, barheight = 6, title="Relative\n Fitness")) +
          geom_point(data=WTresibox, aes(x=x, y=y), color='black', size=0.5) +
          xlab("Position") +
          ylab("Amino acid")
  ggsave('graph/NEP_min_RF_heatmap.png',p_min,width=12, height=2.7, dpi=1200)
  p_average <-  ggplot() +
          geom_tile(data=norm_enrich_table,aes(x=resi,y=aa,fill=average_RF)) +
          labs("Dose (mg)") +
          scale_fill_gradientn(colours=c("blue","blue", "white", "red"),
                limits=c(-2,1.5),
                values=rescale(c(-2,-0.3, 0, 1.5)),
                breaks=c(-2,-1,0,1),
                labels=c('-2','-1','0','1'),
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
                                       barwidth = 0.5, barheight = 6, title="Relative\n Fitness")) +
          geom_point(data=WTresibox, aes(x=x, y=y), color='black', size=0.5) +
          xlab("Position") +
          ylab("Amino acid")
  ggsave('graph/NEP_average_RF_heatmap.png',p_average,width=12, height=2.7, dpi=1200)
  }

aa_level <- rev(c('E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','_'))

Analyzed_onemut_df  <- read_csv('results/NEP_onemut_with_low_input.csv')

arranged_onemut_df <- Analyzed_onemut_df %>%
                       filter(NEP_Mutation!='WT') %>%
                       mutate(Pos=str_sub(NEP_Mutation,2,-2)) %>%
                       mutate(resi=str_sub(NEP_Mutation,1,-2)) %>%
                       mutate(aa=str_sub(NEP_Mutation,-1,-1)) %>%
                       filter(aa %in% aa_level) %>%
                       mutate(aa=factor(aa,levels=aa_level)) %>%
                       mutate(Pos=factor(Pos,levels=as.character(seq(1,121)))) %>%
                       arrange(Pos) %>%
                       mutate(resi=factor(resi,levels=unique(resi))) %>%
                       mutate(Pos=as.numeric(as.character(Pos))) %>%
                       select(NEP_Mutation, resi, Pos, aa, output1_RF,output2_RF,average_RF,min_RF)

WTresibox  <- arranged_onemut_df %>%
  select(resi,Pos) %>%
  unique() %>%
  mutate(WT_resi=str_sub(resi,1,1)) %>%
  mutate(x=seq(1,121)) %>%
  mutate(y=match(WT_resi,aa_level)) %>%
  select(resi,WT_resi,x, y)

plot_enrich_heatmap(arranged_onemut_df, WTresibox)
