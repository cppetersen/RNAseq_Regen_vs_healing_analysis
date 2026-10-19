library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(data.table)
library(gplots)
library(stringr)
library(splitstackshape)
library(forcats)

##*##
#setwd("~/Documents/StudentData/Elias/RNAseq_analysis/analysis032126/")
head<-read.table("~/Documents/StudentData/Elias/RNAseq_analysis/analysis032126/ddv6_HeadR_031225.tsv",sep="\t", header=FALSE)
tail<-read.table("~/Documents/StudentData/Elias/RNAseq_analysis/analysis032126/ddv6_TailR_031225.tsv",sep="\t", header=FALSE)
heal<-read.table("~/Documents/StudentData/Elias/RNAseq_analysis/analysis032126/ddv6_Heal_031225.tsv",sep="\t", header=FALSE)

colnames(head)<-c("Contig","experiment","Experiment.2","time","log2FC")
colnames(tail)<-c("Contig","experiment","Experiment.2","time","log2FC")
colnames(heal)<-c("Contig","experiment","Experiment.2","time","log2FC")

setwd("~/Documents/StudentData/Elias/RNAseq_analysis/analysis101426_RNAseq_github_test/")

head.rename<-separate_wider_delim(head,time,delim=" ",names=c("hours","tag"))
tail.rename<-separate_wider_delim(tail,time,delim=" ",names=c("hours","tag"))
heal.rename<-separate_wider_delim(heal,time,delim=" ",names=c("hours","tag"))


#head.rename$hours<-as.numeric(head.rename$hours)

head.wide<-pivot_wider(head.rename,names_from=hours,values_from=c(log2FC))
head.wide<-pivot_wider(head,names_from=time,values_from=c(log2FC))
head.wide<-select(head.wide,1,2,4,5,8,11,13,14,6,7,9,10,12)


tail.wide<-pivot_wider(tail.rename,names_from=hours,values_from=c(log2FC))
tail.wide<-pivot_wider(tail,names_from=time,values_from=c(log2FC))
tail.wide<-select(tail.wide,1,2,4,5,8,11,13,14,6,7,9,10,12)


heal.wide<-pivot_wider(heal.rename,names_from=hours,values_from=c(log2FC))
heal.wide<-pivot_wider(heal,names_from=time,values_from=c(log2FC))
heal.wide<-select(heal.wide,1,2,4,5,8,11,13,14,6,7,9,10,12)

head.tail<-left_join(head.wide,tail.wide, by="Contig")

head.tail.heal<-left_join(head.tail,heal.wide,by="Contig")

head.tail.heal.analysis<-head.tail.heal %>% mutate(diff.1h = `1 h.x` - `1 h`) %>% mutate(diff.2h = `2 h.x` - `2 h`) %>% mutate(diff.4h = `4 h.x` - `4 h`) %>% mutate(diff.6h = `6 h.x` - `6 h`) %>% mutate(diff.8h = `8 h.x` - `8 h`) %>% mutate(diff.12h = `12 h.x` - `12 h`) %>% mutate(diff.16h = `16 h.x` - `16 h`) %>% mutate(diff.24h = `24 h.x` - `24 h`) %>% mutate(diff.36h = `36 h.x` - `36 h`) %>% mutate(diff.48h = `48 h.x` - `48 h`)

head.tail.heal.analysis<- head.tail.heal.analysis %>% mutate(WI.head.late = (`12 h.x`+`16 h.x`+`24 h.x`)/3) %>% mutate(WI.tail.late = (`12 h.y`+`16 h.y`+`24 h.y`)/3) %>% mutate(WI.head.early = (`1 h.x`+`2 h.x`+`4 h.x`+`6 h.x`)/4)  %>% mutate(WI.tail.early = (`1 h.y`+`2 h.y`+`4 h.y`+`6 h.y`)/4) %>% mutate(diff.early = (`diff.1h`+`diff.2h`+`diff.4h`+`diff.6h`)/4) %>% mutate(diff.late = (`diff.12h`+`diff.16h`+`diff.24h`)/3) %>% mutate(diff.total = (`diff.2h`+`diff.4h`+`diff.6h`+`diff.8h`+`diff.12h`+`diff.16h`+`diff.24h`)/7)

all.data.long<-rbind(head,tail,heal)

all.data.long.delim<-separate_wider_delim(all.data.long,time,delim=" ",names=c("hours","tag"))


#Plots figure 1E
facet_titles <- c('dd_Smed_v6_28398_0_1'='wnt1','dd_Smed_v6_24180_0_1'='notum','dd_Smed_v6_14068_0_1'='nlg-1','dd_Smed_v6_20318_0_1'='equinox','dd_Smed_v6_9584_0_1'='follistatin','dd_Smed_v6_4154_0_1'='glypican-1')
select.genes2<-c('dd_Smed_v6_28398_0_1','dd_Smed_v6_24180_0_1','dd_Smed_v6_14068_0_1','dd_Smed_v6_20318_0_1','dd_Smed_v6_9584_0_1','dd_Smed_v6_4154_0_1')
all.data.long.delim %>% filter(Contig %in% select.genes2) %>% mutate(hours = as.numeric(hours),Contig = factor(Contig, levels = select.genes2)) %>% ggplot(aes(x=hours,y=log2FC,color=experiment, fill=experiment, group=experiment))+
  geom_point(size=2)+
  geom_line(linewidth=1)+
  facet_wrap(~Contig, scales="free_y",labeller = as_labeller(facet_titles))+
  theme_classic()+theme(legend.position = "none")+theme(axis.text= element_text(face = "bold"))
ggsave("plot2.gene.png")

facet_titles <- c('dd_Smed_v6_14711_0_1'='egr-3','dd_Smed_v6_7607_0_1'='inhibin-1','dd_Smed_v6_6053_0_1'='mex3-2','dd_Smed_v6_4392_0_1'='TRAF-1','dd_Smed_v6_5999_0_1'='CPO','dd_Smed_v6_9519_0_1'='DAT')
select.genes2<-c('dd_Smed_v6_14711_0_1','dd_Smed_v6_7607_0_1','dd_Smed_v6_6053_0_1','dd_Smed_v6_4392_0_1','dd_Smed_v6_5999_0_1','dd_Smed_v6_9519_0_1')
all.data.long.delim %>% filter(Contig %in% select.genes2) %>% mutate(hours = as.numeric(hours),Contig = factor(Contig, levels = select.genes2)) %>% ggplot(aes(x=hours,y=log2FC,color=experiment, fill=experiment, group=experiment))+
  geom_point(size=2)+
  geom_line(linewidth=1)+
  facet_wrap(~Contig, scales="free_y",labeller = as_labeller(facet_titles))+
  theme_classic()+  theme(legend.position = "none")+ theme(axis.text= element_text(face = "bold"))
ggsave("plot3.gene.png")


facet_titles <- c('dd_Smed_v6_28398_0_1'='wnt1','dd_Smed_v6_24180_0_1'='notum','dd_Smed_v6_14068_0_1'='nlg-1','dd_Smed_v6_20318_0_1'='equinox','dd_Smed_v6_9584_0_1'='follistatin','dd_Smed_v6_4154_0_1'='glypican-1')
select.genes2<-c('dd_Smed_v6_14068_0_1','dd_Smed_v6_28398_0_1','dd_Smed_v6_24180_0_1','dd_Smed_v6_4154_0_1','dd_Smed_v6_9584_0_1','dd_Smed_v6_20318_0_1')
all.data.long.delim %>% filter(Contig %in% select.genes2) %>% mutate(hours = as.numeric(hours),Contig = factor(Contig, levels = select.genes2)) %>% ggplot(aes(x=hours,y=log2FC,color=experiment, fill=experiment, group=experiment))+
  geom_point(size=2)+
  geom_line(linewidth=1)+
  facet_wrap(~Contig, scales="free_y",labeller = as_labeller(facet_titles))+
  theme_classic()+theme(legend.position = "none")+theme(strip.text = element_text(face = "bold.italic"))+theme(strip.text = element_text(size = 14))+theme(axis.text= element_text(face = "bold"))
ggsave("plot2b.gene.png",height=5.2, width=6)

facet_titles <- c('dd_Smed_v6_14711_0_1'='egr-3','dd_Smed_v6_7607_0_1'='inhibin-1','dd_Smed_v6_6053_0_1'='mex3-2','dd_Smed_v6_4392_0_1'='TRAF-1','dd_Smed_v6_5999_0_1'='CPO','dd_Smed_v6_9519_0_1'='DAT')
select.genes2<-c('dd_Smed_v6_14711_0_1','dd_Smed_v6_6053_0_1','dd_Smed_v6_7607_0_1','dd_Smed_v6_5999_0_1','dd_Smed_v6_4392_0_1','dd_Smed_v6_9519_0_1')
all.data.long.delim %>% filter(Contig %in% select.genes2) %>% mutate(hours = as.numeric(hours),Contig = factor(Contig, levels = select.genes2)) %>% ggplot(aes(x=hours,y=log2FC,color=experiment, fill=experiment, group=experiment))+
  geom_point(size=2)+
  geom_line(linewidth=1)+
  facet_wrap(~Contig, scales="free_y",labeller = as_labeller(facet_titles))+
  theme_classic()+  theme(legend.position = "none")+ theme(strip.text = element_text(face = "bold.italic"))+theme(strip.text = element_text(size = 14))+theme(axis.text= element_text(face = "bold"))
ggsave("plot3b.gene.png",height=5.2, width=6)


#Make a large plot of many different genes in Fig S5
EG.individual.plots.anno.10<-read.table("~/Documents/StudentData/Elias/RNAseq_analysis/analysis032126/EG.individual.plots.anno.082426.2.csv", sep=",", header=TRUE)
contig.label.names<-setNames(EG.individual.plots.anno.10$dd_anno,EG.individual.plots.anno.10$gene_id)

all.data.long.delim.maxdiff<-all.data.long.delim %>% group_by(Contig) %>% mutate(max.head = max(log2FC[experiment=="head"|experiment=="tail"])) %>% ungroup() 
all.data.long.delim.maxdiff<-all.data.long.delim.maxdiff %>% group_by(Contig) %>% mutate(max.heal = max(log2FC[experiment=="heal"])) %>% ungroup() 
all.data.long.delim.maxdiff<-all.data.long.delim.maxdiff %>% mutate(max.diff = max.head - max.heal)
all.data.long.delim.maxdiff<-all.data.long.delim.maxdiff %>% mutate(max.rel.diff = (max.head - max.heal)/max.heal)

all.data.long.delim.maxdiff <- all.data.long.delim.maxdiff %>% 
  mutate(Contig = fct_reorder(Contig, max.diff, .fun = median, .desc = TRUE))

all.data.long.delim.maxdiff %>% filter(Contig %in% EG.individual.plots.anno.10$gene_id) %>% mutate(hours=as.numeric(hours)) %>% ggplot(aes(x=hours,y=log2FC,color=experiment, fill=experiment, group=experiment))+
  geom_point(size=2)+
  geom_line(linewidth=1)+
  #  ggtitle(select.anno16)+
  facet_wrap(~Contig, scales="free_y",labeller = labeller(Contig = contig.label.names))+
  theme_classic()+theme(strip.text = element_text(face = "bold"))+theme(axis.text= element_text(face = "bold"))
ggsave("plot.facet.select.anno16.ww7.png",width=16,height=16)


#Make a large plot of many different genes in Fig S5v2
EG.individual.plots.anno.11<-read.table("~/Documents/StudentData/Elias/RNAseq_analysis/analysis090126b/EG.individual.plots.anno.full.082426.csv", sep=",", header=TRUE)
contig.label.names<-setNames(EG.individual.plots.anno.11$dd_anno,EG.individual.plots.anno.11$gene_id)

all.data.long.delim.maxdiff<-all.data.long.delim %>% group_by(Contig) %>% mutate(max.head = max(log2FC[experiment=="head"|experiment=="tail"])) %>% ungroup() 
all.data.long.delim.maxdiff<-all.data.long.delim.maxdiff %>% group_by(Contig) %>% mutate(max.heal = max(log2FC[experiment=="heal"])) %>% ungroup() 
all.data.long.delim.maxdiff<-all.data.long.delim.maxdiff %>% mutate(max.diff = max.head - max.heal)
all.data.long.delim.maxdiff<-all.data.long.delim.maxdiff %>% mutate(max.rel.diff = (max.head - max.heal)/max.heal)

all.data.long.delim.maxdiff <- all.data.long.delim.maxdiff %>% 
  mutate(Contig = fct_reorder(Contig, max.diff, .fun = median, .desc = TRUE))

Wurtzel_induced_HorT<-read.table("~/Documents/StudentData/Elias/RNAseq_analysis/analysis090126b/Wurtzel_induced_HorT.csv", sep=",", header=TRUE)
Wurtzel_induced_HorTnotIncision<-read.table("~/Documents/StudentData/Elias/RNAseq_analysis/analysis090126b/Wurtzel_induced_HorTnotIncision.csv", sep=",", header=TRUE)

all.data.long.delim.maxdiff %>% filter(Contig %in% Wurtzel_induced_HorT$Wurtzel_induced_HorT) %>% mutate(hours=as.numeric(hours)) %>% ggplot(aes(x=hours,y=log2FC,color=experiment, fill=experiment, group=experiment))+
  geom_point(size=2)+
  geom_line(linewidth=1)+
  #  ggtitle(select.anno16)+
  facet_wrap(~Contig, scales="free_y",labeller = labeller(Contig = contig.label.names),ncol=8)+
  theme_classic()+theme(strip.text = element_text(face = "bold"))+theme(axis.text= element_text(face = "bold"))
ggsave("plot.facet.select.Wurtzel_induced_HorT.png",width=16,height=22)

all.data.long.delim.maxdiff %>% filter(Contig %in% Wurtzel_induced_HorTnotIncision$Wurtzel_induced_HorTnotIncision) %>% mutate(hours=as.numeric(hours)) %>% ggplot(aes(x=hours,y=log2FC,color=experiment, fill=experiment, group=experiment))+
  geom_point(size=2)+
  geom_line(linewidth=1)+
  #  ggtitle(select.anno16)+
  facet_wrap(~Contig, scales="free_y",labeller = labeller(Contig = contig.label.names),ncol=8)+
  theme_classic()+theme(strip.text = element_text(face = "bold"))+theme(axis.text= element_text(face = "bold"))
ggsave("plot.facet.select.Wurtzel_induced_HorTnotIncision.png",width=16,height=12)



##*## Create a spreadsheet of log2FC expression across head, tail, heal datasets (wide format)
Wurtzel.gene.import<-read.table(file="~/Documents/StudentData/Elias/RNAseq_analysis/analysis032126/Wurtzel_WIgene_list.txt",sep="\t",header=FALSE)

Wurtzel.gene.import<-Wurtzel.gene.import[[1]]


#Wurtzel.list<-c('dd_Smed_v6_16222_0_1','dd_Smed_v6_9273_0_1','dd_Smed_v6_36283_0_1','dd_Smed_v6_4944_0_1','dd_Smed_v6_28398_0_1','dd_Smed_v6_15499_0_1','dd_Smed_v6_24180_0_1','dd_Smed_v6_14068_0_1','dd_Smed_v6_585_0_1','dd_Smed_v6_20048_0_1','dd_Smed_v6_4279_0_1','dd_Smed_v6_996_0_1','dd_Smed_v6_5406_0_1','dd_Smed_v6_4392_0_1','dd_Smed_v6_11254_0_1','dd_Smed_v6_8930_0_1','dd_Smed_v6_7607_0_1','dd_Smed_v6_5862_0_1','dd_Smed_v6_18120_0_1','dd_Smed_v6_4619_0_1','dd_Smed_v6_7444_0_1','dd_Smed_v6_13835_0_1','dd_Smed_v6_12587_0_1','dd_Smed_v6_2789_0_1','dd_Smed_v6_14370_0_1','dd_Smed_v6_3606_0_1','dd_Smed_v6_8029_0_1','dd_Smed_v6_8302_0_1','dd_Smed_v6_3632_0_1','dd_Smed_v6_14219_0_1','dd_Smed_v6_678_0_1','dd_Smed_v6_320_0_1','dd_Smed_v6_7293_0_1','dd_Smed_v6_511_0_1','dd_Smed_v6_5781_0_1','dd_Smed_v6_21069_0_1','dd_Smed_v6_11650_0_1','dd_Smed_v6_5749_0_1','dd_Smed_v6_14439_0_1','dd_Smed_v6_2497_0_1','dd_Smed_v6_1921_0_1','dd_Smed_v6_3588_0_1','dd_Smed_v6_10685_0_1','dd_Smed_v6_11561_0_1','dd_Smed_v6_8839_0_1','dd_Smed_v6_10569_0_1','dd_Smed_v6_7731_0_1','dd_Smed_v6_5415_0_1','dd_Smed_v6_20318_0_1','dd_Smed_v6_7188_0_1','dd_Smed_v6_15647_0_1','dd_Smed_v6_1986_0_1','dd_Smed_v6_921_0_1','dd_Smed_v6_12244_0_1','dd_Smed_v6_613_0_1','dd_Smed_v6_1087_0_1','dd_Smed_v6_12472_0_1','dd_Smed_v6_2324_0_1','dd_Smed_v6_3400_0_1','dd_Smed_v6_4154_0_1','dd_Smed_v6_12310_0_1','dd_Smed_v6_3953_0_1','dd_Smed_v6_758_1_1','dd_Smed_v6_9519_0_1','dd_Smed_v6_6808_0_1','dd_Smed_v6_8495_0_1','dd_Smed_v6_9410_0_1','dd_Smed_v6_11629_0_1','dd_Smed_v6_4012_0_1','dd_Smed_v6_8673_0_1','dd_Smed_v6_372_0_1','dd_Smed_v6_5853_0_1','dd_Smed_v6_8634_0_1','dd_Smed_v6_4341_0_1','dd_Smed_v6_3638_0_1','dd_Smed_v6_14289_0_1','dd_Smed_v6_10716_0_1','dd_Smed_v6_1039_0_1','dd_Smed_v6_6668_0_1','dd_Smed_v6_15695_0_1','dd_Smed_v6_3012_0_1','dd_Smed_v6_1778_0_1','dd_Smed_v6_9642_0_1','dd_Smed_v6_13475_0_1','dd_Smed_v6_7994_0_1','dd_Smed_v6_10393_0_1','dd_Smed_v6_6562_0_1','dd_Smed_v6_6364_0_1','dd_Smed_v6_13985_0_1','dd_Smed_v6_7591_0_1','dd_Smed_v6_6053_0_1','dd_Smed_v6_3257_0_1','dd_Smed_v6_10783_0_1','dd_Smed_v6_17481_0_1','dd_Smed_v6_5999_0_1','dd_Smed_v6_2394_0_1','dd_Smed_v6_11501_0_1','dd_Smed_v6_6259_0_1','dd_Smed_v6_3441_0_1','dd_Smed_v6_1531_0_1','dd_Smed_v6_7811_0_1','dd_Smed_v6_4808_0_1','dd_Smed_v6_2921_0_1','dd_Smed_v6_9825_0_1','dd_Smed_v6_12695_0_1','dd_Smed_v6_3939_0_1','dd_Smed_v6_6448_0_1','dd_Smed_v6_102_0_1','dd_Smed_v6_5390_0_1','dd_Smed_v6_9610_0_1','dd_Smed_v6_10259_0_1','dd_Smed_v6_7824_0_1','dd_Smed_v6_12066_0_1','dd_Smed_v6_9204_0_1','dd_Smed_v6_4299_0_1','dd_Smed_v6_9402_0_1','dd_Smed_v6_12455_0_1','dd_Smed_v6_11964_0_1','dd_Smed_v6_9490_0_1','dd_Smed_v6_5786_0_1','dd_Smed_v6_19_0_1','dd_Smed_v6_13188_0_1','dd_Smed_v6_561_0_1','dd_Smed_v6_7295_0_1','dd_Smed_v6_19149_0_1','dd_Smed_v6_13860_0_1','dd_Smed_v6_6995_0_1','dd_Smed_v6_14711_0_1')

wurtzel.head.tail.heal.analysis<-head.tail.heal.analysis %>% filter(Contig %in% Wurtzel.gene.import)

ddv6.anno1<-read.table(file = "~/Documents/StudentData/Elias/RNAseq_analysis/analysis032126/ddv6_RNAseq_anno2.txt",sep="\t", header=TRUE, quote = "")
head.tail.heal.analysis.anno<-left_join(head.tail.heal.analysis,ddv6.anno1,by="Contig")

#WI.head.late, WI.head.early, WI.tail.early, WI.tail.late
head.tail.heal.analysis.anno.early <- head.tail.heal.analysis.anno %>% filter(WI.tail.early>1) %>% filter(WI.head.early>1)
head.tail.heal.analysis.anno.late <- head.tail.heal.analysis.anno %>% filter(WI.tail.late>1) %>% filter(WI.head.late>1)
head.tail.heal.analysis.anno.incision<-head.tail.heal.analysis.anno %>% mutate(inc.late = (`8 h`+`12 h`+`16 h`+`24 h`)/4) %>% filter(inc.late>1)
head.tail.heal.analysis.anno.incision.early<-head.tail.heal.analysis.anno %>% mutate(inc.early = (`1 h`+`2 h`+`4 h`+`6 h`)/4) %>% filter(inc.early>0.8) %>% filter(WI.head.early<0.5)
inc.early.spec<-as.vector(head.tail.heal.analysis.anno.incision.early$Contig)


wurtzel.head.tail.heal.anno.analysis<-head.tail.heal.analysis.anno %>% filter(Contig %in% Wurtzel.gene.import)
write.table(wurtzel.head.tail.heal.anno.analysis,"wurtzel.head.tail.heal.anno.analysis.csv",sep=",")




# Wurtzel.short.anno<-read.table("~/Documents/StudentData/Elias/RNAseq_analysis/analysis032126/Wurtzel_short_anno.txt", sep="\t", header=TRUE)
# wurtzel.head.tail.heal.anno.analysis2<-wurtzel.head.tail.heal.anno.analysis %>% left_join(Wurtzel.short.anno,wurtzel.head.tail.heal.anno.analysis, by="Contig")
# 
# wurtzel.head.tail.heal.anno.analysis2<-wurtzel.head.tail.heal.anno.analysis2 %>% mutate(inc.late = (`12 h`+`16 h`+`24 h`)/3) %>% mutate(inc.early = (`1 h`+`2 h`+`4 h`+`6 h`)/4)
# 
# wurtzel.head.tail.heal.anno.analysis2<-wurtzel.head.tail.heal.anno.analysis2 %>% rowwise() %>% mutate(max.late = max(WI.head.late,WI.tail.late)) %>% mutate(max.early = max(WI.head.early,WI.tail.early))

##*# Heatmap Figure 1A and S1A
library(pheatmap)
WI.anno3<-read.table("/Users/christianpetersen/Documents/StudentData/Elias/RNAseq_analysis/Wurtzel_Wenemoser_lists_summary_anno_031325.txt",sep="\t", header=TRUE)
colnames(WI.anno3)<-c("Contig","gene_anno")
wurtzel.WI3.names<-WI.anno3$Contig
wurtzel.head.tail.heal.analysis3<-head.tail.heal.analysis %>% filter(Contig %in% wurtzel.WI3.names)

wurtzel.head.tail.heal.analysis.heatmap3<- left_join(wurtzel.head.tail.heal.analysis3, WI.anno3, by="Contig")

rownames.w.h.t.h.a.h3<-wurtzel.head.tail.heal.analysis.heatmap3$gene_anno

wurtzel.head.tail.heal.analysis.heatmap3<-wurtzel.head.tail.heal.analysis.heatmap3 %>% select(3:13,15:25,27:37)
wurtzel.head.tail.heal.analysis.heatmap3<-as.matrix(wurtzel.head.tail.heal.analysis.heatmap3)
rownames(wurtzel.head.tail.heal.analysis.heatmap3)<-rownames.w.h.t.h.a.h3

# pheatmap(wurtzel.head.tail.heal.analysis.heatmap3,cluster_cols=FALSE, cluster_rows=TRUE,fontsize=6)
# heatmap_WI3<-pheatmap(wurtzel.head.tail.heal.analysis.heatmap3,cluster_cols=FALSE, cluster_rows=TRUE,fontsize=5)
# ggsave("EG_heatmap1_WI3.pdf",heatmap_WI3,width=7, height=8)

pheatmap(wurtzel.head.tail.heal.analysis.heatmap3,cluster_cols=FALSE, cluster_rows=TRUE,fontsize=5,scale="row")
heatmap_WI3<-pheatmap(wurtzel.head.tail.heal.analysis.heatmap3,cluster_cols=FALSE, cluster_rows=TRUE,fontsize=5,scale="row")
ggsave("EG_heatmap1_WI3b.pdf",heatmap_WI3,width=7, height=8)



#overlay all of the geom_line plots from all the genes at once
all.data.long.delim.overlay<-all.data.long.delim %>% mutate(Contig.experiment=paste0(Contig,".",experiment))
WI.timing<-read.table("/Users/christianpetersen/Documents/StudentData/Elias/RNAseq_analysis/Elias_Wurtzel_analysis050526/Wurtzel_Wenemoser_timing_060226.txt", sep="\t",header=TRUE)
WI.anno3<-read.table("/Users/christianpetersen/Documents/StudentData/Elias/RNAseq_analysis/Wurtzel_Wenemoser_lists_summary_anno_031325.txt",sep="\t", header=TRUE)
colnames(WI.anno3)<-c("Contig","gene_anno")
wurtzel.WI3.names<-WI.anno3$Contig
overlayed.lineplots<- all.data.long.delim.overlay %>% filter(Contig %in% wurtzel.WI3.names) 
overlayed.lineplots.timing<-left_join(overlayed.lineplots,WI.timing,by="Contig")


#normalize within each Contig (gene) by z-score of the log2FC
overlayed.lineplots.timing.zscore<-overlayed.lineplots.timing %>% group_by(Contig) %>%
  mutate(log2FC.zscore = as.numeric(scale(log2FC))) %>%
  ungroup()


overlayed.lineplots.timing.zscore %>% filter(Timing == "Early") %>% mutate(hours=as.numeric(hours)) %>% ggplot(aes(x=hours,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
  #  geom_point(size=2)+
  geom_line(linewidth=1, alpha=0.05)+
  stat_smooth(aes(group=experiment))+
  #  facet_wrap(~Contig, scales="free_y")+
  theme_classic()+theme(axis.text= element_text(face = "bold",size=14))
ggsave("overlay.WI3.early.zscore.pdf",width=4.5, height=6)

overlayed.lineplots.timing.zscore %>% filter(Timing == "Late") %>% mutate(hours=as.numeric(hours)) %>% ggplot(aes(x=hours,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
  #  geom_point(size=2)+
  geom_line(linewidth=1, alpha=0.05)+
  stat_smooth(aes(group=experiment))+
  #  facet_wrap(~Contig, scales="free_y")+
  theme_classic()+theme(axis.text= element_text(face = "bold",size=14))
ggsave("overlay.WI3.late.zscore.pdf",width=4.5, height=6)


#overlayed.lineplots.timing.zscore %>% filter(Timing == "Early") %>% mutate(name = fct_reorder(hours, as.numeric(hours))) %>% ggplot(aes(x=name,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
overlayed.lineplots.timing.zscore %>% filter(Timing == "Early") %>% mutate(hours=as.numeric(hours)) %>% ggplot(aes(x=hours,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
  #  geom_point(size=2)+
  geom_line(linewidth=1, alpha=0.05)+
  stat_smooth(aes(group=experiment))+
  #  facet_wrap(~Contig, scales="free_y")+
  theme_classic()+theme(axis.text= element_text(face = "bold",size=18))
ggsave("overlay.WI3.early.zscore2.pdf",width=8, height=6)


overlayed.lineplots.timing.zscore %>% filter(Timing == "Late") %>% mutate(hours=as.numeric(hours)) %>% ggplot(aes(x=hours,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
  #overlayed.lineplots.timing.zscore %>% filter(Timing == "Late") %>% mutate(name = fct_reorder(hours, as.numeric(hours))) %>% ggplot(aes(x=name,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
  #  geom_point(size=2)+
  geom_line(linewidth=1, alpha=0.05)+
  stat_smooth(aes(group=experiment))+
  #  facet_wrap(~Contig, scales="free_y")+
  theme_classic()+theme(axis.text= element_text(face = "bold",size=18))
ggsave("overlay.WI3.late.zscore2.pdf",width=8, height=6)

#overlayed.lineplots.timing.zscore %>% filter(Timing == "Early") %>% mutate(name = fct_reorder(hours, as.numeric(hours))) %>% ggplot(aes(x=name,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
overlayed.lineplots.timing.zscore %>% filter(Timing == "Early") %>% mutate(hours=as.numeric(hours)) %>% ggplot(aes(x=hours,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
  #  geom_point(size=2)+
  geom_line(linewidth=1, alpha=0.05)+
  stat_smooth(aes(group=experiment))+
  coord_cartesian(ylim = c(-2, 2.5))+
  #  facet_wrap(~Contig, scales="free_y")+
  theme_classic()+theme(axis.text= element_text(face = "bold",size=18))
ggsave("overlay.WI3.early.zscore3.png",width=8, height=6)


overlayed.lineplots.timing.zscore %>% filter(Timing == "Late") %>% mutate(hours=as.numeric(hours)) %>% ggplot(aes(x=hours,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
  #overlayed.lineplots.timing.zscore %>% filter(Timing == "Late") %>% mutate(name = fct_reorder(hours, as.numeric(hours))) %>% ggplot(aes(x=name,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
  #  geom_point(size=2)+
  geom_line(linewidth=1, alpha=0.05)+
  stat_smooth(aes(group=experiment))+
  coord_cartesian(ylim = c(-2, 2.5))+
  #  facet_wrap(~Contig, scales="free_y")+
  theme_classic()+theme(axis.text= element_text(face = "bold",size=18))
ggsave("overlay.WI3.late.zscore3.png",width=8, height=6)


#overlayed.lineplots.timing.zscore %>% filter(Timing == "Early") %>% mutate(name = fct_reorder(hours, as.numeric(hours))) %>% ggplot(aes(x=name,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
overlayed.lineplots.timing.zscore %>% filter(Timing == "Early") %>% mutate(hours=as.numeric(hours)) %>% ggplot(aes(x=hours,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
  #  geom_point(size=2)+
  geom_line(linewidth=1, alpha=0.05)+
  stat_smooth(aes(group=experiment))+
  #  facet_wrap(~Contig, scales="free_y")+
  theme_classic()+theme(axis.text= element_text(face = "bold",size=18))
ggsave("overlay.WI3.early.zscore2.pdf",width=8, height=6)



#overlayed.lineplots.timing.zscore %>% filter(Timing == "Late") %>% mutate(name = fct_reorder(hours, as.numeric(hours))) %>% ggplot(aes(x=name,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
overlayed.lineplots.timing.zscore %>% filter(Timing == "Late") %>% mutate(hours=as.numeric(hours)) %>% ggplot(aes(x=hours,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
  
  #  geom_point(size=2)+
  geom_line(linewidth=1, alpha=0.10)+
  stat_smooth(aes(group=experiment))+
  #  facet_wrap(~Contig, scales="free_y")+
  theme_classic()+theme(axis.text= element_text(face = "bold",size=18))+
  ylim(-2, 2)+
  xlab("time (hours)")+
  theme(
    # Changes the outer background area
    plot.background = element_rect(fill = "black", color = "black"),
    # Changes the inner plotting panel background
    panel.background = element_rect(fill = "black", color = "black"),
    # Formats all text characters to white
    text = element_text(color = "white"),
    # Formats axis labels and titles specifically to white
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.ticks = element_line(color = "white"),
    axis.line = element_line(color = "white")
    # Optional: Lighten grid lines so they are visible but subtle
    #    panel.grid.major = element_line(color = "gray30"),
    #    panel.grid.minor = element_line(color = "gray15")
  )
ggsave("overlay.WI3.late.zscore2_bb.pdf",width=8, height=6)


overlayed.lineplots.timing.zscore %>% filter(Timing == "Early") %>% mutate(hours=as.numeric(hours)) %>% ggplot(aes(x=hours,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
  #overlayed.lineplots.timing.zscore %>% filter(Timing == "Early") %>% mutate(name = fct_reorder(hours, as.numeric(hours))) %>% ggplot(aes(x=name,y=log2FC.zscore, color=experiment, group=Contig.experiment))+
  #  geom_point(size=2)+
  geom_line(linewidth=1, alpha=0.10)+
  stat_smooth(aes(group=experiment))+
  #  facet_wrap(~Contig, scales="free_y")+
  theme_classic()+theme(axis.text= element_text(face = "bold",size=18))+
  xlab("time (hours)")+
  ylim(-2, 2)+
  theme(
    # Changes the outer background area
    plot.background = element_rect(fill = "black", color = "black"),
    # Changes the inner plotting panel background
    panel.background = element_rect(fill = "black", color = "black"),
    # Formats all text characters to white
    text = element_text(color = "white"),
    # Formats axis labels and titles specifically to white
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.ticks = element_line(color = "white"),
    axis.line = element_line(color = "white")
    # Optional: Lighten grid lines so they are visible but subtle
    #    panel.grid.major = element_line(color = "gray30"),
    #    panel.grid.minor = element_line(color = "gray15")
  )
ggsave("overlay.WI3.early.zscore2_bb.pdf",width=8, height=6)




##########################################################################
#Fitting data to impulse function

#setwd("~/Documents/StudentData/Elias/RNAseq_analysis/analysis032126/")

library(minpack.lm)
library(ggplot2)
library(dplyr)

# ---------------------------
# Define simplified impulse function
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

impulse <- function(x, A, B, t1, t2) {
  A * sigmoid(B * (x - t1)) * (1 - sigmoid(B * (x - t2)))
}

# ---------------------------
# Main function to fit all genes
fit_all_genes <- function(df, plot_dir = NULL) {
  timepoints <- 0:10  # 11 timepoints
  results <- list()
  
  if (!is.null(plot_dir)) {
    dir.create(plot_dir, showWarnings = FALSE)
  }
  
  for (i in 1:nrow(df)) {
    gene <- df[i, 1]
    y_vals <- as.numeric(df[i, 2:12])
    gene_df <- data.frame(x = timepoints, y = y_vals)
    
    start_vals <- list(A = max(y_vals, na.rm = TRUE), B = 1, t1 = 2, t2 = 8)
    
    fit <- tryCatch({
      nlsLM(
        y ~ impulse(x, A, B, t1, t2),
        data = gene_df,
        start = start_vals,
        lower = c(0, 0.01, 0, 0),
        upper = c(10, 10, 10, 10),
        control = nls.lm.control(maxiter = 500)
      )
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      coefs <- coef(fit)
      results[[i]] <- data.frame(
        Gene = gene,
        A = coefs["A"],
        B = coefs["B"],
        t1 = coefs["t1"],
        t2 = coefs["t2"],
        Fit = TRUE
      )
      
      # Plot
      pred_df <- data.frame(x = seq(0, 10, length.out = 200))
      pred_df$y <- impulse(pred_df$x, coefs["A"], coefs["B"], coefs["t1"], coefs["t2"])
      
      p <- ggplot(gene_df, aes(x = x, y = y)) +
        geom_point(size = 2) +
        geom_line(data = pred_df, aes(x = x, y = y), color = "blue", size = 1) +
        ylim(0,5)+
        ggtitle(paste("Impulse Fit -", gene)) +
        theme_minimal()
      
      if (!is.null(plot_dir)) {
        ggsave(filename = paste0(plot_dir, "/", gene, ".png"), plot = p, width = 5, height = 4)
      } else {
        print(p)  # Just show it if not saving
      }
    } else {
      results[[i]] <- data.frame(
        Gene = gene,
        A = NA, B = NA, t1 = NA, t2 = NA,
        Fit = FALSE
      )
    }
  }
  
  return(bind_rows(results))
}



# Example: assume your data frame is called `gene_data`
# Column 1: Gene names
# Columns 2-12: Expression at timepoints 0 to 10

#results <- fit_all_genes(gene_data, plot_dir = "impulse_plots")

# View summary of fit results
#head(results)

Wurtzel.short.anno<-read.table("~/Documents/StudentData/Elias/RNAseq_analysis/analysis032126/Wurtzel_short_anno.txt", sep="\t", header=TRUE)
wurtzel.head.tail.heal.anno.analysis2<-wurtzel.head.tail.heal.anno.analysis %>% left_join(Wurtzel.short.anno,wurtzel.head.tail.heal.anno.analysis, by="Contig")
wurtzel.head.tail.heal.anno.analysis2<-wurtzel.head.tail.heal.anno.analysis2 %>% mutate(inc.late = (`12 h`+`16 h`+`24 h`)/3) %>% mutate(inc.early = (`1 h`+`2 h`+`4 h`+`6 h`)/4)
wurtzel.head.tail.heal.anno.analysis2<-wurtzel.head.tail.heal.anno.analysis2 %>% rowwise() %>% mutate(max.late = max(WI.head.late,WI.tail.late)) %>% mutate(max.early = max(WI.head.early,WI.tail.early))

wurtzel.head.timeseries<-wurtzel.head.tail.heal.anno.analysis2 %>% select(1,3:13)

results.head <- fit_all_genes(wurtzel.head.timeseries, plot_dir = "impulse_plots_head2")
head(results.head)


wurtzel.tail.timeseries<-wurtzel.head.tail.heal.anno.analysis2 %>% select(1,15:25)
results.tail <- fit_all_genes(wurtzel.tail.timeseries, plot_dir = "impulse_plots_tail2")
head(results.tail)

wurtzel.heal.timeseries<-wurtzel.head.tail.heal.anno.analysis2 %>% select(1,27:37)
results.heal <- fit_all_genes(wurtzel.heal.timeseries, plot_dir = "impulse_plots_heal2")
head(results.heal)


#results.total<-left_join()


results.tail <-results.tail %>% mutate (type="tail")
results.head <- results.head %>% mutate (type="head")
results.heal <- results.heal %>% mutate (type="heal")

results.tail<-results.tail %>% mutate(tmax=(t1+t2)/2) %>% mutate(ymax = impulse((t1+t2)/2,A, B, t1, t2))
results.head<-results.head %>% mutate(tmax=(t1+t2)/2) %>% mutate(ymax = impulse((t1+t2)/2,A, B, t1, t2))
results.heal<-results.heal %>% mutate(tmax=(t1+t2)/2) %>% mutate(ymax = impulse((t1+t2)/2,A, B, t1, t2))

results.total<-bind_rows(results.head,results.tail,results.heal)

#sigmoid = 1/(1+e^-x)
#A * sigmoid(B * (x - t1)) * (1 - sigmoid(B * (x - t2)))
#results.total<-results.total %>% mutate(tmax=(t1+t2)/2) %>% mutate(ymax = impulse((t1+t2)/2,A, B, t1, t2))

# results.total %>% ggplot(aes(x=type, y=A))+
#   geom_jitter()+
#   theme_classic()
# ggsave("impulse.jitter.head.heal.tail.A.png")
# 
# results.total %>% ggplot(aes(x=type, y=B))+
#   geom_jitter()+
#   theme_classic()
# ggsave("impulse.jitter.head.heal.tail.B.png")
# 
# 
# results.total %>% ggplot(aes(x=type, y=t1))+
#   geom_jitter()+
#   theme_classic()
# ggsave("impulse.jitter.head.heal.tail.t1.png")
# 
# 
# results.total %>% ggplot(aes(x=type, y=t2))+
#   geom_jitter()+
#   theme_classic()
# ggsave("impulse.jitter.head.heal.tail.t2.png")
# 
# results.total %>% ggplot(aes(x=type, y=tmax))+
#   geom_jitter()+
#   theme_classic()
# ggsave("impulse.jitter.head.heal.tail.tmax.png")
# 
# results.total %>% ggplot(aes(x=type, y=ymax))+
#   geom_jitter()+
#   theme_classic()
# ggsave("impulse.jitter.head.heal.tail.ymax.png")


#install.packages("ggpubr")
library(ggpubr)

results.headheal.wide<-left_join(results.head,results.heal,by="Contig")

results.tailheal.wide<-left_join(results.tail,results.heal,by="Contig")

results.headheal.wide <- results.headheal.wide %>% mutate(dA = A.y - A.x) %>% mutate(dB = B.y-B.x) %>% mutate(dt1 = t1.y-t1.x) %>% mutate(dt2 = t2.y-t2.x)
results.headheal.wide <- results.headheal.wide %>% mutate(dA.perc = 100*dA / A.x) %>% mutate(dB.perc = 100*dB/B.x) %>% mutate(dt1.perc = 100*dt1/t1.x) %>% mutate(dt2.perc = 100*dt2/t2.x)


results.tailheal.wide <- results.tailheal.wide %>% mutate(dA = A.y - A.x) %>% mutate(dB = B.y-B.x) %>% mutate(dt1 = t1.y-t1.x) %>% mutate(dt2 = t2.y-t2.x)
results.tailheal.wide <- results.tailheal.wide %>% mutate(dA.perc = 100*dA / A.x) %>% mutate(dB.perc = 100*dB/B.x) %>% mutate(dt1.perc = 100*dt1/t1.x) %>% mutate(dt2.perc = 100*dt2/t2.x)

results.headheal.wide<-left_join(results.headheal.wide,Wurtzel.short.anno,by="Contig")
results.tailheal.wide<-left_join(results.tailheal.wide,Wurtzel.short.anno,by="Contig")

#save parameters from fit
write.table(results.headheal.wide,"results.headheal.wide.csv")
write.table(results.tailheal.wide,"results.tailheal.wide.csv")

#plot parameters Figure S3
ggpaired(results.headheal.wide, "A.x","A.y",line.color = "gray", line.size = 0.4, palette = "npg", title="head.vs.heal.A")
ggsave("impulse.headheal.paired.A.png")

ggpaired(results.headheal.wide, "B.x","B.y",line.color = "gray", line.size = 0.4, palette = "npg", title="head.vs.heal.B")
ggsave("impulse.headheal.paired.B.png")

ggpaired(results.headheal.wide, "t1.x","t1.y",line.color = "gray", line.size = 0.4, palette = "npg", title="head.vs.heal.t1")
ggsave("impulse.headheal.paired.t1.png")

ggpaired(results.headheal.wide, "t2.x","t2.y",line.color = "gray", line.size = 0.4, palette = "npg", title="head.vs.heal.t2")
ggsave("impulse.headheal.paired.t2.png")

ggpaired(results.headheal.wide, "tmax.x","tmax.y",line.color = "gray", line.size = 0.4, palette = "npg", title="head.vs.heal.tmax")
ggsave("impulse.headheal.paired.tmax.png")

ggpaired(results.headheal.wide, "ymax.x","ymax.y",line.color = "gray", line.size = 0.4, palette = "npg", title="head.vs.heal.ymax")
ggsave("impulse.headheal.paired.ymax.png")


ggpaired(results.tailheal.wide, "A.x","A.y",line.color = "gray", line.size = 0.4, palette = "npg", title="tail.vs.heal.A")
ggsave("impulse.tailheal.paired.A.png")

ggpaired(results.tailheal.wide, "B.x","B.y",line.color = "gray", line.size = 0.4, palette = "npg", title="tail.vs.heal.B")
ggsave("impulse.tailheal.paired.B.png")

ggpaired(results.tailheal.wide, "t1.x","t1.y",line.color = "gray", line.size = 0.4, palette = "npg", title="tail.vs.heal.t1")
ggsave("impulse.tailheal.paired.t1.png")

ggpaired(results.tailheal.wide, "t2.x","t2.y",line.color = "gray", line.size = 0.4, palette = "npg", title="tail.vs.heal.t2")
ggsave("impulse.tailheal.paired.t2.png")

ggpaired(results.tailheal.wide, "tmax.x","tmax.y",line.color = "gray", line.size = 0.4, palette = "npg", title="tail.vs.heal.tmax")
ggsave("impulse.tailheal.paired.tmax.png")

ggpaired(results.tailheal.wide, "ymax.x","ymax.y",line.color = "gray", line.size = 0.4, palette = "npg", title="tail.vs.heal.ymax")
ggsave("impulse.tailheal.paired.ymax.png")

#plot parameters Fig 1C
ggpaired(results.headheal.wide, "ymax.x","ymax.y",
         #        color="black",
         line.color = "darkgray", 
         line.size = 0.4, 
         #         palette = "npg", 
         title="head.vs.heal.t2",
         #         fill="condition", 
         #         palette = c("#E7B800", "#2E9FDF"),
         #         palette= c("blue","green"),
         point.size = 4)+
  #  scale_fill_manual(values = alpha(c("t2.y" = "#4bb254", "t2.x" = "#1E90FF"),0.75))+
  theme(element_text(face = "bold"))+
  theme(axis.text.x = element_text(size = 24))+
  theme(axis.text.y = element_text(size = 24))
ggsave("impulse.headheal.paired.ymax_2.png",width=6, height=8)


ggpaired(results.headheal.wide, "t2.x","t2.y",
         #        color="black",
         line.color = "darkgray", 
         line.size = 0.4, 
         #         palette = "npg", 
         title="head.vs.heal.t2",
         #         fill="condition", 
         #         palette = c("#E7B800", "#2E9FDF"),
         #         palette= c("blue","green"),
         point.size = 4)+
  #  scale_fill_manual(values = alpha(c("t2.y" = "#4bb254", "t2.x" = "#1E90FF"),0.75))+
  theme(element_text(face = "bold"))+
  theme(axis.text.x = element_text(size = 24))+
  theme(axis.text.y = element_text(size = 24))
ggsave("impulse.headheal.paired.t2_2.png",width=6, height=8)


#statistical tests comparing impulse parameter fits between head vs heal and between tail vs heal conditions 

#stat tests on the impulse model fitting
#results.headheal.wide,results.tailheal.wide
#"A.x","A.y"; "B.x","B.y"; "t1.x","t1.y" ; "t2.x","t2.y" ; "ymax.x","ymax.y"

ks.test(results.headheal.wide$A.x,results.headheal.wide$A.y,
        alternative=c("two.sided")
)
#D = 0.23803, p-value = 0.002723


ks.test(results.headheal.wide$B.x,results.headheal.wide$B.y,
        alternative=c("two.sided")
)
#D = 0.18703, p-value = 0.034

ks.test(results.headheal.wide$t1.x,results.headheal.wide$t1.y,
        alternative=c("two.sided")
)
#D = 0.22491, p-value = 0.005523

ks.test(results.headheal.wide$t2.x,results.headheal.wide$t2.y,
        alternative=c("two.sided")
)
#D = 0.41864, p-value = 2.723e-09

ks.test(results.headheal.wide$tmax.x,results.headheal.wide$tmax.y,
        alternative=c("two.sided")
)
#D = 0.33324, p-value = 4.826e-06

ks.test(results.headheal.wide$ymax.x,results.headheal.wide$ymax.y,
        alternative=c("two.sided")
)
#D = 0.3, p-value = 5.6e-05

###
###
ks.test(results.tailheal.wide$A.x,results.tailheal.wide$A.y,
        alternative=c("two.sided")
)
#D = 0.20624, p-value = 0.01474


ks.test(results.tailheal.wide$B.x,results.tailheal.wide$B.y,
        alternative=c("two.sided")
)
#D = 0.084521, p-value = 0.8041


ks.test(results.tailheal.wide$t1.x,results.tailheal.wide$t1.y,
        alternative=c("two.sided")
)
#D = 0.13042, p-value = 0.28

ks.test(results.tailheal.wide$t2.x,results.tailheal.wide$t2.y,
        alternative=c("two.sided")
)
#D = 0.25994, p-value = 0.0008192

ks.test(results.tailheal.wide$tmax.x,results.tailheal.wide$tmax.y,
        alternative=c("two.sided")
)
#D = 0.25994, p-value = 0.03378


ks.test(results.headheal.wide$ymax.x,results.headheal.wide$ymax.y,
        alternative=c("two.sided")
)
#D = 0.3, p-value = 5.6e-05

#091226 wilcoxon paired-sample tests of the impulse function parameters
#(A, B, t1, t2, tmax, ymax) x (head vs heal, tail vs heal) 

library(dplyr)
library(rstatix)



results.total.wilcox.headvsheal.A <-results.total %>% filter(type %in% c("head","heal")) %>%
  pairwise_wilcox_test(A ~ type, p.adjust.method = "holm")
results.total.wilcox.headvsheal.A$p.adj
#p=4.365596e-05

results.total.wilcox.headvsheal.B <-results.total %>% filter(type %in% c("head","heal")) %>%
  pairwise_wilcox_test(B ~ type, p.adjust.method = "holm")
results.total.wilcox.headvsheal.B$p.adj
#p=0.005003838

results.total.wilcox.headvsheal.t1 <-results.total %>% filter(type %in% c("head","heal")) %>%
  pairwise_wilcox_test(t1 ~ type, p.adjust.method = "holm")
results.total.wilcox.headvsheal.t1$p.adj
#p=0.02715112

results.total.wilcox.headvsheal.t2 <-results.total %>% filter(type %in% c("head","heal")) %>%
  pairwise_wilcox_test(t2 ~ type, p.adjust.method = "holm")
results.total.wilcox.headvsheal.t2$p.adj
#p=2.390736e-09

results.total.wilcox.headvsheal.tmax <-results.total %>% filter(type %in% c("head","heal")) %>%
  pairwise_wilcox_test(tmax ~ type, p.adjust.method = "holm")
results.total.wilcox.headvsheal.tmax$p.adj
#p=2.569947e-07

results.total.wilcox.headvsheal.ymax <-results.total %>% filter(type %in% c("head","heal")) %>%
  pairwise_wilcox_test(ymax ~ type, p.adjust.method = "holm")
results.total.wilcox.headvsheal.ymax$p.adj
#p=2.909718e-07

#
results.total.wilcox.tailvsheal.A <-results.total %>% filter(type %in% c("tail","heal")) %>%
  pairwise_wilcox_test(A ~ type, p.adjust.method = "holm")
results.total.wilcox.tailvsheal.A$p.adj
#p=0.001686925

results.total.wilcox.tailvsheal.B <-results.total %>% filter(type %in% c("tail","heal")) %>%
  pairwise_wilcox_test(B ~ type, p.adjust.method = "holm")
results.total.wilcox.tailvsheal.B$p.adj
#p=0.8615165

results.total.wilcox.tailvsheal.t1 <-results.total %>% filter(type %in% c("tail","heal")) %>%
  pairwise_wilcox_test(t1 ~ type, p.adjust.method = "holm")
results.total.wilcox.tailvsheal.t1$p.adj
#p=0.2951852

results.total.wilcox.tailvsheal.t2 <-results.total %>% filter(type %in% c("tail","heal")) %>%
  pairwise_wilcox_test(t2 ~ type, p.adjust.method = "holm")
results.total.wilcox.tailvsheal.t2$p.adj
#p=0.000328818

results.total.wilcox.tailvsheal.tmax <-results.total %>% filter(type %in% c("tail","heal")) %>%
  pairwise_wilcox_test(tmax ~ type, p.adjust.method = "holm")
results.total.wilcox.tailvsheal.tmax$p.adj
#p=0.02704853

results.total.wilcox.tailvsheal.ymax <-results.total %>% filter(type %in% c("tail","heal")) %>%
  pairwise_wilcox_test(ymax ~ type, p.adjust.method = "holm")
results.total.wilcox.tailvsheal.ymax$p.adj
#p=1.390402e-06

wilcoxon.summary.impulse.labels<-c("A","B","t1","t2","tmax","ymax")
wilcoxon.summary.impulse.head<-c(results.total.wilcox.headvsheal.A$p.adj,results.total.wilcox.headvsheal.B$p.adj,results.total.wilcox.headvsheal.t1$p.adj,results.total.wilcox.headvsheal.t2$p.adj,results.total.wilcox.headvsheal.tmax$p.adj,results.total.wilcox.headvsheal.ymax$p.adj)
wilcoxon.summary.impulse.tail<-c(results.total.wilcox.tailvsheal.A$p.adj,results.total.wilcox.tailvsheal.B$p.adj,results.total.wilcox.tailvsheal.t1$p.adj,results.total.wilcox.tailvsheal.t2$p.adj,results.total.wilcox.tailvsheal.tmax$p.adj,results.total.wilcox.tailvsheal.ymax$p.adj)


####
#plot given gene of interest overlay log2FC datapoints with impulse model, eg Fig 1C
#here is another gene for showing the method:egal1
gene.of.interest<-"dd_Smed_v6_8302_0_1" 

#this is the actual wnt1 head data
goi.head.series<-wurtzel.head.timeseries %>% filter(Contig==gene.of.interest) %>% select(2:11)

#this is the impulse parameters fitting wnt1 head data
goi.results.head<-results.head %>% filter(Contig==gene.of.interest)

#results.head <- fit_all_genes(wurtzel.head.timeseries, plot_dir = "impulse_plots_head2")

#the impulse parameters fitting wnt1 heal data
goi.results.heal<-results.heal %>% filter(Contig==gene.of.interest)

#the actual wnt1 heal data
goi.heal.series<-wurtzel.heal.timeseries %>% filter(Contig==gene.of.interest) %>% select(2:11)

goi.tail.series<-wurtzel.tail.timeseries %>% filter(Contig==gene.of.interest) %>% select(2:11)
goi.tail.series<-as.data.frame(goi.tail.series)

goi.head.series_df<-as.data.frame(t(goi.head.series))
goi.heal.series_df<-as.data.frame(t(goi.heal.series))

gene_df <- as.data.frame(t(goi.head.series))



#colnames(gene_df)<-c("y.head","y.heal")
#gene_df<-cbind(gene_df,c("1","2","3","4","5","6","7","8","9","10"))
gene_df<-cbind(goi.head.series_df,goi.heal.series_df,1:10)

colnames(gene_df)<-c("y.head","y.heal","x")

pred_df <- data.frame(x = seq(0, 10, length.out = 200))
#pred_df$y <- impulse(pred_df$x, goi.results.head[2], goi.results.head[3], goi.results.head[4], goi.results.head[5])

A.head<-goi.results.head[2]
B.head<-goi.results.head[3]
t1.head<-goi.results.head[4]
t2.head<-goi.results.head[5]

colnames(pred_df)<-c("x.val")

#pred_df <- pred_df %>% mutate(y=A * (1/(1+exp(-(B * ('x.val' - t1))))) * (1-(1/(1+exp(-(B * ('x.val' - t2)))))))

pred_df <- pred_df %>%
  rowwise() %>%
  mutate(y.head = impulse(x.val-1,A.head,B.head,t1.head,t2.head)) %>%
  ungroup() # Always a good practice to ungroup afterwards

pred_df_heal <- data.frame(x = seq(0, 10, length.out = 200))
A.heal<-goi.results.heal[2]
B.heal<-goi.results.heal[3]
t1.heal<-goi.results.heal[4]
t2.heal<-goi.results.heal[5]
colnames(pred_df_heal)<-c("x.val")

pred_df_heal <- pred_df_heal %>%
  rowwise() %>%
  mutate(y.heal = impulse(x.val-1,A.heal,B.heal,t1.heal,t2.heal)) %>%
  ungroup() # Always a good practice to ungroup afterwards



ggplot(gene_df, aes(x = x, y = y.head)) +
  geom_point(size = 3,shape=21, stroke=2,fill="white",color="red") +
  geom_point(aes(x=x,y=y.heal),size=3,shape=21,stroke=2,color="#4bb254")+
  geom_line(data = pred_df, aes(x = x.val, y = y.head$A), color = "red", size = 2) +
  geom_line(data = pred_df_heal, aes(x = x.val, y = y.heal$A), color = "#4bb254", size = 2) +
  xlim(1,10)+
  ggtitle(paste("egal1")) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 18, face = "bold"))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.line = element_line(linewidth = 1.5))+
  #  scale_x_axis_discrete=c(1="0",2="1",3="2",4="4",5="6",6="8",7="12",8="16",9="24",10="36")
  #  scale_x_axis_discrete=c("1"="0","2"="1","3"="2","4"="4","5"="6","6"="8","7"="12","8"="16","9"="24","10"="36")
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))

#  scale_x_continuous(breaks = c(0,1,2,4,6,8,12,16,24,36))
ggsave("egal1_dd8302_impulsefit_HeadvsHeal_061626.png",width=4, height=5)



##
#####032126 save this######
#manually find max of each series

library(ggrepel)
library(ggpmisc)
wurtzel.head.tail.heal.anno.analysis3<-wurtzel.head.tail.heal.anno.analysis2 %>% rowwise() %>% mutate(max.head = max(c_across(3:13))) %>% mutate(max.tail = max(c_across(15:25))) %>% mutate(max.heal=max(c_across(27:37)))

library(ggpubr)

#max_diff paired dotplots Fig S1C
ggpaired(wurtzel.head.tail.heal.anno.analysis3, "max.head","max.heal")
ggpaired(wurtzel.head.tail.heal.anno.analysis3, "max.head","max.heal", color="condition", line.color="black", line.size=0.1)
ggsave("head.vs.heal.paired.pdf")

ggpaired(wurtzel.head.tail.heal.anno.analysis3, "max.tail","max.heal")
ggpaired(wurtzel.head.tail.heal.anno.analysis3, "max.tail","max.heal", color="condition", line.color="black", line.size=0.1)
ggsave("tail.vs.heal.paired.pdf")


#max_diff histograms Fig S1D
wurtzel.head.tail.heal.anno.analysis3 <- wurtzel.head.tail.heal.anno.analysis3 %>% mutate(max.diff.head = max.head-max.heal) %>% mutate(max.diff.tail = max.tail-max.heal)

wurtzel.head.tail.heal.anno.analysis3 %>% ggplot(aes(x=max.diff.head))+
  geom_histogram(bins=45)+
  theme_classic()+xlim(-2,2)+ylim(0,20)+
  theme(axis.text= element_text(face = "bold"),text = element_text(size = 20))
ggsave("max.diff.head.histogram.pdf")

wurtzel.head.tail.heal.anno.analysis3 %>% ggplot(aes(x=max.diff.tail))+
  geom_histogram(bins=45)+
  theme_classic()+xlim(-2,2)+ylim(0,20)+
  theme(axis.text= element_text(face = "bold"),text = element_text(size = 20))
ggsave("max.diff.tail.histogram.pdf")


#### histograms of max expression differences (regen vs healing) in Fig 1D and S4
wurtzel.head.tail.heal.anno.analysis3<-wurtzel.head.tail.heal.anno.analysis3 %>% mutate(highlight=wurtzel.name %in% c("wnt1","notum"))

wurtzel.head.tail.heal.anno.analysis3 %>% ggplot(aes(y=fct_reorder(wurtzel.name,max.diff.head),x=max.diff.head,fill=highlight))+
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_manual( values = c( "TRUE"="tomato", "FALSE"="gray" ), guide = FALSE )+
  #  xlim(0,5)+ylim(0.5,2)+
  geom_abline(slope=10000000000, color="blue", linetype="dashed")+
  theme(axis.text.x=element_text(angle = 90, hjust = 0))+
  scale_x_continuous(position = "top")+
  theme(axis.text= element_text(face = "bold"),text = element_text(size = 6))
ggsave("maxdiffhead.barplot.pdf",width=2, height=8)


wurtzel.head.tail.heal.anno.analysis3 %>% ggplot(aes(y=fct_reorder(wurtzel.name,max.diff.tail),x=max.diff.tail,fill=highlight))+
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_manual( values = c( "TRUE"="tomato", "FALSE"="gray" ), guide = FALSE )+
  #  xlim(0,5)+ylim(0.5,2)+
  #  geom_line(aes(y=lwr), color = "", linetype = "dashed")+
  geom_abline(slope=10000000000, color="blue", linetype="dashed")+
  #  geom_abline(a=0,color = "blue")+
  theme(axis.text.x=element_text(angle = 90, hjust = 0))+
  scale_x_continuous(position = "top")+
  theme(axis.text= element_text(face = "bold"),text = element_text(size = 6))
ggsave("maxdifftail.barplot.pdf",width=2, height=8)


wurtzel.head.tail.heal.anno.analysis3 %>% ggplot(aes(x=fct_reorder(wurtzel.name,max.diff.head),y=max.diff.head,fill=highlight))+
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_manual( values = c( "TRUE"="tomato", "FALSE"="gray" ), guide = FALSE )+
  #  xlim(0,5)+ylim(0.5,2)+
  geom_abline(slope=0, color="blue", linetype="dashed")+
  theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.1))+
  theme(axis.text= element_text(face = "bold"),text = element_text(size = 6))
ggsave("maxdiffhead.barplot2.pdf",width=8, height=2)

wurtzel.head.tail.heal.anno.analysis3 %>% ggplot(aes(x=fct_reorder(wurtzel.name,max.diff.tail),y=max.diff.tail,fill=highlight))+
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_manual( values = c( "TRUE"="tomato", "FALSE"="gray" ), guide = FALSE )+
  geom_abline(slope=0, color="blue", linetype="dashed")+
  #  xlim(0,5)+ylim(0.5,2)+
  theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.1))+
  theme(axis.text= element_text(face = "bold"),text = element_text(size = 6))
ggsave("maxdifftail.barplot2.pdf",width=8, height=2)


wurtzel.head.tail.heal.anno.analysis3<-wurtzel.head.tail.heal.anno.analysis3 %>% mutate(norm.maxdiff.head=(max.diff.head)/max.head) %>% mutate(norm.maxdiff.tail=(max.diff.tail)/max.tail)

wurtzel.head.tail.heal.anno.analysis3 %>% ggplot(aes(x=fct_reorder(wurtzel.name,norm.maxdiff.head),y=norm.maxdiff.head,fill=highlight))+
  geom_bar(stat="identity")+
  theme_classic()+
  #  xlim(0,5)+ylim(0.5,2)+
  scale_fill_manual( values = c( "TRUE"="tomato", "FALSE"="gray" ), guide = FALSE )+
  theme(axis.text.x=element_text(angle = 90, hjust = 0))+
  theme(axis.text= element_text(face = "bold"),text = element_text(size = 8))
ggsave("norm.maxdiff.head.barplot.pdf",width=8, height=4)




