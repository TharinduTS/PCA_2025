# PCA_2025

I am using the thinned vcf from Ben

trop_only_allchrs_concat_maxmissingcount_0_biallelic_genoqual30.vcf.gz

copied this file to my mac and use R to calculate and plot PI

```
#https://github.com/zhengxwen/SNPRelate/issues/13
#http://corearray.sourceforge.net/tutorials/SNPRelate/
#https://github.com/zhengxwen/SNPRelate/wiki/Preparing-Data

# ***************   These are the things you should edit  **********************

input_file_name<-"trop_WGS_no_cal_mello_niger_all_chrs_thinned_5000.vcf.recode.vcf"
plot_title_to_use<-"PCA_no_Niger"
saving_file_name<-"No_cal_no_mello_no_niger_whole_genome.pdf"
sample_data_file_name<-"sample_list_no_cal_mello_niger.txt"
summary_table_name<-"No_cal_no_mello_no_niger_whole_genome.txt"

# if you added or changed colors in excel sheet with new colors, you have to edit line for colors below  as well
# ******************************************************************************
# ******install all the packages needed if you are running it on computecanada*********

# install.packages("devtools")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("gdsfmt")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("SNPRelate")
# install.packages("plyr")
# install.packages("ggplot2")
# install.packages("ggrepel")
# ************

library("devtools")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("gdsfmt")
library(gdsfmt)
library(SNPRelate)
#install.packages("ggpp")
#install.packages("ggrepel")
library(ggpp)

# set working directory to current script directory 
# **uncomment this if you use this in local computer. KEEP COMMENTED OUT IF YOU ARE ON COMPUTECANADA) ***
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#enter input vcf filename here
vcf.fn <- input_file_name


# ************* if you retry this, you may have to restart r session before trying ************

snpgdsVCF2GDS(vcf.fn, "test.gds", method="copy.num.of.ref",ignore.chr.prefix = "chr")

# *********************************************************************************************

#snpgdsSummary("test.gds")

#genofile = openfn.gds("test.gds", readonly=FALSE)
genofile = snpgdsOpen("test.gds", readonly=FALSE)

#samp.annot<-data.frame(pop.group = c("brunnescens","hecki","hecki","hecki","hecki","hecki","hecki","maura","maura","maura","maura","maura","maura","Borneo","Sumatra","Malay","Sumatra","Borneo","Borneo","Borneo","siberu","nigra","nigra","nigrescens","nigra","nigrescens","ochreata","ochreata","ochreata","togeanus","togeanus","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana"))

#trying to read sample names and populations from a tsv file instead typing all here
#name the tsv file as samples.list with sample name in the first column , population in second column and colour in third column
sample_data<-read.table(sample_data_file_name)

#create a sample name array in the same order
sample_name_order<-as.array(sample_data$V1)
population_order<-as.array(sample_data$V2)
color_order<-as.array(sample_data$V3)

# get simpler names  for those complex sample names
library(plyr)
#first renaming the sample names which does not follow the specific format manually**DO NOT USE"_" HERE AS SPLIT WILL USE THIS IN NEXT LINE

#copy sample name list to edit
full_sample_names<-sample_name_order

#rename the samples that does not follow that format manually
simple_sample_list_changing<-mapvalues(full_sample_names, from = c("Vred_8_Vred_GCTGTGGA_cuttrim_sorted.bam","JM_no_label1_Draken_CCACGT_cuttrim_sorted.bam","JM_no_label2_Draken_TTCAGA_cuttrim_sorted.bam","946_Draken_TCGTT_cuttrim_sorted.bam","993_Draken_GGTTGT_cuttrim_sorted.bam","2014_Inhaca_10_Inhaca_ATATGT_cuttrim_sorted.bam","2014_Inhaca_150_Inhaca_ATCGTA_cuttrim_sorted.bam","2014_Inhaca_152_Inhaca_CATCGT_cuttrim_sorted.bam","2014_Inhaca_24_Inhaca_CGCGGT_cuttrim_sorted.bam","2014_Inhaca_38_Inhaca_CTATTA_cuttrim_sorted.bam","2014_Inhaca_52_Inhaca_GCCAGT_cuttrim_sorted.bam","2014_Inhaca_65_Inhaca_GGAAGA_cuttrim_sorted.bam"), 
                                       to = c("Vred8","JM1","JM2","Draken946","Draken993","Inhaca10","Inhaca150","Inhaca152","Inhaca24","Inhaca38","Inhaca52","Inhaca65"))

#convert facter list into chrs to rename
s_list_chr<-as.character(simple_sample_list_changing)

#extract sample IDs from names
shortened_sample_list<-sapply(strsplit(s_list_chr,split = "_"),`[`, 1)

#samp.annot<-data.frame(pop.group = c("brunnescens","hecki","hecki","hecki","hecki","hecki","hecki","maura","maura","maura","maura","maura","maura","Borneo","Sumatra","Malay","Sumatra","Borneo","Borneo","Borneo","siberu","nigra","nigra","nigrescens","nigra","nigrescens","ochreata","ochreata","ochreata","togeanus","togeanus","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana"))
samp.annot<-data.frame(pop.group = population_order)



#and then load
add.gdsn(genofile, "sample.annot", samp.annot)

snpgdsSummary("test.gds")

# LD prinnung
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,  method = c("composite"),missing.rate=0, verbose = TRUE)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))


                    # or without LD prunning
#pca <- snpgdsPCA(genofile, num.thread=2)
#pc.percent <- pca$varprop*100
#head(round(pc.percent, 2))


tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the third eigenvector
                  stringsAsFactors = FALSE)
head(tab)

#plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
#text(tab$EV2, tab$EV1,labels=tab$sample.id, cex= 0.4)

library(ggplot2)
#ggplot(...)+...+ theme(axis.text.x = element_text(angle=60, hjust=1))
#devtools::install_github("slowkow/ggrepel")
library(ggrepel)
library(tidyverse)

pdf(saving_file_name,w=8, h=8, version="1.4", bg="transparent")
tab$Species <- population_order
tab$samp.color <- color_order
tab$samp.fieldid <- shortened_sample_list

# get unique colours and pop names as labels to use in next line
pop_colors<-unique(color_order)

# you will have to use a manual array if you are going to apply same color to more than one pop
pop_labels<-unique(population_order)

# to apply correct colors for the values
#getting same colors for the values
#colors_to_apply<-pop_colors

#creating the array for values
#col_vals<-as.array(paste("'",pop_colors,"'","=","'",colors_to_apply,"'",sep=""))

#pass percentage to a variable
pc1_perc<-format(round(pc.percent[1],2),nsmall = 2)
pc2_perc<-format(round(pc.percent[2],2),nsmall = 2)
pc3_perc<-format(round(pc.percent[3],2),nsmall = 2)

#adding a column with shortened pop names
tab$shortened_pops<-tab$samp.color
#shorten names 
tab$shortened_pops[tab$shortened_pops=="Ghana_East"]<-"GE"
tab$shortened_pops[tab$shortened_pops=="Ghana_West"]<-"GW"
#tab$shortened_pops[tab$shortened_pops=="Lab_tads"]<-""
tab$shortened_pops[tab$shortened_pops=="Sierra_Leone"]<-"SL"
tab$shortened_pops[tab$shortened_pops=="Ivory_coast"]<-"IC"
tab$shortened_pops[tab$shortened_pops=="Liberia"]<-"LB"
tab$shortened_pops[tab$shortened_pops=="Nigeria"]<-"NG"

#remove extra labels from each pop
one_val_per_pop<-tab %>% distinct(shortened_pops, .keep_all = TRUE)
one_val_per_pop<-one_val_per_pop%>%select(sample.id,shortened_pops)

updated_tab<-merge(tab, one_val_per_pop, by = 1, all = TRUE)
updated_tab$shortened_pops.y %>% replace_na('none')

h<-ggplot(data=updated_tab, aes(x=-EV1,y=-EV2, label = sample.id, color = samp.color,shape = samp.color)) +
  # label axis 
  xlab(paste("PC 1 -","  ",pc1_perc,"%",sep = ""))+
  ylab(paste("PC 2 -","  ",pc2_perc,"%",sep = ""))+
  #ggtitle(plot_title_to_use)+
  #title(plot_title_to_use,cex=1)+
  #labs(x=expression(paste("-Eigenvector 1","  ",pc1_perc,"%",sep = "")), y=expression("-Eigenvector 2"),title =plot_title_to_use,cex=1) +
  # legend details
  scale_shape_manual(name="Population",values = c("Ghana_East"=16,"Ghana_West"=16,"Ivory_coast"=16,"Nigeria"=16,"Sierra_Leone"=16,"Nigeria"=16,"Liberia"=16),breaks = tab$samp.color,labels = tab$samp.color)+
  # ************
  
  # if you add more colors in excel sheet, you have to add them here as well
  scale_colour_manual(name="Population", values = c("Ghana_East"="red","Ivory_coast"="lightblue","Ghana_West"="pink","Sierra_Leone"="green","Nigeria"="darkblue","Liberia"="orange"),breaks = tab$samp.color,labels = tab$samp.color)+
  
  # ************
  # add points and fieldID labels
  #geom_text_repel(aes(-EV1,-EV2, label=(shortened_pops.y)), size=6, point.padding = unit(2, "lines"),max.overlaps = 22,position = position_nudge_to(x = c(0.25, -0.15))) + 
  geom_jitter(size=6,stroke=1.5) + 
  #geom_label(aes(fill = factor(samp.color)), colour = "white", fontface = "bold") +
  # change to cleaner theme
  theme_classic(base_size = 16) +
  # make it clean
  theme_bw()+ theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  
  #this is for theme inside plot
  # theme(legend.key.size = unit(0.15, 'cm'), #change legend key size
  #      legend.key.height = unit(0.15, 'cm'), #change legend key height
  #      legend.key.width = unit(0.15, 'cm'), #change legend key width
  #      legend.title = element_text(size=11), #change legend title font size
  #      legend.text = element_text(size=11,face = "italic"), #change legend text font size
  #      legend.position = c(0.2, .15))+
# italicize species names
#theme(legend.text = element_text(face="italic"))+hel
# make the text bigger
#theme(text = element_text(size=20)) +
# move the legend
#theme(legend.position = c(.18, .15)) +
# add space around axis
theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt"),size = 20),
      axis.text.y = element_text(margin=margin(10,10,10,10,"pt"),size = 20)) +
  # remove boxes around legend symbols
  theme(legend.key = element_blank())+
  #these are for annotations. not used here
  annotate(geom = "text", x = 0.03, y = .23, label = "", color = "black", angle = 75, size=7)+
  annotate(geom = "text", x = -.12, y = -.10, label = "", color = "black", angle = 0, size=7)+
  annotate(geom = "text", x = -.25, y = .07, label = "", color = "black", angle = -15, size=7)+
  theme(text = element_text(size = 16))
h

ggsave("PC1and2.pdf",h,height = 6,width = 8)

#. doing same for ev1 and ev3

i<-ggplot(data=tab, aes(x=-EV1,y=-EV3, label = sample.id, color = samp.color,shape = samp.color)) +
  # label axis 
  xlab(paste("PC 1 - ","  ",pc1_perc,"%",sep = ""))+
  ylab(paste("PC 3 - ","  ",pc3_perc,"%",sep = ""))+
  ggtitle(plot_title_to_use)+
  #title(plot_title_to_use,cex=1)+
  #labs(x=expression(paste("-Eigenvector 1","  ",pc1_perc,"%",sep = "")), y=expression("-Eigenvector 2"),title =plot_title_to_use,cex=1) +
  # legend details
  scale_shape_manual(name="Population",values = c("Ghana_East"=0,"Ivory_coast"=16,"Ghana_West"=17,"Sierra_Leone"=1,"Nigeria"=19,"Liberia"=20),breaks = tab$samp.color,labels = tab$samp.color)+
  # ************
  
  # if you add more colors in excel sheet, you have to add them here as well
  scale_colour_manual(name="Population", values = c("Ghana_East"="red","Ivory_coast"="lightblue","Ghana_west"="pink","Sierra_Leone"="green","Nigeria"="darkblue","Liberia"="orange"),breaks = tab$samp.color,labels = tab$samp.color)+
  
  # ************
  
  # add points and fieldID labels
  geom_text_repel(aes(-EV1,-EV3, label=(samp.fieldid)), size=6, point.padding = unit(0.5, "lines"),max.overlaps = 0) + geom_point(size=6) + 
  #geom_label(aes(fill = factor(samp.color)), colour = "white", fontface = "bold") +
  # change to cleaner theme
  theme_classic(base_size = 16) +
  # make it clean
  theme_bw()+ theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  
  #this is for theme inside plot
  # theme(legend.key.size = unit(0.15, 'cm'), #change legend key size
  #      legend.key.height = unit(0.15, 'cm'), #change legend key height
  #      legend.key.width = unit(0.15, 'cm'), #change legend key width
  #      legend.title = element_text(size=11), #change legend title font size
  #      legend.text = element_text(size=11,face = "italic"), #change legend text font size
  #      legend.position = c(0.2, .15))+
# italicize species names
#theme(legend.text = element_text(face="italic"))+hel
# make the text bigger
#theme(text = element_text(size=20)) +
# move the legend
#theme(legend.position = c(.18, .15)) +
# add space around axis
theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt"),size = 11),
      axis.text.y = element_text(margin=margin(10,10,10,10,"pt"),size = 11)) +
  # remove boxes around legend symbols
  theme(legend.key = element_blank())+
  #these are for annotations. not used here
  annotate(geom = "text", x = 0.03, y = .23, label = "", color = "black", angle = 75, size=7)+
  annotate(geom = "text", x = -.12, y = -.10, label = "", color = "black", angle = 0, size=7)+
  annotate(geom = "text", x = -.25, y = .07, label = "", color = "black", angle = -15, size=7)
i
# add PC% to the summary tables

tab$PC1_perc<-pc1_perc
tab$PC2_perc<-pc2_perc
tab$PC3_perc<-pc3_perc

#save the dataframe to be used in the all samples together plot
write.table(tab, file=summary_table_name, quote=FALSE, sep='\t', col.names = NA)

dev.off()




closefn.gds(genofile)
```
