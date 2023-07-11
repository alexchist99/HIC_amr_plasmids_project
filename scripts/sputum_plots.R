
setwd("/home/alexclear/HIC_amr_plasmids_project/sputum_analysis/")


#YOU NEED TO INSTALL "WEBR" PACKAGE SEPARETELY!
#I used this link for that goal: https://rpubs.com/cardiomoon/398623; otherwise it does not work ;(
#But you can try if you are brave and patient enough

library(dplyr)
library(ggplot2)
library(webr)
library(ggpubr)
library(optparse)
options(scipen = 999)


#parse arguments
option_list = list(
  make_option(c("-f", "--sample1"), type = "character", default = NA,
              help = "input sample1 name to retrieve all necessary files"),
  make_option(c("-s", "--sample2"), type = "character", default = NA,
              help = "input sample2 name to retrieve all necessary files"),
  make_option(c("-m", "--output_mapped"), type = "character", default = NA,
              help = "output percentage of mapped reads"),
  make_option(c("-p", "--output_pie1"), type = "character", default = NA,
              help = "output pie figure one"),
  make_option(c("-o", "--output_pie2"), type = "character", default = NA,
            help = "output pie figure two"))


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#assign
SD2 = opt$sample1
SD8 = opt$sample2
mapping = opt$output_mapped
pieSD2 = opt$output_pie1
pieSD8 = opt$output_pie2



#-------------------General function for parsing and merging files ---------------------------------------
merge_fun <- function(sample){ 
  class = read.table(paste0(sample,"_classif_intersected.txt"))
  amr = read.csv(paste0(sample,"_amr_intersected.txt"))
  intersected = read.table(paste0(sample,"_intersected_nodes_list.txt"))
  all_cont = read.table(paste0(sample,"_all_nodes_list.txt"))
  
  #preprocess amr info
  colnames(amr) <- c("V1","Status","Gene","Class","Action","Protein")
  amr$V1 <- sub('_[^_]+$','',amr$V1)
  


  class_intersected <- merge(intersected,class)
 
  
  
  #delete singletons and get % of mapped reads from the rest taxa 
  threshold <- round(nrow(intersected)*0.05) #5% 
  dt <- class_intersected %>% count(V2) %>% filter(n>=threshold) %>% 
    as.data.frame() %>% 
    mutate(percent_mapped = round(n/nrow(all_cont)*100,2))
  dt$sample = sample
  out_taxa <- dt$V2 %>% unlist() %>% as.vector()
  
  
  
  dt_amr <-  merge(amr,class_intersected,all.x = T)
  dt_amr <- dt_amr[dt_amr$V2 %in% out_taxa,]
 
  dt_amr$sample = sample
  
  return(list(dt,dt_amr))
  }

dt_mapping <- rbind(merge_fun(SD2)[[1]],merge_fun(SD8)[[1]])


dt_mapping$sample <- factor(dt_mapping$sample)

levels(dt_mapping$sample) <- c("B-9782","B-9817")

#-------------------Plot for mapping sputum reads % ---------------------------------------
pdf(mapping)
ggplot(dt_mapping, aes(y=sample,x = percent_mapped,fill = V2))+
  geom_bar(stat="identity")+ coord_flip()+ xlab("Mapped reads,%")+ 
  ylab("Sample")+
  theme_classic()+scale_fill_discrete(name = "Taxa")+
  theme(legend.position="bottom")+
  ggtitle("Percentage of mapped sputume reads to stool plasmidome for two samples")
dev.off()

#-------------------Amr plots----------------------------------------------------------

dt_amr <- rbind(merge_fun(SD2)[[2]],merge_fun(SD8)[[2]])

dt_amr$sample <- factor(dt_amr$sample)
levels(dt_amr$sample) <- c("B-9782","B-9817")
colnames(dt_amr)[7] <- "Plasmid"

sd2 = dt_amr %>%  filter(sample=="B-9782")
sd2$Gene = factor(sd2$Gene)
sd2$Plasmid = factor(sd2$Plasmid)

pdf(pieSD2)
plot2 <- PieDonut(sd2, aes(pies = Plasmid, donuts = Gene), 
         title = "AMR genes shared by sputume and stool plasmidome (B-9782/SD2)",color="white")
dev.off()


sd8 = dt_amr %>% filter(sample=="B-9817") 
sd8$Gene = factor(sd8$Gene)
sd8$Plasmid = factor(sd8$Plasmid)

pdf(pieSD8)
plot3 <-PieDonut(sd8, aes(pies = Plasmid, donuts = Gene), 
         title = "AMR genes shared by sputume and stool plasmidome (B-9817/SD8)")+theme(panel.border = element_blank())
dev.off()


# pdf("hzz.pdf")
# ggarrange(plot1,plot2,plot3)
# dev.off()



