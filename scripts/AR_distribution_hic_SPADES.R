library(ggplot2)
library(tidyr)
library(dplyr)
library(ape)
library(gplots)
library(optparse)
library(stringr)
options(scipen = 999)


#parse arguments
option_list = list(
  make_option(c("-m", "--matrix"), type = "character", default = NA,
              help = "input adjaency matrix of metagenomic MAG-MAG contacts"),
  make_option(c("-t", "--tree"), type = "character", default = NA,
              help = "input tree of taxa"),
  make_option(c("-a", "--outputh"), type = "character", default = NA,
              help = "output heatmap with all samples"),
  make_option(c("-r", "--outputr"), type = "character", default = NA,
              help = "output arg distribution"))


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


adj_data = opt$matrix
sample_tree = opt$tree
output_file_name1 = opt$outputh
output_file_name2 = opt$outputr


# #download data 

dt = read.csv(adj_data)

# #downloading metagenomic tree
rep_tree <- read.tree(sample_tree)


# # Preprocessing for mag-mag heatmap

#first take that contacts where plasmid contig (according to ncbi annotation) has its host (MAG) - in case of avoiding spurious contatcs
#and find false positive amr genes
#but keep the initial number of hosts as in the tree:
taxa_list = table(c(dt$from,dt$to))


rare = names(taxa_list[taxa_list<250])

dt$prevalent = ifelse((dt$from %in% rare | dt$to %in% rare),"no","yes")

dt1 <- dt %>% filter(dt$compliment=="yes") 

dt2 <- dt %>% filter(dt$prevalent=="no" & dt$compliment=="no") 

dt <- rbind(dt1,dt2)



dt$ARG <- ifelse(startsWith("NODE",dt$Contig),"-","*")
aa1  <- dt %>% dplyr::count(from,to,sample)  
arg1 <- dt %>% dplyr::filter(ARG == "*") %>% dplyr::select(from,to,ARG,sample) %>% unique()
final_dt1 <- merge(aa1,arg1, by.y = c("from","to","sample"),all.x = TRUE) 

final_dt1[is.na(final_dt1)] <- " "
final_dt1$from <- sub("^[^_]*__", "",final_dt1$from)
final_dt1$to <- sub("^[^_]*__", "",final_dt1$to)
final_dt1 <- final_dt1[-1,]



#MATRIX MAKE for mag-mag heatmap

from_to <- unique(c(final_dt1$from,final_dt1$to))

#mag&mag contacts
m <- matrix(data = rep(0,length(from_to)*length(from_to)),nrow = length(from_to),
            ncol= length(from_to))

# #ARG presence
m1 <- matrix(data = rep(0,length(from_to)*length(from_to)),nrow = length(from_to),
            ncol= length(from_to))


# #Give names
colnames(m) <- from_to 
rownames(m) <- from_to 



colnames(m1) <- from_to 
rownames(m1) <- from_to 


#fun for avoiding errors while filling the matrixes
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

#fill the matrixes

for(j in unique(c(final_dt1$to))){
  for(i in unique(c(final_dt1$from))){
    m[i,j] <- ifelse(
      is.integer0(final_dt1[final_dt1$to==i & final_dt1$from==j,]$n),
      0,
      final_dt1[final_dt1$to==i & final_dt1$from==j,]$n)
    #print(is.integer0(final_dt1[final_dt1$to==i & final_dt1$from==j,]$n))
  }
}



# # for log transformation (not to confuse 1 as 1 contact and 1 as 0+1 for log)
m <- m+1

for(j in unique(c(final_dt1$to))){
  for(i in unique(c(final_dt1$from))){
    m1[i,j] <- ifelse(
      is.integer0(final_dt1[final_dt1$to==i & final_dt1$from==j,]$ARG),
      " ",
      final_dt1[final_dt1$to==i & final_dt1$from==j,]$ARG)
    #print(is.integer0(final_dt1[final_dt1$to==i & final_dt1$from==j,]$n))
  }
}


#TREE 
#Some preprocessing...

rep_tree$tip.label <- gsub("\\..*","",rep_tree$tip.label)
tree_tips <- rep_tree$tip.label

#root
rep_tree<- root(rep_tree,
                   resolve.root = T,
                   outgroup="Methanobrevibacter_A")


#can’t have any branch lengths of zero or downstream commands will collapse those nodes…
rep_tree$edge.length[which(rep_tree$edge.length == 0)] <- 0.00001
rep_tree_um <- chronopl(rep_tree,
                        lambda = 0.1,
                        tol = 0)
rep_tree_d <- as.dendrogram(as.hclust.phylo(rep_tree_um))


##HEATMAPS WITH TREE
#make the heatmap symmetric regarding to the tree topology
#some tips have not the same order as in the tree - change it:
tree_tips1 = c(tree_tips[3:length(tree_tips)],rev(tree_tips[1:2]))
print(tree_tips1)

m <- m[tree_tips,tree_tips1]
m1 <- m1[tree_tips,tree_tips1]



#Final heatmap

# #Select pallete
color <- colorRampPalette(c('white','blue','red'))(100)


pdf(output_file_name1)
heatmap.2(log(m),Colv = F ,Rowv = rep_tree_d,dendrogram="row",
          margins=c(10,10),
          trace='none',
          keysize=1.2,
          key.title="Number of contacts (log)",
          col=color,
          cellnote=m1,
          notecol="black",
          main="InterMag plasmid contacts",
          xlab="MAG",
          ylab="MAG",
          density.info="density")
dev.off() 



#ARG:

dt$Best_Hit_ARO <-sub("^\\d+|\\d+$","",
                      gsub("'","",
                                            sub("-.*","",
                                                gsub("\\(|\\)", "",  dt$Best_Hit_ARO))))

dt$Best_Hit_ARO <- str_split(dt$Best_Hit_ARO,pattern = " ",simplify = TRUE)[,1]
#dt$Drug.Class <- str_split(dt$Drug.Class,pattern = " ",simplify = TRUE)[,1]

#did not handle regex here -did manually()
dt$Best_Hit_ARO <- sub("van[A-Z]*","van",dt$Best_Hit_ARO)
dt$Best_Hit_ARO <- sub("tet[A-Z]*","tet",dt$Best_Hit_ARO)
dt$Best_Hit_ARO <- sub("evg[A-Z]*","evg",dt$Best_Hit_ARO)
dt$Best_Hit_ARO <- sub("Qnr[A-Z]*","Qnr",dt$Best_Hit_ARO)

dt$from <- sub("g__","",dt$from )
dt$to <- sub("g__","",dt$to)
dt$fromto <- factor(paste(dt$from,dt$to,sep="\n"))

d <- dt %>% dplyr::filter(ARG=="*")  %>% count(fromto,Best_Hit_ARO,sort = T)
#%>% filter(n>1) 


#%>% filter(n>1)
#forcats::fct_rev(y)
pdf(output_file_name2)
ggplot(d, aes(forcats::fct_reorder(fromto,n,.desc = T), Best_Hit_ARO,size = n)) +
  geom_point(shape = 21, stroke = 0,fill="black") +
  scale_x_discrete(position = "top") +
  scale_radius(range = c(1, 10)) +
  #scale_fill_discrete(low = "orange", high = "blue") +
  theme_minimal() +
  theme(legend.position = "bottom", 
        text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.1),
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  guides(size = guide_legend(override.aes = list(fill = "black", color = "black", stroke = .25), 
                             label.position = "bottom",
                             title.position = "right", 
                             order = 1),
         fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
  labs(size = "Area = Number of shared ARG genes" ,fill = "black",x = NULL, y = NULL)
  dev.off() 