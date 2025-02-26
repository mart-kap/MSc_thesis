library("seqinr")
library("ape")
library("ggtree")
library('ggnewscale')
library('ggtreeExtra')
library("ggplot2")
library("stringr")

setwd("C:/Users/martk/Documents/wur/masterthesis/bioinformatica/bioinfR/Output/Tree")

#Read taxonomy files here for colorcoding
DF_Taxonomy <- read.csv("C:/Users/martk/Documents/wur/masterthesis/bioinformatica/Slots_Taxonomy_modifiedlist.csv",sep=";")
DF_Taxonomy$Kingdom[DF_Taxonomy$Kingdom == ""] <- "Protists"
rownames(DF_Taxonomy) <- DF_Taxonomy$ID

#Read alignment here (for current trees we first use alignment with ALL
#HMMERs found, then second the one returned from the Python file + Clustal
# Omega website)
Alignment <- read.alignment("C:/Users/martk/Downloads/clustalo-I20241118-135337-0780-58647974-p1m.fa", "fasta")

#Calculate distanceMatrix (possibly your distance matrix can also be used)
#Edit; yes that is possible. Read your DisMat (the output as given by the browser)
# This output will be the "Distances_Matrix" -> so skip assigning "Distances"
# You might have to add colnames(Distances_Matrix) <- rownames(Distances_Matrix)
# Also skip the first row
Distances <- dist.alignment(Alignment, matrix = "identity",gap=TRUE)
Distances_Matrix <- as.matrix(Distances)
#heatmap(Distances_Matrix)

#Calculate tree (also here, the ClustalO output might work, but I )
my_nj <- ape::nj(Distances)
nodedf <- data.frame(ID=my_nj$tip.label)
nodedf <- merge(nodedf,DF_Taxonomy,by.x="ID",by.y=0)

#Graphical part of treemaking + make a "TreeDF" dataframe for selective graphs
Tree <- ggtree(my_nj, ladderize = TRUE, branch.length = "none", layout="circular")
TreeDF <- as.data.frame(Tree$data)
TreeDF$Taxid <- str_split_fixed(TreeDF$label,"_", 2)[,1]
TreeDF <- merge(TreeDF,DF_Taxonomy,by="Taxid",all=T)

# the TreeDF is altered according to our needs, so create TreeDF_temp
TreeDF_temp <- TreeDF[,c("label","Kingdom")]
TreeDF_temp <- TreeDF_temp[!is.na(TreeDF_temp$Kingdom),]


#Set colors for colorcoding
Colors <- unique(TreeDF_temp$Kingdom)
Colors = setNames(c("#0072b2","#009e73","#cc79a7","#d55e00", "#56b4e9"), Colors) 

#Here the kingdom is colorcoded (you need more colors due to the Asgard)
#When running the first round, run the code below first to filter the actins.
Tree <- Tree + new_scale_fill() + geom_fruit(data=TreeDF_temp, offset=0.055, width=5,geom=geom_tile, mapping=aes(y=label, fill=Kingdom)) + scale_fill_manual(values=Colors)
ggsave("Tree1_General.png", Tree, width=15, height=15, bg = "transparent")
TreeDF <- TreeDF[TreeDF$isTip,]
write.csv(TreeDF, "treedf.csv")


#ACTIN_DF <- TreeDF[TreeDF$angle < 196.2925,]
#write.csv(ACTIN_DF, "actindf.csv")

#ACTINlike_DF <- TreeDF[TreeDF$angle >= 196.2925,]
#write.csv(ACTINlike_DF, "actinlikedf.csv")

#When running for the first time, run this first before making the tree.
#It removes all proteins we thought are NOT actins from the tree.
#They were clustered in certain regions (hence the filter by angle).
#TreeDF['isACTIN'] <- FALSE
#TreeDF$isACTIN[TreeDF$angle < 196.2925] <- TRUE
#TreeDF$isACTIN[TreeDF$angle < 156.9298 & TreeDF$angle > 153.1499] <- FALSE
#TreeDF <- TreeDF[TreeDF$isTip,]
write.csv(TreeDF, "treedf.csv")

#Below are 3 examples to create custom trees to highlight certain aspects

#Here I labelled all species + uniprot accession

#So first I created a new column with the name (paste species+uniprot)
TreeDF$comb <- paste(TreeDF$Species,TreeDF$label,sep="_")

#Then I created a copy of the main tree created above
Tree2 <- Tree

#copy the DataFrame as well
TreeDF_temp <- TreeDF
TreeDF_temp <- TreeDF_temp[!is.na(TreeDF_temp$label),]

#Label all items with the values of the "comb" column
Tree2 <- Tree2 + geom_tiplab(data=TreeDF_temp, mapping=aes(node=node, label=comb), size=1.5)
ggsave("Tree2_Species.png", Tree2, width=60, height=60, bg = "transparent", limitsize = FALSE)


