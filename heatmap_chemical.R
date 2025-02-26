setwd("C:/Users/martk/Documents/wur/masterthesis/bioinformatica/bioinfR/Output/Chemical")

#Read all packages required for script
library(ggfortify)
library(pheatmap)
library(tidyr)
library(dplyr)

#DF <- read.csv("C:/Users/martk/Documents/wur/masterthesis/bioinformatica/bioinfR/Output/Tree/new-soks-unalignedgaps.csv")

#Load results calculated in calculateGlobalProps.R
DF <- results


for(property in c("windowCharge","WindowHydropathy")){
DFpivot <- DF %>% pivot_wider(id_cols=AlignmentPosition, names_from=Sequence_Name, values_from=property)
DFpivot <- as.data.frame(DFpivot)
rownames(DFpivot) <- DFpivot$AlignmentPosition
DFpivot$AlignmentPosition <- NULL
#DFpivot[is.na(DFpivot)] <- 0
DFpivot <- t(DFpivot)

#SETTINGS AND COLORS OF HEATMAP
fraction_reduced_black = 0
number_of_colors <- 200
max_cutoff <- max(DFpivot[!is.na(DFpivot)])
min_cutoff <- min(DFpivot[!is.na(DFpivot)])
if(min_cutoff > 0){min_cutoff <- 0}

#fraction_pos = max(DFpivot)/(max(DFpivot)-min(DFpivot))
#fraction_neg = abs(min(DFpivot)/(max(DFpivot)-min(DFpivot)))

fraction_pos = max_cutoff/(max_cutoff-min_cutoff)
fraction_neg = abs(min_cutoff/(max_cutoff-min_cutoff))

colors <- c("#56b4e9", "#FFFFFF", "#d55e00", "#FFFFFF")
ramp_neg <- colorRampPalette(colors[1:2])(number_of_colors*fraction_neg * (1+fraction_reduced_black))[1:(number_of_colors*fraction_neg)]
ramp_pos <- colorRampPalette(colors[3:4])(number_of_colors*fraction_pos * (1+fraction_reduced_black))[1:(number_of_colors*fraction_pos)]
cols <- c(ramp_neg, rev(ramp_pos))

x_interval <- 50
x <- as.numeric(colnames(DFpivot))
x_labels <- ifelse(x%%x_interval == 0, x, "")
#y_labels <- select(newDF, "Kingdom")

heatmap <- pheatmap(DFpivot, cluster_rows = T ,cluster_cols = F,cellwidth = 1, cellheight = 2,color = cols, cex=0.4, labels_col=x_labels, angle_col=0)
ggsave(paste(property,"_heatmap.png",sep=""), plot = heatmap, width = 12, height = 39.5, dpi=400, limitsize = FALSE, bg = "white")
dev.off()}
