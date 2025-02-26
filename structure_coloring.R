#Run the tree script before this, TreeDF is needed.

# Load required libraries
# install.packages("dplyr")
library(dplyr)
library(Biostrings)
library(idpr)
library(ggfortify)
library(pheatmap)
library(tidyr)
library(dplyr)
# BiocManager::install("idpr")
# https://bioconductor.org/packages/devel/bioc/vignettes/idpr/inst/doc/chargeHydropathy-vignette.html

setwd("C:/Users/martk/Documents/wur/masterthesis/bioinformatica/bioinfR/Output/Structures")

window_size = 9
fasta_file_path = "C:/Users/martk/Documents/wur/masterthesis/bioinformatica/ClustalO_outputs/241101DoubleCuratedFastas/DoubleCuratedAlignment.fa"
#results_file = "new-soks.unalignedgaps.Rdata"

#FILL THESE IN:

results_file = "new-soks-unalignedgaps.csv"
referenceprotein = "3197_tr|A0A2R6W7K7|A0A2R6W7K7_MARPO"
property = "windowCharge" #"windowCharge" OR "WindowHydropathy"
kingdomshort = "plants"        #OR "fungi" OR "animals"
kingdomlong = "Viridiplantae"  #OR "Fungi" OR "Metazoa"


#CALCULATE GLOBAL PROPS

# this function takes a gappy sequence and a data.frame as input, and matches
# the rows of the data.frame to the position of the corresponding amino acids in
# the sequence. The name of the column containing the amino acid position is
# given by the id variable

align_properties = function (gapped_sequence, properties_data, id)
{
  # Initialize a counter for the properties data index
  properties_index <- 1
  
  gapped_df = data.frame(Gapped.sequence = strsplit(as.character(gapped_sequence), "")[[1]], ..ResidueNumber = NA)
  
  for (i in 1:nrow(gapped_df)) {
    if (gapped_df$Gapped.sequence[i] != "-") { 
      gapped_df$..ResidueNumber[i] = properties_index
      properties_index = properties_index + 1
    }
  }
  
  gapped_df <- gapped_df %>% 
    rename_at('..ResidueNumber', ~id) %>% 
    full_join(properties_data, by = id)
  
  return(gapped_df)
  
}


# Read the input FASTA file
sequences <- readAAStringSet(fasta_file_path)

# Initialize results data frame
results <- data.frame(Gapped_sequence = character(),
                      Sequence_Name = character(),
                      Position = integer(),
                      CenterResidue = character(),
                      Window = character(),
                      windowCharge = numeric(),
                      WindowHydropathy = numeric(),
                      stringsAsFactors = FALSE)
i = 1


# Create a Vector with Columns
columns = c("name","length") 

#Create a Empty DataFrame with 0 rows and n columns
length_df = data.frame(matrix(nrow = 0, ncol = length(columns))) 

# Assign column names
colnames(length_df) = columns


# Calculate local charges and hydropathy for each sequence
for (seq_name in names(sequences)) {
  print (paste0("Processing sequence # ",i, " with name ", seq_name ))
  
  # get the sequence, remove gaps, replace X with glycine
  gapped_sequence <- sequences[[seq_name]]
  sequence = gsub("-", "", gapped_sequence)
  sequence = gsub("X", "G", sequence)
  
  lengthv <- c(seq_name, nchar(sequence))
  length_df <- rbind(length_df, lengthv)
  
  properties_df <- chargeCalculationLocal(sequence, window = window_size)
  hydropathy_df = scaledHydropathyLocal (sequence, window = window_size, plotResults = FALSE)
  hydropathy_df$WindowHydropathy = (hydropathy_df$WindowHydropathy * 9) - 4.5
  
  # Add the sequence hydropathy to the charge data frame
  properties_df$WindowHydropathy = hydropathy_df$WindowHydropathy
  
  # Match the properties to the alignment
  properties_df = align_properties(gapped_sequence, properties_df, "Position")
  
  # Add the sequence name and Alignment position to the charge data frame
  properties_df$Sequence_Name <- seq_name
  properties_df$AlignmentPosition = 1:nrow (properties_df)
  
  # Append to results data frame
  results <- rbind(results, properties_df)
  i=i+1
}

# View the results
head(results)

#SUMMARISE BY GROUP, CALCULATE STANDARD DEVIATION

DF <- results


DFpivot <- DF %>% pivot_wider(id_cols=AlignmentPosition, names_from=Sequence_Name, values_from=property)
DFpivot <- as.data.frame(DFpivot)
rownames(DFpivot) <- DFpivot$AlignmentPosition
DFpivot$AlignmentPosition <- NULL
#DFpivot[is.na(DFpivot)] <- 0
DFpivot <- t(DFpivot)
#write.csv(DFpivot, "chargeDF.csv", row.names = TRUE)
chemicalDF <- DFpivot
phyloDF <- TreeDF
#phyloDF$comb <- paste(phyloDF$Kingdom,phyloDF$label)
  
newDF <- merge(phyloDF, chemicalDF, by.x = "label", by.y = 0)
onlyKDDF <- newDF[-c(1:18, 20)]
groupedSDDF <- onlyKDDF %>% mutate_all(~replace(., .=='NA', NA)) %>% 
  dplyr::group_by(Kingdom) %>% 
  summarise(across(everything(), ~ sd(., na.rm=T)))
groupedSDDF <- as.data.frame(groupedSDDF)
rownames(groupedSDDF) <- groupedSDDF$Kingdom
groupedSDDF$Kingdom <- NULL
  
df_transposed <- as.data.frame(t(groupedSDDF))
df_transposed$rowid <- as.numeric(row.names(df_transposed))


#dev.off()

structureSD <- df_transposed[,kingdomlong]
reference <- DFpivot[referenceprotein,]

both <- rbind(structureSD, reference)
both <- as.data.frame(both)

transposeboth <- as.data.frame(t(both))
write.csv(transposeboth, paste("strucSD_", property, "_", kingdomshort, ".csv", sep = ""))

# Now proceed with the python file for PyMOL coloring.
