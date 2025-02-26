
# Load required libraries
# install.packages("dplyr")
library(dplyr)
library(Biostrings)
library(idpr)
# BiocManager::install("idpr")
# https://bioconductor.org/packages/devel/bioc/vignettes/idpr/inst/doc/chargeHydropathy-vignette.html

setwd("C:/Users/martk/Documents/wur/masterthesis/bioinformatica/bioinfR/Output/Tree")

window_size = 9
fasta_file_path = "C:/Users/martk/Documents/wur/masterthesis/bioinformatica/ClustalO_outputs/241101DoubleCuratedFastas/DoubleCuratedAlignment.fa"
#results_file = "new-soks.unalignedgaps.Rdata"
results_file = "new-soks-unalignedgaps.csv"



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

# Save the data
save(results, file = results_file)
write.csv(results, results_file)
