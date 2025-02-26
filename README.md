What follows here is an overview of the scripts and what they were used for.

-	tree.R is used to create the large phylogenetic tree.

-	properties.R calculates the values for the hydropathy and charge analysis of the actins.

-	heatmap_chemical.R uses these values and visualizes them.

-	heatmap_grouped_by_kingdom.R makes the heatmaps per kingdom that are shown in this thesis.

-	structure_coloring.R calculates the standard deviation values used to color the actin structures. tree.R has to be run first.

-	After running structure_coloring.R, PyMOL_coloring.py can be run in PyMOL to color the structures.

-	thallus_area.py determines the thallus area from the TIFF files taken by the scanner during phenotyping.

-	ablation_analysis.py is used for analyzing the ablation data. It takes the TIFF data from the ablation and masks them for a given threshold, then calculates the average change over all pixel between every frame in the TIFF file. It gives these values as a .csv file.
