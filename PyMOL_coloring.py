# Run the structure coloring R file before this file.

import numpy as np
import pandas as pd

def coloring(property, kingdom):                # Give the property (charge/hydropathy) and the kingdom (plants/animals/fungi)
    df = pd.read_csv(f"C:/Users/martk/Documents/wur/masterthesis/bioinformatica/bioinfR/Output/Structures/strucSD_{property}_{kingdom}.csv")
    df = df[["structureSD", "reference"]]       # Change the DF to only include the SD (or other values
                                                # to be mapped on the structure) and the reference.
    selected = 5                                # Window = 9, so first 4 amino acids are not colored.
    Color = np.linspace([255, 255, 255], [255, 0, 0], num=256)       # Sets the colorscale
    maxval = np.nanmax(df["structureSD"])       # For normalising the colorscale.
    for i in range(len(Color)):
        color_rgb = Color[i].astype(int)        # This loop sets the colors in Pymol.
        color_name = f'color{i}'
        colorline = f'set_color {color_name}, [{color_rgb[0]}, {color_rgb[1]}, {color_rgb[2]}]'
        cmd.do (colorline)
    for index,row in df.iterrows():
        # Check for every value in column structureSD if there is also a value at the same row in column reference.
        if pd.notna(row['reference']):
            value_for_color = row['structureSD']                # If there is, take the value from the structure SD
                                                                # and use it to find the right color.
            colorindex = round(value_for_color * 256 / maxval)  # This 0.4 can also be changed to maxval (now changed
                                                                # to lower value because of one very high value).
            if colorindex > 255:
                colorindex = 255                                # This is for when not using maxval above, sets the
                                                                # color to the 'highest'.
            line = f'color color{colorindex}, resi {selected}'
            cmd.do (line)                                       # Colors the amino acid.
            selected = selected + 1                             # And move to the next amino acid.
cmd.extend("coloring",coloring)
