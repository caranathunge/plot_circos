import pandas as pd
from pycirclize import Circos
from pycirclize.utils import calc_group_spaces, ColorCycler
import numpy as np

def plot_circos(file_path, 
                palette: 'str' = 'Set3',
                track_1_data: 'str'= 'Cancer_type',
                track_2_data: 'str' = 'Tumors',
                bar_pal_start: 'int' = 2,
                hmap_pal_start: 'int' = -1,
                plot_title: 'str' = '',
                track_2_label: 'str' = '',
                file_name: 'str' = 'plot.pdf'):
    """
    Function to create circos plots from data containing cancer driver gene data.
    First track shows cancer types, second track is a bar plot showing frequency 
    data (i.e. tumor frequency), and the third track shows the driver gene status
    of the fusion genes in heatmap colored based on absence/presence

    Args:
    - file_path (str): Path to the file containing data to be plotted.
    - palette (str): Set color palette for the two outer tracks (Default: `Set3`)
    - track_2_data (str): Column name for track 2 plot data (Default: `Tumors`)
    - bar_pal_start (int): Palette counter start for the bar plot and the outer ring (Default: 2)
    - hmap_pal_start (int): Palette counter start for the inner heatmap (Default: -1)
    - plot_title (str): Title for the plot (Default = '')
    - track_2_label (str): Label for track 2 (Default = '')
    - track_1_data (str): Column name for track 1 plot data (Default: 'Cancer_type')
    - file_name (str): Name for the plot file in pdf format (Default: 'plot')
    
    Returns:
    - None
    """
 
    # Create a data frame from the txt file
    df = pd.read_csv(file_path, sep='\t', header=0)
    
    # Build sectors (which is a dictionary) using unique items in the cancer_type column and their occurences
    # sectors = {name: 30 for name in unique_cancer_types_list}
    sectors = df[track_1_data].value_counts().to_dict()
    
    # Assign spaces among sectors
    circos = Circos(sectors, space=7)
    
    # Set a color map. We are using 'Set3' as default
    ColorCycler.set_cmap(palette)
    
    ## CANCER TYPE TRACK
    
    # To reset the colors for each track so that each sector (cancer_type) has the same color, we are using 2 counters. Otherwise,
    # color picking will continue on after the first track. We are starting the counter at color 3 (python series starts at 0) of
    # 'Set3' palette
    col_1 = bar_pal_start
    
    for sector in circos.sectors:
        # Now we are on color 4
        col_1 = col_1 + 1
    
        # Add the cancer_type track between radius 94 and 100 for each sector
        track1 = sector.add_track(r_lim=(94, 100))
    
        # Fill the track with our color and set the edge color to 'lightgrey'
        track1.axis(fc=ColorCycler(col_1), ec='lightgrey')
    
        # Add the cancer_type text to the first track
        track1.text(sector.name, size=10, fontweight='book')
    
    # Outside of the for loop, we print 'A' for the first track, so that it only appears once, right after the last sector
    track1.text(
        text='A',
        orientation='horizontal',
        x=track1.end + 0.8,
        size=7,
        ignore_range_error=True,
        r=96,
        color='#636363',
        fontweight='medium',
    )
    
    ## TUMOR FREQUENCY TRACK - BAR PLOT

    # Find the maximum value from Tumors column and round it up to nearest ten. We will use this as the max val in the plots
    max_value = np.max(df[track_2_data])
    rounded_max_value = np.ceil(max_value / 10) * 10


    # Starting counter number 2 for track two at the same color
    col_2 = bar_pal_start
    
    for sector in circos.sectors:
        # Now we are on color 4
        col_2 = col_2 + 1
    
        # Add the bar plot track between radius 48 and 66 for each sector
        track2 = sector.add_track((66, 48), r_pad_ratio=0.1)
    
        # Filter DataFrame based on the desired cancer type to plot in each sector
        filtered_df1 = df[df[track_1_data] == sector.name]
        
        # Get x-axis and y-axis data for the filtered DataFrame.
        # x-axis is sector size.
        # y-axis is tumor frequency
        x = np.arange(sector.size)
        y = filtered_df1[track_2_data].values
        
        # We assign the rounded_max_value to vmax in the sector object
        vmax2 = rounded_max_value
    
        # Add the fusion gene names from 'query' column as x-axis tick labels on the outside
        track2.xticks(
            x+0.4,
            labels=filtered_df1['query'].values,
            label_orientation='vertical',
            label_size=5,
            show_bottom_line=True,
            line_kws=dict(ec=ColorCycler(col_2), lw=2),
            text_kws=dict(fontweight='book'),
        )

        # Set ticks and labels on the y-axis using the calculated vmax
        track2.yticks(
            [0, vmax2 / 2, vmax2],
            labels=[str(int(value)) for value in np.arange(0, vmax2 + 1, vmax2 / 2)],
            vmin=0,
            vmax=vmax2,
            label_size=4,
            tick_length=0,
            side='left',
        )
    
        # Show grid
        track2.grid(y_grid_num=5)
    
        # Set axis edge colors and line weight
        track2.axis(ec='None', lw=1)
    
        # Plot the bars to represent the Tumor frequency
        track2.bar(x, 
                   y, 
                   vmin=0,
                   vmax=vmax2,
                   align='edge', 
                   width=0.75, 
                   color=ColorCycler(col_2))
    
    # Outside of the for loop, we print 'B' for the second track, so that it only appears once, right after the last sector
    track2.text(
        text='B',
        orientation='horizontal',
        x=track2.end + 0.8,
        size=7,
        ignore_range_error=True,
        r=68,
        color='#636363',
        fontweight='medium',
    )
    
    
    ## DRIVER STATUS TRACK - HEATMAP
    
    # Set vmin2 and vmax2 as 0 and 1 to represent presence and absence of driver status
    vmin3, vmax3 = 0, 1
    
    # We will be using a heatmap to plot this track, so we can only use a sequential color palette. We will try to be as close
    # as possible to the Set3 color palette colors
    pal_list = ['Reds', 'Greens','YlOrBr', 'Blues', 'Reds', 'Purples', 'Greys']
    
    # Initiate a counter
    col_pal = hmap_pal_start
    
    for sector in circos.sectors:
        # STarting on color palette one from pal_list. Python series starts at 0.
        col_pal = col_pal + 1
    
        # Plot heatmap between radius 44 and 46
        track3 = sector.add_track((44, 46))
    
        # Filter data based on cancer type
        filtered_df2 = df[df[track_1_data] == sector.name]
    
        # Get the drive status column from the filtered df
        data = filtered_df2['driver'].values
    
        # Plot the heatmap
        track3.heatmap(
            data,
            vmin=vmin3,
            vmax=vmax3,
            cmap=pal_list[col_pal],
            rect_kws=dict(ec='lightgrey', lw=1),
        )
    
    # Outside of the for loop, we print 'C' for the third track, so that it only appears once, right after the last sector
    track3.text(
        text='C',
        orientation='horizontal',
        x=track3.end + 0.8,
        size=7,
        ignore_range_error=True,
        r=45,
        color='#636363',
        fontweight='medium',
    )
    
    # We put the track labels/legends outside of the circle
    circos.text(
        'A: Cancer type',
        size=8,
        r=106,
        deg=180,
        ha='left',
        va='bottom',
        color='#636363',
        fontweight='medium',
    )
    circos.text(
        track_2_label,
        size=8,
        r=108,
        deg=180,
        ha='left',
        va='center',
        color='#636363',
        fontweight='medium',
    )
    circos.text(
        'C: Driver gene status',
        size=8,
        r=110,
        deg=180,
        ha='left',
        va='top',
        color='#636363',
        fontweight='medium',
    )
    circos.text(
        plot_title,
        size=14,
        r=110,
        deg=0,
        ha='center',
        fontweight='medium',
    )
    # Show plot
    fig = circos.plotfig()
    
    # Save plot as a pdf
    fig.savefig(file_name, format='pdf', dpi=300)