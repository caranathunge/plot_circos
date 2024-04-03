import pandas as pd
from pycirclize import Circos
from pycirclize.utils import calc_group_spaces, ColorCycler
import numpy as np


def plot_circos_hla(
    file_path,
    palette: 'str' = 'Set3',
    track_1_data: 'str' = '',
    track_2_data: 'str' = '',
    track_3_data: 'str' = '',
    track_3_labels: 'str' = '',
    pal_start: 'int' = 0,
    plot_title: 'str' = '',
    file_name: 'str' = 'plot.pdf',
):
    """
    Function to create circos plots from data containing cancer driver gene data.

    Args:
    - file_path (str): Path to the file containing data to be plotted.
    - palette (str): Set color palette for the two outer tracks
    - track_1_data (str): Column name for track 1 plot data
    - track_2_data (str): Column name for track 2 plot data
    - track_3_data (str): Column name for track 3 plot data
    - pal_start (int): Palette counter start for the tracks (Default: 0)
    - plot_title (str): Title for the plot (Default = '')
    - track_3_labels (str): Labels for track 2 (Default = '')
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
    circos = Circos(sectors, space=7, start=4, end=360)

    # Set a color map. We are using 'Set3' as default
    ColorCycler.set_cmap(palette)

    ## CANCER TYPE TRACK

    # To reset the colors for each track so that each sector (cancer_type) has the same color, we are using 2 counters. Otherwise,
    # color picking will continue on after the first track. We are starting the counter at color 3 (python series starts at 0) of
    # 'Set3' palette
    col_1 = pal_start

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
        x=track1.end + 0.65,
        size=7,
        ignore_range_error=True,
        r=96,
        color='#636363',
        fontweight='medium',
    )

    ## PEPTIDE LENGTH TRACK - RECTANGLES
    col_3 = pal_start
    for sector in circos.sectors:
        track1 = sector.add_track(r_lim=(91, 84))
        filtered_df2 = df[df[track_1_data] == sector.name]
        length_list = list(filtered_df2[track_2_data].values)
        col_3 = col_3 + 1
        for i in range(int(track1.size)):
            start, end = i, i + 1
            track1.rect(
                start=start, end=end, fc=ColorCycler(col_3), ec='white', lw=1
            )
            track1.text(
                text=length_list[i],
                x=(end + start) / 2,
                adjust_rotation=True,
                size=6,
                color='black',
                orientation='vertical',
            )
    track1.text(
        text='B',
        orientation='horizontal',
        x=track1.end + 0.65,
        size=7,
        ignore_range_error=True,
        r=86,
        color='#636363',
        fontweight='medium',
    )

    ## HLA AFFINITY TRACK - BAR PLOT - NEGATIVE SCALE

    min_value = np.min(df[track_3_data])
    rounded_min_value = np.floor(min_value / 10) * 10
    vmin = rounded_min_value
    vmax = 0

    col_2 = pal_start
    for sector in circos.sectors:
        track = sector.add_track((66, 48), r_pad_ratio=0.1)
        col_2 = col_2 + 1
        # Filter DataFrame based on the desired cancer type
        filtered_df = df[df[track_1_data] == sector.name]
        # Get x-axis and y-axis data for the filtered DataFrame
        x = np.arange(sector.size)
        y = filtered_df[
            track_3_data
        ].values  # Get values from the 'tumor' column as y-axis

        # xticks = filtered_df['query'].values  # Get values from the 'query' column as x-axis

        track.xticks(
            x + 0.4,
            labels=filtered_df[track_3_labels].values,
            label_orientation='vertical',
            label_size=5,
            show_bottom_line=True,
            line_kws=dict(ec=ColorCycler(col_2), lw=2),
            text_kws=dict(fontweight='book'),
        )

        track.yticks(
            [vmax, vmin / 2, vmin],
            labels=[
                str(int(value)) for value in np.arange(0, vmin - 1, vmin / 2)
            ],
            vmin=vmin,
            vmax=vmax,
            label_size=4,
            tick_length=0,
            side='left',
        )
        track.grid(y_grid_num=5)
        track.axis(ec='None', lw=1)

        track.bar(
            x,
            y,
            align='edge',
            width=0.75,
            color=ColorCycler(col_2),
            vmin=vmin,
            vmax=vmax,
        )

    track.text(
        text='C',
        orientation='horizontal',
        x=track.end + 0.65,
        size=7,
        ignore_range_error=True,
        r=72,
        color='#636363',
        fontweight='medium',
    )
    track.text(
        text='D',
        orientation='horizontal',
        x=track.end + 0.65,
        size=7,
        ignore_range_error=True,
        r=58,
        color='#636363',
        fontweight='medium',
    )
    circos.text(
        'A: HLA',
        size=8,
        r=105,
        deg=180,
        ha='left',
        # va='bottom',
        color='#636363',
        fontweight='medium',
    )
    circos.text(
        'B: Length',
        size=8,
        r=109,
        deg=180,
        ha='left',
        # va='center',
        color='#636363',
        fontweight='medium',
    )
    circos.text(
        'C: Peptide',
        size=8,
        r=113,
        deg=180,
        ha='left',
        # va='top',
        color='#636363',
        fontweight='medium',
    )
    circos.text(
        'D: HLA arena binding energy',
        size=8,
        r=117,
        deg=180,
        ha='left',
        # va='bottom',
        color='#636363',
        fontweight='medium',
    )
    # Show plot
    fig = circos.plotfig()

    # Save plot as a pdf
    fig.savefig(file_name, format='pdf', dpi=300)
