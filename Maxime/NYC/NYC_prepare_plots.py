import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
plot_dataframe = gpd.plotting.plot_dataframe
import numpy as np
import itertools
plt.interactive(False)

EPSG=2263 # projection
RAWDATA_DIR = "NYC_data/raw_data"
nycdistrs=gpd.read_file("NYC_data/nysd_16c/nysd.shp").to_crs(epsg=EPSG)

def schdistr_labels(ax, color="black", alpha=0.5, scaleup=1.0, **kwargs):
    xmin,xmax = plt.xlim()
    ymin,ymax = plt.ylim()
    # Get dimensions of y-axis in pixels
    y1, y2 = ax.get_window_extent().get_points()[:, 1]
    # Get unit scale
    yscale = (y2-y1)/(ymax-ymin)
    ax.set_autoscale_on(False)
    for irow in range(nycdistrs.shape[0]):
        row = nycdistrs.loc[irow, :]
        geom = row.geometry
        schdistr = row.SchoolDist
        square_side = np.sqrt(geom.area)
        center = geom.representative_point()
        if not (xmin < center.x < xmax):
            continue
        if not (ymin < center.y < ymax):
            continue
        plt.annotate(schdistr,
            xy=(center.x, center.y),
            xycoords="data",
            xytext=(0,0),
            textcoords="offset points",
            horizontalalignment="center",
            verticalalignment="center",
            fontweight="black",
            color=color,
            alpha=alpha,
            fontsize=max(square_side * yscale / 3 * scaleup, 10),
            **kwargs
            )
    ax.set_autoscale_on(True)

def background_schdistrs(ax, **kwargs):
    plot_dataframe(nycdistrs, ax=ax, **kwargs)
    #  plt.axes().set_aspect('equal', 'datalim')
    #  plt.axis("off")
    return None
