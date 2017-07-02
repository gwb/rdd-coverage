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

def background_schdistrs(ax, **kwargs):
    plot_dataframe(nycdistrs, ax=ax, **kwargs)
    #plt.axes().set_aspect('equal', 'datalim')
    #plt.axis("off")
    return None
