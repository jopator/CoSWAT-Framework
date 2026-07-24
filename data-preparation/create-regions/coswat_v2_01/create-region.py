"""
File: create-region.py
Author: Jose P. Teran
Github: jopator
Date: 2026-07-20
Description: Generate region basic files used for model creation based on reference geopackage

This doesn't need to be run on the container. It just needs python >=3.11
"""

import os
import sys
import pandas as pd
import geopandas as gpd
import numpy as np
from pathlib import Path
from datetime import datetime
from region_dict import coswat_reg_dict,continent_dict


os.chdir("/data/brussel/vo/000/bvo00033/coswat_model/coswat_framework")

refereceDataDir = Path("../reference_data/COSWAT/regions_v2")
dstRegionDir    = Path("resources/regions_v2")
continentsDir   = Path("resources/continents")

for mainRegion, item in coswat_reg_dict.items():
    os.makedirs(dstRegionDir/mainRegion,exist_ok=True)
    subRegions_gdf = gpd.GeoDataFrame(columns=['name','class','upstream','area_km2','geometry'])
    for subRegion in item:
        os.makedirs(dstRegionDir/mainRegion/subRegion,exist_ok=True)
        subRegion_FN = refereceDataDir/f"{subRegion}.gpkg"

        subRegion_gdf = gpd.read_file(subRegion_FN)
        subRegion_gdf = subRegion_gdf.to_crs("ESRI:54003") #Convert to ESRI-54003

        # Clip to land mass
        continentName = continent_dict[subRegion.split("-")[0]]
        landMass_gdf  = gpd.read_file(continentsDir/f"{continentName}-ESRI-54003.gpkg")

        out_gdf = subRegion_gdf.clip(landMass_gdf.to_crs(landMass_gdf.crs))

        out_gdf['name']     = subRegion
        out_gdf['upstream'] = str(item[subRegion])
        out_gdf['area_km2'] = (out_gdf.geometry.area/1000000).astype(int)
        # Define class
        out_gdf['class']  = list(subRegion.split("-")[2])[0]

        out_gdf = out_gdf[['name','class','upstream','area_km2','geometry']].copy()
        out_gdf['geometry'] = out_gdf['geometry'].copy().buffer(5*1000)
        out_gdf.to_file(dstRegionDir/mainRegion/subRegion/"region_mask-ESRI-54003.gpkg")
        subRegions_gdf = pd.concat([subRegions_gdf,out_gdf])

    subRegions_gdf = subRegions_gdf.reset_index(drop=True)
    subRegions_gdf.to_file(dstRegionDir/mainRegion/"region_structure.gpkg")
        

