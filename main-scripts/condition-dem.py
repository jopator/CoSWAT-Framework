

import sys, os
from cjfx import *
import argparse
import numpy as np
import pandas as pd
import geopandas as gpd
from pathlib import Path
import warnings
import datavariables as variables
import rasterio
import rioxarray as rxr
from rasterio.features import rasterize
from shapely.geometry import Polygon, MultiPolygon
from shapely.geometry import shape as shapely_shape
from rasterio.features import shapes
from shapely.ops import unary_union
from collections import defaultdict, deque

import datavariables as variables

warnings.filterwarnings("ignore") 



def cells_along_segment(x0, y0, x1, y1, transform):
    # convert real-world coords to fractional row/col space
    col0 = (x0 - transform.c) / transform.a
    row0 = (y0 - transform.f) / transform.e
    col1 = (x1 - transform.c) / transform.a
    row1 = (y1 - transform.f) / transform.e

    dx = col1 - col0
    dy = row1 - row0

    if dx == 0 and dy == 0:
        return [(int(np.floor(row0)), int(np.floor(col0)))]

    cur_col, cur_row = int(np.floor(col0)), int(np.floor(row0))
    end_col, end_row = int(np.floor(col1)), int(np.floor(row1))

    step_col = 1 if dx > 0 else -1
    step_row = 1 if dy > 0 else -1

    t_delta_col = abs(1 / dx) if dx != 0 else np.inf
    t_delta_row = abs(1 / dy) if dy != 0 else np.inf

    next_col_boundary = cur_col + (1 if step_col > 0 else 0)
    next_row_boundary = cur_row + (1 if step_row > 0 else 0)

    t_max_col = (next_col_boundary - col0) / dx if dx != 0 else np.inf
    t_max_row = (next_row_boundary - row0) / dy if dy != 0 else np.inf

    cells = [(cur_row, cur_col)]

    while (cur_col, cur_row) != (end_col, end_row):
        if t_max_col < t_max_row:
            t_max_col += t_delta_col
            cur_col += step_col
        else:
            t_max_row += t_delta_row
            cur_row += step_row
        cells.append((cur_row, cur_col))

    return cells


# change working directory
me = os.path.realpath(__file__)
os.chdir(os.path.dirname(me))

if __name__ == '__main__':

    # create argument parser
    parser = argparse.ArgumentParser(description="a script to initialise a SWAT+ project for a given region")

    parser.add_argument("r", help="the name of the region to initialise the model for. If not specified, all regions will be processed.", nargs='*', default=[])
    parser.add_argument("--v", help="the version of the model setup to use. If not specified, the datavariables value will be used.", nargs='?', default=None)

    args = parser.parse_args()
    version = variables.version

    if args.v:
        version = args.v

    if len(args.r) > 0: 
        regions = args.r
        if len(regions) == 1 and regions[0] == 'all': regions = list_folders('../data-preparation/resources/regions/')
    else: regions = list_folders('../data-preparation/resources/regions/')

    details = {
        'auth': variables.final_proj_auth,
        'code': variables.final_proj_code,
        'coswat_data_version':variables.coswat_data_version,
        'region_source':variables.region_source
    }

    for region in regions:
        print(f"\t> Conditioning DEM for {region} - {version}                \n")

        #dirs
        modelData  = Path('../model-data')  / f'{variables.coswat_data_version}' / f'{region}'
        modelSetup = Path('../model-setup') / f'CoSWATv{version}' / f'{region}'

        # Prepare burn - in (from hydroRivers)
        burnInRiv_FN = modelData / 'shapes/hydroRivers_burn_in-ESRI-54003.gpkg'

        if variables.region_source == "regions_v2":
            baseRegion = region.split("-")[0]+"-"+region.split("-")[1]
            details['baseRegion'] = baseRegion
        details['region'] = region
        region_FN  = variables.cutline.format(**details)
        region_gdf = gpd.read_file(region_FN)

        if not os.path.isfile(burnInRiv_FN):
            print("\n \t\t > Preparing reference network (HydroRIVERS Flow Order 5)")
            hydroRivers_FN = "../resources/hydro-rivers/hydroRIVERS_ordclas5_v10.gpkg" #"../resources/hydro-rivers/hydroRIVERS_ordclas5_v10.gpkg"
            region_gdf = region_gdf.to_crs("epsg:4326")
            hydroRivers_global_gdf = gpd.read_file(hydroRivers_FN)
            region_rivers_gdf = hydroRivers_global_gdf.clip(region_gdf)
            region_rivers_gdf = region_rivers_gdf.to_crs("ESRI:54003")

            burn_shp = modelData / 'shapes/hydroRivers_burn_in-ESRI-54003.shp'
            region_rivers_gdf.to_file(burnInRiv_FN)
            region_rivers_gdf.to_file(burn_shp)

        else:
            region_rivers_gdf = gpd.read_file(burnInRiv_FN)

        # Generate copy of original if not existing
        dem_FN_copy    = modelData / 'raster/dem-aster-{auth}-{code}_original.tif'.format(**details)
        dem_FN         = modelData / 'raster/dem-aster-{auth}-{code}.tif'.format(**details)
        if os.path.isfile(dem_FN_copy):
            dem_FN = dem_FN_copy
        else:
            with rasterio.open(dem_FN) as src:
                dem = src.read(1)
                profile = src.profile
            with rasterio.open(dem_FN_copy, 'w', **profile) as dst:
                dst.write(dem, 1)


        # recondition DEM (Burn)
        print("\t\t > Reconditioning DEM")
        with rasterio.open(dem_FN) as src:
            dem_bounds = src.bounds
            shape = src.shape
            crs = src.crs
            transform = src.transform
            nodata = src.nodata
            dem_arr = src.read(1)

        dem_da = rxr.open_rasterio(dem_FN).squeeze()
        dem_da = dem_da.where(dem_da != dem_da.rio.nodata)

        # ==Rasterize river network==
        valid = (dem_arr!=nodata).astype('uint8')
        shapes_gen = shapes(valid, mask=valid.astype(bool), transform=transform)
        valid_polygons = [shapely_shape(geom) for geom, val in shapes_gen]
        dem_valid_extent = unary_union(valid_polygons)


        if region_rivers_gdf.crs != crs:
            region_rivers_gdf = region_rivers_gdf.to_crs(crs)

        rivers_filtered = region_rivers_gdf[region_rivers_gdf.geometry.within(dem_valid_extent)].copy()

        print(f"\t\t\t Kept {len(rivers_filtered)} of {len(region_rivers_gdf)} reaches")

        rivers = rivers_filtered.copy()
        # Keep UPLAND_SKM from hydroRivers
        rivers_sorted = rivers.sort_values("UPLAND_SKM")
        shapes = zip(rivers_sorted.geometry, rivers_sorted["UPLAND_SKM"])
        upland_raster = rasterize(shapes=shapes,out_shape=shape,transform=transform,fill=np.nan,all_touched=True,dtype="float64",)



        base_depth   = variables.base_depth
        upland_scale = variables.upland_scale
        max_depth    = variables.max_depth
        burn_depth = base_depth + upland_scale * np.log10(upland_raster)
        burn_depth = np.clip(burn_depth, base_depth, max_depth)

        # ==Get elevation to each river cell==
        print(f"\t\t\t Getting elevation profile...")
        reach_profiles = {}

        for _, row in rivers.iterrows():
            geom = row.geometry
            if geom is None or geom.is_empty:
                continue

            lines = geom.geoms if geom.geom_type == 'MultiLineString' else [geom]

            ordered_cells = []
            for line in lines:
                coords = list(line.coords)
                for (x0, y0), (x1, y1) in zip(coords[:-1], coords[1:]):
                    ordered_cells.extend(cells_along_segment(x0, y0, x1, y1, transform))

            de_dupli = [cell for i, cell in enumerate(ordered_cells) if i == 0 or cell != ordered_cells[i-1]]
            de_dupli = [(r, c) for r, c in de_dupli if 0 <= r < dem_arr.shape[0] and 0 <= c < dem_arr.shape[1]]

            if not de_dupli:
                continue

            elevations = [dem_arr[r, c] for r, c in de_dupli]
            reach_profiles[row['HYRIV_ID']] = {'cells': de_dupli, 'elev': elevations}  
        
        # ==Monotonic minimum==
        print(f"\t\t\t Monotonic minimum...")
        # Upstream contributor per reach
        id_nextDown = dict(zip(rivers['HYRIV_ID'],rivers['NEXT_DOWN']))
        upstream_of = defaultdict(list)

        inlet_nr = {hyriv_id:0 for hyriv_id in reach_profiles}

        for hyriv_id, next_down in id_nextDown.items():
            if hyriv_id in reach_profiles and next_down in reach_profiles:
                upstream_of[next_down].append(hyriv_id)
                inlet_nr[next_down] +=1      
                                    
        # Topological order
        queue = deque([r for r in reach_profiles if inlet_nr[r]==0])
        topo_order = [] # To process all reaches after their upstream dependencies are finished
        inlet_nr_copy = inlet_nr.copy()

        while queue:
            current = queue.popleft()
            topo_order.append(current)
            next_down = id_nextDown.get(current)
            if next_down in inlet_nr_copy:
                inlet_nr_copy[next_down] -= 1
                if inlet_nr_copy[next_down]== 0:
                    queue.append(next_down)
        # Cummulative minimum
        end_value = {}
        for hyriv_id in topo_order:
            elev = reach_profiles[hyriv_id]['elev']
            upstream_ids = upstream_of.get(hyriv_id,[])
            start_ref = min(end_value[u] for u in upstream_ids) if upstream_ids else np.inf

            corrected = []

            running_min = start_ref
            for e in elev:
                if e == nodata:
                    corrected.append(running_min)
                    continue

                running_min = min(running_min,e)
                corrected.append(running_min)
            reach_profiles[hyriv_id]['corrected_elev'] = corrected
            end_value[hyriv_id] = running_min

        all_inf_ids = [hyriv_id for hyriv_id, data in reach_profiles.items() if all(np.isinf(v) for v in data['corrected_elev'])]             # Get rid of nodata cells
        lengths = [len(reach_profiles[h]['cells']) for h in all_inf_ids]
        corrected_raster = np.full(dem_arr.shape, np.nan, dtype='float64')
        for hyriv_id, data in reach_profiles.items():
            for (r, c), elev in zip(data['cells'], data['corrected_elev']):
                if np.isfinite(elev):
                    corrected_raster[r, c] = elev

        # ==Burn in on corrected profile== 
        print(f"\t\t\t Creating corrected profile...")
        for hyriv_id, data in reach_profiles.items():
            cells = data['cells']
            corrected = data['corrected_elev']

            burned = []
            for (r, c), elev in zip(cells, corrected):
                depth = burn_depth[r, c]
                if np.isfinite(elev) and np.isfinite(depth):
                    burned.append(elev - depth)
                else:
                    burned.append(elev)

            reach_profiles[hyriv_id]['burned_elev'] = burned

        # ==Apply to DEM== 
        print(f"\t\t\t DEM correction...")
        dem_burned = dem_arr.copy().astype('float64')

        for hyriv_id, data in reach_profiles.items():
            for (r, c), elev in zip(data['cells'], data['burned_elev']):
                if np.isfinite(elev):
                    dem_burned[r, c] = elev

        dem_burned_da = dem_da.copy(data=dem_burned)
        dem_burned_da = dem_burned_da.where(dem_burned_da != dem_da.rio.nodata)

        import gc

        del dem_arr, dem_burned, corrected_raster
        gc.collect()

        # ==Burn lakes==
        print(f"\t\t\t Lake burn...")                        
        if variables.burn_lakes:
            channel_mask = np.zeros(shape, dtype=bool)
            for hyriv_id, data in reach_profiles.items():
                for (r, c) in data['cells']:
                    channel_mask[r, c] = True
            res_rast_to_burn_gdf = gpd.read_file(modelData / 'shapes/lakes-grand-{auth}-{code}-demAligned.shp'.format(**details))
            interior = res_rast_to_burn_gdf.geometry.buffer(variables.data_resolution*(-1),cap_style="square", join_style="mitre")
            border   = res_rast_to_burn_gdf.geometry.difference(interior)
            border_gdf = res_rast_to_burn_gdf.copy()
            border_gdf['geometry'] = border.geometry
            interior_gdf = res_rast_to_burn_gdf.copy()
            interior_gdf['geometry'] = interior.geometry
            interior_gdf_valid = interior_gdf[~interior_gdf.geometry.is_empty].copy()
            interior_mask = rasterize(interior_gdf_valid.geometry, out_shape=shape, transform=transform)
            border_mask   = rasterize(border_gdf.geometry, out_shape=shape, transform=transform)
            dem_burned_values = dem_burned_da.values
            combined_interior_mask = (interior_mask == 1) & (valid == 1) & (~channel_mask)
            combined_border_mask = (border_mask == 1) & (valid == 1) & (~channel_mask)
            dem_burned_values[combined_interior_mask] -= 14
            dem_burned_values[combined_border_mask] -= 10
            dem_burned_da = dem_burned_da.copy(data=dem_burned_values)

        dem_burned_da = dem_burned_da.fillna(nodata) 
        dem_burned_da = dem_burned_da.rio.write_nodata(nodata)

        out_dem_FN = modelData / 'raster/dem-aster-{auth}-{code}.tif'.format(**details)
        dem_burned_da.rio.to_raster(out_dem_FN) # On model - data
        out_dem_FN = modelSetup / 'Watershed/Rasters/DEM/dem-aster-{auth}-{code}.tif'.format(**details)
        dem_burned_da.rio.to_raster(out_dem_FN) # On model - setup

# =============================================================================================================================================================
