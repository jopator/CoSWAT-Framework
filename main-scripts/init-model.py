#!/bin/python3

import sys, os
from cjfx import *
import argparse

import numpy as np
import numpy as np
import pandas as pd
import geopandas as gpd
import numpy as np
import sys
import os

import rasterio
from rasterio import features
from rasterio.features import rasterize
from rasterio.enums import MergeAlg
from rasterio.crs import CRS
from shapely.geometry import Polygon, MultiPolygon, shape
from shapely.validation import make_valid
from simplification.cutil import simplify_coords_vw

ignore_warnings()

def remove_holes(geometry):
    if geometry.geom_type == 'Polygon':
        return Polygon(geometry.exterior)
    elif geometry.geom_type == 'MultiPolygon':
        return MultiPolygon([Polygon(part.exterior) for part in geometry.geoms])
    else:
        return geometry
    

def simplify_geometry_vw(geom, tolerance):
    if geom.is_empty:
        return geom
    
    if geom.geom_type == "Polygon":
        exterior = simplify_coords_vw(list(geom.exterior.coords), tolerance)
        interiors = [simplify_coords_vw(list(ring.coords), tolerance) for ring in geom.interiors]
        return shape({"type": "Polygon", "coordinates": [exterior] + interiors})

    elif geom.geom_type == "MultiPolygon":
        polygons = []
        for poly in geom.geoms:
            exterior = simplify_coords_vw(list(poly.exterior.coords), tolerance)
            interiors = [simplify_coords_vw(list(ring.coords), tolerance) for ring in poly.interiors]
            polygons.append({"type": "Polygon", "coordinates": [exterior] + interiors})
        return shape({"type": "MultiPolygon", "coordinates": [p["coordinates"] for p in polygons]})

    else:
        # For lines or points (if any)
        return geom
    
def remove_overlaps_by_area(orig_gdf, area_col, buffer_dist=0):
    gdf = orig_gdf.copy()

    gdf["to_remove"] = "no"

    # Ensure parea is float
    gdf[area_col] = gdf[area_col].astype(float)

    if buffer_dist != 0:
        gdf["geometry"] = gdf.buffer(buffer_dist)

    changed = True
    
    while changed:
        changed = False

        for i, row in gdf[gdf["to_remove"] == "no"].iterrows():
            overlaps = gdf[
                (gdf["to_remove"] == "no") &
                (gdf.index != i) &
                (gdf.geometry.intersects(row.geometry))
            ]

            if not overlaps.empty:
                cluster = gdf.loc[[i] + overlaps.index.tolist()]
                largest_idx = cluster[area_col].idxmax()

                for idx in cluster.index:
                    if idx != largest_idx:
                        print(f"\t\t ** Reservoir {gdf.loc[idx, 'LakeId']} was marked for removal because it overlaps with reservoir {gdf.loc[largest_idx, 'LakeId']}")
                        gdf.at[idx, "to_remove"] = "yes"
                        changed = True


    check_df = gdf[['LakeId','to_remove']]

    return_gdf = orig_gdf.merge(check_df,how='left',on='LakeId')

    return return_gdf[return_gdf["to_remove"] == "no"].drop(columns="to_remove").reset_index(drop=True)

def rasterize_lakes(dem_fp, lakes_gdf, out_fp, attr="LakeId", all_touched=True, nodata=0):
    with rasterio.open(dem_fp) as dem:
        profile   = dem.profile.copy()
        transform = dem.transform
        crs       = dem.crs
        height    = dem.height
        width     = dem.width

    gdf = lakes_gdf.copy()
    
    gdf['geometry'] = gdf['geometry'].apply(remove_holes)
    
    if gdf.crs != crs:
        gdf = gdf.to_crs(crs)
    if attr not in gdf.columns:
        raise ValueError(f"Attribute '{attr}' not found.")

    gdf = gdf[gdf[attr].notna() & gdf.geometry.notna() & (~gdf.geometry.is_empty)]
    vals = gdf[attr].astype(np.uint32)

    arr = np.full((height, width), nodata, dtype=np.uint32)
    shapes = zip(gdf.geometry, vals.astype(int))  # ints are fine for uint32 target
    
    rasterize(
        shapes=shapes,
        out=arr,
        transform=transform,
        all_touched=all_touched,
    )

    profile.update(dtype="uint32", nodata=nodata, count=1, compress="lzw")
    with rasterio.open(out_fp, "w", **profile) as dst:
        dst.write(arr, 1)

def polygonize_lakes(raster_fp, dissolve=True, min_area=None, connectivity=8):
    with rasterio.open(raster_fp) as src:
        img = src.read(1)
        transform = src.transform
        crs = src.crs
        nodata = src.nodata if src.nodata is not None else 0

    if img.dtype not in (np.int16, np.int32, np.uint8, np.uint16, np.float32):
        img = img.astype(np.int32)

    mask = img != nodata
    recs = []
    for geom, val in features.shapes(img, mask=mask, transform=transform, connectivity=connectivity):
        v = int(val)
        if v == nodata:
            continue
        g = shape(geom)
        if min_area is not None and g.area < min_area:
            continue
        recs.append({"geometry": g, "LakeId": v})

    gdf = gpd.GeoDataFrame(recs, crs=crs)
    if dissolve and not gdf.empty:
        gdf = gdf.dissolve(by="LakeId", as_index=False)
    return gdf

def update_dem(dem_fp, lakes_fp, out_fp, force_crs=None):

    with rasterio.open(dem_fp) as dem_ds:
        dem = dem_ds.read(1)
        prof = dem_ds.profile.copy()
        transform = dem_ds.transform
        crs = dem_ds.crs
        dem_nodata = dem_ds.nodata
        dem_dtype = dem_ds.dtypes[0]

    gdf = gpd.read_file(lakes_fp)

    if gdf.crs != crs:
        gdf = gdf.to_crs(crs)

    gdf = gdf[(gdf["elev_masl"].notna()) & (~gdf.geometry.is_empty)]

    gdf = gdf[gdf.geometry.notna()]
    gdf['geometry'] = gdf.geometry.apply(make_valid)
    gdf = gdf[~gdf.geometry.is_empty]

    out_dtype = prof["dtype"]
    if np.issubdtype(np.dtype(out_dtype), np.integer):
        out_dtype = "float32"

    out_arr = dem.astype(out_dtype, copy=True)

    shapes = ((geom, float(val)) for geom, val in zip(gdf.geometry, gdf["elev_masl"]))

    rasterize(
        shapes=shapes,
        out=out_arr,
        transform=transform,
        all_touched=False,
        merge_alg=MergeAlg.replace,
    )

    # Keep driver, transform, etc.; fix nodata dtype; allow CRS override
    prof.update(
        driver="GTiff",
        dtype=out_dtype,
        count=1,
        compress="lzw",
        crs=(CRS.from_user_input(force_crs) if force_crs is not None else crs),
        transform=transform,
        nodata=(np.nan if out_dtype.startswith("float") and dem_nodata is None
                else (np.array(dem_nodata, dtype=np.float32).item() if out_dtype.startswith("float") and dem_nodata is not None
                      else dem_nodata))
    )

    with rasterio.open(out_fp, "w", **prof) as dst:
        dst.write(out_arr, 1)


# change working directory
me = os.path.realpath(__file__)
os.chdir(os.path.dirname(me))

import datavariables as variables
from resources.template_proj import template_string

if __name__ == '__main__':

    # create argument parser
    parser = argparse.ArgumentParser(description="a script to initialise a SWAT+ project for a given region")

    parser.add_argument("r", help="the name of the region to initialise the model for. If not specified, all regions will be processed.", nargs='*', default=[])
    parser.add_argument("--v", help="the version of the model setup to use. If not specified, the datavariables value will be used.", nargs='?', default=None)

    args = parser.parse_args()


    print('\n# initialising SWAT+ project')
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
    }

    for region in regions:
        report(f"\t> initializing {region}.qgs                ")

        continent = region.split('-')[0]
        zone = region.split('-')[1]

        dst_dir = create_path(f'../model-setup/CoSWATv{version}/')

        if exists(f'{dst_dir}/{region}/{region}.qgs'):
            # remove the directory path before continuing
            print()
            delete_path(f'{dst_dir}/{region}/')
            print("\t> creating a new project...")

        proj_name   = f"{region}"
        proj_dir    = f'{dst_dir}/{proj_name}'

        data_dir    = f'../model-data/{proj_name}'


# =============================================================================================================================================================
        # Preprocess DEM if new reservoir methods are used
        if variables.new_res_methods:
            res_shp     = os.path.join(data_dir,'shapes/lakes-grand-{auth}-{code}.shp'.format(**details))
            dem_path    = os.path.join(data_dir,'raster/dem-aster-{auth}-{code}.tif'.format(**details))
            # Read GeoDataFrames
            res_gdf     = gpd.read_file(res_shp)

            # DEM Resolution (to check for short channels)
            dem_resolution  = variables.data_resolution  
            dem_diagonal    = (dem_resolution**2+dem_resolution**2)**0.5


            # Reservoir preprocessing
            res_gdf = res_gdf[res_gdf['Hylak_id'].notna()].copy()

            #Simplify lake geometry if necessary this will now be done inside the loop if the geometry cannot be solved
            if variables.simplify_geometry:

                print('\t Pre-processing reservoirs: Burning into DEM to adjust network \n')

                if variables.simplify_method    == 'DP':
                    print("\t \t > Simplifying geometries with DP algorithm")
                    res_gdf["geometry"] = res_gdf["geometry"].simplify(tolerance=dem_resolution, preserve_topology=True)
                
                elif variables.simplify_method  == "VW":
                    print("\t \t > Simplifying geometries with VW algorithm")
                    res_gdf["geometry"] = res_gdf["geometry"].apply(lambda g: simplify_geometry_vw(g, dem_resolution))

                elif variables.simplify_method  == 'ConV':
                    res_gdf['geometry']  = res_gdf.geometry.convex_hull
                    print("\t \t > Simplifying geometries with VW algorithm")

                elif variables.simplify_method  == 'ConC':
                    res_gdf['geometry']  = res_gdf.geometry.concave_hull(ratio=0.20, allow_holes=False)
                    print("\t \t > Simplifying geometries with Concave Hull")


            # Remove lakes that overlay or are too close (based on buffer)
            res_gdf_clean = remove_overlaps_by_area(res_gdf, area_col='parea', buffer_dist=0)
            res_gdf_final = remove_overlaps_by_area(res_gdf_clean, area_col='parea', buffer_dist=variables.lake_buffer_step*3)
            res_gdf = res_gdf_final.copy()
            
            


            # Rasterized lakes
            res_gdf_for_raster = res_gdf.copy()
            res_gdf_for_raster['geometry'] = res_gdf['geometry'].buffer(dem_resolution) 

            out_rasterized = os.path.join(data_dir,'shapes/lakes-grand-{auth}-{code}-rasterized.tif'.format(**details))
            rasterize_lakes(dem_path,res_gdf_for_raster,out_rasterized)
            res_rast_gdf = polygonize_lakes(out_rasterized)


            res_data_df = res_gdf.copy()[['LakeId', 'Hylak_id', 'smax', 'pvol', 'evol', 'parea', 'earea', 'br1',
                                            'br2', 'RES', 'calcAreas', 'calcVol', 'elev_masl', 'dis_avg',
                                            'Lake_name', 'Country', 'Continent', 'Lake_type', 'Depth_avg',
                                            'Res_time', 'Grand_id', 'RES_NAME', 'DAM_NAME', 'RIVER', 'MAIN_BASIN',
                                            'SUB_BASIN', 'YEAR', 'REM_YEAR', 'DAM_HGT_M', 'DAM_LEN_M', 'DEPTH_M',
                                            'DOR_PC', 'MAIN_USE', 'LAKE_CTRL', 'TIMELINE', 'a', 'b', 'c', 'd',
                                            'maxDepth', 'meanDepth', 'MaxArea_km']]

            res_rast_gdf = res_rast_gdf.merge(res_data_df,on='LakeId',how='left')


            res_rast_gdf = gpd.GeoDataFrame(data=res_rast_gdf,geometry=res_rast_gdf.geometry)

            out_simplified = os.path.join(data_dir,'shapes/lakes-grand-{auth}-{code}-demAligned.shp'.format(**details))
            res_rast_gdf['geometry']    = res_rast_gdf['geometry'].apply(remove_holes)
            res_rast_to_burn_gdf        = res_rast_gdf.copy()
            res_rast_to_burn_gdf['geometry'] = res_rast_to_burn_gdf.buffer(
                0,
                join_style=2,
                mitre_limit=1e6)

            res_rast_to_burn_gdf.to_file(out_simplified)

            # Buffer inwards half dem so we can snap directly
            import shapely
            def keep_polygons(g):
                g = shapely.make_valid(g)
                parts = [p for p in shapely.get_parts(g) if isinstance(p, Polygon)]
                return shapely.union_all(parts) if parts else None


            res_rast_buff_gdf = res_rast_gdf.copy()
            res_rast_buff_gdf['geometry'] = res_rast_buff_gdf.buffer(
                -dem_resolution/2 - 0.1,
                join_style=2,
                mitre_limit=1e6)

            res_rast_buff_gdf['geometry'] = res_rast_buff_gdf['geometry'].apply(remove_holes)
            res_rast_buff_gdf['geometry'] = res_rast_buff_gdf['geometry'].apply(keep_polygons)

            res_rast_buff_gdf['geometry'] = res_rast_buff_gdf.buffer(0.1, join_style=2, mitre_limit=1e6)
            
            res_rast_buff_gdf = res_rast_buff_gdf[res_rast_buff_gdf.geometry.type.isin(["Polygon","MultiPolygon"])]
            

            out_simplified_buffered = os.path.join(data_dir,'shapes/lakes-grand-{auth}-{code}-demAligned-forSnap.shp'.format(**details))
            res_rast_buff_gdf.to_file(out_simplified_buffered)

            # Burn DEM
            dem_path_new    = os.path.join(data_dir,'raster/dem-aster-{auth}-{code}-lakeBurnt.tif'.format(**details))
            update_dem(dem_path,out_simplified,dem_path_new)
# =============================================================================================================================================================


        # data source paths
        if variables.new_res_methods:
            dem_fn          = f"{data_dir}/raster/dem-aster-{variables.final_proj_auth}-{variables.final_proj_code}-lakeBurnt.tif"
        else:
            dem_fn          = f"{data_dir}/raster/dem-aster-{variables.final_proj_auth}-{variables.final_proj_code}.tif"


        landuse_fn      = f"{data_dir}/raster/landuse-esa-{variables.esa_landuse_year}-{variables.final_proj_auth}-{variables.final_proj_code}.tif"
        soils_fn        = f"{data_dir}/raster/soils-fao-{variables.final_proj_auth}-{variables.final_proj_code}.tif"

        lakes_fn        = f"{data_dir}/shapes/lakes-grand-{variables.final_proj_auth}-{variables.final_proj_code}.shp"
        burn_shape_fn   = f"{data_dir}/shapes/burn-shape-{variables.final_proj_auth}-{variables.final_proj_code}.shp"

        # create project structure
        create_path(f"{proj_dir}/")
        dir_DEM         = create_path(f"{proj_dir}/Watershed/Rasters/DEM/")
        dir_Landscape   = create_path(f"{proj_dir}/Watershed/Rasters/Landscape/")
        dir_Landuse     = create_path(f"{proj_dir}/Watershed/Rasters/Landuse/")
        dir_Soil        = create_path(f"{proj_dir}/Watershed/Rasters/Soil/")

        dir_Shapes      = create_path(f"{proj_dir}/Watershed/Shapes/")

        copy_file(dem_fn, f"{dir_DEM}/{file_name(dem_fn)}")
        copy_file(landuse_fn, f"{dir_Landuse}/{file_name(landuse_fn)}")
        copy_file(soils_fn, f"{dir_Soil}/{file_name(soils_fn)}")
        
        
        with zipfile.ZipFile("../data-preparation/resources/shapes.dat", 'r') as zip_ref:
            zip_ref.extractall(dir_Shapes)
        
        shapes_files = list_files(f'{dir_Shapes}')
        for shapes_file in shapes_files:
            if "[dem]" in shapes_file:
                copy_file(shapes_file, shapes_file.replace('[dem]', f'{file_name(dem_fn, extension=False)}'), delete_source=True)

        geopandas.read_file(burn_shape_fn).to_file(f"{dir_Shapes}/{file_name(burn_shape_fn)}")
        geopandas.read_file(lakes_fn).to_file(f"{dir_Shapes}/{file_name(lakes_fn)}")

        # prepare qgs project
        project_string = template_string.format(
            project_name        = proj_name,
            authid              = '{auth}:{code}'.format(**details),

            rivs_1_id           = f'{rand_apha_num(8)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(12)}',
            channel_shape_id    = f'{rand_apha_num(8)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(12)}',
            dem_id              = f'{rand_apha_num(8)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(12)}',
            lsus_shape_id       = f'{rand_apha_num(8)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(12)}',
            hillshade_id        = f'{rand_apha_num(8)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(12)}',
            outlets_id          = f'{rand_apha_num(8)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(12)}',
            landuse_id          = f'{rand_apha_num(8)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(12)}',
            reservoir_shape_id  = f'{rand_apha_num(8)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(12)}',
            se_outlets_shape_id = f'{rand_apha_num(8)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(12)}',
            soils_id            = f'{rand_apha_num(8)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(12)}',
            burn_shape_id       = f'{rand_apha_num(8)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(12)}',
            stream_shape_id     = f'{rand_apha_num(8)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(12)}',
            subbasins_id        = f'{rand_apha_num(8)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(12)}',
            lakes_id            = f'{rand_apha_num(8)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(4)}_{rand_apha_num(12)}',
            
            thresholdCh         = variables.thresholdCh,
            thresholdSt         = variables.thresholdSt,

            dem_file_name       = file_name(dem_fn, extension=False),
            land_use_file_name  = file_name(landuse_fn, extension=False),
            soils_file_name     = file_name(soils_fn, extension=False),
            burn_file_name      = file_name(burn_shape_fn, extension=False),
            lakes_file_name     = file_name(lakes_fn, extension=False),
            
            dem_file_name_underscore_hyphens        = file_name(dem_fn, extension=False).replace('-', '_'),
            land_use_file_name_underscore_hyphens   = file_name(landuse_fn, extension=False).replace('-', '_'),
            soils_file_name_underscore_hyphens      = file_name(soils_fn, extension=False).replace('-', '_'),
            burn_file_name_underscore_hyphens       = file_name(burn_shape_fn, extension=False).replace('-', '_'),
            lakes_file_name_underscore_hyphens      = file_name(lakes_fn, extension=False).replace('-', '_'),

        )

        write_to(f'{proj_dir}/{proj_name}.qgs', project_string)
        print(f'\n\t> initialised {proj_name}.qgs\n')

print()
