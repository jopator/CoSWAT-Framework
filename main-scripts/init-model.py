#!/bin/python3

import sys, os
from cjfx import *
import argparse
import numpy as np
import pandas as pd
import geopandas as gpd

import rasterio
from rasterio import features
from rasterio.features import rasterize
from shapely.geometry import Polygon, MultiPolygon, shape
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
    
def remove_overlaps_by_area(orig_gdf, area_col, buffer_dist=0,print_out = False):
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
                        if print_out:
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
        'coswat_data_version':variables.coswat_data_version,
        'region_source':variables.region_source
    }

    for region in regions:
        report(f"\t> initializing {region}.qgs                \n")

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

        data_dir    = f'../model-data/{variables.coswat_data_version}/{proj_name}'

        # Align lakes to DEM
        if variables.include_lakes:
            res_shp     = os.path.join(data_dir,'shapes/lakes-grand-{auth}-{code}.shp'.format(**details))
            dem_path    = os.path.join(data_dir,'raster/dem-aster-{auth}-{code}.tif'.format(**details))

            out_simplified = os.path.join(data_dir,'shapes/lakes-grand-{auth}-{code}-demAligned.shp'.format(**details))
            out_simplified_snap = os.path.join(data_dir,'shapes/lakes-grand-{auth}-{code}-demAligned-forSnap.shp'.format(**details))

            if os.path.isfile(res_shp):
                print("\t\t > Aligning lakes to DEM")
                res_gdf     = gpd.read_file(res_shp)
                dem_resolution  = variables.data_resolution  
                dem_diagonal    = (dem_resolution**2+dem_resolution**2)**0.5
                res_gdf = res_gdf[res_gdf['Hylak_id'].notna()].copy()
                
                #Simplify lake geometry if necessary this will now be done inside the loop if the geometry cannot be solved
                if variables.simplify_geometry:
                    if variables.simplify_method    == 'DP':
                        res_gdf["geometry"] = res_gdf["geometry"].simplify(tolerance=dem_resolution, preserve_topology=True)                
                    elif variables.simplify_method  == "VW":
                        res_gdf["geometry"] = res_gdf["geometry"].apply(lambda g: simplify_geometry_vw(g, dem_resolution))
                    elif variables.simplify_method  == 'ConV':
                        res_gdf['geometry']  = res_gdf.geometry.convex_hull
                    elif variables.simplify_method  == 'ConC':
                        res_gdf['geometry']  = res_gdf.geometry.concave_hull(ratio=0.20, allow_holes=False)

                # Remove lakes that overlay (take largest)
                res_gdf_clean = remove_overlaps_by_area(res_gdf, area_col='parea', buffer_dist=0)
                # res_gdf_final = remove_overlaps_by_area(res_gdf_clean, area_col='parea', buffer_dist=variables.lake_buffer_step*1)
                res_gdf = res_gdf_clean.copy()
                

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

                res_rast_gdf['geometry']    = res_rast_gdf['geometry'].apply(remove_holes)
                res_rast_to_burn_gdf        = res_rast_gdf.copy()
                res_rast_to_burn_gdf['geometry'] = res_rast_to_burn_gdf.buffer(0,join_style=2,mitre_limit=1e6)
                res_rast_to_burn_gdf.to_file(out_simplified)
                lakesExsist = True    
            else:
                print("\t\t No lake shapefile for region, they wont be included...")
                lakesExsist = False    

        # data source paths
        dem_fn          = f"{data_dir}/raster/dem-aster-{variables.final_proj_auth}-{variables.final_proj_code}.tif"
        landuse_fn      = f"{data_dir}/raster/landuse-esa-{variables.esa_landuse_year}-{variables.final_proj_auth}-{variables.final_proj_code}.tif"
        soils_fn        = f"{data_dir}/raster/soils-fao-{variables.final_proj_auth}-{variables.final_proj_code}.tif"
        lakes_fn        = f"{data_dir}/shapes/lakes-grand-{variables.final_proj_auth}-{variables.final_proj_code}.shp"
        burn_shape_fn   = f"{data_dir}/shapes/hydroRivers_burn_in-ESRI-54003.shp"  #f"{data_dir}/shapes/burn-shape-{variables.final_proj_auth}-{variables.final_proj_code}.shp"

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
        
        
        with zipfile.ZipFile("../resources/shapes.dat", 'r') as zip_ref:
            zip_ref.extractall(dir_Shapes)
        
        shapes_files = list_files(f'{dir_Shapes}')
        for shapes_file in shapes_files:
            if "[dem]" in shapes_file:
                copy_file(shapes_file, shapes_file.replace('[dem]', f'{file_name(dem_fn, extension=False)}'), delete_source=True)

        tmp_burn = gpd.GeoDataFrame(geometry=[], crs=f"{variables.final_proj_auth}:{variables.final_proj_code}")
        tmp_burn.to_file(f"{dir_Shapes}/{file_name(burn_shape_fn)}")

        # gpd.read_file(burn_shape_fn).to_file(f"{dir_Shapes}/{file_name(burn_shape_fn)}")

        if not lakesExsist:
                empty_gdf = gpd.GeoDataFrame(geometry=[], crs=f"{variables.final_proj_auth}:{variables.final_proj_code}")
                empty_gdf.to_file(f"{dir_Shapes}/{file_name(lakes_fn)}")
        
        else:
            gpd.read_file(lakes_fn).to_file(f"{dir_Shapes}/{file_name(lakes_fn)}")

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

        os.remove(f"{dir_Shapes}/{file_name(burn_shape_fn)}")

print()
