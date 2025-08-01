'''
This script deals with issues related to lakes and reservoir delination on the river network.

It 'Resolves' these water bodies by:
    - Simplifying the geometry by applying algorithms and buffering the lake polygon until
        it is certain that the connections on the SWAT+ files wont break.

    - Creates snapping points for the QSWAT+ delineation and topology algorithms so that the same logic
        used here is replicated on the delination process.

    - Creates a reference flagged channels shape file in the Watershed folder of the basin
        To verify if they are adequately defined (Flags: LakeIn, LakeOut, LakeWithin)

> This script is run when setting up the model on the run-qswatplus.py file

It can be also executed alone to test different buffer thresholds and simplification algorithms


Author  : Jose Pablo teran
Date    : May 2025
Contact : jose.pablo.teran.orsini@vub.be
GitHub  : github.com/celray - github.com/jopator
'''

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, LineString, MultiLineString, Polygon, MultiPolygon, mapping, shape, box
import numpy as np
import datavariables as variables
import sys
import os
from simplification.cutil import simplify_coords_vw
from cjfx import exists, goto_dir


# Functions
def get_endpoints(geom):
    if geom.is_empty:
        return None, None

    if isinstance(geom, LineString):
        coords = list(geom.coords)
        if not coords:
            return None, None
        return Point(coords[-1]), Point(coords[0])

    elif isinstance(geom, MultiLineString):
        parts = list(geom.geoms)
        if not parts:
            return None, None
        start_coords = list(parts[-1].coords)
        end_coords = list(parts[0].coords)
        if not start_coords or not end_coords:
            return None, None
        return Point(start_coords[-1]), Point(end_coords[0])

    return None, None

def classify_clipped(row, lake_geometry, tol=1e-6):
    start, end  = get_endpoints(row.geometry)
    if start is None or end is None:
        return None
    start_in    = lake_geometry.distance(start) < tol
    end_in      = lake_geometry.distance(end) < tol

    if start_in and end_in:
        return 'LakeWithin'
    elif not start_in and end_in:
        return 'LakeIn'
    elif start_in and not end_in:
        return 'LakeOut'
    return None

def adjust_flag_in(idx, row,clipped):
    if row['lake_flag'] != 'LakeIn':
        return row['lake_flag']
    
    for j, other in clipped.iterrows():
        if j == idx:
            continue
        if row['start'].equals(other['start']) or row['start'].equals(other['end']):
            if other['lake_flag'] in {'LakeWithin', 'LakeOut'}:
                return 'LakeWithin'
    return row['lake_flag']

def adjust_flag_out(idx, row,clipped):
    if row['lake_flag'] != 'LakeOut':
        return row['lake_flag']
    
    for j, other in clipped.iterrows():
        if j == idx:
            continue
        if row['end'].equals(other['start']) or row['end'].equals(other['end']):
            if other['lake_flag'] in {'LakeWithin', 'LakeIn'}:
                return 'LakeWithin'
    return row['lake_flag']

def adjust_flag_nan(idx, row,clipped):
    if pd.isna(row['lake_flag']):
        for j, other in clipped.iterrows():
            if j == idx:
                continue
            # Case if end of NaN river connects with start of LakeIn river
            if row['end'].equals(other['start']):
                if other['lake_flag'] in {'LakeIn'}:
                    return 'Additional_OK'
            
            # Case if start of NaN rivers connect with end of LakeOut river
            if row['start'].equals(other['end']):
                if other['lake_flag'] in {'LakeOut'}:
                    return 'Additional_OK'
    else:
        return row['lake_flag']

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
    
def remove_overlaps_by_area(gdf, area_col, buffer_dist=0):
    gdf = gdf.copy()
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

    return gdf[gdf["to_remove"] == "no"].drop(columns="to_remove").reset_index(drop=True)
    

if __name__ == '__main__':
    
    # change working directory
    goto_dir(__file__)
    
    # get model setup version
    version = sys.argv[1]
    
    if not exists(f"../model-setup/CoSWATv{version}"):
        print(f'\t! the version CoSWATv{version} does not exist')
        sys.exit(1)
        
    # get regions 
    if len(sys.argv) >= 3: 
        region = sys.argv[2]
    else: 
        print("Did not select a region")
        quit()
    
    # Framework details
    details = {
        'auth'      : variables.final_proj_auth,
        'code'      : variables.final_proj_code,
        'version'   : version,
        'region'    : region
    }
    
    # Paths and files
    projDir     = f'../model-setup/CoSWATv{version}/{region}'
    dataDir     = f'../model-data/{region}'
    
    chan_shp    = os.path.join(projDir,'Watershed/Shapes/dem-aster-{auth}-{code}channel.shp'.format(**details))   
    res_shp     = os.path.join(dataDir,'shapes/lakes-grand-{auth}-{code}.shp'.format(**details))                                                       
    
    # Read GeoDataFrames
    chan_gdf    = gpd.read_file(chan_shp)
    res_gdf     = gpd.read_file(res_shp)
    
    # DEM Resolution (to check for short channels)
    dem_resolution = variables.data_resolution 
    


    '''
    Simplification of lake geometries
    '''

    dem_diagonal = (dem_resolution**2+dem_resolution**2)**0.5

    # Original source/outlets and direction vectors
    orig_geom_info = {}
    warning_printed = False
    
    for idx, row in chan_gdf.iterrows():
        geom = row.geometry
        coords = list(geom.coords) if isinstance(geom, LineString) else list(geom.geoms[0].coords)


        if len(coords) < 2:
            continue

        # Source direction = last segment (moving away from source)
        src_vec = np.array([coords[-2][0] - coords[-1][0], coords[-2][1] - coords[-1][1]])
        src_norm = np.linalg.norm(src_vec)

        
        if src_norm != 0:
            src_vec = src_vec / src_norm
        else:
            if not warning_printed:
                print("Warning: Some channels have Null source or outlet vectors. This normally shouldn't be an issue.")
                warning_printed = True

            src_vec = np.array([0, 0])

        # Outlet direction = first segment (moving away from outlet)
        out_vec = np.array([coords[0][0] - coords[1][0], coords[0][1] - coords[1][1]])
        out_norm = np.linalg.norm(out_vec)

        if out_norm != 0:
            out_vec = out_vec / out_norm
        else:
            out_vec = np.array([0, 0])

        orig_geom_info[row['LINKNO']] = {
            'source': Point(coords[-1]),
            'outlet': Point(coords[0]),
            'geometry' : geom,
            'src_vec': src_vec,
            'out_vec': out_vec
        }
    
    # Empty lists of geometries to create DataFrames
    flagged_chans = []
    res_gdf.loc[:,"resolved"] = False
    res_gdf.loc[:,'bbox'] = None
    res_gdf = res_gdf[res_gdf['Hylak_id'].notna()].copy()
    
    
    network_union = chan_gdf.geometry.union_all()

    # Filter res_gdf by checking if each geometry intersects with the network
    res_gdf = res_gdf[res_gdf.geometry.intersects(network_union)]
    
    # Simplify lake geometry if necessary
    if variables.simplify_geometry:
        if variables.simplify_method    == 'DP':
            print("\t \t > Simplifying geometries with DP algorithm")
            res_gdf["geometry"] = res_gdf["geometry"].simplify(tolerance=dem_resolution, preserve_topology=True)
        
        elif variables.simplify_method  == "VW":
            print("\t \t > Simplifying geometries with VW algorithm")
            res_gdf["geometry"] = res_gdf["geometry"].apply(lambda g: simplify_geometry_vw(g, dem_resolution))

        elif variables.simplify_method  == 'ConV':
            res_gdf['geometry']  = res_gdf.geometry.convex_hull
            print("\t \t > Simplifying geometries with VW algorithm")

    
    # Remove lakes that overlay or are too close (based on buffer)
    res_gdf_clean = remove_overlaps_by_area(res_gdf, area_col='parea', buffer_dist=0)
    res_gdf_final = remove_overlaps_by_area(res_gdf_clean, area_col='parea', buffer_dist=variables.lake_buffer_thres)

    res_gdf = res_gdf_final.copy()


    for index,row in res_gdf.iterrows():                                                                        # Iterating through all reservoirs/lakes
        
        threshold           = variables.lake_buffer_thres
        step                = variables.lake_buffer_step
        attempts            = 0
        lake_id             = row['LakeId']
        lake_geom           = row['geometry']                                                                   # Get the geometry of the current lake / reservoir
        lake_boundary       = lake_geom.boundary
        
        while attempts * step <= threshold:
            buffered = lake_geom.buffer(attempts * step, resolution=8)
            buffered = remove_holes(buffered)                                                                   # Remove Holes just in case

            # Generate bounding box and intersect channels
            minx, miny, maxx, maxy = buffered.bounds
            bbox                = box(minx, miny, maxx, maxy)
            bbox_gdf            = gpd.GeoDataFrame({'LakeId': [lake_id]}, geometry=[bbox], crs=res_gdf.crs)
            
            chan_bbox_intersect = chan_gdf[chan_gdf.intersects(bbox_gdf.geometry.union_all())].copy()           # Get channels that intersect the bounding box
            chan_bbox_intersect = chan_bbox_intersect.copy()                                                    # For the typical pandas shenannigans >:(
            
            chan_bbox_intersect = gpd.sjoin(
                chan_bbox_intersect, 
                bbox_gdf[['LakeId','geometry']], 
                how='left', 
                predicate='intersects').drop(columns='index_right')
            
            chan_bbox_intersect['connected'] = False 
            
            # Check intersection        
            tmp                 = chan_bbox_intersect[chan_bbox_intersect['LakeId'] == lake_id].copy()
            tmp['intersects']   = tmp.geometry.intersects(buffered)                                             # Define if channels that are inside the Bounding Box intersect the geometry
            tmp['within']       = tmp.geometry.within(buffered)
            
            intersecting        = tmp[tmp['intersects'] | tmp['within']]
            intersect_ids       = set(intersecting.index)                                                       # Get index of those that intersect
            
            tmp['endpoints']    = tmp.geometry.map(lambda g: list(get_endpoints(g)))                            # Assign tuple of start and end point for each channel
            connected_ids       = set(intersect_ids)
            new_pending         = set()
            
            for j in tmp.index.difference(connected_ids):
                endpoints = tmp.at[j, 'endpoints']
                for pid in intersect_ids:
                    p_endpoints = tmp.at[pid, 'endpoints']
                    if any(e1.equals(e2) for e1 in endpoints for e2 in p_endpoints):
                        new_pending.add(j)
                        break
                    
            connected_ids.update(new_pending)
            chan_bbox_intersect.loc[chan_bbox_intersect.index.isin(connected_ids), 'connected'] = True
            
            curr_chans = chan_bbox_intersect[(chan_bbox_intersect['LakeId'] == lake_id) & ((chan_bbox_intersect['connected'] == True))]
            
            if curr_chans.empty:
                print(f"\t\t Lake {lake_id} cannot be resolved into the current channel structure within buffer threshold (no channels found) - will be skiped.")
                break  
            
            clipped = gpd.overlay(  curr_chans, 
                                    gpd.GeoDataFrame(geometry=[buffered], crs=res_gdf.crs), 
                                    how='identity')
            
            clipped = clipped.explode(index_parts=False).reset_index(drop=True)
            clipped = clipped[~clipped.geometry.is_empty]
            
            clipped['lake_flag'] = clipped.apply(lambda row: classify_clipped(row, lake_geometry=buffered), axis=1)
            clipped['start'], clipped['end'] = zip(*clipped.geometry.map(get_endpoints))
            
            clipped['lake_flag'] = [adjust_flag_in(i, r,clipped=clipped)  for i, r in clipped.iterrows()]
            clipped['lake_flag'] = [adjust_flag_out(i, r,clipped=clipped) for i, r in clipped.iterrows()]
            clipped['lake_flag'] = [adjust_flag_nan(i, r,clipped=clipped) for i, r in clipped.iterrows()]
            
            
            # Check for short channels
            clipped['length_n']  = clipped.geometry.length
            clipped.loc[(clipped['lake_flag'].isin(['LakeIn', 'LakeOut'])) & (clipped['length_n'] < dem_resolution), 'lake_flag'] = np.nan
            
            if clipped['lake_flag'].isna().sum() == 0:
                break
            
            attempts += 1
        
        # If there are NaN channels
        if 'clipped' not in locals() or clipped.empty or clipped['lake_flag'].isna().any():
            print(f"\t\t Lake {lake_id} cannot be resolved into the current channel structure within buffer threshold.")
            continue
        
        else:
            print(f"\t\t Lake {lake_id} was resolved after {attempts} attempts.")
            res_gdf.loc[index,'resolved'] = True
        
        res_gdf.loc[index,'geometry'] = buffered        
        flagged_chans.append(clipped)
        
    combined_gdf = gpd.GeoDataFrame(pd.concat(flagged_chans, ignore_index=True), crs=clipped.crs)
    res_gdf = res_gdf[res_gdf['resolved']].copy()                                                                                                   # Filter only resolved lakes
    
    
    # Save resolved lakes    
    chan_fn    = os.path.join(projDir,'Watershed/Shapes/lakes-flagged-channels-{auth}-{code}.shp'.format(**details))   
    res_fn     = os.path.join(projDir,'Watershed/Shapes/lakes-grand-{auth}-{code}.shp'.format(**details))
    
    # Format to export
    combined_gdf['Length']      = combined_gdf['Length'].map(lambda x: f"{x:.2f}")
    combined_gdf['DSContArea']  = combined_gdf['DSContArea'].map(lambda x: f"{x:.1f}")
    combined_gdf['strmDrop']    = combined_gdf['strmDrop'].map(lambda x: f"{x:.2f}")
    combined_gdf['StraightL']   = combined_gdf['StraightL'].map(lambda x: f"{x:.1f}")
    combined_gdf['USContArea']  = combined_gdf['USContArea'].map(lambda x: f"{x:.1f}")
    combined_gdf['DOUTEND']     = combined_gdf['DOUTEND'].map(lambda x: f"{x:.1f}")
    combined_gdf['DOUTSTART']   = combined_gdf['DOUTSTART'].map(lambda x: f"{x:.1f}")
    combined_gdf['DOUTMID']     = combined_gdf['DOUTMID'].map(lambda x: f"{x:.1f}")
    
    try:
        res_gdf['calcVol']          = res_gdf['calcVol'].map(lambda x: f"{x:.1f}")
    
    except:
        print('calcVol not Formated')

    res_gdf.to_file(res_fn)
    combined_gdf.to_file(chan_fn)

    '''
    Creation of snapping points

    '''
    snap_points = []
    tolerance = 1e-6

    for idx, row in combined_gdf.iterrows():
        flag    = row['lake_flag']
        lakeId  = row['LakeId']
        linkno  = row['LINKNO']

        if linkno not in orig_geom_info:
            continue

        # Get original geometry info
        info = orig_geom_info[linkno]
        full_coords = list(info['geometry'].coords)  

        start_point, end_point = get_endpoints(row.geometry)

        if flag == "LakeOut":
            res       = 1
            candidate = start_point
            reference = info['source']
            vec       = info['src_vec']
            allowed_vertices = full_coords          # all vertices upstream from source
            ref_index = len(full_coords) - 1        # source is last

        elif flag == "LakeIn":
            res       = 0
            candidate = end_point
            reference = info['outlet']
            vec       = -info['out_vec']
            allowed_vertices = full_coords          # all vertices downstream from outlet
            ref_index = 0                           # outlet is first

        else:
            continue
        
        

        if candidate.distance(reference) < 0.5 * dem_diagonal:                      # If the 'candidate' snapping point is too close to a source, we desplace it one pixel
                                                                                    # It is always half the diagonal because all sources, outlets and joints are on pixel centroids
            
            moved_point = Point(candidate.x + vec[0] * dem_resolution,
                                candidate.y + vec[1] * dem_resolution)

            if not info['geometry'].buffer(tolerance).contains(moved_point):

                # Snap to next vertex
                if flag == "LakeOut":
                    candidates = full_coords[ref_index-1:] if ref_index > 0 else [full_coords[-1]]                              # Only vertices upstream of source (towards the end of coords)
                else:
                    candidates = full_coords[:ref_index+2] if ref_index < len(full_coords)-1 else [full_coords[0]]              # Only vertices downstream of outlet (towards the start of coords)

                closest_vertex = min(candidates, key=lambda c: moved_point.distance(Point(c)))
                point = Point(closest_vertex)
            else:
                point = moved_point
        else:
            point = candidate

        snap_points.append({'geometry': point, 'RES': res, 'LakeId': lakeId})


    # Most important output: Snapping points to be read on QSWAT+
    snap_gdf = gpd.GeoDataFrame(snap_points,crs=combined_gdf.crs)
    snap_gdf.to_file(variables.resolved_snaps.format(**details))


