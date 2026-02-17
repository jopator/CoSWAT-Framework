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
from math import hypot
import numpy as np
import datavariables as variables
import sys
import os
from simplification.cutil import simplify_coords_vw
from shapely.ops import linemerge
from cjfx import exists, goto_dir

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


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

def adjust_flag_out_wrongs(idx, row,clipped):
    if row['lake_flag'] == 'LakeOut':
        for j, other in clipped.iterrows():
            if j == idx:
                continue

            if row['end'].equals(other['end']):
                if other['lake_flag'] in {'Additional_OK','LakeOut'}:
                    return 'Additional_OK'
            
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
            
            # Case if end of NaN rivers connect with end of LakeOut river
            if row['end'].equals(other['end']):
                if other['lake_flag'] in {'LakeOut'}:
                    return 'Additional_OK'

            # Case if end of NaN rivers connect with start of lake Within
            if row['end'].equals(other['start']):
                if other['lake_flag'] in {'LakeWithin'}:
                    return 'LakeIn'
                
            # Case if start of NaN rivers connect with end of lake Within
            if row['start'].equals(other['end']):
                if other['lake_flag'] in {'LakeWithin'}:
                    return 'LakeOut'
                
    else:
        return row['lake_flag']

def adjust_aditional_flag(idx,row,clipped):
    if row['lake_flag'] != 'Additional_OK':
        return row['lake_flag']
    
    for j, other in clipped.iterrows():
        if j == idx:
            continue
        if row['end'].equals(other['start']):
            if other['lake_flag'] in {'LakeWithin'}:
                return 'LakeIn'
            
        if row['start'].equals(other['end']):
            if other['lake_flag'] in {'LakeWithin'}:
                return 'LakeOut'

    return row['lake_flag']


def remove_holes(geometry):
    if geometry.geom_type == 'Polygon':
        return Polygon(geometry.exterior)
    elif geometry.geom_type == 'MultiPolygon':
        return MultiPolygon([Polygon(part.exterior) for part in geometry.geoms])
    else:
        return geometry
    
def other_links_count(point, self_index, lake_subset):
    count = 0
    for j, other in lake_subset.iterrows():
        if j == self_index:
            continue
        other_start, other_end = get_endpoints(other.geometry)
        if other_start.distance(point) <= tolerance or other_end.distance(point) <= tolerance:
            count += 1
    return count

def other_links_count(point, self_index, lake_subset,tolerance = 1e-6):
    count = 0
    for j, other in lake_subset.iterrows():
        if j == self_index:
            continue
        other_start, other_end = get_endpoints(other.geometry)
        if other_start.distance(point) <= tolerance or other_end.distance(point) <= tolerance:
            count += 1
    return count

def lines_to_segments(gdf):
    recs = []
    for _, row in gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        parts = [geom] if geom.geom_type == "LineString" else list(geom.geoms) if geom.geom_type == "MultiLineString" else []
        base = row.drop(labels="geometry").to_dict()
        for part_id, ls in enumerate(parts):
            coords = list(ls.coords)
            for i in range(len(coords) - 1):
                seg = LineString([coords[i], coords[i+1]])  # keeps Z if present
                rec = {**base, "part_id": part_id, "segment_id": i, "geometry": seg}
                recs.append(rec)
    return gpd.GeoDataFrame(recs, crs=gdf.crs)



def remerge_segments(seg_gdf):
    groups = []
    for (linkno, lake_flag), group in seg_gdf.groupby(["LINKNO", "lake_flag"]):
        merged = linemerge(group.geometry.tolist())
        rec = {
            "LINKNO": linkno,
            "lake_flag": lake_flag,
            "geometry": merged
        }
        groups.append(rec)
    return gpd.GeoDataFrame(groups, crs=seg_gdf.crs)
    
# Union buffered pieces into their assigned lake
def _union(series):
    try:
        return union_all(series.tolist())
    except Exception:
        from shapely.ops import unary_union
        return unary_union(series.tolist())
        
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

    simplify_geometry = variables.simplify_geometry
    
    # Paths and files
    projDir     = f'../model-setup/CoSWATv{version}/{region}'
    dataDir     = f'../model-data/{region}'
    
    chan_shp    = os.path.join(projDir,'Watershed/Shapes/dem-aster-{auth}-{code}-lakeBurntchannel.shp'.format(**details))  
    res_shp     = os.path.join(dataDir,'shapes/lakes-grand-{auth}-{code}-demAligned-forSnap.shp'.format(**details))
    subs_shp    = os.path.join(projDir,'Watershed/Shapes/dem-aster-{auth}-{code}-lakeBurntsubbasins.shp'.format(**details))  

    # res_shp     = os.path.join(dataDir,'shapes/lakes-grand-{auth}-{code}.shp'.format(**details))
            
    
    # Read GeoDataFrames
    chan_gdf    = gpd.read_file(chan_shp)
    res_gdf     = gpd.read_file(res_shp)
    subs_gdf    = gpd.read_file(subs_shp)

    # DEM Resolution (to check for short channels)
    dem_resolution  = variables.data_resolution 
    dem_diagonal    = (dem_resolution**2+dem_resolution**2)**0.5


    # Empty lists of geometries to create DataFrames
    flagged_chans               = []
    unresolved_channels         = []
    connected_channels          = []
    res_gdf.loc[:,"resolved"]   = False
    res_gdf.loc[:,'bbox']       = None
    res_gdf.loc[:,'attempts']   = 0
    
    res_gdf = res_gdf[res_gdf['Hylak_id'].notna()].copy()
    
    network_union = chan_gdf.geometry.union_all()
    


    # Filter res_gdf by checking if each geometry intersects with the network
    res_gdf = res_gdf[res_gdf.geometry.intersects(network_union)]


    for index,row in res_gdf.iterrows():                                                                        # Iterating through all reservoirs/lakes
        lake_name           = row['Lake_name']
        threshold           = variables.lake_buffer_thres
        step                = variables.lake_buffer_step
        attempts            = 0
        lake_id             = row['LakeId']
        lake_geom           = row['geometry']                                                                   # Get the geometry of the current lake / reservoir
        lake_boundary       = lake_geom.boundary


        
        while attempts * step <= threshold:
            if attempts > 0:
                res_gdf.loc[index,'attempts'] = attempts

            buffered = lake_geom.buffer((attempts * step), join_style=2, mitre_limit=1e6)
            buffered = remove_holes(buffered)                                                                               # Remove Holes just in case

            buffered_gdf = gpd.GeoDataFrame({'LakeId': [lake_id]}, geometry=[buffered], crs=res_gdf.crs)

            # Generate bounding box and intersect channels and subbasins
            minx, miny, maxx, maxy = buffered.bounds

            bbox                = box(minx-threshold, miny-threshold, maxx+threshold, maxy+threshold)
            bbox_gdf            = gpd.GeoDataFrame({'LakeId': [lake_id]}, geometry=[bbox], crs=res_gdf.crs)
            
            chan_bbox_intersect = chan_gdf[chan_gdf.intersects(bbox_gdf.geometry.union_all())].copy()           # Get channels that intersect the bounding box
            chan_bbox_intersect = chan_bbox_intersect.copy()                                                    # For the typical pandas shenannigans >:(
            


            chan_bbox_intersect = gpd.sjoin(
                chan_bbox_intersect, 
                bbox_gdf[['LakeId','geometry']], 
                how='left', 
                predicate='intersects').drop(columns='index_right')

            
            chan_bbox_intersect['connected'] = False 
            

            subs_intersect = subs_gdf[subs_gdf.intersects(buffered_gdf.geometry.union_all())].copy()                # Get subbasins that intersect the bounding box


            noDataBasins = True

            while noDataBasins:
                buffered_res = buffered_gdf.copy()
                buffered_res['geometry'] = buffered_res['geometry'].buffer(variables.data_resolution/2, join_style=2, mitre_limit=1e6)

                expected_hru_gdf = gpd.overlay(subs_intersect, buffered_res, how='difference', keep_geom_type = True, make_valid = True)
                subs_overlay_gdf = gpd.overlay(subs_intersect, buffered_gdf, how='difference', keep_geom_type = True, make_valid = True)
                expected_hru_gdf = expected_hru_gdf.dissolve()
                
                inter = gpd.sjoin(subs_overlay_gdf, expected_hru_gdf[['geometry']], how="inner", predicate="intersects")
                touch = gpd.sjoin(subs_overlay_gdf, expected_hru_gdf[['geometry']], how="inner", predicate="touches")

                # keep subs with zero positive-area overlap
                hit_idx = inter.index.difference(touch.index)
                no_data_subs = subs_overlay_gdf.loc[~subs_overlay_gdf.index.isin(hit_idx)]


                if len(no_data_subs)>0:
                    # crossed = gpd.sjoin(no_data_subs, chan_gdf, how="inner", predicate="crosses")
                    # result  = no_data_subs.loc[no_data_subs.index.isin(crossed.index)]                                                          # Only those that cross channels matter
                    # no_data_subs = result.copy()
                
                    buf = no_data_subs.copy()
                    buf["geometry"] = buf.geometry.buffer(variables.data_resolution/2, join_style=2, mitre_limit=1e6)
                    buf = buf[~buf.geometry.is_empty & buf.geometry.notna()]


                    #Map each buffered sub to nearest lake
                    nearest = gpd.sjoin_nearest(
                        buf[["geometry"]],
                        buffered_gdf[["geometry"]],
                        how="left",
                        distance_col="dist").rename(columns={"index_right":"lake_idx"})



                    for lake_idx, grp in nearest.groupby("lake_idx"):
                        if lake_idx is None or len(grp)==0:
                            continue
                        to_add = _union(grp.geometry)
                        buffered_gdf.at[lake_idx, "geometry"] = buffered_gdf.geometry.at[lake_idx].union(to_add)
                        
                    
                
                else:
                    noDataBasins = False


            buffered = buffered_gdf.loc[0,'geometry']

            # Check intersection        
            tmp                 = chan_bbox_intersect[chan_bbox_intersect['LakeId'] == lake_id].copy()
            tmp['crosses']      = tmp.geometry.crosses(buffered)                                                # Define if channels that are inside the Bounding Box intersect the geometry
            tmp['within']       = tmp.geometry.within(buffered)
            
            intersecting        = tmp[tmp['crosses'] | tmp['within']]
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


            # inset lake boundary by 0.5 m to avoid boundary-coincident degeneracy
            eps = -0.001
            lake_inset = buffered.buffer(eps, join_style=2, mitre_limit=1e6)
            if lake_inset.is_empty:
                lake_inset = buffered

            clipped = gpd.overlay(
            curr_chans,
            gpd.GeoDataFrame(geometry=[lake_inset], crs=res_gdf.crs),
            how="identity",
            keep_geom_type=False,           # keep all types
            )
                      
            # We will also add the first level ds and us linknos in case there are bifurcations    
            current_linknos = clipped['LINKNO'].to_list()
            new_linknos     = []
            
            for idx, row in clipped.iterrows():
                dslinkno    = row['DSLINKNO']
                uslinkno1   = row['USLINKNO1']
                uslinkno2   = row['USLINKNO2']
                
                if not dslinkno in current_linknos and not dslinkno in new_linknos:
                    new_linknos.append(dslinkno)
                  
                elif not uslinkno1 in current_linknos and not uslinkno1 in new_linknos:
                    new_linknos.append(uslinkno1)
                    
                elif not uslinkno2 in current_linknos and not uslinkno2 in new_linknos:
                    new_linknos.append(uslinkno2)
                    
            
            extra_chans = chan_gdf[chan_gdf['LINKNO'].isin(new_linknos)]
            extra_chans = extra_chans.copy()
            extra_chans['connected']    = True
            extra_chans['LakeId']       = lake_id
            
            clipped_with_extra = pd.concat([clipped,extra_chans])
            clipped = clipped_with_extra.reset_index(drop=True)
            
            # Preprocess
            clipped = clipped[clipped.geometry.type.isin(["LineString", "MultiLineString"])]
            clipped = clipped.explode(index_parts=False).reset_index(drop=True)
            clipped = clipped[~clipped.geometry.is_empty]
            clipped['start'], clipped['end'] = zip(*clipped.geometry.map(get_endpoints))
            
            
            # First round of classification
            clipped['lake_flag'] = clipped.apply(lambda row: classify_clipped(row, lake_geometry=buffered), axis=1)            
            # First round of adjust flags
            clipped['lake_flag'] = [adjust_flag_in(i, r,clipped=clipped)  for i, r in clipped.iterrows()]
            clipped['lake_flag'] = [adjust_flag_out(i, r,clipped=clipped) for i, r in clipped.iterrows()]
            clipped['lake_flag'] = [adjust_flag_nan(i, r,clipped=clipped) for i, r in clipped.iterrows()]
            
            
            # Check for short channels
            clipped['length_n']  = clipped.geometry.length
            clipped.loc[(clipped['lake_flag'].isin(['LakeIn', 'LakeOut'])) & (clipped['length_n'] < dem_resolution), 'lake_flag'] = np.nan

            # Second pass for NaN
            clipped['lake_flag'] = [adjust_flag_nan(i, r,clipped=clipped) for i, r in clipped.iterrows()]

            # Third pass for NaN and fake outs and ins
            clipped['lake_flag'] = [adjust_flag_nan(i, r,clipped=clipped) for i, r in clipped.iterrows()]
            clipped['lake_flag'] = [adjust_flag_out_wrongs(i, r,clipped=clipped) for i, r in clipped.iterrows()]
            clipped['lake_flag'] = [adjust_aditional_flag(i, r,clipped=clipped) for i, r in clipped.iterrows()]


            # Check to avoid infinite loops
            tol = 1e-9
            to_nan = []

            for idx, row in clipped.iterrows():
                if row['lake_flag'] != 'Additional_OK':
                    continue

                connections_in, connections_out = [], []
                s0, s1 = row['start'], row['end']

                for j, other in clipped.iterrows():
                    if j == idx:
                        continue
                    o0, o1 = other['start'], other['end']

                    # end-node connections
                    if s1.equals(o0) or s1.equals(o1) or s1.distance(o0) <= tol or s1.distance(o1) <= tol:
                        connections_out.append(other['lake_flag'])

                    # start-node connections
                    if s0.equals(o0) or s0.equals(o1) or s0.distance(o0) <= tol or s0.distance(o1) <= tol:
                        connections_in.append(other['lake_flag'])

                if ('LakeOut' in connections_in) and ('LakeIn' in connections_out):
                    to_nan.append(idx)

            if to_nan:
                clipped.loc[to_nan, 'lake_flag'] = np.nan

            clipped = clipped[clipped['length_n']>0]
            
            # Correct bifurcation channels
            clipped['count_end'] = 1
            clipped['count_start'] = 1
            cols = ['LINKNO', 'DSLINKNO', 'USLINKNO1', 'USLINKNO2', 'DSNODEID', 'strmOrder',
                    'Length', 'Magnitude', 'DSContArea', 'strmDrop', 'Slope', 'StraightL',
                    'USContArea', 'WSNO', 'DOUTEND', 'DOUTSTART', 'DOUTMID', 'BasinNo',
                    'LakeId', 'connected', 'length_n']

            for idx, row in clipped.iterrows():
                if row['lake_flag'] == 'LakeOut':
                    start = row['start']
                    count = other_links_count(start,idx,clipped)
                    if count>1:
                        clipped.loc[idx,'count_start'] = count

            for idx, row in clipped.iterrows():
                if row['lake_flag'] == 'LakeIn':
                    start = row['end']
                    count = other_links_count(start,idx,clipped)
                    if count>1:
                        clipped.loc[idx,'count_end'] = count


            # Separate inlet channels that end in bifurcation
            inlets_bif_gdf = clipped[clipped['count_end']>1]

            # Separate outlets channels that end in bifurcation
            outlets_bif_gdf = clipped[clipped['count_start']>1]

            if len(inlets_bif_gdf) > 0:
                segments_gdf = lines_to_segments(inlets_bif_gdf)
                segments_gdf.loc[segments_gdf['segment_id'] == 0,'lake_flag'] = 'LakeWithin'
                remerged_segments = remerge_segments(segments_gdf)
                cols = ['LINKNO', 'DSLINKNO', 'USLINKNO1', 'USLINKNO2', 'DSNODEID', 'strmOrder',
                        'Length', 'Magnitude', 'DSContArea', 'strmDrop', 'Slope', 'StraightL',
                        'USContArea', 'WSNO', 'DOUTEND', 'DOUTSTART', 'DOUTMID', 'BasinNo',
                        'LakeId', 'connected', 'length_n']

                new_fixed = pd.merge(remerged_segments,inlets_bif_gdf[cols],how='left',on='LINKNO')
                new_fixed = gpd.GeoDataFrame(data=new_fixed,geometry=new_fixed['geometry'])

                new_fixed['count_end'] = 1
                new_fixed['count_start'] = 1
                fixed_clipped = pd.concat([clipped[clipped['count_end']==1].reset_index(drop=True),new_fixed.reset_index(drop=True)]).reset_index(drop=True)
                clipped = fixed_clipped.copy()


            
            if len(outlets_bif_gdf) > 0:
                segments_gdf    = lines_to_segments(outlets_bif_gdf)

                if len(segments_gdf) == 1:
                    # We get the one downstream
                    outlets_bif_gdf = pd.merge(outlets_bif_gdf[['LakeId','LINKNO','lake_flag','connected','geometry']],chan_gdf[chan_gdf.columns.drop('geometry')],how='left',on='LINKNO')
                    outlets_bif_gdf = gpd.GeoDataFrame(data = outlets_bif_gdf, geometry = outlets_bif_gdf.geometry)

                    linkno_problem  = outlets_bif_gdf['LINKNO'].iloc[0]
                    dslinkno        = float(outlets_bif_gdf['DSLINKNO'].iloc[0])

                    ds_gdf = clipped[clipped['LINKNO']==dslinkno].reset_index(drop=True)

                    # First correct if the upstream channel to this new one is lake out >:(
                    ds_us_linknos = [ds_gdf.loc[0,'USLINKNO1'],ds_gdf.loc[0,'USLINKNO2']]

                    if ds_us_linknos[0] == linkno_problem:
                        ds_us_linkno = ds_us_linknos[1]
                    
                    else:
                        ds_us_linkno = ds_us_linknos[0]
                    
                    
                    ds_us_gdf = clipped[clipped['LINKNO'] == ds_us_linkno].reset_index(drop=True)


                    # Reclassify them

                    #  1: The problematic one
                    outlets_bif_gdf.loc[0,'lake_flag'] = 'LakeWithin'

                    #  2: The one downstream
                    ds_gdf['lake_flag'] = 'LakeOut'
                    seg_ds_gdf = lines_to_segments(ds_gdf)

                    last_mask = seg_ds_gdf['segment_id'].eq(
                    seg_ds_gdf.groupby('LINKNO')['segment_id'].transform('max'))
                    seg_ds_gdf.loc[last_mask, 'lake_flag'] = 'LakeWithin'

                    seg_remerg_ds_gdf = remerge_segments(seg_ds_gdf)


                    # 3: Adjust other upstream of new outlet
                    for idx, row in ds_us_gdf.iterrows():
                        if row['lake_flag']!='Additional_OK':
                            ds_us_gdf.loc[idx,'lake_flag'] = 'LakeWithin'


                    outlets_bif_gdf = pd.concat([outlets_bif_gdf,seg_remerg_ds_gdf,ds_us_gdf])

                    outlets_bif_gdf['LakeId']       = lake_id
                    outlets_bif_gdf['connected']    = True

                    new_fixed = outlets_bif_gdf.copy() 
                    new_fixed['count_end'] = 1
                    new_fixed['count_start'] = 1

                else:
                    last_mask       = segments_gdf['segment_id'].eq(segments_gdf.groupby('LINKNO')['segment_id'].transform('max'))
                    segments_gdf.loc[last_mask, 'lake_flag'] = 'LakeWithin'
                    remerged_segments = remerge_segments(segments_gdf)

                    new_fixed = pd.merge(remerged_segments,outlets_bif_gdf[cols],how='left',on='LINKNO')
                    new_fixed = gpd.GeoDataFrame(data=new_fixed,geometry=new_fixed['geometry'])

                    new_fixed['count_end'] = 1
                    new_fixed['count_start'] = 1

                fixed_clipped = pd.concat([clipped[clipped['count_start']==1].reset_index(drop=True),new_fixed.reset_index(drop=True)]).reset_index(drop=True)
                clipped = fixed_clipped.copy()
            
            clipped['start'], clipped['end'] = zip(*clipped.geometry.map(get_endpoints))
            
            # One last check for bifurcations (no fix, allows up to 2 consecutive bifurcations, if more, buffer the lake)
            biffurcations_ok = True

            for idx, row in clipped.iterrows():
                if row['lake_flag'] == 'LakeOut':
                    start = row['start']
                    count = other_links_count(start,idx,clipped)
                    if count>1:
                        biffurcations_ok = False

            for idx, row in clipped.iterrows():
                if row['lake_flag'] == 'LakeIn':
                    start = row['end']
                    count = other_links_count(start,idx,clipped)
                    if count>1:
                        biffurcations_ok = False          


            # Final check: If there are lake-within but outside reservoir/lake, then flag them as NaN
            # Basically, if the lake within has to be connected on both ends, if it is not connected, there is an issue -> Flag it as Nan
            # But also, this is only the case if this channel is outside the lake geometry (lake geometry DOES NOT cover the channel)

            buffer_for_check = buffered.buffer(20, join_style=2, mitre_limit=1e6)

            for idx, row in clipped.iterrows():
                if row['lake_flag'] != 'LakeWithin':
                    continue

                start, end = row['start'], row['end']
                start_in, end_in = buffer_for_check.contains(start), buffer_for_check.contains(end)


                # Check if both endpoints are connected to at least one other LakeWithin segment
                connected_start = False
                connected_end = False
                for j, other in clipped.iterrows():
                    if j == idx:
                        continue
                    
                    o_start, o_end = other['start'], other['end']

                    if start.equals(o_start) or start.equals(o_end):
                        connected_start = True
                    
                    if end.equals(o_start) or end.equals(o_end):
                        connected_end = True
                    
                    if connected_start and connected_end:
                        break

                # If either endpoint is not connected to another LakeWithin, flag as NaN
                if not (connected_start and connected_end) and not (start_in and end_in):
                    linkno = row['LINKNO']
                    print(' The stupid last check is failing')
                    print(f'linkno problem: {linkno}')
                    clipped.loc[idx, 'lake_flag'] = np.nan



            # If all ok, break while
            if clipped['lake_flag'].isna().sum() == 0 and biffurcations_ok:

                clipped = clipped[['lake_flag','LINKNO', 'DSLINKNO', 'USLINKNO1', 'USLINKNO2', 'DSNODEID', 'strmOrder',
                        'Length', 'Magnitude', 'DSContArea', 'strmDrop', 'Slope', 'StraightL',
                        'USContArea', 'WSNO', 'DOUTEND', 'DOUTSTART', 'DOUTMID', 'BasinNo',
                        'LakeId', 'connected', 'length_n','geometry']]

                break


            attempts += 1

            # connected_channels.append(clipped)


            
        res_gdf.loc[index,'geometry'] = buffered

        


        # If there are NaN 
    
        if 'clipped' not in locals() or clipped.empty or clipped['lake_flag'].isna().any():
            print(f"\t\t Lake: id = {lake_id} ; name {lake_name} >  cannot be resolved into the current channel structure within buffer threshold.")
            

            if 'clipped' in locals() and not clipped.empty:
                problematic = clipped.copy()
                # print(f"\t\t Problematic segments for Lake {lake_id}:")
                
                # for i, row in problematic.iterrows():
                #     print(f"\t\t   - LINKNO: {row['LINKNO']}")
                unresolved_channels.append(problematic)


            continue


        elif not ((clipped['lake_flag'] == 'LakeIn').any() and (clipped['lake_flag'] == 'LakeOut').any()):
            print(f"\t\t Lake {lake_id} cannot be resolved into the current channel structure - It either has no outlet or no inlets.")

        elif not biffurcations_ok:
            print(f"\t\t Lake {lake_id} cannot be resolved into the current channel structure - Issue with channel biffurcations on boundaries")

        else:
            print(f"\t\t Lake: id = {lake_id} ; name = {lake_name} >   was resolved after {attempts} attempts.")
            res_gdf.loc[index,'resolved'] = True
            flagged_chans.append(clipped)


    
    combined_gdf    = gpd.GeoDataFrame(pd.concat(flagged_chans, ignore_index=True), crs=clipped.crs)
    # connected_gdf   = gpd.GeoDataFrame(pd.concat(connected_channels, ignore_index=True), crs=clipped.crs)
    unresolved_gdf  = res_gdf[~res_gdf['resolved']].copy()
    res_gdf         = res_gdf[res_gdf['resolved']].copy()

    # Make small buffer in the end to avoid weird intersections not being taken into account
    buffered = res_gdf['geometry'].buffer(20, join_style=2, mitre_limit=1e6)

    res_gdf['geometry'] = buffered

    # if unresolved_channels:
    #     problem_gdf = gpd.GeoDataFrame(pd.concat(unresolved_channels, ignore_index=True), crs=clipped.crs)


    # Save resolved lakes    
    chan_fn         = os.path.join(projDir,'Watershed/Shapes/lakes-flagged-channels-{auth}-{code}.shp'.format(**details))   
    res_fn          = os.path.join(projDir,'Watershed/Shapes/lakes-grand-{auth}-{code}.shp'.format(**details))
    unres_shp       = os.path.join(projDir,'Watershed/Shapes/unresolved-{auth}-{code}.shp'.format(**details))
    unflag_shp      = os.path.join(projDir,'Watershed/Shapes/unflagged-channels-{auth}-{code}.shp'.format(**details))
    conn_shp        = os.path.join(projDir,'Watershed/Shapes/connected-lake-channels-{auth}-{code}.shp'.format(**details))

    # Format to export
    combined_gdf['Length']      = combined_gdf['Length'].map(lambda x: f"{x:.2f}")
    combined_gdf['DSContArea']  = combined_gdf['DSContArea'].map(lambda x: f"{x:.1f}")
    combined_gdf['strmDrop']    = combined_gdf['strmDrop'].map(lambda x: f"{x:.2f}")
    combined_gdf['StraightL']   = combined_gdf['StraightL'].map(lambda x: f"{x:.1f}")
    combined_gdf['USContArea']  = combined_gdf['USContArea'].map(lambda x: f"{x:.1f}")
    combined_gdf['DOUTEND']     = combined_gdf['DOUTEND'].map(lambda x: f"{x:.1f}")
    combined_gdf['DOUTSTART']   = combined_gdf['DOUTSTART'].map(lambda x: f"{x:.1f}")
    combined_gdf['DOUTMID']     = combined_gdf['DOUTMID'].map(lambda x: f"{x:.1f}")

    # if unresolved_channels:
    #     problem_gdf['Length']      = problem_gdf['Length'].map(lambda x: f"{x:.2f}")
    #     problem_gdf['DSContArea']  = problem_gdf['DSContArea'].map(lambda x: f"{x:.1f}")
    #     problem_gdf['strmDrop']    = problem_gdf['strmDrop'].map(lambda x: f"{x:.2f}")
    #     problem_gdf['StraightL']   = problem_gdf['StraightL'].map(lambda x: f"{x:.1f}")
    #     problem_gdf['USContArea']  = problem_gdf['USContArea'].map(lambda x: f"{x:.1f}")
    #     problem_gdf['DOUTEND']     = problem_gdf['DOUTEND'].map(lambda x: f"{x:.1f}")
    #     problem_gdf['DOUTSTART']   = problem_gdf['DOUTSTART'].map(lambda x: f"{x:.1f}")
    #     problem_gdf['DOUTMID']     = problem_gdf['DOUTMID'].map(lambda x: f"{x:.1f}")

    try:
        res_gdf['calcVol']          = res_gdf['calcVol'].map(lambda x: f"{x:.1f}")
        unresolved_gdf['calcVol']   = unresolved_gdf['calcVol'].map(lambda x: f"{x:.1f}")

    except:
        print('calcVol not Formated')


    res_gdf.to_file(res_fn)
    combined_gdf.to_file(chan_fn)
    unresolved_gdf.to_file(unres_shp)
    # connected_gdf.to_file(conn_shp)

    # if unresolved_channels:
    #     problem_gdf.to_file(unflag_shp)



    '''
    Creation of snapping points
    '''

    snap_points = []
    tolerance = 1e-6

    for idx, row in combined_gdf.iterrows():
        flag    = row['lake_flag']
        lakeId  = row['LakeId']
        linkno  = row['LINKNO']

        start_point, end_point = get_endpoints(row.geometry)

        if flag == "LakeOut":
            res = 1
            point = start_point
        elif flag == "LakeIn":
            res = 0
            point = end_point
        else:
            continue
        
        snap_points.append({'geometry': point, 'RES': res, 'LakeId': lakeId})

    snap_gdf = gpd.GeoDataFrame(snap_points, geometry='geometry', crs=combined_gdf.crs)

    # Most important output: Snapping points to be read on QSWAT+
    snap_gdf = gpd.GeoDataFrame(snap_points,crs=combined_gdf.crs)
    snap_gdf = snap_gdf.drop_duplicates(subset='geometry')
    snap_gdf.to_file(variables.resolved_snaps.format(**details))



    



    # def solve_noDataBasins(res_gdf, subs_gdf, chan_gdf):
    #     noDataBasins = True
    #     print("\t \t > There could be some no data subbasins adjacent to lakes, they will be added as part of the lake to avoid issues")
    #     while noDataBasins:
    #         buffered_res = res_gdf.copy()
    #         buffered_res['geometry'] = buffered_res['geometry'].buffer(variables.data_resolution/2, join_style=2, mitre_limit=1e6)

    #         expected_hru_gdf = gpd.overlay(subs_gdf, buffered_res, how='difference', keep_geom_type = True, make_valid = True)
    #         subs_overlay_gdf = gpd.overlay(subs_gdf, res_gdf, how='difference', keep_geom_type = True, make_valid = True)
    #         expected_hru_gdf = expected_hru_gdf.dissolve()
            
    #         inter = gpd.sjoin(subs_overlay_gdf, expected_hru_gdf[['geometry']], how="inner", predicate="intersects")
    #         touch = gpd.sjoin(subs_overlay_gdf, expected_hru_gdf[['geometry']], how="inner", predicate="touches")

    #         # keep subs with zero positive-area overlap
    #         hit_idx = inter.index.difference(touch.index)
    #         no_data_subs = subs_overlay_gdf.loc[~subs_overlay_gdf.index.isin(hit_idx)]

    #         if len(no_data_subs)>0:
    #             # crossed = gpd.sjoin(no_data_subs, chan_gdf, how="inner", predicate="crosses")
    #             # result  = no_data_subs.loc[no_data_subs.index.isin(crossed.index)]                                                          # Only those that cross channels matter
    #             # no_data_subs = result.copy()
            
    #             buf = no_data_subs.copy()
    #             buf["geometry"] = buf.geometry.buffer(variables.data_resolution/2, join_style=2, mitre_limit=1e6)
    #             buf = buf[~buf.geometry.is_empty & buf.geometry.notna()]


    #             #Map each buffered sub to nearest lake
    #             nearest = gpd.sjoin_nearest(
    #                 buf[["geometry"]],
    #                 res_gdf[["geometry"]],
    #                 how="left",
    #                 distance_col="dist").rename(columns={"index_right":"lake_idx"})



    #             for lake_idx, grp in nearest.groupby("lake_idx"):
    #                 if lake_idx is None or len(grp)==0:
    #                     continue
    #                 to_add = _union(grp.geometry)
    #                 res_gdf.at[lake_idx, "geometry"] = res_gdf.geometry.at[lake_idx].union(to_add)
            
    #         else:
    #             print("\t \t \t ... No potential no data sub basins anymore")
    #             noDataBasins = False
    #             solved = True

    #     return res_gdf,solved
