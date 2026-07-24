"""
Author: Jose P. Teran
Date: 2026-03-31

Description:
Adjust lake geometry to allow integration into network
"""

# Basic info

import os
import sys
import datavariables as variables
import argparse
import geopandas as gpd
from shapely import BufferCapStyle, BufferJoinStyle, Point
from shapely.ops import unary_union, polygonize
from shapely.strtree import STRtree
import pandas as pd
import numpy as np
import warnings
import rasterio
from rasterio import features
from rasterio.features import rasterize
from shapely.geometry import Polygon, MultiPolygon, shape, box, MultiLineString, LineString
from rasterio import features
from shapely.strtree import STRtree

def remove_holes(geometry):
    if geometry.geom_type == 'Polygon':
        return Polygon(geometry.exterior)
    elif geometry.geom_type == 'MultiPolygon':
        return MultiPolygon([Polygon(part.exterior) for part in geometry.geoms])
    else:
        return geometry
    
def rasterize_lakes(dem_fp, lakes_gdf, out_fp, attr="Hylak_id", all_touched=True, nodata=0, scale=1):
    with rasterio.open(dem_fp) as dem:
        profile   = dem.profile.copy()
        transform = dem.transform
        crs       = dem.crs
        height    = dem.height
        width     = dem.width

    fine_transform = transform * transform.scale(1/scale, 1/scale)
    fine_height = height * scale
    fine_width = width * scale

    gdf = lakes_gdf.copy()
    gdf['geometry'] = gdf['geometry'].apply(remove_holes)
    if gdf.crs != crs:
        gdf = gdf.to_crs(crs)
    if attr not in gdf.columns:
        raise ValueError(f"Attribute '{attr}' not found.")

    gdf = gdf[gdf[attr].notna() & gdf.geometry.notna() & (~gdf.geometry.is_empty)]
    vals = gdf[attr].astype(np.uint32)

    arr = np.full((fine_height, fine_width), nodata, dtype=np.uint32)
    shapes = zip(gdf.geometry, vals.astype(int))
    rasterize(
        shapes=shapes,
        out=arr,
        transform=fine_transform,
        all_touched=all_touched,
    )

    profile.update(dtype="uint32", nodata=nodata, count=1, compress="lzw",
                   transform=fine_transform, height=fine_height, width=fine_width)
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
        recs.append({"geometry": g, "Hylak_id": v})

    gdf = gpd.GeoDataFrame(recs, crs=crs)
    if dissolve and not gdf.empty:
        gdf = gdf.dissolve(by="Hylak_id", as_index=False)
    return gdf


warnings.filterwarnings('ignore', module='pyogrio')
warnings.filterwarnings('ignore', message='Column names longer than 10 characters')


# Functions
def classify_cat0(lakes_snap_gdf, chans_gdf):
    case_0_list = []
    lake_boundary = lakes_snap_gdf.copy()
    lake_boundary['geometry'] = lakes_snap_gdf.geometry.boundary
    for idxLake, lake in lakes_snap_gdf.iterrows():
        lakebox = lake.geometry.bounds
        lake_perim = lake_boundary.loc[lake.name, 'geometry']
        candidates = chans_gdf.cx[lakebox[0]:lakebox[2], lakebox[1]:lakebox[3]]
        for _, chan in candidates.iterrows():
            coords = list(chan.geometry.coords)
            source = Point(coords[-1])
            outlet = Point(coords[0])
            source_inside = lake.geometry.contains(source)
            outlet_inside = lake.geometry.contains(outlet)
            if source_inside and outlet_inside:
                crosses = chan.geometry.intersects(lake_perim)
                case_0_list.append({'chanId': chan.name,'lakeId': lake['LakeId'],'category': 0,'flagged': crosses,'geometry': chan.geometry})
    return gpd.GeoDataFrame(case_0_list, crs=chans_gdf.crs)

def fix_cat0(lakes_snap_gdf, cat0_gdf, dem_res):
    cat_0_flagged = cat0_gdf[cat0_gdf['flagged']==True]
    for lake_id, group in cat_0_flagged.groupby('lakeId'):
        lake_row    = lakes_snap_gdf[lakes_snap_gdf['LakeId'] == lake_id].iloc[0]
        lake_geom   = lake_row.geometry
        for _, row in group.iterrows():
            chan_geom   = row['geometry']
            lake_perim  = lake_geom.boundary
            outside_segment = chan_geom.difference(lake_geom)
            combined    = unary_union([outside_segment, lake_perim])
            polygons    = list(polygonize(combined))
            notch       = min(polygons, key=lambda p: p.area)
            lake_geom   = unary_union([lake_geom, notch]).buffer(dem_res * 0.05,cap_style=BufferCapStyle.square,join_style=BufferJoinStyle.mitre)
        lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'geometry'] = lake_geom
    return lakes_snap_gdf

def classify_cat1(lakes_snap_gdf, chans_gdf, dem_res):
    lake_boundary = lakes_snap_gdf.copy()
    lake_boundary['geometry'] = lakes_snap_gdf.geometry.boundary
    results = []
    for _, lake in lakes_snap_gdf.iterrows():
        lake_box = lake.geometry.bounds
        lake_perim = lake_boundary.loc[lake.name, 'geometry']
        candidates = chans_gdf.cx[lake_box[0]:lake_box[2], lake_box[1]:lake_box[3]]
        for _, chan in candidates.iterrows():
            coords = list(chan.geometry.coords)
            source = Point(coords[-1])
            outlet = Point(coords[0])
            source_inside = lake.geometry.contains(source)
            outlet_inside = lake.geometry.contains(outlet)
            if not source_inside and outlet_inside:
                inside_segment  = chan.geometry.intersection(lake.geometry)
                outside_segment = chan.geometry.difference(lake.geometry)
                inside_length   = inside_segment.length  if not inside_segment.is_empty  else 0
                outside_length  = outside_segment.length if not outside_segment.is_empty else 0
                intersection    = chan.geometry.intersection(lake_perim)
                num_crossings   = len(intersection.geoms) if hasattr(intersection, 'geoms') else 1
                flagged_short_inside  = inside_length  < dem_res
                flagged_short_outside = outside_length < dem_res*0.95
                flagged_crossing      = num_crossings  > 1
                results.append({'chanId': chan.name, 'lakeId': lake['LakeId'], 'category': 1, 'inside_length': inside_length, 'outside_length': outside_length, 
                                'num_crossings': num_crossings, 'flagged_short_inside': flagged_short_inside, 'flagged_short_outside': flagged_short_outside, 
                                'flagged_crossing': flagged_crossing, 'geometry': chan.geometry})
    return gpd.GeoDataFrame(results, crs=chans_gdf.crs)


def fix_cat1(lakes_snap_gdf, cat1_gdf, chans_gdf, dem_res):
    cat1_flagged = cat1_gdf[cat1_gdf['flagged_crossing'] | cat1_gdf['flagged_short_inside'] | cat1_gdf['flagged_short_outside']]
    for lake_id, group in cat1_flagged.groupby('lakeId'):
        lake_row = lakes_snap_gdf[lakes_snap_gdf['LakeId'] == lake_id].iloc[0]
        lake_geom = lake_row.geometry
        for _, row in group.iterrows():
            chan_geom = row['geometry']
            if row['flagged_crossing']:
                lake_perim      = lake_geom.boundary
                outside_segment = chan_geom.difference(lake_geom)
                if outside_segment.is_empty:
                    continue
                combined = unary_union([outside_segment, lake_perim])
                polygons = list(polygonize(combined))
                if not polygons:
                    continue
                notch       = min(polygons, key=lambda p: p.area)
                lake_geom   = unary_union([lake_geom, notch]).buffer(dem_res * 0.05, cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)

            if row['flagged_short_inside']:
                coords    = list(chan_geom.coords)
                outlet    = Point(coords[0])
                half      = dem_res * 1.05
                square    = box(outlet.x - half, outlet.y - half, outlet.x + half, outlet.y + half)
                lake_geom = unary_union([lake_geom, square])

            if row['flagged_short_outside']:
                outside_segment = chan_geom.difference(lake_geom)
                if not outside_segment.is_empty:
                    # buffer the outside segment to cover diagonal cases
                    outside_buffered = outside_segment.buffer(dem_res * 0.5, cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)
                    coords      = list(chan_geom.coords)
                    source      = Point(coords[-1]) # add a square to source
                    half        = dem_res * 1.05
                    square      = box(source.x - half, source.y - half, source.x + half, source.y + half)
                    lake_geom   = unary_union([lake_geom, outside_buffered, square])
        lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'geometry'] = lake_geom
    return lakes_snap_gdf

def runFix_cat1(lakes_snap_gdf, chans_gdf, dem_res):
    cat1_fixed = False
    cat1_iter  = 10
    cat1_iteration = 0
    while not cat1_fixed and cat1_iteration < cat1_iter:
        cat1_iteration += 1
        cat1_gdf = classify_cat1(lakes_snap_gdf, chans_gdf, dem_res)
        cat1_flagged = cat1_gdf[cat1_gdf['flagged_crossing'] | cat1_gdf['flagged_short_inside'] | cat1_gdf['flagged_short_outside']]
        if cat1_flagged.empty:
            cat1_fixed = True
        else:
            lakes_snap_gdf = fix_cat1(lakes_snap_gdf, cat1_gdf, chans_gdf, dem_res)
    return lakes_snap_gdf,cat1_fixed

def classify_cat2(lakes_snap_gdf, chans_gdf, dem_res):
    lake_boundary = lakes_snap_gdf.copy()
    lake_boundary['geometry'] = lakes_snap_gdf.geometry.boundary
    dem_diagonal = (dem_res**2 + dem_res**2)**0.5
    tol = dem_res * 0.1
    results = []
    for _, lake in lakes_snap_gdf.iterrows():
        lake_box   = lake.geometry.bounds
        lake_perim = lake_boundary.loc[lake.name, 'geometry']
        candidates = chans_gdf.cx[lake_box[0]:lake_box[2], lake_box[1]:lake_box[3]]
        for _, chan in candidates.iterrows():
            coords = list(chan.geometry.coords)
            source = Point(coords[-1])
            outlet = Point(coords[0])
            source_inside = lake.geometry.contains(source)
            outlet_inside = lake.geometry.contains(outlet)
            if source_inside and not outlet_inside:
                inside_segment  = chan.geometry.intersection(lake.geometry)
                outside_segment = chan.geometry.difference(lake.geometry)
                inside_length   = inside_segment.length  if not inside_segment.is_empty else 0
                outside_length  = outside_segment.length if not outside_segment.is_empty else 0
                intersection    = chan.geometry.intersection(lake_perim)
                num_crossings   = len(intersection.geoms) if hasattr(intersection, 'geoms') else 1

                flagged_short_inside  = inside_length < dem_res
                if flagged_short_inside:
                    flagged_short_outside = False
                else:
                    flagged_short_outside = outside_length < dem_res*0.95
                flagged_crossing      = num_crossings > 1

                # detect special case: channel is dem_res*2 or dem_diagonal*2
                dx = abs(outlet.x - source.x)
                dy = abs(outlet.y - source.y)
                angle = np.degrees(np.arctan2(dy, dx))
                is_diagonal = 15 < angle < 75
                expected_length = dem_diagonal * 2 if is_diagonal else dem_res * 2
                chan_length = chan.geometry.length
                flagged_special = flagged_short_inside and abs(chan_length - expected_length) < tol

                # if special, unflag short_inside to avoid double fixing
                if flagged_special:
                    flagged_short_inside = False

                results.append({'chanId': chan.name, 'lakeId': lake['LakeId'], 'category': 2, 'inside_length': inside_length, 'outside_length': outside_length, 'num_crossings': num_crossings, 'flagged_short_inside': flagged_short_inside, 'flagged_short_outside': flagged_short_outside, 'flagged_crossing': flagged_crossing, 'flagged_special': flagged_special, 'geometry': chan.geometry})
    return gpd.GeoDataFrame(results, crs=chans_gdf.crs)


def fix_cat2(lakes_snap_gdf, cat2_gdf, chans_gdf, dem_res):
    cat2_flagged = cat2_gdf[cat2_gdf['flagged_crossing'] | cat2_gdf['flagged_short_inside'] | cat2_gdf['flagged_short_outside'] | cat2_gdf['flagged_special']]
    for lake_id, group in cat2_flagged.groupby('lakeId'):
        lake_row  = lakes_snap_gdf[lakes_snap_gdf['LakeId'] == lake_id].iloc[0]
        lake_geom = lake_row.geometry
        for _, row in group.iterrows():
            chan_geom = row['geometry']

            if row['flagged_crossing']:
                lake_perim      = lake_geom.boundary
                outside_segment = chan_geom.difference(lake_geom)
                if outside_segment.is_empty:
                    continue
                combined = unary_union([outside_segment, lake_perim])
                polygons = list(polygonize(combined))
                if not polygons:
                    continue
                notch     = min(polygons, key=lambda p: p.area)
                lake_geom = unary_union([lake_geom, notch]).buffer(dem_res * 0.05, cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)

            if row['flagged_short_inside']:
                coords = list(chan_geom.coords)
                source = Point(coords[-1])
                half   = dem_res * 1.05
                square = box(source.x - half, source.y - half, source.x + half, source.y + half)
                lake_geom = unary_union([lake_geom, square])

            if row['flagged_short_outside']:
                coords = list(chan_geom.coords)
                source = Point(coords[0])
                half   = dem_res * 1.05
                square = box(source.x - half, source.y - half, source.x + half, source.y + half)
                lake_geom = unary_union([lake_geom, square])

            if row['flagged_special']:
                midpoint = chan_geom.interpolate(0.75, normalized=True)
                source   = Point(list(chan_geom.coords)[-1])
                half_line = LineString([source, midpoint])
                half_buffered = half_line.buffer(dem_res / 2, cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)
                lake_geom = unary_union([lake_geom, half_buffered])

        lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'geometry'] = lake_geom
    return lakes_snap_gdf


def runFix_cat2(lakes_snap_gdf, chans_gdf, dem_res):
    cat2_fixed     = False
    cat2_iter      = 10
    cat2_iteration = 0
    fixed_special  = set()
    while not cat2_fixed and cat2_iteration < cat2_iter:
        cat2_iteration += 1
        cat2_gdf = classify_cat2(lakes_snap_gdf, chans_gdf, dem_res)
        # collect special channels before unflagging
        new_special = set(cat2_gdf[cat2_gdf['flagged_special']]['chanId'].tolist())
        # unflag already processed special channels
        cat2_gdf.loc[cat2_gdf['chanId'].isin(fixed_special), 'flagged_special'] = False
        cat2_gdf.loc[cat2_gdf['chanId'].isin(fixed_special), 'flagged_short_inside'] = False
        cat2_flagged = cat2_gdf[cat2_gdf['flagged_crossing'] | cat2_gdf['flagged_short_inside'] | cat2_gdf['flagged_short_outside'] | cat2_gdf['flagged_special']]
        
        if cat2_flagged.empty:
            cat2_fixed = True
        else:
            fixed_special.update(new_special)
            lakes_snap_gdf = fix_cat2(lakes_snap_gdf, cat2_gdf, chans_gdf, dem_res)
    
    return lakes_snap_gdf, cat2_fixed

def classify_cat3(lakes_snap_gdf, chans_gdf, dem_res):
    lake_boundary = lakes_snap_gdf.copy()
    lake_boundary['geometry'] = lakes_snap_gdf.geometry.boundary
    dem_diagonal = (dem_res**2 + dem_res**2)**0.5
    results = []
    for _, lake in lakes_snap_gdf.iterrows():
        lake_box  = lake.geometry.bounds
        lake_perim = lake_boundary.loc[lake.name, 'geometry']
        candidates = chans_gdf.cx[lake_box[0]:lake_box[2], lake_box[1]:lake_box[3]]
        for _, chan in candidates.iterrows():
            coords = list(chan.geometry.coords)
            source = Point(coords[-1])
            outlet = Point(coords[0])
            source_inside = lake.geometry.contains(source)
            outlet_inside = lake.geometry.contains(outlet)
            if not source_inside and not outlet_inside:
                intersection  = chan.geometry.intersection(lake_perim)
                num_crossings = len(intersection.geoms) if hasattr(intersection, 'geoms') else 1
                if num_crossings >= 2:
                    flagged_multi = num_crossings > 2
                    # for exactly 2 crossings, check if internal segment is smaller than diagonal
                    flagged_short_internal = False
                    if num_crossings == 2:
                        internal_segment = chan.geometry.intersection(lake.geometry)
                        if not internal_segment.is_empty and internal_segment.length < dem_diagonal:
                            flagged_short_internal = True
                    results.append({'chanId': chan.name, 'lakeId': lake['LakeId'], 'category': 3, 'num_crossings': num_crossings, 'flagged': flagged_multi, 'flagged_short_internal': flagged_short_internal, 'geometry': chan.geometry})
    if not results:
        return gpd.GeoDataFrame(columns=['chanId', 'lakeId', 'category', 'num_crossings', 'flagged', 'flagged_short_internal', 'geometry'], geometry='geometry', crs=chans_gdf.crs)
    return gpd.GeoDataFrame(results, crs=chans_gdf.crs)



def fix_cat3(lakes_snap_gdf, cat3_gdf, chans_gdf, dem_res):
    cat3_flagged = cat3_gdf[cat3_gdf['flagged'] | cat3_gdf['flagged_short_internal']]
    for lake_id, group in cat3_flagged.groupby('lakeId'):
        lake_row  = lakes_snap_gdf[lakes_snap_gdf['LakeId'] == lake_id].iloc[0]
        lake_geom = lake_row.geometry
        for _, row in group.iterrows():
            chan_geom = row['geometry']

            if row['flagged']:
                outside_segment = chan_geom.difference(lake_geom)
                if outside_segment.is_empty:
                    continue
                lake_perim = lake_geom.boundary
                combined  = unary_union([outside_segment, lake_perim])
                polygons  = list(polygonize(combined))
                if not polygons:
                    continue
                notches = [p for p in polygons if p.area < lake_geom.area]
                if not notches:
                    continue
                lake_geom = unary_union([lake_geom] + notches).buffer(dem_res * 0.05, cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)


            if row['flagged_short_internal']:
                internal_segment = chan_geom.intersection(lake_geom)
                if internal_segment.is_empty:
                    continue
                carve     = internal_segment.buffer(dem_res * 0.05, cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)
                lake_geom = lake_geom.difference(carve)
                # keep only the largest geometry
                if hasattr(lake_geom, 'geoms'):
                    lake_geom = max(lake_geom.geoms, key=lambda g: g.area)

        lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'geometry'] = lake_geom
    return lakes_snap_gdf


def runFix_cat3(lakes_snap_gdf, chans_gdf, dem_res):
    cat3_fixed     = False
    cat3_iter      = 10
    cat3_iteration = 0
    while not cat3_fixed and cat3_iteration < cat3_iter:
        cat3_iteration += 1
        cat3_gdf     = classify_cat3(lakes_snap_gdf, chans_gdf, dem_res)
        cat3_flagged = cat3_gdf[cat3_gdf['flagged'] | cat3_gdf['flagged_short_internal']]
        if cat3_flagged.empty:
            cat3_fixed = True
        else:
            lakes_snap_gdf = fix_cat3(lakes_snap_gdf, cat3_gdf, chans_gdf, dem_res)
    return lakes_snap_gdf, cat3_fixed

def remove_holes(geometry):
    if geometry.geom_type == 'Polygon':
        return Polygon(geometry.exterior)
    elif geometry.geom_type == 'MultiPolygon':
        return MultiPolygon([Polygon(part.exterior) for part in geometry.geoms])
    else:
        return geometry

def carve_short_interlake_channels(lakes_snap_gdf, chans_gdf, dem_res):
    half = dem_res
    tree = STRtree(lakes_snap_gdf.geometry.values)

    chans_outside = gpd.overlay(chans_gdf, lakes_snap_gdf, how='difference', keep_geom_type=True)
    chans_outside = chans_outside.explode(index_parts=False).reset_index(drop=True)

    def get_endpoints(geom):
        if isinstance(geom, MultiLineString):
            geoms = list(geom.geoms)
        else:
            geoms = [geom]
        source = Point(list(geoms[-1].coords)[-1])
        outlet = Point(list(geoms[0].coords)[0])
        return source, outlet

    def both_endpoints_touch_lake(chan_geom):
        source, outlet = get_endpoints(chan_geom)
        source_touches = any(lakes_snap_gdf.iloc[cidx].geometry.touches(source) for cidx in tree.query(source))
        outlet_touches = any(lakes_snap_gdf.iloc[cidx].geometry.touches(outlet) for cidx in tree.query(outlet))
        return source_touches and outlet_touches

    chans_outside['both_touch'] = chans_outside.geometry.apply(both_endpoints_touch_lake)
    min_len = (dem_res**2 + dem_res**2)**0.5
    chans_interlake = chans_outside[chans_outside['both_touch'] & (chans_outside.geometry.length < min_len)]

    for _, chan in chans_interlake.iterrows():
        source, outlet = get_endpoints(chan.geometry)

        source_lake_idx = None
        outlet_lake_idx = None
        for cidx in tree.query(chan.geometry):
            lake = lakes_snap_gdf.iloc[cidx]
            if source_lake_idx is None and lake.geometry.touches(source):
                source_lake_idx = lakes_snap_gdf.index[cidx]
            if outlet_lake_idx is None and lake.geometry.touches(outlet):
                outlet_lake_idx = lakes_snap_gdf.index[cidx]
            if source_lake_idx is not None and outlet_lake_idx is not None:
                break

        if source_lake_idx is None or outlet_lake_idx is None or source_lake_idx == outlet_lake_idx:
            continue

        carve_source = box(source.x - half, source.y - half, source.x + half, source.y + half)
        carve_outlet = box(outlet.x - half, outlet.y - half, outlet.x + half, outlet.y + half)
        lakes_snap_gdf.loc[source_lake_idx, 'geometry'] = lakes_snap_gdf.loc[source_lake_idx, 'geometry'].difference(carve_source)
        lakes_snap_gdf.loc[outlet_lake_idx, 'geometry'] = lakes_snap_gdf.loc[outlet_lake_idx, 'geometry'].difference(carve_outlet)

    return lakes_snap_gdf

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

def classify_channels(lakes_snap_gdf, chans_gdf, dem_res, tolerance=None):
    if tolerance is None:
        tolerance = dem_res * 0.1
    
    tree = STRtree(lakes_snap_gdf.geometry.values)
    
    chans_outside = gpd.overlay(chans_gdf, lakes_snap_gdf, how='difference', keep_geom_type=True)
    chans_outside = chans_outside.drop(columns=[c for c in lakes_snap_gdf.columns if c != 'geometry' and c in chans_outside.columns], errors='ignore')
    chans_outside = chans_outside.explode(index_parts=False).reset_index(drop=True)

    chans_outside['LakeIn']  = None
    chans_outside['LakeOut'] = None

    for idx, chan in chans_outside.iterrows():
        source, outlet = get_endpoints(chan.geometry)
        
        for cidx in tree.query(chan.geometry):
            lake = lakes_snap_gdf.iloc[cidx]
            lake_id = lake['LakeId']
            
            if lake.geometry.distance(source) < tolerance:
                chans_outside.at[idx, 'LakeOut'] = lake_id
            if lake.geometry.distance(outlet) < tolerance:
                chans_outside.at[idx, 'LakeIn'] = lake_id
    
    return chans_outside

def flag_infinite_loop_lakes(lakes_snap_gdf, combined_gdf):
    lakes_snap_gdf['infinite_loop'] = False
    loop_lines = []
    for idx, lake in lakes_snap_gdf.iterrows():
        lake_id = lake['LakeId']
        outlets = set(combined_gdf[combined_gdf['LakeOut'] == lake_id]['LINKNO'].tolist())
        inlets  = combined_gdf[combined_gdf['LakeIn']  == lake_id]['LINKNO'].tolist()
        if not outlets or not inlets:
            continue
        
        for inlet in inlets: 
            upstream = set()
            current_level = {inlet}
            for _ in range(3):           # for each inlet, go upstream 3 levels and check if we hit an outlet
                next_level = set()
                for linkno in current_level:
                    chan = combined_gdf[combined_gdf['LINKNO'] == linkno]
                    if chan.empty:
                        continue
                    for us_col in ['USLINKNO1', 'USLINKNO2']:
                        us_vals = chan[us_col].dropna().tolist()
                        for uslinkno in us_vals:
                            if int(uslinkno) >= 0:
                                next_level.add(int(uslinkno))
                upstream.update(next_level)
                current_level = next_level
            loop_outlets = upstream.intersection(outlets)
            if not loop_outlets:
                continue
            lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'infinite_loop'] = True

            # Trace path of infinite loop
            for loop_outlet in loop_outlets:
                path_links = [loop_outlet]
                current = loop_outlet

                for _ in range(3):
                    chan = combined_gdf[combined_gdf['LINKNO']==current]
                    dslinkno = int(chan['DSLINKNO'].iloc[0])
                    if dslinkno == inlet:
                        path_links.append(inlet)
                        break
                    path_links.append(dslinkno)
                    current = dslinkno

                # Create geometry
                outlet_segs = combined_gdf[(combined_gdf['LINKNO'] == loop_outlet) & (combined_gdf['LakeOut'] == lake_id)]
                inlet_segs  = combined_gdf[(combined_gdf['LINKNO'] == inlet) & (combined_gdf['LakeIn'] == lake_id)]
                if outlet_segs.empty or inlet_segs.empty:
                    continue
                path_geoms = []
                current_geom   = outlet_segs.iloc[0].geometry
                path_geoms.append(current_geom)
                source, outlet_pt = get_endpoints(current_geom)
                current_end    = outlet_pt
                current_linkno = loop_outlet
                tol = 300
                for _ in range(10):
                    chan = combined_gdf[combined_gdf['LINKNO'] == current_linkno]
                    ds_linkno = int(chan['DSLINKNO'].iloc[0])
                    if ds_linkno < 0:
                        break
                    candidates = combined_gdf[combined_gdf['LINKNO'] == ds_linkno]
                    next_seg = candidates[candidates.geometry.apply(lambda g: get_endpoints(g)[0].distance(current_end) < tol)]
                    if next_seg.empty:
                        break
                    next_geom = next_seg.iloc[0].geometry
                    path_geoms.append(next_geom)
                    current_linkno = ds_linkno
                    _, current_end = get_endpoints(next_geom)
                    if current_linkno == inlet:
                        break
                if path_geoms:
                    loop_lines.append({'lakeId':lake_id,'inlet':inlet,'outlet':loop_outlet,'geometry':unary_union(path_geoms)})
    loop_lines_gdf = gpd.GeoDataFrame(loop_lines, crs=combined_gdf.crs) if loop_lines else gpd.GeoDataFrame(columns=['lakeId', 'inlet', 'outlet', 'geometry'])
    return lakes_snap_gdf,loop_lines_gdf


def fix_infinite_loop_lakes(lakes_snap_gdf, loop_lines_gdf, dem_res):
    for _, loop in loop_lines_gdf.iterrows():
        lake_id   = loop['lakeId']
        lake_geom = lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'geometry'].iloc[0]
        area_to_add = loop.geometry.buffer(dem_res * 0.05, cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)
        lake_geom   = remove_holes(unary_union([lake_geom, area_to_add]))
        lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'geometry'] = lake_geom
    return lakes_snap_gdf

def runFix_inf(lakes_snap_gdf, combined_gdf, chans_gdf, dem_res):
    loop_fixed = False
    max_iter   = 10
    iteration  = 0
    while not loop_fixed and iteration < max_iter:
        iteration += 1
        lakes_snap_gdf, loop_lines_gdf = flag_infinite_loop_lakes(lakes_snap_gdf, combined_gdf)
        if lakes_snap_gdf[lakes_snap_gdf['infinite_loop'] == True].empty:
            loop_fixed = True
        else:
            lakes_snap_gdf = fix_infinite_loop_lakes(lakes_snap_gdf, loop_lines_gdf, dem_res)
            combined_gdf   = classify_channels(lakes_snap_gdf, chans_gdf, dem_res, tolerance=50)
            combined_gdf   = combined_gdf.explode(index_parts=False).reset_index(drop=True)
    return lakes_snap_gdf, combined_gdf


def checkSmallSubs(subs_gdf,lakes_snap_gdf):
    subs_clipped = gpd.overlay(subs_gdf, lakes_snap_gdf, how ='difference', keep_geom_type=True)
    subs_clipped['area_new'] = subs_clipped.geometry.area
    subs_clipped['flagged_small'] = subs_clipped['area_new'] < (dem_res * dem_res)*0.99999
    flagged_subs = subs_clipped[subs_clipped['flagged_small']]
    return flagged_subs

def checkNoHRUSubs(subs_gdf,lakes_snap_gdf):
    subs_clipped = gpd.overlay(subs_gdf, lakes_snap_gdf, how ='difference', keep_geom_type=True)
    subs_clipped = subs_clipped.explode(index_parts=False).reset_index(drop=True)
    lakes_outer = lakes_snap_gdf.copy()
    buff_dist = dem_res/2*0.95 - dem_res*0.05
    lakes_outer['geometry'] = lakes_snap_gdf.geometry.buffer(buff_dist, cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)
    lakes_border = lakes_snap_gdf.copy()
    lakes_border['geometry'] = lakes_outer.geometry.difference(lakes_border.geometry)
    lakes_border['geometry'] = lakes_border.geometry.buffer(dem_res*0.15, cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)

    border_union = lakes_border.geometry.union_all()
    lakes_border.to_file(f'{out_dir}/lakeBorder.gpkg')

    def mostly_inside(geom, threshold=0.98):
        if geom.is_empty:
            return False
        intersection = geom.intersection(border_union)
        return (intersection.area / geom.area) >= threshold

    subs_clipped['nodataflag'] = subs_clipped.geometry.apply(mostly_inside)

    # subs_clipped['nodataflag'] = subs_clipped.geometry.apply(lambda g: lakes_border.geometry.contains(g).any())
    flagged_subs = subs_clipped[subs_clipped['nodataflag']]
    flagged_subs.to_file(f'{out_dir}/flagged_subs.gpkg')
    return flagged_subs

def fixNoDataBasin(flagged_subsSmall, flagged_subsNoData, lakes_snap_gdf):
    tree = STRtree(lakes_snap_gdf.geometry.values)
    
    flagged_subs = gpd.GeoDataFrame(
        pd.concat([flagged_subsSmall, flagged_subsNoData], ignore_index=True),crs=flagged_subsSmall.crs).drop_duplicates(subset='geometry')
    
    for idx, sub in flagged_subs.iterrows():
        nearest_idx = tree.nearest(sub.geometry)
        lake_id  = lakes_snap_gdf.iloc[nearest_idx]['LakeId']
        lake_geom = lakes_snap_gdf.iloc[nearest_idx].geometry
        new_geom = unary_union([lake_geom, sub.geometry])
        lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'geometry'] = new_geom
    
    return lakes_snap_gdf

def fix_short_outlets(lakes_snap_gdf, combined_gdf, dem_res):
    half = dem_res / 2
    short_outlets = combined_gdf[combined_gdf['LakeOut'].notna() & (combined_gdf.geometry.length < dem_res * 0.95)]

    for lake_id, group in short_outlets.groupby('LakeOut'):
        lake_geom = lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'geometry'].iloc[0]
        for _, chan in group.iterrows():
            source, outlet = get_endpoints(chan.geometry)
            square    = box(outlet.x - half, outlet.y - half, outlet.x + half, outlet.y + half)
            lake_geom = unary_union([lake_geom, square])
        lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'geometry'] = lake_geom

    return lakes_snap_gdf

def fix_short_inlets(lakes_snap_gdf, combined_gdf, dem_res):
    short_outlets = combined_gdf[combined_gdf['LakeIn'].notna() & (combined_gdf.geometry.length < dem_res * 0.95)]

    for lake_id, group in short_outlets.groupby('LakeIn'):
        lake_geom = lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'geometry'].iloc[0]
        for _, chan in group.iterrows():
            source, outlet = get_endpoints(chan.geometry)
            missing = 1.05 * dem_res - chan.geometry.length
            half = missing / 2
            square = box(source.x - half, source.y - half, source.x + half, source.y + half)
            lake_geom = lake_geom.difference(square)

        lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'geometry'] = lake_geom

    return lakes_snap_gdf


def fix_inlet_subbasin_boundary(lakes_snap_gdf, combined_gdf, subs_original, dem_res):
    half = dem_res * 1.05
    tol  = dem_res * 0.1
    subs_boundary = subs_original.geometry.boundary.union_all()

    inlets  = combined_gdf[combined_gdf['LakeIn'].notna()]
    outlets = combined_gdf[combined_gdf['LakeOut'].notna()]

    for lake_id, group in inlets.groupby('LakeIn'):
        lake_geom = lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'geometry'].iloc[0]
        for _, chan in group.iterrows():
            source, outlet = get_endpoints(chan.geometry)
            if subs_boundary.distance(outlet) < tol:
                square    = box(outlet.x - half, outlet.y - half, outlet.x + half, outlet.y + half)
                lake_geom = unary_union([lake_geom, square])
        lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'geometry'] = lake_geom

    for lake_id, group in outlets.groupby('LakeOut'):
        lake_geom = lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'geometry'].iloc[0]
        for _, chan in group.iterrows():
            source, outlet = get_endpoints(chan.geometry)
            if subs_boundary.distance(source) < tol:
                square    = box(source.x - half, source.y - half, source.x + half, source.y + half)
                lake_geom = unary_union([lake_geom, square])
        lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'] == lake_id, 'geometry'] = lake_geom

    return lakes_snap_gdf

def check_cat3(lakes_snap_gdf, chans_gdf, dem_res, label):
    cat3 = classify_cat3(lakes_snap_gdf, chans_gdf, dem_res)
    flagged = cat3[cat3['flagged'] | cat3['flagged_short_internal']]
    if not flagged.empty:
        print(f"  ! cat3 flagged after {label}: {len(flagged)}")


def remove_holes(geometry):
    if geometry.geom_type == 'Polygon':
        return Polygon(geometry.exterior)
    elif geometry.geom_type == 'MultiPolygon':
        return MultiPolygon([Polygon(part.exterior) for part in geometry.geoms])
    else:
        return geometry
    
def rasterize_lakes(dem_fp, lakes_gdf, out_fp, attr="Hylak_id", all_touched=True, nodata=0, scale=1):
    with rasterio.open(dem_fp) as dem:
        profile   = dem.profile.copy()
        transform = dem.transform
        crs       = dem.crs
        height    = dem.height
        width     = dem.width

    fine_transform = transform * transform.scale(1/scale, 1/scale)
    fine_height = height * scale
    fine_width = width * scale

    gdf = lakes_gdf.copy()
    gdf['geometry'] = gdf['geometry'].apply(remove_holes)
    if gdf.crs != crs:
        gdf = gdf.to_crs(crs)
    if attr not in gdf.columns:
        raise ValueError(f"Attribute '{attr}' not found.")

    gdf = gdf[gdf[attr].notna() & gdf.geometry.notna() & (~gdf.geometry.is_empty)]
    vals = gdf[attr].astype(np.uint32)

    arr = np.full((fine_height, fine_width), nodata, dtype=np.uint32)
    shapes = zip(gdf.geometry, vals.astype(int))
    rasterize(
        shapes=shapes,
        out=arr,
        transform=fine_transform,
        all_touched=all_touched,
    )

    profile.update(dtype="uint32", nodata=nodata, count=1, compress="lzw",
                   transform=fine_transform, height=fine_height, width=fine_width)
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
        recs.append({"geometry": g, "Hylak_id": v})

    gdf = gpd.GeoDataFrame(recs, crs=crs)
    if dissolve and not gdf.empty:
        gdf = gdf.dissolve(by="Hylak_id", as_index=False)
    return gdf

if __name__ == '__main__':

    print("Integrating lakes/reservoirs into model network...")


    # change working directory
    os.chdir(os.path.dirname(__file__))

    # get model setup version
    parser = argparse.ArgumentParser(description="a terminal script for running the model setup and delineation.")

    parser.add_argument("r", help="the name of the region to run the model for. If not specified, all regions will be processed.", nargs='*', default=[])
    parser.add_argument("--v", help="the version of the model setup to use. If not specified, the datavariables value will be used.", nargs='?', default=None)

    args = parser.parse_args()

    if args.v is None: version = variables.version
    else: version = args.v  
    
    if len(args.r) > 0: regions = args.r
    else: regions = os.listdir(f"..data/model-setup/CoSWATv{version}/")


    details = {
    'auth': variables.final_proj_auth,
    'code': variables.final_proj_code,
    'coswat_data_version':variables.coswat_data_version,
    'region_source':variables.region_source}

    for region in regions:
        details['region']   = region
        details['version']  = version
        
        setup_dir  = f'../model-setup/CoSWATv{version}/{region}'
        data_dir  = f'../model-data/{variables.coswat_data_version}/{region}'
        
        out_dir    = f'../model-setup/CoSWATv{version}/{region}/Watershed/Shapes/integrateLakes'

        os.makedirs(out_dir, exist_ok=True)

        # Shapes
        lakes_path = 'shapes/lakes-grand-ESRI-54003.shp'
        lake_processed = 'shapes/lakes-grand-ESRI-54003-demAligned.shp'
        chans_path = 'Watershed/Shapes/dem-aster-ESRI-54003channel.shp'
        subs_path  = 'Watershed/Shapes/dem-aster-ESRI-54003subbasins.shp'
        dem_path   = 'Watershed/Rasters/DEM/dem-aster-ESRI-54003.tif'

        # PROCESS
        #===========================
        dem_res = variables.data_resolution

        lakes_gdf = gpd.read_file(f"{data_dir}/{lakes_path}")
        chans_gdf = gpd.read_file(f"{setup_dir}/{chans_path}")
        subs_gdf = gpd.read_file(f"{setup_dir}/{subs_path}")
        watershed_gdf = gpd.GeoDataFrame(geometry=[unary_union(subs_gdf.geometry)])

        subs_original = subs_gdf.copy()
        lakes_snap_gdf = lakes_gdf.copy()

        # Create simplified geom
        lakes_simp = lakes_gdf.copy()
        lakes_simp['geometry'] = lakes_simp.geometry.simplify(tolerance=dem_res*0.5, preserve_topology=True)

        out_rasterized = f'{setup_dir}/Watershed/Rasters/lake_out-ESRI-54003_rast.tif'

        rasterize_lakes(f'{setup_dir}/{dem_path}',lakes_simp,out_rasterized,scale=2)
        res_rast_gdf = polygonize_lakes(out_rasterized)
        res_rast_gdf['geometry'] = res_rast_gdf['geometry'].apply(remove_holes)

        # Align for snapping
        lakes_snap_gdf = res_rast_gdf.copy()
        lakes_snap_gdf = lakes_snap_gdf[lakes_snap_gdf.geometry.apply(lambda g: chans_gdf.geometry.intersects(g).any())]
        lakes_snap_gdf['geometry'] = lakes_snap_gdf.geometry.buffer(dem_res*0.05, cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)   
        lakes_snap_gdf['geometry'] = lakes_snap_gdf.geometry.apply(lambda g: max(g.geoms, key=lambda p: p.area) if g.geom_type == 'MultiPolygon' else g)

        lakes_snap_gdf.to_file(f"{data_dir}/shapes/lakes-grand-ESRI-54003_snap.gpkg")

        # Remove lakes that are not contained in watershed and do not intersect with any channel
        lakes_snap_gdf = lakes_snap_gdf[lakes_snap_gdf.geometry.apply(lambda g: watershed_gdf.geometry.iloc[0].contains(g))]
        lakes_snap_gdf = lakes_snap_gdf[lakes_snap_gdf.geometry.apply(lambda g: chans_gdf.geometry.intersects(g).any())]

        lakes_snap_gdf['LakeId'] = lakes_snap_gdf['Hylak_id'].copy()
        lakes_snap_gdf['RES'] = int(1)
        lakes_snap_gdf.to_file(f"{out_dir}/lakes_snap.gpkg")
        lakes_perim_gdf = lakes_snap_gdf.copy()
        lakes_perim_gdf['geometry'] = lakes_snap_gdf.geometry.boundary
        
        # Integrated approach
        max_iter    = 10
        iteration   = 0
        fixed       = False

        cat_0_fixed = False
        cat_1_fixed = False
        cat_2_fixed = False
        cat_3_fixed = False
        noDataSubs_fixed = False
        infiniteLooops_fixed = False


        while not fixed and iteration < max_iter:
            iteration += 1
            # REMOVE OVERLAPS
            geoms = lakes_snap_gdf.geometry.values
            idx = lakes_snap_gdf.index.values

            new_geoms = {i: g for i, g in zip(idx, geoms)}
            overlap_recs = []

            for i in range(len(geoms)):
                for j in range(i + 1, len(geoms)):
                    gi, gj = geoms[i], geoms[j]
                    if gi.intersects(gj):
                        overlap = gi.intersection(gj)
                        if not overlap.is_empty and overlap.area > 0:
                            overlap_recs.append({'geometry': overlap, 'Hylak_id': len(overlap_recs) + 1})

            if overlap_recs:
                overlaps_gdf = gpd.GeoDataFrame(overlap_recs, crs=lakes_snap_gdf.crs)
                out_rasterized = f'{out_dir}/lake_overlap-ESRI-54003_rast.tif'
                rasterize_lakes(f'{setup_dir}/{dem_path}', overlaps_gdf, out_rasterized, attr='Hylak_id', scale=2)
                overlap_rast_gdf = polygonize_lakes(out_rasterized)
                overlap_union = unary_union(overlap_rast_gdf.geometry).buffer(dem_res*0.05, cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)

                lakes_snap_gdf['geometry'] = lakes_snap_gdf.geometry.difference(overlap_union)

            print(f"  > iteration {iteration}...".ljust(50), end='\r')
            if iteration == 1: # Save initial flagged chans
                cat0_gdf = classify_cat0(lakes_snap_gdf,chans_gdf)
                cat1_gdf = classify_cat1(lakes_snap_gdf, chans_gdf, dem_res)
                cat2_gdf = classify_cat2(lakes_snap_gdf, chans_gdf, dem_res)
                cat3_gdf = classify_cat3(lakes_snap_gdf, chans_gdf, dem_res)

                cat0_gdf.to_file(f"{out_dir}/chans_cat0_init.gpkg")
                cat1_gdf.to_file(f"{out_dir}/chans_cat1_init.gpkg")
                cat2_gdf.to_file(f"{out_dir}/chans_cat2_init.gpkg")
                cat3_gdf.to_file(f"{out_dir}/chans_cat3_init.gpkg")

            # CATEGORY 0: Chans completely within
            print(f"  > Fixing category 0 (Within)...".ljust(50), end='\r')
            cat0_gdf       = classify_cat0(lakes_snap_gdf,chans_gdf) # Flag those that intersect the boundary
            lakes_snap_gdf = fix_cat0(lakes_snap_gdf,cat0_gdf,dem_res)
            cat0_fixed_gdf = classify_cat0(lakes_snap_gdf,chans_gdf)

            if cat0_fixed_gdf[cat0_fixed_gdf['flagged'] == True].empty:
                cat_0_fixed = True

            # CATEGORY 1: Inlets
            print(f"  > Fixing category 1 (In)...".ljust(50), end='\r')
            lakes_snap_gdf, cat1_fixed = runFix_cat1(lakes_snap_gdf, chans_gdf, dem_res)

            # CATEGORY 2: Outlets
            print(f"  > Fixing category 2 (Out)...".ljust(50), end='\r')
            lakes_snap_gdf, cat2_fixed = runFix_cat2(lakes_snap_gdf, chans_gdf, dem_res)

            # CATEGORY 3: Single chan
            print(f"  > Fixing category 3  (Single Channel)...".ljust(50), end='\r')
            lakes_snap_gdf, cat3_fixed = runFix_cat3(lakes_snap_gdf, chans_gdf, dem_res)

            # No data basin
            print(f"  > Fixing no data basins...".ljust(50), end='\r')
            flagged_subs_small   = checkSmallSubs(subs_gdf,lakes_snap_gdf)
            flagged_subs_nodata  = checkNoHRUSubs(subs_gdf,lakes_snap_gdf)
            lakes_snap_gdf       = fixNoDataBasin(flagged_subs_small,flagged_subs_nodata,lakes_snap_gdf)
            # INFINITE LOOPS
            print(f"  > Classifying channels ...".ljust(50), end='\r')
            combined_gdf = classify_channels(lakes_snap_gdf, chans_gdf, dem_res,tolerance = 50)
            print(f"  > Fixing infinite loops ...".ljust(50), end='\r')
            lakes_snap_gdf,loops_gdf = flag_infinite_loop_lakes(lakes_snap_gdf,combined_gdf)
            if iteration == 1:
                loops_gdf.to_file(f'{out_dir}/loops_initial.gpkg')
            lakes_snap_gdf,combined_gdf = runFix_inf(lakes_snap_gdf,combined_gdf,chans_gdf,dem_res)
            print(f"  > Fixing short outlets and other shenannigans...".ljust(50), end='\r')
            lakes_snap_gdf = fix_short_outlets(lakes_snap_gdf, combined_gdf, dem_res)
            
            lakes_snap_gdf = fix_inlet_subbasin_boundary(lakes_snap_gdf, combined_gdf, subs_original, dem_res)
            
            print(f"  > Carving on chains ...".ljust(50), end='\r')
            lakes_snap_gdf = carve_short_interlake_channels(lakes_snap_gdf,chans_gdf,dem_res)
            lakes_snap_gdf = fix_short_inlets(lakes_snap_gdf, combined_gdf, dem_res)
            
            # Check if they remain ok after fix
            cat0_gdf_final = classify_cat0(lakes_snap_gdf,chans_gdf)
            cat1_gdf_final = classify_cat1(lakes_snap_gdf, chans_gdf, dem_res)
            cat1_flagged_final   = cat1_gdf_final[cat1_gdf_final['flagged_crossing'] | cat1_gdf_final['flagged_short_inside'] | cat1_gdf_final['flagged_short_outside']]
            cat2_gdf_final = classify_cat2(lakes_snap_gdf, chans_gdf, dem_res)
            cat2_flagged_final = cat2_gdf_final[cat2_gdf_final['flagged_crossing'] | cat2_gdf_final['flagged_short_inside'] | cat2_gdf_final['flagged_short_outside']]
            cat3_gdf_final = classify_cat3(lakes_snap_gdf, chans_gdf, dem_res)
            cat3_flagged_final = cat3_gdf_final[cat3_gdf_final['flagged'] | cat3_gdf_final['flagged_short_internal']]
            lakes_snap_gdf_final,loop_gdf_final = flag_infinite_loop_lakes(lakes_snap_gdf,combined_gdf)
            lakes_snap_flagged_final   = lakes_snap_gdf_final[lakes_snap_gdf_final['infinite_loop']]
            flagged_subs_final_small   = checkSmallSubs(subs_gdf,lakes_snap_gdf)
            flagged_subs_nodata_final  = checkNoHRUSubs(subs_gdf,lakes_snap_gdf)

            if cat0_gdf_final[cat0_gdf_final['flagged'] == True].empty:
                cat_0_fixed = True
            if cat1_flagged_final.empty:
                cat_1_fixed = True
            if cat2_flagged_final.empty:
                cat_2_fixed = True
            if cat3_flagged_final.empty:
                cat_3_fixed = True
            if flagged_subs_final_small.empty and flagged_subs_nodata_final.empty:
                noDataSubs_fixed = True
            if lakes_snap_flagged_final.empty:
                infiniteLooops_fixed = True

            # Reclassify channels
            if cat_0_fixed and cat_1_fixed and cat_2_fixed and cat_3_fixed and noDataSubs_fixed and infiniteLooops_fixed and not overlap_recs:

                fixed = True

                lakes_snap_gdf['geometry'] = lakes_snap_gdf.geometry.buffer(-dem_res/8, cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)\
                    .buffer(dem_res/8, cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)
                lakes_snap_gdf['geometry'] = lakes_snap_gdf['geometry'].buffer(50, cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)\
                    .buffer(-50,cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)# Second buffer to avoid small annoying notches
                
                lakes_snap_gdf, _ = runFix_cat3(lakes_snap_gdf, chans_gdf, dem_res)
                lakes_snap_gdf['geometry'] = lakes_snap_gdf['geometry'].apply(remove_holes)

                lakes_snap_gdf = lakes_snap_gdf[lakes_snap_gdf.geometry.apply(lambda g: chans_gdf.geometry.intersects(g).any())]
                
                lakes_snap_gdf = lakes_snap_gdf.drop(columns=['infinite_loop'], errors='ignore')
                lakes_snap_gdf = lakes_snap_gdf.sort_values('LakeId').reset_index(drop=True)

                lake_data    = gpd.read_file(f'{data_dir}/{lake_processed}')[['Hylak_id', 'smax', 'pvol', 'evol', 'parea', 'earea', 'br1', 'br2']]
                lake_data    = lake_data[lake_data['Hylak_id'].isin(list(lakes_snap_gdf['Hylak_id']))]

                lakes_snap_gdf = pd.merge(lakes_snap_gdf,lake_data,on='Hylak_id')
                ordered_cols = ['LakeId', 'Hylak_id', 'smax', 'pvol', 'evol', 'parea', 'earea', 'br1', 'br2', 'RES']
                
                lakes_snap_gdf = lakes_snap_gdf[[c for c in ordered_cols if c in lakes_snap_gdf.columns] + ['geometry']]
                lakes_snap_gdf['RES'] = lakes_snap_gdf['RES'].fillna(1).astype(int)
                lakeFN_out     = "lakes-grand-{auth}-{code}.gpkg".format(**details)
                lakeFN_out_shp = "lakes-grand-{auth}-{code}.shp".format(**details)
                lakes_snap_gdf.to_file(f"{out_dir}/{lakeFN_out}")
                lakes_snap_gdf.to_file(f"{setup_dir}/Watershed/Shapes/{lakeFN_out_shp}")
                combined_gdf = classify_channels(lakes_snap_gdf, chans_gdf, dem_res, tolerance=50)
                print(f"Fixed after {iteration} attempts...")
                total_lakes  = sum(len(x) if isinstance(x, list) else 1 for x in lakes_snap_gdf['Hylak_id'])

            elif iteration ==  max_iter-1:
                print("Some lakes could not be integrated... they should be removed to integrate them. The run will proceed without them")                
                # Export flaggged that remains
                excluded_lakes = []
                if not cat0_gdf_final[cat0_gdf_final['flagged'] == True].empty:
                    cat0_gdf_final[cat0_gdf_final['flagged'] == True].to_file(f"{out_dir}/unresolved_cat_0.gpkg")
                
                    excluded_lakes.extend(cat0_gdf_final[cat0_gdf_final['flagged'] == True]['LakeId'].to_list())
                
                if not cat1_flagged_final.empty:
                    cat1_flagged_final.to_file(f"{out_dir}/unresolved_cat_1.gpkg")
                    excluded_lakes.extend(cat0_gdf_final[cat0_gdf_final['flagged'] == True]['LakeId'].to_list())

                if not cat2_flagged_final.empty:
                    cat2_flagged_final.to_file(f"{out_dir}/unresolved_cat_2.gpkg")
                    excluded_lakes.extend(cat2_flagged_final['LakeId'].to_list())

                if not cat3_flagged_final.empty:
                    cat3_flagged_final.to_file(f"{out_dir}/unresolved_cat_3.gpkg")
                    excluded_lakes.extend(cat3_flagged_final['LakeId'].to_list())

                if not lakes_snap_flagged_final.empty:
                    lakes_snap_flagged_final.to_file(f"{out_dir}/unresolved_infinite_loop.gpkg")
                    excluded_lakes.extend(lakes_snap_flagged_final['LakeId'].to_list())

                lakes_snap_gdf = lakes_snap_gdf[~lakes_snap_gdf['LakeId'].isin(excluded_lakes)]    

        lakes_perim_gdf = lakes_snap_gdf.copy()
        lakes_perim_gdf['geometry'] = lakes_snap_gdf.geometry.boundary

        print(f"  > A total of {total_lakes} lakes/reservoirs have been integrated")

        # SNAPPING POINTS
        snap_points = []
        print(f"\t Creating snapping points for QSWAT+ ...")
        for idx, row in combined_gdf.iterrows():
            lakeIN  = row['LakeIn']
            lakeOUT = row['LakeOut']
            start_point, end_point = get_endpoints(row.geometry)

            if lakeOUT is not None:
                snap_points.append({'geometry': start_point, 'RES': 1, 'LakeId': lakeOUT})

            if lakeIN is not None:
                snap_points.append({'geometry': end_point, 'RES': 0, 'LakeId': lakeIN})

        snap_gdf = gpd.GeoDataFrame(snap_points, geometry='geometry', crs=combined_gdf.crs)
        snap_gdf = snap_gdf.drop_duplicates(subset='geometry')
        snap_gdf.to_file(variables.resolved_snaps.format(**details))




# Cemetery
#=======================================================================================================================================================
# def merge_overlapping_lakes(lakes_snap_gdf):
#     visited = set()
#     merged_rows = []
#     new_lake_id = int(lakes_snap_gdf['LakeId'].max()) + 1

#     def get_cluster(idx, gdf):
#         cluster = set()
#         stack = [idx]
#         while stack:
#             current = stack.pop()
#             if current in cluster:
#                 continue
#             cluster.add(current)
#             overlaps = gdf.loc[current, 'overlap_with']
#             for lid in overlaps:
#                 next_idx = gdf[gdf['LakeId'] == lid].index
#                 if len(next_idx) > 0:
#                     stack.append(next_idx[0])
#         return cluster

#     for idx in lakes_snap_gdf.index:
#         if idx in visited:
#             continue
#         if not lakes_snap_gdf.loc[idx, 'overlap_flag']:
#             row = lakes_snap_gdf.loc[idx].to_dict()
#             row['merged_flag'] = False
#             row['merged_ids'] = []
#             merged_rows.append(row)
#             visited.add(idx)
#             continue

#         cluster_idx = get_cluster(idx, lakes_snap_gdf)
#         visited.update(cluster_idx)
#         group = lakes_snap_gdf.loc[list(cluster_idx)]

#         group = group.copy()
#         group['smax'] = pd.to_numeric(group['smax'], errors='coerce')
        
#         largest = group.loc[group['smax'].idxmax()]

#         def nansum(col):   return group[col].replace('', np.nan).apply(pd.to_numeric, errors='coerce').sum(skipna=True)
#         def nanavg(col):   return group[col].replace('', np.nan).apply(pd.to_numeric, errors='coerce').mean(skipna=True)
#         def tolist(col):   return group[col].dropna().tolist()
#         def fromlargest(col): return largest[col]

#         new_row = {
#             'geometry':     unary_union(group.geometry.values),
#             'LakeId':       new_lake_id,
#             'merged_flag':  True,
#             'merged_ids':   group['LakeId'].tolist(),
#             'Hylak_id':     tolist('Hylak_id'),
#             'smax':         nansum('smax'),
#             'pvol':         nansum('pvol'),
#             'evol':         nansum('evol'),
#             'parea':        nansum('parea'),
#             'earea':        nansum('earea'),
#             'br1':          fromlargest('br1'),
#             'br2':          fromlargest('br2'),
#             'RES':          fromlargest('RES'),
#             'calcAreas':    nansum('calcAreas'),
#             'calcVol':      nansum('calcVol'),
#             'elev_masl':    nanavg('elev_masl'),
#             'dis_avg':      fromlargest('dis_avg'),
#             'Lake_name':    np.nan,
#             'Country':      fromlargest('Country'),
#             'Continent':    fromlargest('Continent'),
#             'Lake_type':    fromlargest('Lake_type'),
#             'Depth_avg':    fromlargest('Depth_avg'),
#             'Res_time':     fromlargest('Res_time'),
#             'Grand_id':     tolist('Grand_id'),
#             'RES_NAME':     tolist('RES_NAME'),
#             'DAM_NAME':     tolist('DAM_NAME'),
#             'RIVER':        tolist('RIVER'),
#             'MAIN_BASIN':   tolist('MAIN_BASIN'),
#             'SUB_BASIN':    tolist('SUB_BASIN'),
#             'YEAR':         tolist('YEAR'),
#             'REM_YEAR':     tolist('REM_YEAR'),
#             'DAM_HGT_M':    tolist('DAM_HGT_M'),
#             'DAM_LEN_M':    tolist('DAM_LEN_M'),
#             'DEPTH_M':      fromlargest('DEPTH_M'),
#             'DOR_PC':       fromlargest('DOR_PC'),
#             'MAIN_USE':     fromlargest('MAIN_USE'),
#             'LAKE_CTRL':    fromlargest('LAKE_CTRL'),
#             'TIMELINE':     tolist('TIMELINE'),
#             'a':            fromlargest('a'),
#             'b':            fromlargest('b'),
#             'c':            fromlargest('c'),
#             'd':            fromlargest('d'),
#             'maxDepth':     fromlargest('maxDepth'),
#             'meanDepth':    fromlargest('meanDepth'),
#             'MaxArea_km':   fromlargest('MaxArea_km'),
#         }
#         merged_rows.append(new_row)
#         new_lake_id += 1

#     result_gdf = gpd.GeoDataFrame(merged_rows, crs=lakes_snap_gdf.crs)
#     result_gdf = result_gdf.drop(columns=['overlap_flag', 'overlap_with'], errors='ignore')
#     return result_gdf

            # # OVERLAYING Lakes
            # print(f"  > Checking overlaps ...".ljust(50), end='\r')
            # tree = STRtree(lakes_snap_gdf.geometry.values)
            # lakes_snap_gdf['overlap_flag'] = False
            # lakes_snap_gdf['overlap_with'] = [[] for _ in range(len(lakes_snap_gdf))]
            # for idx, lake in lakes_snap_gdf.iterrows():
            #     candidates_idx = tree.query(lake.geometry)
            #     overlapping = []
            #     for cidx in candidates_idx:
            #         other = lakes_snap_gdf.iloc[cidx]
            #         if other.name == idx:
            #             continue
            #         if lake.geometry.overlaps(other.geometry) or lake.geometry.contains(other.geometry) or other.geometry.contains(lake.geometry) or lake.geometry.touches(other.geometry):
            #             overlapping.append(other['LakeId'])
            #     if overlapping:
            #         lakes_snap_gdf.at[idx, 'overlap_flag'] = True
            #         lakes_snap_gdf.at[idx, 'overlap_with'] = overlapping
            
            # previously_merged = set(lakes_snap_gdf[lakes_snap_gdf['merged_flag'] == True]['LakeId'].tolist())
            # lakes_snap_gdf = merge_overlapping_lakes(lakes_snap_gdf)
            # lakes_snap_gdf.loc[lakes_snap_gdf['LakeId'].isin(previously_merged), 'merged_flag'] = True

        # # BUFFER: Not perfectly align to the DEM to avoid weird line intersection, smooth out to remove protrusions
        # lakes_snap_gdf['geometry'] = lakes_snap_gdf.geometry.buffer(-dem_res/2 * 0.95, cap_style=BufferCapStyle.square, join_style=BufferJoinStyle.mitre)\
        # .buffer(-dem_res * 0.5, cap_style="square", join_style="mitre")\
        # .buffer(dem_res * 0.5, cap_style="square", join_style="mitre")