'''
This script identifies irrigated HRUS and their source (reservoir/river/groundwater) and generates decision tables:

    - Based on FAO Global Irrigated Area maps
    - Applying Irrigation Topology logic from Vanderkelen et al. (2022)


========================================================================================================================
    ! This will only work if using swatplus-61.0.2.11-31-coswatn5-ifx-lin_x86_64 executable
========================================================================================================================

Author  : Jose Pablo teran
Date    : May 2025
Contact : jose.pablo.teran.orsini@vub.be
GitHub  : github.com/celray - github.com/jopator
'''

import pandas as pd
import geopandas as gpd
from read_swat import *
import matplotlib.pyplot as plt
import os
import numpy as np
import rioxarray as rxr
from cjfx import list_folders, goto_dir,exists
from rasterio.features import shapes
from shapely.geometry import shape
import datavariables as variables
import ast
import argparse
import sys
import os

def cells_polys_above(ds, var='pctg', thr=50):
    da = ds[var].squeeze()
    arr = da.values
    mask = (~np.isnan(arr)) & (arr > thr)
    geoms = []
    for geom, val in shapes(arr.astype('float32'), mask=mask, transform=da.rio.transform()):
        geoms.append(shape(geom))
    return gpd.GeoDataFrame(geometry=geoms, crs=da.rio.crs)

def assign_irrig_flag(lsus_gdf, poly_gdf, out_col, landuse_col='Landuse', allowed=('AGRL','AGRR')):
    gdf = lsus_gdf.copy()
    gdf[out_col] = 0
    # only AGRL/AGRR can be 1
    mask_ag = gdf[landuse_col].isin(allowed)
    if not poly_gdf.empty:
        lsu_ag = gdf.loc[mask_ag].to_crs(poly_gdf.crs)
        base = lsu_ag.reset_index().rename(columns={'index':'_idx'})[['_idx','geometry']]
        # intersection check
        hit = gpd.sjoin(poly_gdf, base, how='inner', predicate='intersects')['_idx'].unique()
        gdf.loc[hit, out_col] = 1
    return gdf


if __name__ == '__main__':

    # change working directory
    goto_dir(__file__)

    # get model setup version
    parser = argparse.ArgumentParser(description="a terminal script for running the model setup and delineation.")

    parser.add_argument("r", help="the name of the region to run the model for. If not specified, all regions will be processed.", nargs='*', default=[])
    parser.add_argument("--v", help="the version of the model setup to use. If not specified, the datavariables value will be used.", nargs='?', default=None)

    args = parser.parse_args()

    # get model setup version
    if args.v is None: version = variables.version
    else: version = args.v  

    # get regions
    if len(args.r) > 0: regions = args.r
    else: regions = list_folders(f"../model-setup/CoSWATv{version}/")

    if not exists(f"../model-setup/CoSWATv{version}"):
        print(f'\t! the version, CoSWATv{version}, does not exist, the following versions are available:')
        for v in list_folders('../model-setup/'):
            if v.startswith('CoSWATv'):
                print(f'\t\t- {v}')
        print(f'\t> please specify a valid version using the --v argument')
        sys.exit(1)



    details = {
    'auth': variables.final_proj_auth,
    'code': variables.final_proj_code,
    }


    for region in regions:
        details['region']   = region
        details['version']  = version


        model_dir  = f'../model-setup/CoSWATv{version}/{region}'
        
        # Shapes
        lakes_path = 'Watershed/Shapes/lakes_coswat-ESRI-54003.gpkg'
        chans_path = 'Watershed/Shapes/dem-aster-ESRI-54003channel.shp'
        rivs1_path = 'Watershed/Shapes/rivs1.shp'
        lsus_path  = 'Watershed/Shapes/lsus2.shp'
        hrus_path  = 'Watershed/Shapes/hrus2.shp'

        # FAO irrigation maps
        surface_water_fn = '../data-preparation/resources/area_irrigated_surface_water/gmia_v5_aeisw_pct_aei.asc'
        ground_water_fn  = '../data-preparation/resources/area_irrigated_groundwater/gmia_v5_aeigw_pct_aei.asc'

        # Threshold for Topology
        main_threshold              = variables.main_threshold
        tributary_order             = variables.tributary_order
        fao_irrg_area_pctg_thres    = variables.fao_irrg_area_pctg_thres

        # Model files
        hru_data_fn         = f'{model_dir}/Scenarios/Default/TxtInOut/hru-data.hru'
        landuse_lum_fn      = f'{model_dir}/Scenarios/Default/TxtInOut/landuse.lum'
        lum_dtl_fn          = f'{model_dir}/Scenarios/Default/TxtInOut/lum.dtl'
        management_sch_fn   = f'{model_dir}/Scenarios/Default/TxtInOut/management.sch'


        # Read and process
        chan_gdf    = gpd.read_file(f'{model_dir}/{chans_path}')
        lsus_gdf    = gpd.read_file(f'{model_dir}/{lsus_path}')
        rivs1_gdf   = gpd.read_file(f'{model_dir}/{rivs1_path}')
        hrus_gdf    = gpd.read_file(f'{model_dir}/{hrus_path}')

        lakeFN = f'{model_dir}/{lakes_path}'
        
        if os.path.isfile(lakeFN):
            lakes_gdf   = gpd.read_file(lakeFN)
            noLakes = False
        
        else:
            noLakes = True
            
        sfw_irr_da  = rxr.open_rasterio(surface_water_fn).squeeze("band",drop=True).rio.write_crs("EPSG:4326")
        grw_irr_da  = rxr.open_rasterio(ground_water_fn).squeeze("band",drop=True).rio.write_crs("EPSG:4326")

        # bounds
        wshed_gdf   = lsus_gdf.dissolve()
        bounds      = wshed_gdf.to_crs('epsg:4326').bounds

        sfw_irr_da  = sfw_irr_da.rio.clip_box(bounds['minx'].iloc[0],bounds['miny'].iloc[0],bounds['maxx'].iloc[0],bounds['maxy'].iloc[0])
        grw_irr_da  = grw_irr_da.rio.clip_box(bounds['minx'].iloc[0],bounds['miny'].iloc[0],bounds['maxx'].iloc[0],bounds['maxy'].iloc[0])

        # To dataset
        sfw_irr_ds  = sfw_irr_da.to_dataset(name='pctg').rename({'x':'lon','y':'lat'})
        grw_irr_ds  = grw_irr_da.to_dataset(name='pctg').rename({'x':'lon','y':'lat'})
        
        
        # Check
        landuse_lum    = swatModelFile(landuse_lum_fn)
        landuse_lum_df = landuse_lum.dframe

        # add new rows
        if landuse_lum_df[landuse_lum_df['name']=='agrl_lum'].empty:
            print("No agricultural HRUs in this region. Skipping...")
            quit()


        """
        ========================================================================================================================
        Detection of irrigated HRUS
        ========================================================================================================================
        """

        # Intersect LSUS to the irrigated areas
        hrus_gdf['irr_sfw'] = 0
        hrus_gdf['irr_grw'] = 0

        # create polygons from raster where pctg > fao_irrg_area_pctg_thres (50 by default)
        sfw_polys = cells_polys_above(sfw_irr_ds, var='pctg', thr=fao_irrg_area_pctg_thres)
        grw_polys = cells_polys_above(grw_irr_ds, var='pctg', thr=fao_irrg_area_pctg_thres)

        # assign binary flags
        hrus_gdf = assign_irrig_flag(hrus_gdf, sfw_polys, out_col='irr_sfw', landuse_col='Landuse')
        hrus_gdf = assign_irrig_flag(hrus_gdf, grw_polys, out_col='irr_grw', landuse_col='Landuse')

        hrus_gdf.loc[hrus_gdf['irr_sfw'] == 1, 'irr_grw'] = 0                   # If surface, then no GW (in case there was a weird intersection)
        irr_sfw_hrus_gdf = hrus_gdf[hrus_gdf['irr_sfw']==1]
        irr_sfw_hrus_gdf_orig = irr_sfw_hrus_gdf.copy()




        """
        ========================================================================================================================
        Reservoir Irrigation Topology
        ========================================================================================================================
        """

        links_dict = {}

        order = tributary_order
        
        if noLakes:
            print("This region does not have reservoirs...")
            links_dict = {}
        
        else:
            

            for idx, row in lakes_gdf.iterrows():
                smax_k          = float(row["smax"])/1000000
                if row['merged_fla']:
                    merged_hylaks = row['Hylak_id']
                    ratios        = row['ratios']
                    main_use_all  = row['MAIN_USE']
                    smax_biggest  = 0
                    
                    frac_inactive_list = []

                    for hylak_id, ratio,mainUse in zip(ast.literal_eval(merged_hylaks),
                                                        ast.literal_eval(ratios),
                                                        ast.literal_eval(main_use_all)
                        ):
                        ratio,mainUse = float(ratio), str(mainUse)
                        smax_single = smax_k*ratio
                        if smax_single > smax_biggest:
                            mainUse_merged  = mainUse
                            
                    main_use = mainUse_merged
                else:
                    
                    main_use    = str(row['MAIN_USE'])

                if main_use == 'Irrigation':
                    
                    lake_id     = row['LakeId']
                    lake_gdf    = lakes_gdf[lakes_gdf['LakeId']==lake_id]                                         # get lakes
                    outlet_gdf  = chan_gdf[chan_gdf['LakeMain']==lake_id]                                         # get main outlet
                    outlet_gdf  = outlet_gdf.reset_index(drop=True)
                    complete    = False
                    links = []
                    outletLINKNO    = outlet_gdf['LINKNO'].iloc[0]
                    outletDSLINKNO  = outlet_gdf['DSLINKNO'].iloc[0]
                    outletLENGTH    = outlet_gdf['Length'].iloc[0]

                    total_length = outletLENGTH
                    links.append(outletLINKNO)

                    currentLINKNO = outletDSLINKNO
                    
                    if currentLINKNO == -1: # Just one linkno downstream
                        continue

                    while not complete:
                        new_chan_gdf    = chan_gdf[chan_gdf['LINKNO'] == currentLINKNO].reset_index(drop=True)
                        new_chan_lenght = new_chan_gdf['Length'].iloc[0]         
                        total_length += new_chan_lenght

                        currentUSLINKNO1 = new_chan_gdf['USLINKNO1'].iloc[0]
                        currentUSLINKNO2 = new_chan_gdf['USLINKNO2'].iloc[0]

                        links.append(currentLINKNO)

                        if currentUSLINKNO1 not in links:
                            links.append(currentUSLINKNO1)
                            
                            if order > 1 and currentUSLINKNO1 > 0:
                                extra_chan_gdf = chan_gdf[chan_gdf['LINKNO'] == currentUSLINKNO1].reset_index(drop=True)
                                
                                extraUSLINKNO1 = extra_chan_gdf['USLINKNO1'].iloc[0]
                                extraUSLINKNO2 = extra_chan_gdf['USLINKNO2'].iloc[0]
                                
                                if extraUSLINKNO1 > 0:    
                                    links.append(extraUSLINKNO1)

                                if extraUSLINKNO2 >0:
                                    links.append(extraUSLINKNO2)

                        elif currentUSLINKNO2 not in links:
                            links.append(currentUSLINKNO2)

                            if order > 1 and currentUSLINKNO2 > 0:
                                extra_chan_gdf = chan_gdf[chan_gdf['LINKNO'] == currentUSLINKNO2].reset_index(drop=True)

                                extraUSLINKNO1 = extra_chan_gdf['USLINKNO1'].iloc[0]
                                extraUSLINKNO2 = extra_chan_gdf['USLINKNO2'].iloc[0]
                                
                                if extraUSLINKNO1 > 0:    
                                    links.append(extraUSLINKNO1)

                                if extraUSLINKNO2 >0:
                                    links.append(extraUSLINKNO2)


                        currentLINKNO = new_chan_gdf['DSLINKNO'].iloc[0]
                        lakeInFlag    = new_chan_gdf['LakeIn'].iloc[0]

                        if total_length > main_threshold*1000:
                            break

                        if currentLINKNO == -1:
                            break

                        if lakeInFlag != 0:
                            break
                    
                    links_dict[lake_id] = links                 # For each irrigation reservoir, we get the links from the topology logic




        # Transform links dict to rivs dict
        rivs_dict = {}

        for reservoir,links in links_dict.items():
            irrig_chans_gdf      = rivs1_gdf[rivs1_gdf['LINKNO'].isin(links)]
            irrig_chans          = irrig_chans_gdf['Channel'].to_list()
            rivs_dict[reservoir] = irrig_chans

        
        # Adjust LSU shapefile
        # We will add two columns: number of reservoirs it is associated to, and the list of reservoirs
        lsus_gdf['ir_res_nr']   = 0
        lsus_gdf['ir_res']      = None

        for reservoir, channels in rivs_dict.items():
            lsus_subset = lsus_gdf.loc[lsus_gdf['Channel'].isin(channels)].copy()
            idx = lsus_subset.index

            lsus_gdf.loc[idx, 'ir_res_nr'] = lsus_gdf.loc[idx, 'ir_res_nr'] + 1

            curr = lsus_gdf.loc[idx, 'ir_res'].fillna('')
            lsus_gdf.loc[idx, 'ir_res'] = curr.apply(lambda s: reservoir if s == '' else f'{s},{reservoir}')
            
            
        # Export irrigation topology
        lsus_gdf.to_file(f'{model_dir}/Watershed/Shapes/Irrig_topo_lsus.gpkg')
            
        

        """
        ========================================================================================================================
        Classify by source for Irrigation
        ========================================================================================================================
        """

        # Spatial join and identification
        intersect_lsus_gdf = lsus_gdf[['ir_res_nr','ir_res','geometry']]
        intersect_lsus_gdf = intersect_lsus_gdf[intersect_lsus_gdf['ir_res_nr']>0].reset_index(drop=True)

        irr_sfw_hrus_gdf = irr_sfw_hrus_gdf_orig.reset_index(drop=True).sjoin(intersect_lsus_gdf, how='left',predicate='intersects').drop('index_right',axis=1).drop_duplicates()

        # Separate river source and reservoir source
        irr_sfw_hrus_gdf['res_source'] = None
        irr_sfw_hrus_gdf['riv_source'] = None

        # Surface Water: Either Reservoir source or River Source
        irr_sfw_hrus_gdf.loc[irr_sfw_hrus_gdf['ir_res_nr'] >= 1,   'res_source'] = 1 
        irr_sfw_hrus_gdf.loc[irr_sfw_hrus_gdf['ir_res_nr'].isna(), 'riv_source'] = 1

        # Groundwater 
        irr_grw_hrus_gdf = hrus_gdf[(hrus_gdf['irr_grw']==1) & (~hrus_gdf['HRUS'].isin(irr_sfw_hrus_gdf['HRUS']))]


        irr_sfw_hrus_gdf.to_file(f'{model_dir}/Watershed/Shapes/Irrig_topo_hrus_sfw.gpkg')
        irr_grw_hrus_gdf.to_file(f'{model_dir}/Watershed/Shapes/Irrig_topo_hrus_grw.gpkg')
        
        
        '''
        Ok so now we've got:
            - irr_sfw_hrus_gdf where we know reservoir or riv source
            - irr_grw_hrus_gdf where we know ground water

            * For reservoir source, we need to apply irrigation on the dtl and take it from the reservoir on the attributes
            but before that, if it has more than one reservoir, the irrigation demand is weighted based on the reservoir max capacity;

                --> Total capacity for irrigated HRUS: sum of reservoir capacities;
                    weight for each reservoir: reservoir capacity / sum of reservoir capacities      :)

            * For the river just take from the channel that corresponds to the HRU, however we need to see if there is no more water in the river,
            no irrigation can be applied !!

            * For groundwater it will just be taken from the deep aquifer ("infinite" source in the model)

            
            * Later on, we can benefit from the water allocation module for a more realistic approach

        '''

        '''
        1 - Create new land use mgt type per HRU agrl_lum_<HRU> and adjust hru-data.hru
        
        '''


        irrigated_hrus = pd.concat([irr_sfw_hrus_gdf[irr_sfw_hrus_gdf['irr_sfw']==1][['HRUS']],irr_grw_hrus_gdf[['HRUS']]]).reset_index(drop=True)
        irrigated_hrus['lu_mgt'] = 'agrl_lum_' + irrigated_hrus['HRUS'].astype(str)
        irrigated_hrus = irrigated_hrus.rename(columns={'HRUS':'id'}).drop_duplicates().reset_index(drop=True)
        hru_data = swatModelFile(hru_data_fn)
        hru_data_df = hru_data.dframe
        hru_data_df['id']    = hru_data_df['id'].astype(str)
        irrigated_hrus['id'] = irrigated_hrus['id'].astype(str)
        hru_data_df.set_index('id', inplace=True)
        irrigated_hrus.set_index('id', inplace=True)
        hru_data_df.update(irrigated_hrus[['lu_mgt']])
        hru_data_df.reset_index(inplace=True)
        hru_data_df = hru_data_df.fillna('null')
        hru_data.dframe = hru_data_df.copy()

        hru_data.write(name='hru-data.hru updated by CoSWAT+ framework after write-irrigation.py',overwrite=True)


        '''
        2 - Create land use classes in landuse.lum

        reference: copy from agrl lum, just update name and management schedule name (agrl_<HRU>)

        '''
        irrigated_hrus = irrigated_hrus.drop_duplicates()
        irrigated_hrus.reset_index(inplace=True)
        irrigated_hrus['mgt'] = 'agrl_' + irrigated_hrus['id'].astype(str)
        irrigated_hrus = irrigated_hrus.rename(columns={'lu_mgt':'name'})

        landuse_lum    = swatModelFile(landuse_lum_fn)
        landuse_lum_df = landuse_lum.dframe
        
        tmpl_row = (landuse_lum_df.loc[landuse_lum_df['name'] == 'agrl_lum'].iloc[0].to_dict())
        n = len(irrigated_hrus)
        new_rows = pd.DataFrame([tmpl_row] * n)

        new_rows['name'] = irrigated_hrus['name'].to_numpy()
        new_rows['mgt']  = irrigated_hrus['mgt'].to_numpy()

        landuse_lum_df = pd.concat([landuse_lum_df, new_rows], ignore_index=True)
        landuse_lum_df = landuse_lum_df.drop_duplicates(subset=['name'], keep='first').reset_index(drop=True)

        landuse_lum_df = landuse_lum_df.fillna('null')
        landuse_lum.dframe = landuse_lum_df.copy()

        landuse_lum.write(name='landuse.lum updated by CoSWAT+ framework after define-irrigation.py',overwrite=True)



        '''
        3 -  Create land use tables for irrigation
        '''


        def wstress_Conds(wstress):
            # Condition Dataframe --> This is different for each reservoir
            cond_df = pd.DataFrame(columns=['var','obj','obj_num','lim_var','lim_op','lim_const','alt1'])
            cond_df["var"] = ["w_stress"]
            cond_df["obj"] = ["hru"]
            cond_df["obj_num"] = [int(x) for x in [0]]
            cond_df["lim_var"] = ["null"]
            cond_df["lim_op"] = ["-"]
            cond_df["lim_const"] = [float(x) for x in [wstress]]
            cond_df["alt1"] = ["<"]
            return cond_df
        
        def res_irrig_actions(amount,times,channel,nr_of_res,res_nr_list:list,res_frac_list:list):

            amount = float(amount)          
            times  = int(times)
            amount_per_res = (np.array(res_frac_list, dtype=float) * amount).tolist()

            data = {
                'act_typ': ['irrigate'] + ['res_irr_dmd'] * nr_of_res,
                'obj':     ['cha']      + ['res'] * nr_of_res,
                'obj_num': [int(channel)] + [int(x) for x in res_nr_list],
                'name':    ['drip_high'] + ['record_irr'] * nr_of_res,
                'option':  ['drip']      + ['null'] * nr_of_res,
                'const':   [amount]      + amount_per_res,
                'const2':  [times]       + [0] * nr_of_res,
                'fp':      ['null']      + ['null'] * nr_of_res,
                'alt1':    ['y']         + ['y'] * nr_of_res,
            }
            return pd.DataFrame(data, columns=['act_typ','obj','obj_num','name','option','const','const2','fp','alt1'])
        

        def river_irrig_actions(amount,times,channel):

            amount = float(amount)          
            times  = int(times)

            data = {
                'act_typ': ['irrigate'],
                'obj':     ['cha'],
                'obj_num': [int(channel)],
                'name':    ['drip_high'],
                'option':  ['drip'],
                'const':   [amount],
                'const2':  [times],
                'fp':      ['null'],
                'alt1':    ['y']
            }
            return pd.DataFrame(data, columns=['act_typ','obj','obj_num','name','option','const','const2','fp','alt1'])
        
        def groundW_irrig_actions(amount,times,channel):

            amount = float(amount)          
            times  = int(times)

            data = {
                'act_typ': ['irrigate'],
                'obj':     ['aqu'],
                'obj_num': [int(channel)],
                'name':    ['drip_high'],
                'option':  ['drip'],
                'const':   [amount],
                'const2':  [times],
                'fp':      ['null'],
                'alt1':    ['y']
            }
            return pd.DataFrame(data, columns=['act_typ','obj','obj_num','name','option','const','const2','fp','alt1'])
        
        """
        ========================================================================================================================
        Decision table structure for HRUS with reservoir source

        * Water will be taken from channel, but irrigation demand will be recorded to condition release (Hanazaki et al., 2006)
        * Irrigate 25 mm with sprinkler if water stress < 0.8

        it will be a 'unique' dtl per channel/reservoir(s) combination based on irr_sfw_hrus_gdf where res_source is 1


        name                   conds      alts      acts     
        coswat_hru_res_irrig       1         1         2
        var                      obj   obj_num           lim_var            lim_op     lim_const      alt1
        w_stress                 hru         0              null                 -           0.9         <                    !water stress < this value
        act_typ                  obj   obj_num              name            option         const        const2                fp  outcome
        irrigate                 cha         x         drip_high              drip            80            50              null  y
        res_irr_dmd              res         x        record_irr              null            80             0              null  y
        ========================================================================================================================
        """

        res_irr_hrus = irr_sfw_hrus_gdf[irr_sfw_hrus_gdf['res_source']==1][['Channel','HRUS','ir_res_nr','ir_res']]
        res_irr_hrus = res_irr_hrus.drop_duplicates().reset_index(drop=True)

        for index, row in res_irr_hrus.iterrows():
            chan = row['Channel']
            hru  = row['HRUS']
            nr_of_res = row['ir_res_nr']
            res_id    = row['ir_res']

            dtl_name = 'agrl_dtl_'+str(hru)
            # convert reservoir list string to a list
            res_list_str = str(res_id).split(',')
            res_list_gis = [int(float(x)) for x in res_list_str]

            # Map LakeId -> reservoir id (index + 1)
            lakes_gdf = lakes_gdf.sort_values(by='LakeId').reset_index(drop=True)
            res_list = [lakes_gdf.index[lakes_gdf['LakeId'] == gis_id].tolist()[0] + 1 for gis_id in res_list_gis]
            # Get fractions
            source_res_gdf = lakes_gdf[lakes_gdf['LakeId'].isin(res_list_gis)]

            fraction_list = []

            total_vol = source_res_gdf['pvol'].astype(float).sum()

            for idx,reservoir in source_res_gdf.iterrows():
                pvol = float(reservoir['pvol'])
                fract = pvol/total_vol
                fraction_list.append(round(fract,5))

            lum_dtl = swat_Dtl(lum_dtl_fn)                               # read lum.dtl

            actions = 1+int(nr_of_res)

            conds_df = wstress_Conds(0.9)
            acts_df  = res_irrig_actions(80,50,int(chan),int(nr_of_res),res_list,fraction_list)

            lum_dtl.add_dtl(dtl_name,1,1,actions,conds_df,acts_df,overwrite=True,replace=True)


        print('\t\t > Irrigation decision tables have been created for reservoir sourced HRUS')

        """
        ========================================================================================================================
        Decision table structure for HRUS with river source

        * By default, if there is not enough water in the river, irrigation wont be applied, otherwise, it is applied just like when the source is reservoir
        * Irrigate 25 mm with sprinkler if water stress < 0.8

        it will be a 'unique' dtl per channel based on irr_sfw_hrus_gdf where riv_source is 1

        
        name                   conds      alts      acts     
        coswat_hru_rriv_irrig      2         1         2
        var                      obj   obj_num           lim_var            lim_op     lim_const      alt1
        w_stress                 hru         0              null                 -           0.8         <                    !water stress < this value
        act_typ                  obj   obj_num              name            option         const        const2                fp  outcome
        irrigate                 cha         x         sprinkler     sprinkler_ilm            25            20              null  y

        ========================================================================================================================
        """

        riv_irr_hrus = irr_sfw_hrus_gdf[irr_sfw_hrus_gdf['riv_source']==1][['Channel','HRUS']]
        riv_irr_hrus = riv_irr_hrus.drop_duplicates().reset_index(drop=True)
        for index, row in riv_irr_hrus.iterrows():
            chan = row['Channel']
            hru  = row['HRUS']


            dtl_name = 'agrl_dtl_'+str(hru)

            lum_dtl = swat_Dtl(lum_dtl_fn)                               # read lum.dtl


            conds_df = wstress_Conds(0.9)
            acts_df  = river_irrig_actions(80,50,int(chan))

            lum_dtl.add_dtl(dtl_name,1,1,1,conds_df,acts_df,overwrite=True,replace=True)


        print('\t\t > Irrigation decision tables have been created for river sourced HRUS')


        """
        ========================================================================================================================
        Decision table structure for HRUS with groundwater

        * By default, if there is not enough water in the aquifer, irrigation wont be applied, otherwise, it is applied just like when the source is reservoir or channel
        * Irrigate 25 mm with sprinkler if water stress < 0.8

        it will be a unique dtl per channel (will be matched with aquifer of LSU) for all hrus with groundwater source in irr_grw_hrus_gdf

        
        name                   conds      alts      acts     
        coswat_hru_rriv_irrig      2         1         2
        var                      obj   obj_num           lim_var            lim_op     lim_const      alt1
        w_stress                 hru         0              null                 -           0.8         <                    !water stress < this value
        act_typ                  obj   obj_num              name            option         const        const2                fp  outcome
        irrigate                 aqu         x         sprinkler     sprinkler_ilm            25            20              null  y

        ========================================================================================================================
        """

        grw_irr_hrus = irr_grw_hrus_gdf[['Channel','HRUS']]
        grw_irr_hrus = irr_grw_hrus_gdf.drop_duplicates().reset_index(drop=True)

        for index, row in grw_irr_hrus.iterrows():
            chan = row['Channel']
            hru  = row['HRUS']


            dtl_name = 'agrl_dtl_'+str(hru)

            lum_dtl = swat_Dtl(lum_dtl_fn)                               # read lum.dtl


            conds_df = wstress_Conds(0.9)
            acts_df  = groundW_irrig_actions(80,50,int(chan))

            lum_dtl.add_dtl(dtl_name,1,1,1,conds_df,acts_df,overwrite=True,replace=True)


        print('\t\t > Irrigation decision tables have been created for groundwater sourced HRUS')


        '''
        4 - Update management.sch file: for each agrl_<HRU>, an auto action is created

        '''
        def build_management_blocks(landuse_lum_df):
            """
            From landuse_lum_df, create text blocks like:
            <mgt>  0  1
                    <dtl_name>
            """
            # keep rows that have both values
            df = landuse_lum_df[['mgt', 'dtl_name']].dropna(subset=['mgt', 'dtl_name'])

            # crude but safe spacing; tweak widths if your editor expects exact columns
            name_w    = 34   # width for 'name' column
            indent_op = 51   # indent so second line starts under 'op_typ' column in your file

            lines = []
            for mgt, dtl in df.itertuples(index=False):
                lines.append(f"{str(mgt):<{name_w}}  0          2")
                lines.append(f"{'':<{indent_op}}{str(dtl)}")
                lines.append(f"{'':<{indent_op}}pl_hv_summer1   agrl")
            return "\n".join(lines) + "\n"

        irrigated_hrus['dtl_name'] = 'agrl_dtl_' + irrigated_hrus['id'].astype(str)
        text = build_management_blocks(irrigated_hrus)
        
        with open(management_sch_fn, "a", encoding="utf-8") as f:
            f.write(text)