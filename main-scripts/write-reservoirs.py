'''
This script modifies the model files to:

    - Correctly define area and storage properties of the lakes based on global datasets: HydroLakes, GranD, GLOBAthy
    - Define release decision tables for lakes and reservoirs.


========================================================================================================================
Schemes used for release decision tables:

    * Unregulated lakes: Doell et al. (2003)
    * Water Supply, flood control and irrigation reservoirs: Hanazaki et al. (2006) / Vanderkelen et al. (2022) / Gharari et al. (2024)

    > Main use is determined by global datasets - If not on GranD , assumed as natural lake.
========================================================================================================================

Author  : Jose Pablo teran
Date    : May 2025
Contact : jose.pablo.teran.orsini@vub.be
GitHub  : github.com/celray - github.com/jopator
'''

import pandas as pd
import geopandas as gpd
import numpy as np
from read_swat import swatModelFile, swat_Dtl
from cjfx import list_folders, goto_dir,exists
import datavariables as variables
import argparse
import ast
import os

# Functions
def doell_DtlActions(const,const2):
    # Actions Dataframe --> This is different for each reservoir
    act_df = pd.DataFrame(columns=['act_typ','obj','obj_num','name','option','const','const2','fp','alt1','alt2'])
    act_df["act_typ"]   = ["release","release"]
    act_df["obj"]       = ["res","res"]
    act_df["obj_num"]   = [int(x) for x in [0,0]]
    act_df["name"]      = ["no_rel","natlake"]
    act_df["option"]    = ["rate","natlake"]
    act_df["const"]     = [int(0),float(const)]
    act_df["const2"]    = [float(x) for x in [0,const2]]
    act_df["fp"]        = ["null","null"]
    act_df["alt1"]      = ["y","n"]
    act_df["alt2"]      = ["n","y"]
    return act_df

def doell_DtlConditions(const):
    # Condition Dataframe --> This is different for each reservoir
    cond_df = pd.DataFrame(columns=['var','obj','obj_num','lim_var','lim_op','lim_const','alt1','alt2'])
    cond_df["var"]          = ["vol"]
    cond_df["obj"]          = ["res"]
    cond_df["obj_num"]      = [int(x) for x in [0]]
    cond_df["lim_var"]      = ["pvol"]
    cond_df["lim_op"]       = ["*"]
    cond_df["lim_const"]    = [float(x) for x in [const]]
    cond_df["alt1"]         = ["<"]
    cond_df["alt2"]         = [">"]
    return cond_df

def Hana06_DtlActions(hana_type,constdoell1,constdoell2,consthana1,consthana2):
    """
    hana_type   : Either "hanazaki_06_gen" or "hanazaki_06_irr"
    constdoell1 : % of pvol at 5 m depth
    constdoell2 : Release rate
    consthana1  : Alpha coefficient (0.85)
    consthana2  : Beta - for irrigation (0.10) 
    """
    # Actions Dataframe --> This is different for each reservoir
    act_df = pd.DataFrame(columns=['act_typ','obj','obj_num','name','option','const','const2','fp','alt1','alt2','alt3','alt4',"alt5"])
    act_df["act_typ"]   = ["release","release","release","release","release"]
    act_df["obj"]       = ["res","res","res","res","res"]
    act_df["obj_num"]   = [int(x) for x in [0,0,0,0,0]]
    act_df["name"]      = ["no_rel","natlake","h06_sch","inflow","overflow"]
    act_df["option"]    = ["rate","natlake",hana_type,"inflo_frac",'ab_emer']
    act_df["const"]     = [float(0),float(constdoell1),float(consthana1),float(1.000),float(0.000)]
    act_df["const2"]    = [float(0),float(constdoell2),float(consthana2),float(0.000),float(0.000)]
    act_df["fp"]        = ["null","null","null","null","null"]
    act_df["alt1"]      = ["y","n","n","n","n"]
    act_df["alt2"]      = ["n","y","n","n","n"]
    act_df["alt3"]      = ["y","n","n","n","n"]
    act_df["alt4"]      = ["n","n","y","n","n"]
    act_df["alt5"]      = ["n","n","y","y","n"]
    act_df["alt6"]      = ["n","n","n","y","y"]
    return act_df
    

def Res_DtlConditions(const):
    # Condition Dataframe --> This is different for each reservoir that is regulated
    cond_df = pd.DataFrame(columns=['var','obj','obj_num','lim_var','lim_op','lim_const','alt1','alt2','alt3','alt4',"alt5","alt6"])
    cond_df["var"] = ["vol","vol","year_seq","year_seq","vol","vol"]
    cond_df["obj"] = ["res","res","res","res","res","res"]
    cond_df["obj_num"] = [int(x) for x in [0,0,0,0,0,0]]
    cond_df["lim_var"] = ["pvol","pvol","null","null","evol","evol"]
    cond_df["lim_op"] = ["*","*","*","*","*","*"]
    cond_df["lim_const"] = [float(x) for x in [0.15,const,6.0,5.0,0.95,1.05]]
    cond_df["alt1"] = ["<","<","<","-","-","-"]
    cond_df["alt2"] = [">",">","<","-","-","-"]
    cond_df["alt3"] = ["<","-","-",">","-","-"]
    cond_df["alt4"] = [">","-","-",">","<","<"]
    cond_df["alt5"] = [">","-","-",">",">","<"]
    cond_df["alt6"] = [">","-","-",">",">",">"]
    return cond_df


def fixLakegdf(processed_gdf,original_gdf):
    def parse_ids(val):
        if isinstance(val, list):
            result = []
            for item in val:
                result.extend(parse_ids(item))
            return result
        try:
            parsed = ast.literal_eval(val)
            return parse_ids(parsed)
        except:
            return [val]
    merged_gdf = processed_gdf[processed_gdf['merged_fla']]
    single_gdf = processed_gdf[~processed_gdf['merged_fla']]
    new_cols = ['Grand_id','smax', 'pvol', 'evol', 'parea', 'earea','br1', 'br2','RES','Lake_type','Lake_name', 'elev_masl', 'Res_time','DOR_PC', 'MAIN_USE','ratios']
    
    if not merged_gdf.empty:
        merged_lakes_data = merged_gdf[['LakeId','Hylak_id','smax','merged_fla','geometry']]
    

        merged_lakes_data_fixed = merged_lakes_data.copy()

        # Initialize cols
        for col in new_cols:
            merged_lakes_data_fixed[col] = None

        # Fix values
        for idx, row in merged_lakes_data_fixed.iterrows():
            hylak_list = [int(x) for x in parse_ids(row['Hylak_id'])]
            original_data = original_gdf[original_gdf['Hylak_id'].isin(hylak_list)]
            smax_new  = float(original_data['smax'].astype(float).sum())

            ratios          = []
            grand_ids       = []
            a_list          = []
            b_list          = []
            c_list          = []
            d_list          = []
            laketype_list   = []
            name_list       = []
            resname_list    = []
            damname_list    = []
            dor_list        = []
            year_list       = []
            rem_year_list   = []
            dam_hgt_list    = []
            dam_len_list    = []
            dam_depth_list  = []
            main_use_list   = []
            smax_largest = 0
            new_hylak_ids   = []
            for hylak_id in hylak_list:
                subset = original_gdf[original_gdf['Hylak_id']==hylak_id]
                orig_smax       = float(subset['smax'].iloc[0])
                grand_ids.append(int(subset['Grand_id'].iloc[0]))
                a_list.append(float(subset['a'].iloc[0]))
                b_list.append(float(subset['b'].iloc[0]))
                c_list.append(float(subset['c'].iloc[0]))
                d_list.append(float(subset['d'].iloc[0]))
                laketype_list.append(str(subset['Lake_type'].iloc[0]))
                name_list.append(str(subset['Lake_name'].iloc[0]))
                resname_list.append(str(subset['RES_NAME'].iloc[0]))
                damname_list.append(str(subset['DAM_NAME'].iloc[0]))
                year_list.append(str(subset['YEAR'].iloc[0]))
                rem_year_list.append(str(subset['REM_YEAR'].iloc[0]))
                dor_list.append(str(subset['DOR_PC'].iloc[0]))
                main_use_list.append(str(subset['MAIN_USE'].iloc[0]))
                new_hylak_ids.append(hylak_id)
                
                if orig_smax>smax_largest:
                    smax_largest = orig_smax
                    hylak_id_largest = hylak_id
                ratios.append(orig_smax/smax_new)
            
            merged_lakes_data_fixed.at[idx,'pvol']  = float(original_data['pvol'].astype(float).sum())
            merged_lakes_data_fixed.at[idx,'evol']  = float(original_data['evol'].astype(float).sum())
            merged_lakes_data_fixed.at[idx,'parea'] = float(original_data['parea'].astype(float).sum())
            merged_lakes_data_fixed.at[idx,'earea'] = float(original_data['earea'].astype(float).sum())

            merged_lakes_data_fixed.at[idx,'br1'] = float(np.average(original_data['br1'].astype(float), weights=original_data['smax'].astype(float)))
            merged_lakes_data_fixed.at[idx,'br2'] = float(np.average(original_data['br2'].astype(float), weights=original_data['smax'].astype(float)))
            merged_lakes_data_fixed.at[idx,'elev_masl'] = float(np.average(original_data['elev_masl'].astype(float), weights=original_data['smax'].astype(float)))
            merged_lakes_data_fixed.at[idx,'Res_time'] = float(np.average(original_data['Res_time'].astype(float), weights=original_data['smax'].astype(float)))

            merged_lakes_data_fixed.at[idx,'smax'] = smax_new

            merged_lakes_data_fixed.at[idx,'Grand_id'] = str(grand_ids)
            merged_lakes_data_fixed.at[idx,'a'] = str(a_list)
            merged_lakes_data_fixed.at[idx,'b'] = str(b_list)
            merged_lakes_data_fixed.at[idx,'c'] = str(c_list)
            merged_lakes_data_fixed.at[idx,'d'] = str(d_list)
            merged_lakes_data_fixed.at[idx,'Lake_type'] = str(laketype_list)
            merged_lakes_data_fixed.at[idx,'Lake_name'] = str(name_list)
            merged_lakes_data_fixed.at[idx,'YEAR'] = str(year_list)
            merged_lakes_data_fixed.at[idx,'REM_YEAR'] = str(rem_year_list)
            merged_lakes_data_fixed.at[idx,'DOR_PC'] = str(dor_list)
            merged_lakes_data_fixed.at[idx,'MAIN_USE'] = str(main_use_list)
            merged_lakes_data_fixed.at[idx,'RES'] = int(1)
            merged_lakes_data_fixed.at[idx,'ratios'] = str(ratios)
            merged_lakes_data_fixed.at[idx,'Hylak_id'] = str(new_hylak_ids)

        columns = ['LakeId', 'Hylak_id','Grand_id', 'smax', 'pvol', 'evol', 'parea', 'earea', 'br1','br2', 'RES','Lake_type', 'elev_masl','Lake_name','DOR_PC', 
                'MAIN_USE','ratios', 'a', 'b', 'c', 'd','merged_fla','geometry']
        merged_lakes_data_fixed = merged_lakes_data_fixed[columns]
        single_gdf = single_gdf.copy()
        single_gdf['ratios'] = '[1]'
        single_gdf = single_gdf[columns].copy()
        return pd.concat([single_gdf,merged_lakes_data_fixed])

    else:
        columns = ['LakeId', 'Hylak_id','Grand_id', 'smax', 'pvol', 'evol', 'parea', 'earea', 'br1','br2', 'RES','Lake_type', 'elev_masl','Lake_name','DOR_PC', 
                'MAIN_USE','ratios', 'a', 'b', 'c', 'd','merged_fla','geometry']
        single_gdf = single_gdf.copy()
        single_gdf['ratios'] = '[1]'
        single_gdf = single_gdf[columns].copy()
        return single_gdf



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

        '''
        Adjust basic model files
        '''
        # Lake and reservoir file
        lakes_grand_FN      = variables.grand_lake_final_shp.format(**details)                                              # Processed reservoir and lake file path
        lakes_original_FN   = f"../model-data/{region}/shapes/lakes-grand-ESRI-54003-demAligned.shp"
        lakes_df            = gpd.read_file(lakes_grand_FN)

        if lakes_df.empty:
            print("No lakes were found, skipping ...")
            quit()

        lakes_df.sort_values(by="LakeId").reset_index(drop=True)
        lakes_original = gpd.read_file(lakes_original_FN).sort_values(by="LakeId").reset_index(drop=True)
        
        # Model files
        txt_in_out_dir    = f"../model-setup/CoSWATv{version}/{region}/Scenarios/Default/TxtInOut"
        hydrology_res     = swatModelFile(f"{txt_in_out_dir}/hydrology.res")
        reservoir_con     = swatModelFile(f"{txt_in_out_dir}/reservoir.con")
        om_water_ini      = swatModelFile(f"{txt_in_out_dir}/om_water.ini")
        
        # Fix lake geodataframe
        print('Preparing lake data ...')
        lakes_df = fixLakegdf(lakes_df,lakes_original)
        lakes_df['LakeId'] = lakes_df['LakeId'].astype(int)
        lakes_df = lakes_df.sort_values(by='LakeId').reset_index(drop=True)
        lakes_df.to_file(f"../model-setup/CoSWATv{version}/{region}/Watershed/Shapes/lakes_coswat-ESRI-54003.gpkg")

        #NOTE: Problems: It's converting emtpy or null to NaN in reservoir.con

        # Overwrite model files
        #===================================================================
        #hydrology.res
        #===================================================================
        hydrology_res.dframe['area_ps'] = lakes_df['parea'].astype(float)
        hydrology_res.dframe['area_es'] = lakes_df['earea'].astype(float)
        hydrology_res.dframe['vol_ps']  = lakes_df['pvol'].astype(float)
        hydrology_res.dframe['vol_es']  = lakes_df['evol'].astype(float)
        
        
        hydrology_res.write(name="hydrology.res updated with adjust-reservoirs.py ",overwrite=True)
        print(f"Updated hydrology.res for {region} \n")
        
        
        #===================================================================
        #reservoir.con
        #===================================================================
        reservoir_con.dframe["area"]    = lakes_df['parea'].astype(float)
        reservoir_con.write(name="reservoir.con updated with adjust-reservoirs.py ",overwrite=True)
        print(f"Updated reservoir.con for {region} \n")

        #===================================================================
        #om_water.ini
        #===================================================================
        om_water_ini.dframe["flo"]    = 1.0000
        
        om_water_ini.write(name="om_water.ini updated with adjust-reservoirs.py ",overwrite=True)
        print(f"Updated om_water.ini for {region}")
        
        '''
        Set up decision tables
        '''

        print("Setting up decision tables ...")
        
        # Model files
        reservoir_res     = swatModelFile(f"{txt_in_out_dir}/reservoir.res")
        res_rel_dtl       = f"{txt_in_out_dir}/res_rel.dtl"
        
        dtl_names = []
        
        
        for index, row in lakes_df.iterrows():
            # Iterate across lakes and reservoirs
            lakeId      = int(row['LakeId'])
            smax_k          = float(row["smax"])/1e9
            dtl_name    = f"res_{lakeId}"
            
            
            if row['merged_fla']:
                merged_hylaks = row['Hylak_id']
                ratios        = row['ratios']
                a_all = row['a'] 
                b_all = row['b']
                c_all = row['c']
                d_all = row['d']
                lakeType_all = row['Lake_type']
                main_use_all = row['MAIN_USE']
                smax_biggest = 0
                
                frac_inactive_list = []

                for hylak_id, ratio,a,b,c,d,lakeType,mainUse in zip(ast.literal_eval(merged_hylaks),
                                                                    ast.literal_eval(ratios),
                                                                    ast.literal_eval(a_all),
                                                                    ast.literal_eval(b_all),
                                                                    ast.literal_eval(c_all),
                                                                    ast.literal_eval(d_all),
                                                                    ast.literal_eval(lakeType_all),
                                                                    ast.literal_eval(main_use_all)
                    ):
                    
                    ratio,a,b,c,d,lakeType,mainUse = float(ratio),float(a),float(b),float(c),float(d),int(lakeType), str(mainUse)
                    smax_single = smax_k*ratio                    
                    try:
                        h_max      = (smax_single/c)**(1/d)
                        h_inactive = h_max - 5
                        if h_inactive<0:
                            h_inactive = 0.1
                        v_inactive = c*(h_inactive**d)
                    except:
                        h_max      = float(lakes_original[lakes_original['Hylak_id']==int(hylak_id)]['maxDepth'].iloc[0])
                        h_inactive = h_max - 5
                        
                        if h_inactive<0:
                            h_inactive = 0.1
                            
                        v_inactive = c*(h_inactive**d)
                    
                    
                    frac_inactive    = round(v_inactive/smax_single,5)
                    frac_inactive_list.append(frac_inactive*ratio)
                    
                    if smax_single > smax_biggest:
                        lakeType_merged = lakeType
                        mainUse_merged  = mainUse
                
                frac_inactive = sum(frac_inactive_list)
                lakeType = lakeType_merged
                mainUse = mainUse_merged
                

            else:
                lakeType    = int(row['Lake_type'])
                hylak_id      = int(row['Hylak_id'])
                a = float(row['a']) 
                b = float(row['b'])
                c = float(row['c'])
                d = float(row['d'])
                mainUse = str(row['MAIN_USE'])
                # We will get the corresponding volume for a depth of 5 m below smax depth
                try:
                    h_max      = (smax_k/c)**(1/d)
                    h_inactive = h_max - 5
                    if h_inactive<0:
                        h_inactive = 0.1
                    v_inactive = c*(h_inactive**d)
                except:
                    h_max      = float(lakes_original[lakes_original['Hylak_id']==int(hylak_id)]['maxDepth'].iloc[0])
                    h_inactive = h_max - 5
                    if h_inactive<0:
                        h_inactive = 0.1
                    v_inactive = c*(h_inactive**d)
                    
                frac_inactive    = round(v_inactive/smax_k,5)
            print(f"  # Setting up dtl for lake {lakeId} ...".ljust(50), end='\r')
                    
            if lakeType == 1:           # If lake is natural (not regulated)
                #===================================================================
                # Decision table for Doell release method & adjustment of model files
                #=================================================================== 
                const1      = frac_inactive
                const2      = 0.01000
                
                # Create action and condition table
                act_df      = doell_DtlActions(const1,const2)
                cond_df     = doell_DtlConditions(const1) 
                
                # Read decision table
                res_dtl = swat_Dtl(res_rel_dtl)
                
                # Add decision table
                dtl_names.append(dtl_name)
                res_dtl.add_dtl(dtl_name,1,2,2,cond_df,act_df,overwrite=True,replace=True)
            
            else:                       # Otherwise is regulated and we need to get the main use

                
                
                if mainUse == 'Irrigation':
                    doellconst1      = frac_inactive
                    doellconst2      = 0.01000
                    
                    cond_df    = Res_DtlConditions(doellconst1)
                    acts_df    = Hana06_DtlActions("irr-h06",doellconst1,doellconst2,0.85000,0.10000)
                    
                    # Read decision table
                    res_dtl = swat_Dtl(res_rel_dtl)
                    
                    # Add decision table
                    dtl_names.append(dtl_name)
                    res_dtl.add_dtl(dtl_name,6,6,5,cond_df,acts_df,overwrite=True,replace=True)

                else:
                    doellconst1      = frac_inactive
                    doellconst2      = 0.01000
                    
                    cond_df    = Res_DtlConditions(doellconst1)
                    acts_df    = Hana06_DtlActions("nonirr-h06",doellconst1,doellconst2,0.85000,0.00000)
                    
                    # Read decision table
                    res_dtl = swat_Dtl(res_rel_dtl)
                    
                    # Add decision table
                    dtl_names.append(dtl_name)
                    res_dtl.add_dtl(dtl_name,6,6,5,cond_df,acts_df,overwrite=True,replace=True)

        #=========================================================
        #Assign decision tables to reservoir.res

        reservoir_res.dframe['rel'] = dtl_names
        # print(dtl_names)
        reservoir_res.write(name="reservoir.res updated write-reservoirs.py",overwrite=True)

        print("\t\t # All done setting up lakes & reservoirs for {region} \n".format(**details))

