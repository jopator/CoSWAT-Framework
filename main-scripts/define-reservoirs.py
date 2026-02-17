'''
This script modifies the model files to:

    - Correctly define area and storage properties of the lakes based on global datasets: HydroLakes, GranD, GLOBAthy
    - Define release decision tables for lakes and reservoirs.


========================================================================================================================
Algorithms used for release decision tables:

    * Unregulated lakes: Doell et al. (2003)
    * Water Supply, flood control and irrigation reservoirs: Hanazaki et al. (2006) / Vanderkelen et al. (2022) / Gharari et al. (2024)
    * Hydropower reservoirs: Arheimer et al. (2020) / Gharari et al. (2024)

    > Main use is determined by global datasets - If not on GranD , assumed as natural lake.

    ! This will only work if using swatplus-61.0.2.11-31-coswatn5-ifx-lin_x86_64 executable
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
import sys
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
    cond_df["var"] = ["vol","vol","year_seq","vol","vol"]
    cond_df["obj"] = ["res","res","res","res","res"]
    cond_df["obj_num"] = [int(x) for x in [0,0,0,0,0]]
    cond_df["lim_var"] = ["pvol","pvol","null","evol","evol"]
    cond_df["lim_op"] = ["*","*","*","*","*"]
    cond_df["lim_const"] = [float(x) for x in [0.15,const,5.0,0.95,1.05]]
    cond_df["alt1"] = ["<","<","<","-","-"]
    cond_df["alt2"] = [">",">","<","-","-"]
    cond_df["alt3"] = ["<","-",">","-","-"]
    cond_df["alt4"] = [">","-",">","<","<"]
    cond_df["alt5"] = [">","-",">",">","<"]
    cond_df["alt6"] = [">","-",">",">",">"]
    return cond_df


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

        # Model files
        txt_in_out_dir    = f"../model-setup/CoSWATv{version}/{region}/Scenarios/Default/TxtInOut"
        hydrology_res     = swatModelFile(f"{txt_in_out_dir}/hydrology.res")
        reservoir_con     = swatModelFile(f"{txt_in_out_dir}/reservoir.con")
        om_water_ini      = swatModelFile(f"{txt_in_out_dir}/om_water.ini")
        
        # Lake and reservoir file
        lakes_grand_path = variables.grand_lake_final_shp.format(**details)                                              # Processed reservoir and lake file path
        lakes_df         = gpd.read_file(lakes_grand_path).sort_values(by="LakeId").reset_index(drop=True)               
        
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

        # Model files
        reservoir_res     = swatModelFile(f"{txt_in_out_dir}/reservoir.res")
        res_rel_dtl       = f"{txt_in_out_dir}/res_rel.dtl"
        
        dtl_names = []
        for index, row in lakes_df.iterrows():                                                                      # Iterate across lakes and reservoirs
            lakeType    = row['Lake_type']
            lakeId      = int(row['LakeId'])
            dtl_name    = f"res_{lakeId}"
            a = float(row['a']) 
            b = float(row['b'])
            c = float(row['c'])
            d = float(row['d'])
            
            # We will get the corresponding volume for a depth of 5 m below smax depth
            smax_k          = float(row["smax"])/1000000
            h_max           = (smax_k/c)**(1/d)
            h_inactive      = h_max - 5

            v_inactive      = c*(h_inactive**d)
            vol_5m_depth    = round(v_inactive/smax_k,5)
            print(f"# Setting up decision table for lake {lakeId} ...")	 
                    
            if lakeType == 1:           # If lake is natural (not regulated)
                #===================================================================
                # Decision table for Doell release method & adjustment of model files
                #=================================================================== 
                const1      = vol_5m_depth
                const2      = 0.01000
                
                # Create action and condition table
                act_df      = doell_DtlActions(const1,const2)
                cond_df     = doell_DtlConditions(const1) 
                
                # Read decision table
                res_dtl = swat_Dtl(res_rel_dtl)
                
                # Add decision table
                dtl_names.append(dtl_name)
                res_dtl.add_dtl(dtl_name,1,2,2,cond_df,act_df,overwrite=True)
            
            else:                       # Otherwise is regulated and we need to get the main use

                mainUse     = row['MAIN_USE']
                
                if mainUse == 'Irrigation':
                    doellconst1      = vol_5m_depth
                    doellconst2      = 0.01000
                    
                    cond_df    = Res_DtlConditions(doellconst1)
                    acts_df    = Hana06_DtlActions("irr-h06",doellconst1,doellconst2,0.85000,0.10000)
                    
                    # Read decision table
                    res_dtl = swat_Dtl(res_rel_dtl)
                    
                    # Add decision table
                    dtl_names.append(dtl_name)
                    res_dtl.add_dtl(dtl_name,5,6,5,cond_df,acts_df,overwrite=True)

                else:
                    doellconst1      = vol_5m_depth
                    doellconst2      = 0.01000
                    
                    cond_df    = Res_DtlConditions(doellconst1)
                    acts_df    = Hana06_DtlActions("nonirr-h06",doellconst1,doellconst2,0.85000,0.00000)
                    
                    # Read decision table
                    res_dtl = swat_Dtl(res_rel_dtl)
                    
                    # Add decision table
                    dtl_names.append(dtl_name)
                    res_dtl.add_dtl(dtl_name,5,6,5,cond_df,acts_df,overwrite=True)

        #=========================================================
        #Assign decision tables to reservoir.res

        reservoir_res.dframe['rel'] = dtl_names
        reservoir_res.write(name="reservoir.res updated with doell release dtls",overwrite=True)

        print("\n\n \t\t\t # All done setting up lakes & reservoirs for {region} \n".format(**details))

