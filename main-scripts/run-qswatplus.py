#!/bin/python3

'''
This script runs the COmmunity SWAT+ Model
(CoSWAT-Global) one by one.

Author  : Celray James CHAWANDA
Date    : 14/07/2022
Contact : celray@chawanda.com
Licence : MIT
GitHub  : github.com/celray
'''

import os, sys, platform
from cjfx import list_folders, exists, ignore_warnings, ignore_warnings, goto_dir, pandas
from ccfx import createPath, deleteFile, writeFile
import sqlalchemy
import geopandas as gpd

import os.path
import shutil
import sys
import platform
import warnings
import argparse

ignore_warnings()

if platform.system() == "Linux":
    import pyximport  # importing cython needs this on linux
    pyximport.install()
    ignore_warnings()

# skip deprecation warnings when importing PyQt5
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from qgis.core import *
    from qgis.utils import iface
    from PyQt5.QtGui import *
    from PyQt5.QtCore import *

# QgsApplication.setPrefixPath('C:/Program Files/QGIS 3.10/apps/qgis', True)
qgs = QgsApplication([], True)
qgs.initQgis()

goto_dir(__file__)

# Prepare processing framework
if platform.system() == "Windows":
    sys.path.append('{QGIS_Dir}/apps/qgis/python/plugins'.format(
        QGIS_Dir = os.environ['QGIS_Dir'])) # Folder where Processing is located
else:
    sys.path.append('/usr/share/qgis/python/plugins') # Folder where Processing is located

sys.path.append('../data-preparation')
sys.path.append('../main-scripts')

# extract QSWAT+
if not os.path.exists('../data-preparation/resources/QSWATPlus'):
    shutil.unpack_archive('../data-preparation/resources/QSWATPlus.zip', '../data-preparation/resources/')

# skip syntax warnings on linux

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from processing.core.Processing import Processing
    Processing.initialize()

    import processing


from shapely.geometry import Point, LineString, MultiLineString

def count_intersections(line, polygon):
    intersection = line.intersection(polygon)
    if isinstance(intersection, Point):
        return 1
    elif isinstance(intersection, (LineString, MultiLineString)):
        return len(intersection.geoms) if hasattr(intersection, 'geoms') else 1
    return 0


class outFX:
    """
    Class to handle output messages. currently dummy
    """
    def __init__(self, message:str, v = False):
        self.message = message
        self.v = v
        self.messageBank = []
        self.write(message)

    def write(self, message: str) -> None:
        """
        Write message to console and log file.
        """
        
        if self.v: print(message)
        self.messageBank.append(message)
    
    def append(self, message: str) -> None:
        """
        Append message to console and log file.
        """
        if self.v: self.write(message)
        self.messageBank.append(message)

    def moveCursor(self, cursor) -> None:
        """
        Move the cursor to a new position in the output.
        """
        pass

    def __call__(self, message: str) -> None:
        """
        Allow the object to be called like a function.
        """
        self.append(message)


import atexit


from resources.QSWATPlus.QSWATPlusMain import QSWATPlus
from resources.QSWATPlus.delineation import Delineation
from resources.QSWATPlus.floodplain import Floodplain
from resources.QSWATPlus.landscape import Landscape
from resources.QSWATPlus.raster import Raster
from resources.QSWATPlus.hrus import HRUs
from resources.QSWATPlus.QSWATUtils import QSWATUtils
from resources.QSWATPlus.parameters import Parameters

import datavariables as variables

from glob import glob

atexit.register(QgsApplication.exitQgis)


details = {
    'auth': variables.final_proj_auth,
    'code': variables.final_proj_code,
}


class DummyInterface(object):
    """Dummy iface to give access to layers."""

    def __getattr__(self, *args, **kwargs):
        """Dummy function."""
        def dummy(*args, **kwargs):
            return self
        return dummy

    def __iter__(self):
        """Dummy function."""
        return self

    def __next__(self):
        """Dummy function."""
        raise StopIteration

    def layers(self):
        """Simulate iface.legendInterface().layers()."""
        return list(QgsProject.instance().mapLayers().values())





if __name__ == '__main__':

    # change working directory
    goto_dir(__file__)

    print(os.getcwd())
    # create argument parser
    parser = argparse.ArgumentParser(description="a terminal version of QSWAT+ for running the model setup and delineation.")

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

    print(f"\nregions to run: {', '.join(regions)}")
    print(f"CoSWAT version: {version}")

    for region in regions:

        details['version'] = version
        details['region']  = region

        print(f'\n\nrunning QSWAT+ for region: {region} ({version})')
        iface   = DummyInterface()
        plugin  = QSWATPlus(iface)
        dlg     = plugin._odlg  # useful shorthand for later
        
        projDir = f'../model-setup/CoSWATv{version}/{region}'
        data_dir= f'../model-data/{region}'

        if not os.path.exists(projDir):
            QSWATUtils.error('Project directory {0} not found'.format(projDir), True)
            sys.exit(1)

        projFile = f"{projDir}/{region}.qgs"

        proj = QgsProject.instance()
        
        proj.read(projFile)

        plugin.setupProject(proj, True)

        # make connection and load tables
        landuse_table   = f"{data_dir}/tables/worldLanduseLookup.csv"
        soil_table      = f"{data_dir}/tables/worldSoilsLookup.csv"
        user_soil_table = f"{data_dir}/tables/worldSoilsUsersoil.csv"

        landuse_df      = pandas.read_csv(landuse_table, names=["LANDUSE_ID", "SWAT_CODE"], skiprows=1)
        soil_df         = pandas.read_csv(soil_table, names=["SOIL_ID", "NAME"], skiprows=1)
        user_soil_df    = pandas.read_csv(user_soil_table)

        user_soil_df            = user_soil_df.fillna("")
        user_soil_df['SEQN']    = user_soil_df['SEQN'].astype(str)

        db = sqlalchemy.create_engine(f'sqlite:///{projDir}/{region}.sqlite')

        landuse_df.to_sql('landuse_lookup', db, if_exists="replace", index=False)
        soil_df.to_sql('soil_lookup', db, if_exists="replace", index=False)
        user_soil_df.to_sql('usersoil', db, if_exists="replace", index=False, )

        LOGFILE = os.path.abspath(f'../model-setup/CoSWATv{version}/{region}/qswat_errors.log')
        
        
        if variables.new_res_methods:
            lakesShapefn    = os.path.abspath(f'../model-setup/CoSWATv{version}/{region}/Watershed/Shapes/lakes-grand-{variables.final_proj_auth}-{variables.final_proj_code}.shp')
            rivsShapefn     = os.path.abspath(f'../model-setup/CoSWATv{version}/{region}/Watershed/Shapes/dem-aster-{variables.final_proj_auth}-{variables.final_proj_code}-lakeBurntchannel.shp')
        else:
            rivsShapefn     = os.path.abspath(f'../model-setup/CoSWATv{version}/{region}/Watershed/Shapes/dem-aster-{variables.final_proj_auth}-{variables.final_proj_code}channel.shp')
            lakesShapefn    = ''

        def delineateBasin() -> None:
            plugin._gv.db.clearTable('BASINSDATA')
            plugin.setupProject(proj, True)

            if not (os.path.exists(plugin._gv.textDir) and os.path.exists(plugin._gv.landuseDir)):
                QSWATUtils.error('Directories not created', True)
                sys.exit(1)

            if not dlg.delinButton.isEnabled():
                QSWATUtils.error('Delineate button not enabled', True)
                sys.exit(1)

            delin = Delineation(plugin._gv, plugin._demIsProcessed)
            delin.init()
            delin._dlg.numProcesses.setValue(variables.taudemProcesses)

            QSWATUtils.information('DEM: {0}'.format(os.path.split(plugin._gv.demFile)[1]), True)
            delin.addHillshade(plugin._gv.demFile, None, None, None)
            QSWATUtils.information('Inlets/outlets file: {0}'.format(os.path.split(plugin._gv.outletFile)[1]), True)

            outlets_buffer_gpd  = gpd.read_file(f"../data-preparation/resources/regions/{region}/outlets-buffer.gpkg").to_crs('{auth}:{code}'.format(**details))
            
            outletsCreated = False

            if os.path.exists(os.path.abspath(f'../model-setup/CoSWATv{version}/{region}/Watershed/Shapes/outlets.shp')):
                outletsCreated = True



            delin.runTauDEM2(ver = version, reg = region,
                in_outlet_path = os.path.abspath(f'../model-setup/CoSWATv{version}/{region}/Watershed/Shapes/outlets.shp'),
                Mask_gpd    = outlets_buffer_gpd,
                sel_file    = os.path.abspath(f'../model-setup/CoSWATv{version}/{region}/Watershed/Shapes/outlets_sel.shp'),
                outletsExist2 = outletsCreated
            )


            

            # ===============================================================================================================================================
            # Flood plains will be conditioned in datavariables
            # ===============================================================================================================================================

            if variables.runFloodplains:
                print("Running floodplain...")
                createPath(f'../model-setup/CoSWATv{version}/{region}/Watershed/Rasters/Landscape/Flood/')
                writeFile(f'../model-setup/CoSWATv{version}/{region}/Watershed/Rasters/Landscape/Flood/creatingFloodPlain', 'Creating floodplain...\nThis is just an indicator file\nit will be removed when the floodplain is created')
                fxObj           = outFX('Running floodplain...')
                floodPlain      = Floodplain(plugin._gv, fxObj, 1)
                landScape       = Landscape(plugin._gv, fxObj, 1, fxObj)

                landScape.clipperFile = plugin._gv.subbasinsFile
                landScape.calcHillslopes(variables.thresholdCh, landScape.clipperFile, proj.layerTreeRoot())

                landScape.calcFloodplain(True, proj.layerTreeRoot())
                plugin._gv.floodFile = os.path.abspath(f'../model-setup/CoSWATv{version}/{region}/Watershed/Rasters/Landscape/Flood/invflood0_00.tif')
                deleteFile(f'../model-setup/CoSWATv{version}/{region}/Watershed/Rasters/Landscape/Flood/creatingFloodPlain')
            
            else:
                print("Floodplains skipped...")        

            if variables.new_res_methods:
                print("\t Resolving reservoir shapes into channel network ...")
                os.system(f'python3 resolve-lakes-reservoirs.py {version} {region}')
                
            else:
                delin.lakesDone = True
                plugin._gv.lakeFile == ''

            delin.finishDelineation(details = details)

            # Patch QSWATUtils to log it's errors in a file
            _orig_error = QSWATUtils.error

            def _error_to_file(msg: str, isError: bool = True) -> None:
                try:
                    with open(LOGFILE, "a", encoding="utf-8") as f:
                        f.write(f"{msg.strip()}\n")

                except Exception:
                    pass

                _orig_error(msg, isError)
            QSWATUtils.error = staticmethod(_error_to_file)

        delineateBasin()


        def createHRUS() -> None:
            if not dlg.hrusButton.isEnabled():
                QSWATUtils.error('\t ! HRUs button not enabled', True)
                sys.exit(1)

            hrus = HRUs(plugin._gv, dlg.reportsBox)
            hrus.init()
            hrus._gv.useLandscapes = True
            hrus._dlg.generateFullHRUs.setEnabled(True)
            hrus.fullHRUsWanted = True

            if variables.runFloodplains:
                hrus.initFloodplain()
            
            
            hrus.readFiles()



            if not os.path.exists(QSWATUtils.join(plugin._gv.textDir, Parameters._TOPOREPORT)):
                QSWATUtils.error('\t ! Elevation report not created \n\n\t   Have you run Delineation?\n', True)
                sys.exit(1)

            if not os.path.exists(QSWATUtils.join(plugin._gv.textDir, Parameters._BASINREPORT)):
                QSWATUtils.error('\t ! Landuse and soil report not created', True)
                sys.exit(1)
            
            hrus.calcHRUs()
            if not os.path.exists(QSWATUtils.join(plugin._gv.textDir, Parameters._HRUSREPORT)):
                QSWATUtils.error('\t ! HRUs report not created', True)
                sys.exit(1)
            
            

            if not os.path.exists(QSWATUtils.join(projDir, r'Watershed/Shapes/rivs1.shp')):
                QSWATUtils.error('\t ! Streams shapefile not created', True)
                sys.exit(1)

            if not os.path.exists(QSWATUtils.join(projDir, r'Watershed/Shapes/subs1.shp')):
                QSWATUtils.error('\t ! Subbasins shapefile not created', True)
                sys.exit(1)


            QSWATUtils.information('\t - finished creating HRUs\n', True)

        createHRUS()


        # We need to run again if there were errors :( - in this case, correcting the geometry
        # we check the error file log (to be created and send message based on that)
        if variables.new_res_methods:
            errors_exist = os.path.exists(LOGFILE)
            
            while errors_exist:

                if os.path.exists(LOGFILE):

                    with open(LOGFILE) as logfile:
                        
                        lines = logfile.readlines()
                        
                        if len(lines)>0:
                            print('\t Errors exist --> Need to readjust lakes \n')

                            subsNoLakes_fn  = os.path.abspath(f'../model-setup/CoSWATv{version}/{region}/Watershed/Shapes/subsNoLakes.shp')
                            subs1_fn        = os.path.abspath(f'../model-setup/CoSWATv{version}/{region}/Watershed/Shapes/subs1.shp')
                            hrus2_fn        = os.path.abspath(f'../model-setup/CoSWATv{version}/{region}/Watershed/Shapes/hrus2.shp')
                            lsus_fn         = os.path.abspath(f'../model-setup/CoSWATv{version}/{region}/Watershed/Shapes/lsus2.shp')
                            lake_shape_snap = os.path.abspath(f'../model-data/{region}/shapes/lakes-grand-ESRI-54003-demAligned-forSnap.shp')

                            #rivsShapefn name for channels

                            # Read as geoDataFrames
                            lakes_gdf       = gpd.read_file(lakesShapefn)
                            subsNoLakes_gdf = gpd.read_file(subsNoLakes_fn) 
                            lsus_gdf        = gpd.read_file(lsus_fn)
                            chan_gdf        = gpd.read_file(rivsShapefn)
                            hrus_gdf        = gpd.read_file(hrus2_fn)
                            subs_gdf        = gpd.read_file(subs1_fn)

                            lake_orig = lakes_gdf.copy().drop(columns='geometry')


                            # Down here ...
                            """
                            We deal with the error: routing sink category xxx not found as a source which is because the LSU of that element is to small and is removed

                            We add the small LSUs that cause the error in the channel network to the lake shape
                            These are just small LSUs (smaller than a pixel) that QSWAT+ ignores and therefore ignore its channel - but the channel is part of the lake connections
                            """
                            small_lsus = False
                            for line in lines:
                                line_split = line.split()

                                if line_split[2] == "routing":
                                    small_lsus = True
                                    print('\t -> There are small LSUS that are causing issues with the channel network ...')
                            
                            if small_lsus:
                            # Overlays --> difference --> filter out --> buffer --> add --> Fix indices --> Save new shape
                                subs_overlay_gdf = gpd.overlay(subsNoLakes_gdf, lakes_gdf, how='difference', keep_geom_type = True, make_valid = True)

                                small_lsus_gdf = gpd.overlay(subs_overlay_gdf, lsus_gdf, how='difference', keep_geom_type = True, make_valid = True)
                                small_lsus_gdf['new_area'] = small_lsus_gdf['geometry'].area
                                small_lsus_gdf = small_lsus_gdf[small_lsus_gdf.new_area < variables.data_resolution**2]                                         # Filter out small ones

                                crossed = gpd.sjoin(small_lsus_gdf, chan_gdf, how="inner", predicate="crosses")
                                result  = small_lsus_gdf.loc[small_lsus_gdf.index.isin(crossed.index)]                                                          # Only those that cross channels matter
                                small_lsus_gdf = result.copy()

                                small_lsus_gdf['geometry'] = small_lsus_gdf['geometry'].buffer(variables.data_resolution/2+20, join_style=2, mitre_limit=1e6)   # Buffer to match snap

                                lakes_gdf['_idx']   = lakes_gdf.index
                                small_lsus_gdf_idx  = gpd.sjoin_nearest(small_lsus_gdf[['geometry']],lakes_gdf[['_idx','geometry']],how='left')                 # Get index to match small lsus and lake
                                join = gpd.overlay(lakes_gdf[['LakeId','_idx','geometry']], small_lsus_gdf_idx, how ='union')

                                join["_idx_1"] = join["_idx_1"].fillna(join["_idx_2"])                                                                          # Fix indices
                                join = join[['LakeId','_idx_1','geometry']].reset_index(drop=True)
                                merged = join.dissolve(by="_idx_1").reset_index()                                                                               
                                merged['geometry'] = merged['geometry'].buffer(-20, join_style=2, mitre_limit=1e6)                                              # Merge the corrected geometry with the original data
                                merged = pandas.merge(merged,lake_orig)
                                merged = merged.drop('_idx_1',axis=1)
                                merged.to_file(lake_shape_snap)                                                                                                 # Save as new snapping lake reference shape


                            # Delete logfile
                            os.remove(LOGFILE)

                            # After adjusting, re-do
                            delineateBasin()
                            createHRUS()

                        else:
                            errors_exist = False
                            break
                else:
                    errors_exist = False
                    break

        print()
        print(f'done with running qswat+ for region {region}', '\nQSWAT+ run complete')

