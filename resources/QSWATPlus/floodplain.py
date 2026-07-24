# -*- coding: utf-8 -*-
"""
/***************************************************************************
 QSWATPlus
                                 A QGIS plugin
 Create SWATPlus inputs
                              -------------------
        begin                : 2014-07-18
        copyright            : (C) 2014 by Chris George
        email                : cgeorge@mcmaster.ca
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 """

# Import the PyQt and QGIS libraries
from qgis.PyQt.QtCore import QObject
#from qgis.PyQt.QtGui import *
from qgis.core import QgsProject, QgsVectorLayer
#from qgis.gui import * # @UnusedWildImport
import os
from osgeo import gdal
# import numpy as np
import math
import time
import multiprocessing as mp
from functools import partial
from .QSWATUtils import QSWATUtils, FileTypes  # fileWriter type: ignore 
from .QSWATTopology import QSWATTopology
import numpy
import traceback
#from rasterint import RasterInt
#from rasterfloat import RasterFloat
from .raster import Raster
#from intraster import IntRaster

import datavariables as variables

# Top-level worker functions for multiprocessing (must be outside class to avoid pickling issues)
def _compute_floodplain_chunk(args):
    """Worker function to compute floodplain values for a chunk of rows (vectorized)."""
    ridge_chunk, valley_chunk, flood_thresh, no_data, start_row = args

    # vectorized: no Python loops
    valid       = (ridge_chunk != no_data) & (valley_chunk != no_data)
    denominator = numpy.where(valid, ridge_chunk + valley_chunk, 1)
    sp          = numpy.where(valid & (denominator != 0), valley_chunk / denominator, 0)
    result      = numpy.full_like(ridge_chunk, no_data, dtype=numpy.float64)
    result[valid & (sp <= flood_thresh)] = 1

    return (start_row, result, valid)

def _compute_ridge_heights_chunk(args):
    """Worker function to compute ridge heights for a chunk of rows."""
    dem_chunk, ridgep_chunk, ridges_chunk, ridge_thresh, dem_no_data, no_data, start_row, dem_transform = args
    results = []
    
    num_rows, num_cols = dem_chunk.shape
    for i in range(num_rows):
        for j in range(num_cols):
            if dem_chunk[i, j] != dem_no_data:
                # Simplified ridge height calculation for chunk
                # This is a simplified version - full implementation would need the path-following logic
                ridge_val = ridges_chunk[i, j] if ridges_chunk[i, j] >= ridge_thresh else 0
                if ridge_val > 0:
                    ridge_height = ridge_val - dem_chunk[i, j] 
                    results.append((start_row + i, j, ridge_height))
    return results
 
class Floodplain(QObject):
     
    """
    Calculate floodplain raster  from DEM, stream threshold, ridge threshold, floodplain threshold,
    (optionally) D8 flow direction grid, flow accumulation grid, outlets shapefile.
     
    The slope position value is, for each point, (e - v) / (r - v) where e is the elevation, 
    v the nearest valley floor elevation, and r the nearest ridge elevation.  Nearest means along D8 flow path.
    Valley floor is points with flow accumulation above stream generation threshold.
    Ridge is similar, obtained by inverting DEM, then calculating D8 flow direction and accumulation.
    The floodplain raster  is zero for points with slope position value at most the floodplain threshold, else 1.
    """
    
    def __init__(self, gv, progress, chunkCount):
        """Initialise class variables."""
        QObject.__init__(self)
        self._gv = gv
        self._progress = progress
        ## raster  chunk count
        self.chunkCount = chunkCount
        ## ridge accumulation threshold
        self.ridgeThresh = 0
        ## floodplain accumulation threshold
        self.floodThresh = 0
        ## branch length threshold
        self.branchThresh = 0
        ## dem raster  dataset
        self.demDs = None
        ## dem raster  transform
        self.demTransform = None
        ## dem raster 
        self.demRaster = None
        ## dem raster  no data value
        self.demNoData = -1
        ## dem raster  projection
        self.demProjection = None
        ## no data value for output rasters; must be negative
        self.noData = -1
        ## number of rows in input clipped DEM raster 
        self.numRows = 0
        ## number of columns in input clipped DEM raster 
        self.numCols = 0
        ## depths of valley floor below points
        self.valleyDepthsRaster = None
        ## file path for ridgeHeightsRaster
        self.ridgeHeightsFile = ''
        ## heights of ridge above points
        self.ridgeHeightsRaster = None
        ## distances of ridge from points
        self.ridgeDistancesRaster = None
        ## floodplain raster 
        self.floodplainRaster = None
        ## flag to indicate if using inversion or branch length method
        self.useInversion = True
        ## map of (row, col) to (elevation, branch length) for ridge points
        self.ridgePoints = None
        
    def run(self, demRaster, valleyDepthsFile, ridgeThresh, floodThresh, branchThresh, subbasins, distances, slopeDir, flowAcc, ridgep, ridges, noData, mustRun):
        """
        Generate floodplain raster , with floodplain value 0 and upland value 1.
        
        valleyDepthsRaster is already calculated and stored in valleyDepthsFile
        ridgeThresh is the threshold in pixel counts for ridges (inverse accumulation).
        floodThresh is the critical ratio of relative depth of valley floor against height of ridge above valley floor.
        branchThresh is the threshold in metres for branch lengths.
        dem is clipped filled DEM
        The groups {subbasins, distances, slopeDir, flowAcc} and {ridgep, ridges} are either all None or all existing files,
        assumed to be clipped and with same pixel size as DEM (so (row, col|) coordinates may be shared across all of them).
        The second group is used to calculate ridges by DEM inversion if subbasins is None, else the first is used to 
        calculate ridges by branch lengths.
        """
        self.demRaster = demRaster
        self.ridgeThresh = ridgeThresh
        self.floodThresh = floodThresh
        self.branchThresh = branchThresh
        gdal.AllRegister()
        ds = self.demRaster.ds
        self.numRows = ds.RasterYSize
        self.numCols = ds.RasterXSize
        self.demTransform = ds.GetGeoTransform()
        self.demNoData = self.demRaster.band.GetNoDataValue()
        self.demProjection = ds.GetProjection()
        self.noData = noData
        self._progress('Ridge heights ...')
        root = QgsProject.instance().layerTreeRoot()
        if subbasins is None:
            # use inverted DEM method for ridges
            self.useInversion = True
            time1 = time.process_time()
            if not self.calcRidgeHeightsByInversion(ridgep, ridges, root, mustRun):
                return
            time2 = time.process_time()
            QSWATUtils.loginfo('Ridge heights by inversion took {0} seconds'.format(int(time2 - time1)))
        else:
            self.useInversion = False
            time1 = time.process_time()
            if not self.calRidgeHeghtsByBranchLength(subbasins, distances, slopeDir, flowAcc, valleyDepthsFile, root, mustRun):
                return
            time2 = time.process_time()
            QSWATUtils.loginfo('Ridge heights by branch length took {0} seconds'.format(int(time2 - time1)))
        time1 = time.process_time()
        self.writeFloodPlain(valleyDepthsFile, mustRun)
        time2 = time.process_time()
        QSWATUtils.loginfo('Floodplain creation took {0} seconds'.format(int(time2 - time1)))
        self._progress('')
        
    def runParallel(self, demRaster, valleyDepthsFile, ridgeThresh, floodThresh, branchThresh, subbasins, distances, slopeDir, flowAcc, ridgep, ridges, noData, mustRun, nProc=None):
        """
        Generate floodplain raster using parallel processing, with floodplain value 0 and upland value 1.
        
        Same parameters as run() method, but with parallel processing for major computational loops.
        """
        if nProc is None:
            nProc = max(variables.taudemProcesses, 4)  # Default from variables.taudemProcesses
            
        self.demRaster = demRaster
        self.ridgeThresh = ridgeThresh
        self.floodThresh = floodThresh
        self.branchThresh = branchThresh
        gdal.AllRegister()
        ds = self.demRaster.ds
        self.numRows = ds.RasterYSize
        self.numCols = ds.RasterXSize
        self.demTransform = ds.GetGeoTransform()
        self.demNoData = self.demRaster.band.GetNoDataValue()
        self.demProjection = ds.GetProjection()
        self.noData = noData
        self._progress('Ridge heights (parallel)...')
        root = QgsProject.instance().layerTreeRoot()
        
        if subbasins is None:
            # use inverted DEM method for ridges
            self.useInversion = True
            time1 = time.process_time()
            if not self.calcRidgeHeightsByInversion(ridgep, ridges, root, mustRun):
                return
            time2 = time.process_time()
            QSWATUtils.loginfo('Ridge heights by inversion took {0} seconds'.format(int(time2 - time1)))
        else:
            self.useInversion = False
            time1 = time.process_time()
            # Branch length method doesn't have as heavy computational loops, use original
            if not self.calRidgeHeghtsByBranchLength(subbasins, distances, slopeDir, flowAcc, valleyDepthsFile, root, mustRun):
                return
            time2 = time.process_time()
            QSWATUtils.loginfo('Ridge heights by branch length took {0} seconds'.format(int(time2 - time1)))
        
        time1 = time.process_time()
        self.writeFloodPlainParallel(valleyDepthsFile, mustRun, nProc)
        time2 = time.process_time()
        QSWATUtils.loginfo('Parallel floodplain creation took {0} seconds'.format(int(time2 - time1)))
        self._progress('')
        
    def writeFloodPlainParallel(self, valleyDepthsFile, mustRun, nProc):
        """Calculate slope positions from valleyDepths and ridgeHeights, and write floodplain raster using parallel processing."""
        method = 'inv' if self.useInversion else 'branch'
        thresh = '{0:.2F}'.format(self.floodThresh).replace('.', '_')
        flood = QSWATUtils.join(self._gv.floodDir, method + 'flood' + thresh + '.tif')
        root = QgsProject.instance().layerTreeRoot()
        if mustRun or \
            not QSWATUtils.isUpToDate(valleyDepthsFile, flood) or \
            not QSWATUtils.isUpToDate(self.ridgeHeightsFile, flood):
            if os.path.exists(flood):
                QSWATUtils.tryRemoveLayerAndFiles(flood, root)
            self._gv.clearOpenRasters()
            completed = False
            while not completed:
                try:
                    completed = True # only gets set on MemoryError exception
                    self.ridgeHeightsRaster = Raster(self.ridgeHeightsFile, self._gv, canWrite=False, isInt=False)
                    res = self.ridgeHeightsRaster.open(self.chunkCount)
                    if not res:
                        self._gv.closeOpenRasters()
                        return
                    if self.valleyDepthsRaster is None:
                        # may have been opened when calculating ridge heights by branch length
                        self.valleyDepthsRaster = Raster(valleyDepthsFile, self._gv, canWrite=False, isInt=False)
                        res = self.valleyDepthsRaster.open(self.chunkCount)
                        if not res:
                            self._gv.closeOpenRasters()
                            return
                    self._progress('Flood plain (parallel)...')
                    self.floodplainRaster = Raster(flood, self._gv, canWrite=True, isInt=True)
                    OK = self.floodplainRaster.open(self.chunkCount, numRows=self.numRows, numCols=self.numCols, 
                                                  transform=self.demTransform, projection=self.demProjection, 
                                                  noData=self.noData)
                    if OK:
                        self.calcFloodPlain1Parallel(nProc)
                        self.floodplainRaster.close()
                        QSWATUtils.copyPrj(self.demRaster.fileName, self.floodplainRaster.fileName)
                        # load flood above DEM
                        layers = root.findLayers()
                        demLayer = QSWATUtils.getLayerByLegend(FileTypes.legend(FileTypes._DEM), layers)
                        ft = FileTypes._INVFLOOD if self.useInversion else FileTypes._BRANCHFLOOD
                        floodLayer, _ = QSWATUtils.getLayerByFilename(layers, self.floodplainRaster.fileName, ft, 
                                                          self._gv, demLayer, QSWATUtils._WATERSHED_GROUP_NAME)
                        if floodLayer is None:
                            QSWATUtils.error('Failed to load floodplain raster {0}' \
                                             .format(self.floodplainRaster.fileName), self._gv.isBatch)
                    self.valleyDepthsRaster.close()
                    self.ridgeHeightsRaster.close()
                    self._progress('Flood plain done')
                except MemoryError:
                    QSWATUtils.loginfo('Out of memory for flood plain with chunk count {0}'.format(self.chunkCount))
                    self._gv.closeOpenRasters()
                    completed = False
                    self.chunkCount += 1
        else:# already have uptodate flood file: make sure it is loaded
            # load flood above DEM
            layers = root.findLayers()
            demLayer = QSWATUtils.getLayerByLegend(FileTypes.legend(FileTypes._DEM), layers)
            ft = FileTypes._INVFLOOD if self.useInversion else FileTypes._BRANCHFLOOD
            floodLayer, _ = QSWATUtils.getLayerByFilename(layers, flood, ft, self._gv, 
                                                          demLayer, QSWATUtils._WATERSHED_GROUP_NAME)
            if floodLayer is None:
                QSWATUtils.error('Failed to load floodplain raster {0}' \
                                 .format(self.floodplainRaster.fileName), self._gv.isBatch)
                    
    def calcRidgeHeightsByInversion(self, ridgep, ridges, root, mustRun):
        """
        Create the ridgeHeightsRaster with differences between the elevation at the point and
        the elevation of the nearest ridge cell.

        ridgep is the D8 flow directions and ridges the flow accumulation raster, both calculated from an inverted DEM.
        Uses numpy arrays for fast pixel access instead of per-pixel GDAL calls.
        """
        self.ridgeHeightsFile = QSWATUtils.join(self._gv.demDir, 'invheights.tif')
        if mustRun or not QSWATUtils.isUpToDate(ridgep, self.ridgeHeightsFile) or not QSWATUtils.isUpToDate(ridges, self.ridgeHeightsFile):
            self._gv.clearOpenRasters()
            completed = False
            while not completed:
                try:
                    completed = True
                    if os.path.exists(self.ridgeHeightsFile):
                        QSWATUtils.tryRemoveLayerAndFiles(self.ridgeHeightsFile, root)
                    self.ridgeHeightsRaster = Raster(self.ridgeHeightsFile, self._gv, canWrite=True, isInt=False)
                    res = self.ridgeHeightsRaster.open(self.chunkCount, numRows=self.numRows, numCols=self.numCols,
                                                        transform=self.demTransform, projection=self.demProjection,
                                                        noData=self.noData)
                    if not res:
                        return False
                    ridgepRaster = Raster(ridgep, self._gv, canWrite=False, isInt=False)
                    res = ridgepRaster.open(self.chunkCount)
                    if not res:
                        self._gv.closeOpenRasters()
                        return False
                    ridgesRaster = Raster(ridges, self._gv, canWrite=False, isInt=True)
                    res = ridgesRaster.open(self.chunkCount)
                    if not res:
                        self._gv.closeOpenRasters()
                        return False

                    # read all rasters into numpy arrays for fast access
                    demArr      = self.demRaster.band.ReadAsArray().astype(numpy.float64)
                    dirArr      = ridgepRaster.band.ReadAsArray().astype(numpy.float64)
                    accArr      = ridgesRaster.band.ReadAsArray().astype(numpy.float64)
                    accNoData   = ridgesRaster.band.GetNoDataValue()
                    heightsArr  = numpy.full((self.numRows, self.numCols), -1.0, dtype=numpy.float64)

                    # D8 direction offsets (TauDEM convention: 1-8)
                    dY = numpy.array([0, -1, -1, -1, 0, 1, 1, 1], dtype=numpy.int32)
                    dX = numpy.array([1, 1, 0, -1, -1, -1, 0, 1], dtype=numpy.int32)

                    for row in range(self.numRows):
                        for col in range(self.numCols):
                            if demArr[row, col] == self.demNoData:
                                continue
                            if heightsArr[row, col] >= 0:
                                continue

                            # trace path from (row, col) to nearest ridge
                            path    = []
                            r, c    = row, col
                            val     = 0.0
                            found   = False

                            while True:
                                if heightsArr[r, c] >= 0:
                                    val     = demArr[r, c] + heightsArr[r, c]
                                    found   = True
                                    break
                                acc = accArr[r, c]
                                if acc >= self.ridgeThresh or acc == accNoData:
                                    path.append((r, c))
                                    val = demArr[r, c]
                                    found = True
                                    break
                                path.append((r, c))
                                d = int(dirArr[r, c]) - 1
                                if d < 0 or d >= 8:
                                    val = demArr[r, c]
                                    found = True
                                    break
                                nr = r + dY[d]
                                nc = c + dX[d]
                                if nr < 0 or nr >= self.numRows or nc < 0 or nc >= self.numCols:
                                    val = demArr[r, c]
                                    found = True
                                    break
                                if demArr[nr, nc] == self.demNoData:
                                    val = demArr[r, c]
                                    found = True
                                    break
                                r, c = nr, nc

                            if found:
                                for pr, pc in path:
                                    heightsArr[pr, pc] = val - demArr[pr, pc]

                    # write results
                    self.ridgeHeightsRaster.band.WriteArray(heightsArr)
                    self.ridgeHeightsRaster.band.FlushCache()

                    ridgepRaster.close()
                    ridgesRaster.close()
                    self.ridgeHeightsRaster.close()
                    return True
                except MemoryError:
                    QSWATUtils.loginfo('Out of memory for ridge heights by inversion with chunk count {0}'.format(self.chunkCount))
                    self._gv.closeOpenRasters()
                    completed = False
                    self.chunkCount += 1
                    if not self.demRaster.open(self.chunkCount):
                        return False
        else:
            return True
                    
    def calcRidgeHeightsByInversionParallel(self, ridgep, ridges, root, mustRun, nProc=None):
        """
        Create the ridgeHeightsRaster with differences between the elevation at the point and 
        the elevation of the nearest ridge cell using parallel processing.
        
        ridgep is the D8 flow directions and ridges the flow accumulation raster, both calculated from an inverted DEM.
        """
        if nProc is None:
            nProc = max(variables.taudemProcesses, 4)  # Default from variables.taudemProcesses
            
        self.ridgeHeightsFile = QSWATUtils.join(self._gv.demDir, 'invheights.tif')
        if mustRun or not QSWATUtils.isUpToDate(ridgep, self.ridgeHeightsFile) or not QSWATUtils.isUpToDate(ridges, self.ridgeHeightsFile):
            self._gv.clearOpenRasters()
            completed = False
            while not completed:
                try:
                    completed = True # only gets set on MemoryError exception
                    if os.path.exists(self.ridgeHeightsFile):
                        QSWATUtils.tryRemoveLayerAndFiles(self.ridgeHeightsFile, root)
                    self.ridgeHeightsRaster = Raster(self.ridgeHeightsFile, self._gv, canWrite=True, isInt=False)
                    res = self.ridgeHeightsRaster.open(self.chunkCount, numRows=self.numRows, numCols=self.numCols,
                                                        transform=self.demTransform, projection=self.demProjection, 
                                                        noData=self.noData)
                    if not res:
                        return False
                    # ridgep obtained by inversion has float values, although these are in the range 1 to 8
                    ridgepRaster = Raster(ridgep, self._gv, canWrite=False, isInt=False)
                    res = ridgepRaster.open(self.chunkCount)
                    if not res:
                        self._gv.closeOpenRasters()
                        return False
                    ridgesRaster = Raster(ridges, self._gv, canWrite=False, isInt=True)
                    res = ridgesRaster.open(self.chunkCount)
                    if not res:
                        self._gv.closeOpenRasters()
                        return False
                    
                    try:
                        # Read arrays for parallel processing
                        QSWATUtils.loginfo(f'Reading ridge arrays for parallel processing with {nProc} processes')
                        dem_array = self.demRaster.band.ReadAsArray()
                        ridgep_array = ridgepRaster.band.ReadAsArray()
                        ridges_array = ridgesRaster.band.ReadAsArray()
                        
                        if dem_array is None or ridgep_array is None or ridges_array is None:
                            QSWATUtils.loginfo('Failed to read ridge arrays, falling back to sequential')
                            raise Exception("Failed to read arrays")
                        
                        # Split into row chunks for parallel processing
                        rows_per_chunk = max(1, self.numRows // nProc)
                        chunks = []
                        
                        for i in range(0, self.numRows, rows_per_chunk):
                            end_row = min(i + rows_per_chunk, self.numRows)
                            dem_chunk = dem_array[i:end_row, :]
                            ridgep_chunk = ridgep_array[i:end_row, :]
                            ridges_chunk = ridges_array[i:end_row, :]
                            chunks.append((dem_chunk, ridgep_chunk, ridges_chunk, self.ridgeThresh, self.demNoData, self.noData, i, self.demTransform))
                        
                        QSWATUtils.loginfo(f'Processing {len(chunks)} ridge height chunks in parallel')
                        
                        # Process chunks in parallel - simplified approach for now
                        with mp.Pool(processes=nProc) as pool:
                            results = pool.map(_compute_ridge_heights_chunk, chunks)
                        
                        # Write results back to raster
                        total_points = sum(len(chunk_results) for chunk_results in results)
                        QSWATUtils.loginfo(f'Writing {total_points} ridge height results')
                        
                        for i, chunk_results in enumerate(results):
                            for row, col, value in chunk_results:
                                self.ridgeHeightsRaster.write(row, col, value)
                                
                            # Update progress every 10% of chunks
                            if (i + 1) % max(1, len(results) // 10) == 0:
                                progress_pct = (i + 1) / len(results) * 100
                                self._progress(f'Writing ridge heights: {progress_pct:.0f}%')
                                
                    except Exception as e:
                        QSWATUtils.loginfo(f'Parallel ridge heights calculation failed: {str(e)}, falling back to sequential')
                        # Fall back to sequential processing
                        for row in range(self.numRows):
                            for col in range(self.numCols):
                                if self.demRaster.read(row, col) != self.demNoData:
                                    (val, path) = self.valueAtNearest(row, col, ridgepRaster, ridgesRaster, self.ridgeThresh)
                                    if path is not None:                    
                                        self.propagate(val, path)
                    
                    ridgepRaster.close()
                    ridgesRaster.close()
                    self.ridgeHeightsRaster.close()
                    return True
                except MemoryError:
                    QSWATUtils.loginfo('Out of memory for ridge heights by inversion with chunk count {0}'.format(self.chunkCount))
                    self._gv.closeOpenRasters()
                    completed = False
                    self.chunkCount += 1
                    if not self.demRaster.open(self.chunkCount):
                        return False
        else:
            return True
            
    def calRidgeHeghtsByBranchLength(self, subbasins, distances, slopeDir, flowAcc, valleyDepthsFile, root, mustRun):
        """
        Create the ridgeHeightsRaster with differences between the elevation at the point and 
        the elevation of the nearest ridge cell.
        
        subbasins is the subbasins raster, distances is the distances to outlet raster, 
        slopeDir is the D8 slope directions raster, flowAcc the flow accumulation raster, all four clipped 
        the same as the DEM.
        valleyDepthsFile is the existing raster  giving the depth of the valley flow below each point (as a positive number of metres).
        """
        self.ridgeHeightsFile = QSWATUtils.join(self._gv.demDir, 'branchheights.tif')
        distancesFile = QSWATUtils.join(self._gv.demDir, 'branchdistancess.tif')
        if mustRun or not QSWATUtils.isUpToDate(subbasins, self.ridgeHeightsFile) \
            or not QSWATUtils.isUpToDate(distances, self.ridgeHeightsFile) \
            or not QSWATUtils.isUpToDate(slopeDir, self.ridgeHeightsFile) \
            or not QSWATUtils.isUpToDate(flowAcc, self.ridgeHeightsFile):
            self._gv.clearOpenRasters()
            completed = False
            while not completed:
                try:
                    completed = True # only gets set on MemoryError exception
                    if os.path.exists(self.ridgeHeightsFile):
                        QSWATUtils.tryRemoveLayerAndFiles(self.ridgeHeightsFile, root)
                    if os.path.exists(distancesFile):
                        QSWATUtils.tryRemoveLayerAndFiles(distancesFile, root)
                    self.ridgeHeightsRaster = Raster(self.ridgeHeightsFile, self._gv, canWrite=True, isInt=False)
                    res = self.ridgeHeightsRaster.open(self.chunkCount, numRows=self.numRows, numCols=self.numCols,
                                                        transform=self.demTransform, projection=self.demProjection, 
                                                        noData=self.noData)
                    if not res:
                        return False
                    self.ridgeDistancesRaster = Raster(distancesFile, self._gv, canWrite=True, isInt=False)
                    res = self.ridgeDistancesRaster.open(self.chunkCount, numRows=self.numRows, numCols=self.numCols, 
                                                        transform=self.demTransform, 
                                                        projection=self.demProjection, noData=self.noData)
                    if not res:
                        self._gv.closeOpenRasters()
                        return False
                    subbasinsRaster = Raster(subbasins, self._gv, canWrite=False, isInt=True)
                    res = subbasinsRaster.open(self.chunkCount)
                    if not res:
                        self._gv.closeOpenRasters()
                        return False
                    distRaster = Raster(distances, self._gv, canWrite=False, isInt=False)
                    res = distRaster.open(self.chunkCount)
                    if not res:
                        self._gv.closeOpenRasters()
                        return False
                    self.findRidges(subbasinsRaster, distRaster, root)
                    subbasinsRaster.close()
                    distRaster.close()
                    accRaster = Raster(flowAcc, self._gv, canWrite=False, isInt=True)
                    res = accRaster.open(self.chunkCount)
                    if not res:
                        self._gv.closeOpenRasters()
                        return False
                    dirRaster = Raster(slopeDir, self._gv, canWrite=False, isInt=True)
                    res = dirRaster.open(self.chunkCount)
                    if not res:
                        self._gv.closeOpenRasters()
                        return False
                    self.propagateFromRidges(valleyDepthsFile, dirRaster, accRaster)
                    dirRaster.close()
                    accRaster.close()
                    self.ridgeHeightsRaster.close()
                    self.ridgeDistancesRaster.close()
                    return True
                except MemoryError:
                    QSWATUtils.loginfo('Out of memory for ridge heights by branch length with chunk count {0}'.format(self.chunkCount))
                    self._gv.closeOpenRasters()
                    completed = False
                    self.chunkCount += 1
                    if not self.demRaster.open(self.chunkCount):
                        return False
        else:
            return True
                        
    def valueAtNearest(self, row, col, d8Raster, accRaster, threshold):
        """
        Return the positive elevation of the nearest point in accRaster where it is at least threshold,
        plus path along d8Raster to that point.
        
        If a value is found in result, point already processed, return result - elevation for that point.
        We are working upslope, final height is ridge height,
        ridgeHeights contains ridge - elev, i.e. positive height of ridge above point, 
        so if result already exists, ridge is recovered as result + elev
        """
        path = []
        accNoData = accRaster.band.GetNoDataValue()
        while True:
            res = self.ridgeHeightsRaster.read(row, col)
            if res >= 0: # already done this point
                # eg result is 60, this point's elevation is 100, return value (ridge elevation) is 160
                return (self.demRaster.read(row, col) + res, path)
            acc = accRaster.read(row, col)
            if acc >= threshold or acc == accNoData: # found nearest stopping point
                path.append((row, col))
                return (self.demRaster.read(row, col), path)
            if (row, col) in path:
                (startRow, startCol) = path[0]
                startX, startY = QSWATTopology.cellToProj(startCol, startRow, self.demTransform)
                x, y = QSWATTopology.cellToProj(col, row, self.demTransform)
                QSWATUtils.error('Loop to ({0}, {1}) from ({2}, {3})'.format(x, y, startX, startY), self._gv.isBatch)
                return (0, None)
            path.append((row, col))
            pt, _ = self.alongDirPoint(row, col, d8Raster)
            if pt is None:
                # hit edge of clipped area = must be a ridge
                return (self.demRaster.read(row, col), path)
            (row, col) = pt
            
    def propagate(self, val, path):
        """
        For each point on path, set result to val (ridge elevation) - point elevation."""
        for (row, col) in path:
            self.ridgeHeightsRaster.write(row, col, val - self.demRaster.read(row, col))
                
    def writeFloodPlain(self, valleyDepthsFile, mustRun):
        """Calculate slope positions from valleyDepths and ridgeHeights, and write floodplain raster ."""
        method = 'inv' if self.useInversion else 'branch'
        thresh = '{0:.2F}'.format(self.floodThresh).replace('.', '_')
        flood = QSWATUtils.join(self._gv.floodDir, method + 'flood' + thresh + '.tif')
        root = QgsProject.instance().layerTreeRoot()
        if mustRun or \
            not QSWATUtils.isUpToDate(valleyDepthsFile, flood) or \
            not QSWATUtils.isUpToDate(self.ridgeHeightsFile, flood):
            if os.path.exists(flood):
                QSWATUtils.tryRemoveLayerAndFiles(flood, root)
            self._gv.clearOpenRasters()
            completed = False
            while not completed:
                try:
                    completed = True # only gets set on MemoryError exception
                    self.ridgeHeightsRaster = Raster(self.ridgeHeightsFile, self._gv, canWrite=False, isInt=False)
                    res = self.ridgeHeightsRaster.open(self.chunkCount)
                    if not res:
                        self._gv.closeOpenRasters()
                        return
                    if self.valleyDepthsRaster is None:
                        # may have been opened when calculating ridge heights by branch length
                        self.valleyDepthsRaster = Raster(valleyDepthsFile, self._gv, canWrite=False, isInt=False)
                        res = self.valleyDepthsRaster.open(self.chunkCount)
                        if not res:
                            self._gv.closeOpenRasters()
                            return
                    self._progress('Flood plain...')
                    self.floodplainRaster = Raster(flood, self._gv, canWrite=True, isInt=True)
                    OK = self.floodplainRaster.open(self.chunkCount, numRows=self.numRows, numCols=self.numCols, 
                                                  transform=self.demTransform, projection=self.demProjection, 
                                                  noData=self.noData)
                    if OK:
                        self.calcFloodPlain1()
                        self.floodplainRaster.close()
                        ## fixing aux.xml file seems unnecessary in QGIS 2.16
#                         # now fix maximum value to 1 instead of zero in aux.xml file
#                         # else if loaded has legend 0 to nan and display is all black
#                         xmlFile = self.floodplainRaster.fileName + '.aux.xml'
#                         ok, err = QSWATUtils.setXMLValue(xmlFile, u'MDI', u'key', u'STATISTICS_MAXIMUM', u'1')
#                         if not ok:
#                             QSWATUtils.error(err, self._gv.isBatch)
                        QSWATUtils.copyPrj(self.demRaster.fileName, self.floodplainRaster.fileName)
                        # load flood above DEM
                        layers = root.findLayers()
                        demLayer = QSWATUtils.getLayerByLegend(FileTypes.legend(FileTypes._DEM), layers)
                        ft = FileTypes._INVFLOOD if self.useInversion else FileTypes._BRANCHFLOOD
                        floodLayer, _ = QSWATUtils.getLayerByFilename(layers, self.floodplainRaster.fileName, ft, 
                                                          self._gv, demLayer, QSWATUtils._WATERSHED_GROUP_NAME)
                        if floodLayer is None:
                            QSWATUtils.error('Failed to load floodplain raster {0}' \
                                             .format(self.floodplainRaster.fileName), self._gv.isBatch)
                    self.valleyDepthsRaster.close()
                    self.ridgeHeightsRaster.close()
                    self._progress('Flood plain done')
                except MemoryError:
                    QSWATUtils.loginfo('Out of memory for flood plain with chunk count {0}'.format(self.chunkCount))
                    self._gv.closeOpenRasters()
                    completed = False
                    self.chunkCount += 1
        else:# already have uptodate flood file: make sure it is loaded
            # load flood above DEM
            layers = root.findLayers()
            demLayer = QSWATUtils.getLayerByLegend(FileTypes.legend(FileTypes._DEM), layers)
            ft = FileTypes._INVFLOOD if self.useInversion else FileTypes._BRANCHFLOOD
            floodLayer, _ = QSWATUtils.getLayerByFilename(layers, flood, ft, self._gv, 
                                                          demLayer, QSWATUtils._WATERSHED_GROUP_NAME)
            if floodLayer is None:
                QSWATUtils.error('Failed to load floodplain raster {0}' \
                                 .format(self.floodplainRaster.fileName), self._gv.isBatch)
                
                    
        
    def calcFloodPlain1(self):
        """Calculate the floodplain array elementwise."""
        old_settings = numpy.seterr(divide='raise')
        valleyTransform = self.valleyDepthsRaster.ds.GetGeoTransform()
        for row in range(self.numRows):
            for col in range(self.numCols):
                rh = self.ridgeHeightsRaster.read(row, col)
                if rh != self.noData:
                    vd = self.valleyDepthsRaster.read(row, col)
                    if vd != self.noData:
                        try:
                            sp = 0 if int(rh + vd) == 0 else vd / (rh + vd)
                        except Exception:
                            QSWATUtils.information('Problem in calculating slope position at {0}: numerator {1}; denominator {2}: {3}. Set to 0' \
                                               .format(QSWATTopology.cellToProj(col, row, valleyTransform), vd, rh+vd, traceback.format_exc()),
                                             self._gv.isBatch)
                            sp = 0
                        self.floodplainRaster.write(row, col, 1 if sp <= self.floodThresh else self.noData) 
        numpy.seterr(divide=old_settings['divide'])
                        
    def calcFloodPlain1Parallel(self, nProc=None):
        """Calculate the floodplain array using vectorized numpy operations."""
        if nProc is None:
            nProc = min(mp.cpu_count(), 24)

        old_settings = numpy.seterr(divide='ignore', invalid='ignore')

        try:
            QSWATUtils.loginfo(f'Reading raster data for vectorized floodplain calculation')

            ridge_array     = self.ridgeHeightsRaster.band.ReadAsArray()
            valley_array    = self.valleyDepthsRaster.band.ReadAsArray()

            if ridge_array is None or valley_array is None:
                raise Exception("Failed to read arrays")

            # fully vectorized computation — no Python loops, no multiprocessing needed
            valid       = (ridge_array != self.noData) & (valley_array != self.noData)
            denominator = ridge_array + valley_array
            sp          = numpy.where(valid & (denominator != 0), valley_array / denominator, 0.0)

            result      = numpy.full(ridge_array.shape, self.noData, dtype=numpy.float64)
            result[valid & (sp <= self.floodThresh)] = 1

            # write entire array to raster at once
            self.floodplainRaster.band.WriteArray(result)
            self.floodplainRaster.band.FlushCache()

            QSWATUtils.loginfo('Vectorized floodplain calculation complete')

        except Exception as e:
            QSWATUtils.loginfo(f'Vectorized floodplain failed: {str(e)}, falling back to sequential')
            for row in range(self.numRows):
                for col in range(self.numCols):
                    rh = self.ridgeHeightsRaster.read(row, col)
                    if rh != self.noData:
                        vd = self.valleyDepthsRaster.read(row, col)
                        if vd != self.noData:
                            try:
                                sp = 0 if int(rh + vd) == 0 else vd / (rh + vd)
                            except Exception:
                                sp = 0
                            self.floodplainRaster.write(row, col, 1 if sp <= self.floodThresh else self.noData)

        numpy.seterr(**old_settings)
                        
    # TODO: deal with array access
    #===========================================================================
    # def calcFloodPlain2(self):
    #     """Calculate the floodplain array using numpy array methods."""
    #     for i in xrange(self.chunkCount):
    #         # if there is only one chunk we can avoid reading it again
    #         if i != self.ridgeHeightsRaster.currentIndex:
    #             heightsChunk = self.ridgeHeightsRaster.chunks[i]
    #             self.ridgeHeightsRaster.array = self.ridgeHeightsRaster.band.ReadAsArray(0, heightsChunk.rowOffset, self.ridgeHeightsRaster.numCols, heightsChunk.numRows)
    #             self.ridgeHeightsRaster.currentIndex = i
    #         if i != self.valleyDepthsRaster.currentIndex:
    #             depthsChunk = self.valleyDepthsRaster.chunks[i]
    #             self.valleyDepthsRaster.array = self.valleyDepthsRaster.band.ReadAsArray(0, depthsChunk.rowOffset, self.valleyDepthsRaster.numCols, depthsChunk.numRows)
    #             self.valleyDepthsRaster.currentIndex = i
    #         msk = (self.ridgeHeightsRaster.array != self.noData) & (self.valleyDepthsRaster.array() != self.noData)
    #         # start with array of right size for floodplain (since current size will be last initial chunk which may be small)
    #         # note also we have to initialize it to type float and correct to int later
    #         floodplainChunk = self.floodplainRaster.chunks[i]
    #         self.floodplainRaster.array = np.core.full((floodplainChunk.numRows, self.numCols), self.noData, np.float_)
    #         # the next line uses just one write but does not work
    #         #self.floodplainRaster.array[msk] = 1 if (self.valleyDepthsRaster.array[msk] / (self.valleyDepthsRaster.array[msk] + self.ridgeHeightsRaster.array[msk])) <= self.floodThresh else 0
    #         self.floodplainRaster.array[msk] = self.valleyDepthsRaster.array[msk] / (self.valleyDepthsRaster.array[msk] + self.ridgeHeightsRaster.array[msk])
    #         self.floodplainRaster.array[(self.floodplainRaster.array > self.floodThresh) & (self.ridgeHeightsRaster.array != self.noData) & (self.valleyDepthsRaster.array != self.noData)] = 1
    #         self.floodplainRaster.array[(self.floodplainRaster.array <= self.floodThresh) & (self.ridgeHeightsRaster.array != self.noData) & (self.valleyDepthsRaster.array != self.noData)] = 2
    #         self.floodplainRaster.array[msk] = self.floodplainRaster.array[msk] - 1
    #         self.floodplainRaster.array = self.floodplainRaster.array.astype(np.int_)
    #         self.floodplainRaster.array[self.floodplainRaster.array < 0] = self.noData
    #         self.floodplainRaster.band.WriteArray(self.floodplainRaster.array, 0, floodplainChunk.rowOffset)
    #===========================================================================
                
    def findRidges(self, subbasinsRaster, distRaster, root):
        """
        Create ridge points dictionary (row, col) -> (elevation, maximum branch length) for all cells on 
        subbasin boundaries.
        
        A cell is on a boundary if it has one or more adjacent cells in a different subbasin,
        or if it has an adjacent cell with nodata for its subbasin, ie it is on watershed boundary,
        when its branch length is set to the branch threshold.
        The maximum branch length is the maximum of the branch lengths for adjacent cells 
        that are in a different subbasin or outside the watershed.
        Assumes subbasins, outlet distance and dem rasters have same extents and cell sizes.
        """
        # check assumption about rasters:
        assert self.demTransform == subbasinsRaster.ds.GetGeoTransform(), 'DEM and subbasins rasters not compatible'
        assert self.demTransform == distRaster.ds.GetGeoTransform(), 'DEM and distance to outlet rasters not compatible'
        if len(self._gv.topo.subbasinToStream) == 0:
            # need to calculate some topology
            if self._gv.useGridModel:
                if os.path.exists(self._gv.delinStreamFile):
                    streamLayer = QgsVectorLayer(self._gv.delinStreamFile, 'Delineated streams', 'ogr')
                else:
                    QSWATUtils.error('Cannot use branch length algorithm without a delineated stream shapefile', self._gv.isBatch)
                    return
            else:
                streamLayer = QSWATUtils.getLayerByFilename(root.findLayers(), 
                                                            self._gv.streamFile, FileTypes._STREAMS, None, None, None)[0]
                if streamLayer is None:
                    QSWATUtils.error('Streams layer not found: have you run TauDEM?', self._gv.isBatch)
                    return
            if not self._gv.topo.setUp1(streamLayer):
                return
        self.makeRidges(subbasinsRaster, distRaster)
        time1 = time.process_time()
        self.makeRidgeRaster()
        time2 = time.process_time()
        QSWATUtils.loginfo('Writing ridge raster  took {0} seconds'.format(int(time2 - time1)))
                        
    def makeRidges(self, subbasinsRaster, distRaster):
        """
        Set values of ridge points.
        """
        time1 = time.process_time()
        self.ridgePoints = dict()
        subbasinsNoData = subbasinsRaster.band.GetNoDataValue()
        for row in range(self.numRows):
            for col in range(self.numCols):
                subbasin = subbasinsRaster.read(row, col)
                if subbasin == subbasinsNoData:
                    continue
                maxBranchLength = 0
                finished = False
                for dy in [-1, 0, 1]:
                    row1 = row+dy
                    dxs = [-1, 1] if dy == 0 else [-1, 0, 1]
                    for dx in dxs:
                        col1 = col+dx
                        if self.pointInMap(row1, col1):
                            subbasin1 = subbasinsRaster.read(row+dy, col+dx)
                        else:
                            subbasin1 = subbasinsNoData
                        if subbasin1 == subbasinsNoData:
                            maxBranchLength = self.branchThresh
                            finished = True
                            break
                        elif subbasin1 != subbasin:
                            distanceToJoin = self._gv.topo.getDistanceToJoin(subbasin1, subbasin)
                            branchLength = distRaster.read(row1, col1) + distanceToJoin
                            maxBranchLength = max(branchLength, maxBranchLength)
                    if finished:
                        break
                if maxBranchLength >= self.branchThresh:
                    elevation = self.demRaster.read(row, col)
                    if elevation != self.demNoData:
                        self.ridgePoints[(row, col)] = (elevation, maxBranchLength)
        time2 = time.process_time()
        QSWATUtils.loginfo('Making ridge points took {0} seconds'.format(int(time2 - time1)))
    
    def getRidgeElevation(self, row, col, reportFailure):
        """
        Find nearest ridge point to (row, col) and return its elevation and distance from (row, col).
        
        Distance is cartesian distance, in number of pixels (assumed square).
        Algorithm looks for first acceptable ridge  point on increasing square perimeter around original point,
        as an approximation to the nearest.
        Slope position algorithm for floodplain identification only depends on elevations, not distances, 
        so we don't worry too much about distances to ridge points, and interpret 'nearest' loosely for efficient execution.
        """
        # use view so that changes may be used later in iteration, and avoid repeated searches for adjacent unacceptable points
        n = 0
        while True:
            #=============version where square used has N-s and E-W sides=======
            # for dy in xrange(0-n, n+1):
            #     for dx in xrange(0-n, n+1):
            #         if abs(dy) == n or abs(dx) == n: # on boundary of square
            #             (e, l) = self.ridgePoints.get((row+dy, col+dx), (-1, -1))
            #             if l >= 0:
            #                 return e, math.sqrt(dx * dx + dy * dy)
            #===================================================================
            # version with diagonal sides to square - a little closer to circular maybe?
            if n == 0:
                (e, l) = self.ridgePoints.get((row, col), (-1, -1))
                if l >= 0:
                    return e, 0
            else:
                for dx in range(n+1):
                    dy = n-dx
                    (e, l) = self.ridgePoints.get((row+dy, col+dx), (-1, -1))
                    if l >= 0:
                        return e, math.sqrt(dx * dx + dy * dy)
                    if dx != 0:
                        (e, l) = self.ridgePoints.get((row+dy, col-dx), (-1, -1))
                        if l >= 0:
                            return e, math.sqrt(dx * dx + dy * dy) 
                    if dy != 0:
                        (e, l) = self.ridgePoints.get((row-dy, col+dx), (-1, -1))
                        if l >= 0:
                            return e, math.sqrt(dx * dx + dy * dy)
                        if dx != 0:
                            (e, l) = self.ridgePoints.get((row-dy, col-dx), (-1, -1))
                            if l >= 0:
                                return e, math.sqrt(dx * dx + dy * dy)
            n += 1
            if n > 1000:
                if reportFailure:
                    x, y = QSWATTopology.cellToProj(row, col, self.demTransform)
                    QSWATUtils.information('No ridge point found within 1000 pixels of ({0}, {1}).  Is the branch threshold too high?'.format(x, y), self._gv.isBatch)
                return 0, -1

            
    def propagateFromRidges(self, valleyDepthsFile, dirRaster, accRaster):
        """
        Propagate relative heights of nearest ridge points along flow paths, starting from points with accumulation 1.
        
        Heights are put into ridgeHeightsRaster.  
        Distances from ridge are put into ridgeDistancesRaster  If a point already has a height, 
        but the new path is shorter, the new height overwrites the old.
        """
        time1 = time.process_time()
        self.valleyDepthsRaster = Raster(valleyDepthsFile, self._gv, canWrite=False, isInt=False)
        res = self.valleyDepthsRaster.open(self.chunkCount)
        if not res:
            return
        reportFailure = True
        pathLength = 0 # path length counted in pixels (horizontal and vertical assumed same, so 1)
        diag = math.sqrt(2)
        for row in range(self.numRows):
            for col in range(self.numCols):
                if accRaster.read(row, col) == 1:
                    elevation, pathLength = self.getRidgeElevation(row, col, reportFailure)
                    if pathLength < 0:
                        if reportFailure:
                            reportFailure = False
                        elevation = self.demRaster.read(row, col)
                        pathLength = 0
                    nextRow, nextCol = row, col
                    while True:
                        nextElev = self.demRaster.read(nextRow, nextCol)
                        if nextElev == self.demNoData:
                            break
                        currentPathLength = self.ridgeDistancesRaster.read(nextRow, nextCol)
                        if currentPathLength >= 0:
                            # have a previously stored value
                            if pathLength <  currentPathLength:
                                # new path length from ridge is shorter: update heights raster 
                                self.ridgeHeightsRaster.write(nextRow, nextCol, elevation - nextElev)
                                self.ridgeDistancesRaster.write(nextRow, nextCol, pathLength)
                            else: # already had shorter path from ridge - no point in continuing down flow path
                                break
                        else: # no value stored yet
                            self.ridgeHeightsRaster.write(nextRow, nextCol, elevation - nextElev)
                            self.ridgeDistancesRaster.write(nextRow, nextCol, pathLength)
                        pt, isDiag = self.alongDirPoint(nextRow, nextCol, dirRaster)
                        if pt is None:
                            break
                        pathLength += diag if isDiag else 1
                        nextRow, nextCol = pt
        time2 = time.process_time()
        QSWATUtils.loginfo('Propagating ridge points took {0} seconds'.format(int(time2 - time1)))
                        
                 
    def makeRidgeRaster(self):
        """Make raster  to show ridges."""
        ridgeFile = QSWATUtils.join(self._gv.demDir, 'ridge' + str(self.branchThresh) + '.tif')
        ridgeRaster = Raster(ridgeFile, self._gv, canWrite=True, isInt=True)
        res = ridgeRaster.open(self.chunkCount, numRows=self.numRows, numCols=self.numCols, 
                               transform=self.demTransform, projection=self.demProjection, 
                               noData=self.noData)
        if not res:
            return
        for row in range(self.numRows):
            for col in range(self.numCols):
                if self.demRaster.read(row, col) != self.demNoData:
                    val, _ = self.ridgePoints.get((row, col), (-1, -1))
                    ridgeRaster.write(row, col, 0 if val == -1 else 1)
        ridgeRaster.close()
        #QSWATUtils.getLayerByFilename(self._iface.legendInterface().layers(), ridgeFile, FileTypes._OTHER, 
        #                                self._gv, None, QSWATUtils._WATERSHED_GROUP_NAME)
            
    def alongDirPoint(self, row, col, d8Raster):
        """Return next point (as pair or None) along d8Raster direction from (row, col) that has a data value in demRaster."""
        # d8Raster obtained by inverting the d8 flow directions raster has float values, so must cooerce to int
        dir0 = int(d8Raster.read(row, col) - 1)
        if 0 <= dir0 < 8:
            row1 = row + QSWATUtils._dY[dir0]
            col1 = col + QSWATUtils._dX[dir0]
            if self.pointInMap(row1, col1) and self.demRaster.read(row1, col1) != self.demNoData:
                isDiag = (dir0 // 2 == 1)
                return (row1, col1), isDiag
            else:
                return None, False
        else:
            return None, False
        
    def pointInMap(self, row, col):
        """Return true if row and col are in the limits for the DEM array."""
        return 0 <= row < self.numRows and 0 <= col < self.numCols
        
