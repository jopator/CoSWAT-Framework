# Output Mapping Process

## Overview
The mapping process coordinates the visualization of CoSWAT model results by generating both vector (GPKG) and raster outputs for various hydrological variables.

## Process Flow
1. Vector Data Processing
   - Read HRU shapefile data
   - Process water balance outputs
   - Merge spatial and temporal data
   - Generate GeoPackage outputs

2. Raster Generation
   - Convert vector results to raster format
   - Apply regional projections
   - Generate multi-variable outputs

## Output Variables
### Primary Variables
- precipitation (precip)
- snow fall (snofall)
- snow melt (snomlt)
- surface runoff generation (surq_gen)
- lateral flow (latq)
- water yield (wateryld)
- percolation (perc)
- evapotranspiration (et)

### Secondary Variables
- canopy evaporation (ecanopy)
- plant evaporation (eplant)
- soil evaporation (esoil)
- curve number (cn)
- soil water parameters (sw_init, sw_final, sw_ave)
- snow parameters (sno_init, sno_final, snopack)
- potential evapotranspiration (pet)

## File Structure
```
model-outputs/
├── version-{version}/
│   ├── maps/
│   │   ├── shapefiles/
│   │   │   └── map-data.gpkg
│   │   ├── raster-maps/
│   │   │   ├── precip.tif
│   │   │   ├── snofall.tif
│   │   │   └── ...
│   │   └── map.log
```

## Configuration
Key mapping parameters in `datavariables.py`:
```python
output_re_shape = True  # Force reshaping of output
map_columns = ["precip", "snofall", "snomlt", "surq_gen", ...]
```

## Input Requirements
- HRU shapefile (hrus2.shp)
- Annual water balance file (hru_wb_aa.txt)
- Regional projection information
- Base DEM for rasterization

## Error Handling
- Missing file validation
- Region existence checks
- Projection transformation validation
- Data merge verification

## Related Documentation
- [Model Evaluation](model-evaluation.md)
- [Server Generation](server-generation.md)

## Output Usage
The generated maps can be used for:
- Visualization in GIS software
- Web interface display
- Spatial analysis
- Model result interpretation
- Documentation and reporting
