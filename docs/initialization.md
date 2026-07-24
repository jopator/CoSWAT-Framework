# Model Initialization

## Overview
The initialization process sets up the basic SWAT+ project structure and prepares input files.

## Process Steps
1. Create project directory structure
2. Copy and process input data
3. Set up QGIS project file
4. Initialize SWAT+ database
5. Configure project parameters

## Project Structure
```
model-setup/
├── CoSWATv{version}/
│   └── {region}/
│       ├── Watershed/
│       │   ├── Rasters/
│       │   └── Shapes/
│       └── Scenarios/
```

## Key Files
- Project file (.qgs)
- SWAT+ database (.sqlite)
- Configuration files
- Input rasters and shapefiles

## Configuration
Key initialization parameters in `datavariables.py`:
```python
final_proj_auth = "ESRI"
final_proj_code = 54003
data_resolution = 30.91819138974098635 * 30
```

## Related Steps
- [Data Collection](data-collection.md)
- [QSWAT+ Processing](qswat-processing.md)
