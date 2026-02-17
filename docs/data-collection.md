# Data Collection Process

## Overview
The data collection step gathers required input data for the SWAT+ model setup process.

## Required Data Types
- Digital Elevation Model (DEM)
- Land Use Data
- Soil Data
- Weather Data
- GRDC Station Data
- Lake and Reservoir Data

## Data Sources
- DEM: ASTER Global DEM
- Land Use: ESA Land Cover
- Soil: FAO Soil Database
- Weather: GSWP3-EWEMBI Dataset
- Hydrology: GRDC Database
- Lakes: GRAND Database

## Configuration
Data collection settings are controlled in `datavariables.py`:
```python
redownload_dem = False
esa_landuse_year = 2011
weather_redownload = False
```

## Process Flow
1. Check existing data
2. Download missing datasets
3. Process and prepare data
4. Validate data integrity
5. Store in model-data directory

## Output Structure
```
model-data/
├── {region}/
│   ├── raster/
│   ├── shapes/
│   ├── tables/
│   └── weather/
```

## Related Steps
- [Model Initialization](initialization.md)
- [QSWAT+ Processing](qswat-processing.md)
