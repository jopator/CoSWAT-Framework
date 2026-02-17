# Model Editing Process

## Overview
The model editing step configures and prepares the SWAT+ model files, including weather data setup, parameter adjustments, and final TxtInOut directory preparation.

## Process Steps
1. Database Setup
   - Initialize SQLite database
   - Configure project tables
   - Set up weather generator
   - Link to `swatplus_datasets.sqlite`

2. Weather Data Configuration
   - Import weather generator data
   - Process observed weather data
   - Setup climate stations
   - Configure weather lookup tables

3. Model Parameters
   - Land use parameters
   - Soil parameters
   - Channel parameters
   - Basin parameters
   - Management operations

4. TxtInOut Generation
   - Write model input files
   - Configure file.cio
   - Setup output variables
   - Validate file structure

## File Structure
```
model-setup/CoSWATv{version}/{region}/
├── Scenarios/
│   └── Default/
│       └── TxtInOut/
│           ├── file.cio
│           ├── weather files
│           ├── soil files
│           ├── management files
│           └── ...
└── {region}.sqlite
```

## Weather Files
The following weather files are processed:
- Precipitation (.pcp)
- Temperature (.tmp)
- Relative Humidity (.hmd)
- Solar Radiation (.slr)
- Wind Speed (.wnd)

## Configuration Requirements
```python
# Weather database settings
weather_wgn_db = '../data-preparation/resources/swatplus_wgn.sqlite'
datasets_db = '../data-preparation/resources/swatplus_datasets.sqlite'
weather_dir = '../model-data/{region}/weather/swatplus/observed'
```

## Database Tables Updated
- project_config
- weather_stations
- weather_file_paths
- parameters
- calibration
- climate

## Related Documentation
- [Data Collection](data-collection.md) - Input data preparation
- [Model Initialization](initialization.md) - Project setup
- [QSWAT+ Processing](qswat-processing.md) - Watershed processing
- [Model Execution](model-execution.md) - Running the model
- [Model Evaluation](model-evaluation.md) - Performance assessment

## Error Handling
- Weather file validation
- Database integrity checks
- Parameter range validation
- File permission checks

## Output Validation
Before proceeding to model execution:
- Check TxtInOut completeness
- Validate weather data coverage
- Verify parameter ranges
- Test database connections
