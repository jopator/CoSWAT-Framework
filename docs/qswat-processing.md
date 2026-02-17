# QSWAT+ Processing

## Overview
QSWAT+ processing handles watershed delineation, HRU creation, and model parameter setup.

## Process Steps
1. DEM Processing
   - Burn-in streams
   - Fill sinks
   - Flow direction/accumulation
2. Watershed Delineation
   - Stream definition
   - Outlet definition
   - Subbasin delineation
3. HRU Definition
   - Land use overlay
   - Soil overlay
   - Slope classification
   - HRU creation

## Key Parameters
```python
thresholdSt = 150  # Stream definition threshold
thresholdCh = 150  # Channel threshold
channel_snap_thres = 3500  # Snap threshold for channels
```

## Output Files
- Stream network
- Subbasin map
- HRU definition
- Channel layout
- Reservoir locations

## Related Steps
- [Model Initialization](initialization.md)
- [Model Editing](model-editing.md)
