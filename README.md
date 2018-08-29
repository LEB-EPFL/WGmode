# WGmode
2D simulation functions for slab and ridge waveguide mode study

# Folders structure
## 2Dmode_profile
Includes two main functions to simulate:
- Slab waveguide
- Ridge waveguide: "effective index" approximation
### Input parameters
- n1 = 2.038;     % core Si3N4 3.48 2.038 (glass 1.5)
- n0 = 1.5;     % cladding Sample Media 1.33
- ns = 1.5;     % substrate SiO2 1.4
- lambda = 0.647;   % um 1.5
- NAObj = 0.55;
- Fobj = 4; % [mm]
- beamD = 8; % [mm] % beack aperture obj: beam size at the back obj

## dispertionCurve
Includes two main functions:
- plotDispCurveTE
- plotDispCurveTM

## penetrationDepth
Include two main functions:
- penetraitonDepthCore: to show how the penetration depth/effective index change with wg core size
- penetrationDepth: to show how the penetration depth/effective index change with other set of parameters

# Author
Anna Archetti anna.archetti@epfl.ch

