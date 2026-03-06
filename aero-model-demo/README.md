# Aerodynamic model demo
---
This folder contains a matlab model for the computation of aerodynamic and solar radiation pressure torques on the 3cat-8 satellite. This folder has been adapted from an early-stage study work in `old-functions` to the more relevant functions for this simulator in `formatted_functions` as well as some utils in `validate_functions` for ensuring that the calculations match the simplified analytical cases.

> Caution!
> The validation folder is relatively "hard coded" and can be hard to understand and use readily.

## Generating a geometrical model for the main C++ simulator
To generate a geometric model with the right properties to calculate the aerodynamic and solar radiation pressure torques on the satellite with the appropriate format one must:
1. Go to the `formatted_functions` directory.
2. Run the `create_geometry.m` script:
    - Change the `config` parameter to account for variations in the geometry after deployments (1 = completely stowed, 2 = panel solar out, 3 = completely deployed).
3. Run the `geometry_exporter.m` code to generate a `.txt` file with the geometry data formatted in `C`.
4. Copy and paste the contents of the `.txt` file into `satellite_shape_model.hpp` at the correct location (where the current data is already palced).
5. Recompile the code following the steps from the global readme.