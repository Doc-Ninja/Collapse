# Collapse

Program for the numerical integration of a 2D aAdS system

This Program uses netCDF libraries to store data on file, please download the latest version at <https://www.unidata.ucar.edu/software/netcdf/>

CHANGELOG

v 0.3.1: Corrected a mistake in the equation for A prime, added the number of cycles elapsed in the report.txt file (created when an horizon is found)

v 0.3.0: Input update: implemented the option to load a custom input for the fields, simply rename a checkpoint or final file to start.nc and put it into the Input folder, rememeber to set the corresponding option in the Parameters.h file, if the gridsize of the input data does not match the gridsize set for the simulation, the function will silently interpolate linearly or subsample.

v 0.2.0: Output update: added new options for the output, and brought the existing one up to standard.

v 0.1.0: Fixed the numerical instabilities around the border when evolving Pi, further tests for accuracy needed. Fixed the declaration of Arrays Kpi1-4 and Kphi1-4: this avoids allocating unused space in the heap and prevents a crash from trying to free a pointer to the stack memory.
