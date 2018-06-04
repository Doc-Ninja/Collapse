# Collapse

Program for the numerical integration of a 2D aAdS system

This Program uses netCDF libraries to store data on file, please download the latest version at <https://www.unidata.ucar.edu/software/netcdf/>

<<<<<<< HEAD
CHANGELOG

v 0.2.0: Output update: added new options for the output, and brought the existing one up to standard.

v 0.1.0: Fixed the numerical instabilities around the border when evolving Pi, further tests for accuracy needed. Fixed the declaration of Arrays Kpi1-4 and Kphi1-4: this avoids allocating unused space in the heap and prevents a crash from trying to free a pointer to the stack memory.
=======
CHANGELOG:

v 0.1.0: Fixed the numerical instabilities around the border when evolving Pi, further tests for accuracy needed. Changed the allocation in memory of KPhi and KPi arrays to not use mallocs, this seems to fix the crash.
>>>>>>> dc61b7e324c760f07438da9f6df4006978c0023d
