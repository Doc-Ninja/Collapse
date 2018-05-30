# Collapse

Program for the numerical integration of a 2D aAdS system

This Program uses netCDF libraries to store data on file, please download the latest version at <https://www.unidata.ucar.edu/software/netcdf/>

v 0.1.0: Fixed the numerical instabilities around the border when evolving Pi, further tests for accuracy needed. Changed the allocation in memory of KPhi and KPi arrays to not use mallocs, this seems to fix the crash.

Version testing importing of parameters from ini file, using Ini not invented here:
<https://github.com/benhoyt/inih>