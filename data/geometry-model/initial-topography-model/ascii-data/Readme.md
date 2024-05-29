# Provenance of data

The `global_1deg.txt` file describes topography/bathymetry data on a
global lat/long grid with 1 degree spacing. The data was downloaded
from the website
http://research.jisao.washington.edu/data_sets/elevation/, using their
1-degree latitude-longitude resolution elevation dataset, as a netcdf
file. One can open this file in Paraview and save it as a .csv
file. From there, the data was interpolated to a spherical structured
grid with 2 degrees spacing in python.
