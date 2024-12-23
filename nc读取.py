import xarray as xr
import os

p = r"D:\data\wind_2013.nc"
data = xr.open_dataset(p)

print(data.variables)