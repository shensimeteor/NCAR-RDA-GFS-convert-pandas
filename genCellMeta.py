#!/usr/bin/env python

from netCDF4 import Dataset
import pandas as pd
import numpy as np

ncfile="../data/data2016/gfs.0p25.2016010100.f000.grib2.shen327593.nc"
ds=Dataset(ncfile)
lon=np.array(ds.variables["lon_0"])
lat=np.array(ds.variables["lat_0"])

lons,lats=np.meshgrid(lon,lat)
#print(lons)
#print(lats) # test passed

nlon=lon.shape[0]
nlat=lat.shape[0]

lons1d=np.reshape(lons, (nlon*nlat,))
lats1d=np.reshape(lats, (nlon*nlat,))
cellid=np.arange(0, nlon*nlat)

dct={}
dct["cellid"]=cellid
dct["lon"]=lons1d
dct["lat"]=lats1d
df = pd.DataFrame(dct)
df.to_csv("cellMeta.csv")
