#!/usr/bin/env python
from netCDF4 import Dataset
import pandas as pd
import numpy as np
from scipy import interpolate
import glob
import os
import datetime

# return list of numpy arrays
def read_nc_by_features(ncfile, feats, ny=9, nx=13):
    try:
        ds = Dataset(ncfile)
        retlist=[]
        for feat in feats:
            varname = feat[1]
            try:
                varfield = np.squeeze(np.array(ds.variables[varname]))
                if(feat[2] >= 0 and len(varfield.shape)==3):
                    varfield = varfield[feat[2],:,:]
            except:
                print("Error %s FAIL to find var: %s, set to NaN"%(ncfile, feat[1]))
                varfield = np.empty((ny,nx), np.float)
                varfield.fill(np.nan)
            retlist.append(varfield)
    except FileNotFoundError:
        print("Error %s NOT FOUND"%ncfile)
        retlist=[]
        for feat in feats:
            varfield = np.empty((ny,nx), np.float)
            varfield.fill(np.nan)
        retlist.append(varfield)
    return retlist


def date10_addhr(date10, hour):
    dt = datetime.datetime.strptime(date10, "%Y%m%d%H")
    delta = datetime.timedelta(hours = hour)
    newdt = dt + delta
    return newdt.strftime("%Y%m%d%H")


# return [ [ .. ]] , ret[i_initdate][j_fcsthr][list of numpy_arrays]
def read_multiple_ncs(nc_dir, start_initdate, end_initdata, int_initdate, fcsthr1, fcsthr2, int_fcsthr, ncfile_pattern, feats):
    start_dt = datetime.datetime.strptime(start_initdate, "%Y%m%d%H")
    end_dt = datetime.datetime.strptime(end_initdata, "%Y%m%d%H")
    dlt = datetime.timedelta(hours = int_initdate)
    dt = start_dt
    retlist=[]
    dates=[]
    while ( dt <= end_dt):
        hr=fcsthr1
        cycle_data=[]
        while ( hr <= fcsthr2):
            ncfile = nc_dir +"/" + ncfile_pattern %(dt.strftime("%Y%m%d%H"), hr) 
            print(ncfile)
            ncarrays = read_nc_by_features(ncfile, feats)
            cycle_data.append(ncarrays)
            hr += int_fcsthr
        retlist.append(cycle_data)
        dates.append(dt)
        dt += dlt
    return retlist, dates, range(fcsthr1, fcsthr2+1, int_fcsthr)


def ncdata_interp_fcstto1hr(ncdata, fcsthrs, output_hrs):
    fcsthrs=np.array(fcsthrs)
    output_hrs=np.array(output_hrs)
    n_inhour=len(fcsthrs)
    n_outhour=len(output_hrs)
    indices=np.zeros((n_outhour,2), np.int)
    weights=np.zeros((n_outhour,2), np.float)
    for i in range(n_outhour):
        j1=np.nonzero(fcsthrs <= output_hrs[i])[0][-1]
        j2=np.nonzero(fcsthrs >= output_hrs[i])[0][0]
        if(j1==j2):
            indices[i,0]=j1
            indices[i,1]=j2
            weights[i,0]=1
            weights[i,1]=0
        else:
            indices[i,0]=j1
            indices[i,1]=j2
            weights[i,0]=(fcsthrs[j2] - output_hrs[i]) / (fcsthrs[j2] - fcsthrs[j1])
            weights[i,1]=1 - weights[i,0]
    #print(weights)
    #print(indices)  # weights/indices passed
    nvar=len(ncdata[0][0])
    for icyc in range(len(ncdata)):
        cycledata=ncdata[icyc]
        outputdata = []
        for t in range(len(output_hrs)):
            t1=indices[t,0]
            t2=indices[t,1]
            w1=weights[t,0]
            w2=weights[t,1]
            timedata=[]
            for k in range(nvar):
                try:
                    vardata = cycledata[t1][k] * w1 + cycledata[t2][k] * w2
                except:
                    print((t1,t2,k,w1,w2,len(cycledata), len(cycledata[t1]), len(cycledata[t2])))
                timedata.append(vardata)
            outputdata.append(timedata)
        ncdata[icyc] = outputdata
        #print([x[3][5,5] for x in cycledata])
        #print([x[3][5,5] for x in outputdata])  # linear interp test passed
    return 

def ncdata_to_DataFrame(ncdata, feats, initdts, fcsthrs):
    var0=ncdata[0][0][0]
    ny,nx = var0.shape
    ncells=ny*nx
    
    fcst_deltas=[]
    for hr in fcsthrs:
        fcst_deltas.append(datetime.timedelta(hours=hr))
    dates=[]
    for initdt in initdts:
        for delta in fcst_deltas:
            dates.extend( [(initdt+delta).strftime("%Y%m%d%H")] * ncells )
    ndate=len(initdts)*len(fcsthrs)
#    print(dates)
#    print((nx,ny,ncells,ndate,len(initdts),len(fcsthrs)))

    dct={}
    dct["date"]=dates
    dct["cellid"] = list(range(0, ncells)) * ndate
    cols=["cellid", "date"]
    for ifeat, feat in enumerate(feats):
        feat_field=np.zeros((ncells*ndate,), np.float)
        start_idx=0
        for idt in range(len(initdts)):
            for ihr in range(len(fcsthrs)):
                var=ncdata[idt][ihr][ifeat]
                feat_field[start_idx: start_idx+ncells] = np.reshape(var, (ncells,))
                start_idx += ncells
        dct[feat[0]] = list(feat_field)
        cols.append(feat[0])
#    for key in dct.keys(): 
#        print(key)
#        print(len(dct[key]))
    df = pd.DataFrame(dct, columns=cols)
    return df
    
if __name__ == "__main__":
    #parameters: [ (feature name, variable name in nc, variable level: -1 for 2d var, 0-based) ]
    features=[ ("Uwind_100m","UGRD_P0_L103_GLL0", 2), 
               ("Vwind_100m", "VGRD_P0_L103_GLL0", 2),
               ("Temp_50m", "TMP_P0_L104_GLL0", -1),
               ("RH_50m", "RH_P0_L104_GLL0", -1) ]
    
    nc_dir="../data/data2017/"
    #start_initdate="2016010100"
    #end_initdata=  "2016011000"
    #end_initdata=  "2017010112"
    start_initdate="2017012700"
    end_initdata="2018013118"

    int_initdate=6  # 6 hour interval
    fcsthr1=0
    fcsthr2=6
    int_fcsthr=3  # 3 hour interval
    #ncfile_pattern="gfs.0p25.%s.f%0.3d.grib2.shen327593.nc"
    ncfile_pattern="gfs.0p25.%s.f%0.3d.grib2.shen327093.nc"
    (ncdata, init_dts, fcsthrs) = read_multiple_ncs(nc_dir, start_initdate, end_initdata, int_initdate, fcsthr1, fcsthr2, int_fcsthr, ncfile_pattern, features)
    
    ncdata_interp_fcstto1hr(ncdata, fcsthrs, [0,1,2,3,4,5])
    df=ncdata_to_DataFrame(ncdata, features, init_dts, [0,1,2,3,4,5])
    #df.to_csv("test.csv")
    df.to_pickle("test_2017.pkl")
