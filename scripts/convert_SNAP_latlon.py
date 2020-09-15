import pyproj
import xarray

def get_lonlat(data):

    crs=pyproj.CRS(data.crs_wkt)
    p=pyproj.Transformer.from_crs(crs,crs.geodetic_crs,always_xy=True)
    lon,lat=p.transform(data.xc.data,data.yc.data)

    return (lon,lat)

def find_point(lonlat,lonlat_data):
    from numpy import argmin

    lonpoint=argmin(abs(lonlat_data[0]-lonlat[0]))
    latpoint=argmin(abs(lonlat_data[1]-lonlat[1]))

    return (lonpoint,latpoint)

def extract_point(data,lon,lat):
    lonlat_data=get_lonlat(data)
    x,y=find_point((lon,lat),lonlat_data)
    return data.assign(
            long=xarray.DataArray(lonlat_data[0],coords={'xc':data['xc']},dims=('xc',),attrs={'long_name':'longitude','units':'degrees_east'}),
            lati=xarray.DataArray(lonlat_data[1],coords={'yc':data['yc']},dims=('yc',),attrs={'long_name':'latitude','units':'degrees_north'})
            ).isel(xc=x,yc=y) 

def concat_point(filenames,lon,lat):
    if isinstance(filenames,str):
        import glob
        if glob.has_magic(filenames):
            file_list=sorted(glob.glob(filenames))
        else:
            file_list=[filenames]
    else:
        file_list=filenames

    ncs_out=[]
    for f in file_list:
        print(f)
        with xarray.open_dataset(f) as d:
            ncs_out.append(extract_point(d,lon,lat))

    return xarray.concat(ncs_out,'time',data_vars='minimal')

if __name__ == '__main__':
    # Run all files for Kougarok
    Kougarok_coords=(-164.82,65.16)

    # historical
    import glob,os,sys
    dirs=glob.glob('/nfs/data/ccsi/b0u/SNAP_met/hourly/NCAR-CCSM4/historical/*')
    for vdir in dirs:
        var=os.path.basename(vdir)
        concat_point(vdir+'/'+var+'_hourly_wrf_NCAR-CCSM4_historical_*.nc',Kougarok_coords[0],Kougarok_coords[1]).to_netcdf('/home/b0u/driver_data/SNAP_Kougarok/NCAR-CCSM4/%s_historical_1970-2005.nc'%var,format='NETCDF4_CLASSIC')


    dirs=glob.glob('/nfs/data/ccsi/b0u/SNAP_met/hourly/NCAR-CCSM4/rcp85/*')
    for vdir in dirs:
        var=os.path.basename(vdir)
        concat_point(vdir+'/'+var+'_hourly_wrf_NCAR-CCSM4_rcp85_*.nc',Kougarok_coords[0],Kougarok_coords[1]).to_netcdf('/home/b0u/driver_data/SNAP_Kougarok/NCAR-CCSM4/%s_rcp85_2006-2100.nc'%var,format='NETCDF4_CLASSIC')

