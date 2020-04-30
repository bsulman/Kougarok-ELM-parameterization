from pylab import *
import xarray


# Meteorology comparison/correction
snap=xarray.open_dataset('../../output_data/snap_allhourly.nc',decode_times=False).isel(gridcell=0)
daymet=xarray.merge([xarray.open_dataset('../../output_data/GSWP3_TBOT_1980-2010_z16.nc',decode_times=False).isel(n=0),
                    xarray.open_dataset('../../output_data/GSWP3_PRECTmms_1980-2010_z16.nc',decode_times=False).isel(n=0),
                    xarray.open_dataset('../../output_data/GSWP3_FSDS_1980-2010_z16.nc',decode_times=False).isel(n=0),
                    ])
                    
snap_rcp85=xarray.open_dataset('../../params/snap_rcp85.nc',decode_times=False).isel(gridcell=0)

daymet_seasonal=daymet.groupby_bins(daymet['DTIME']%365,52).mean()
daymet_seasonal_std=daymet.groupby_bins(daymet['DTIME']%365,52).std()
t_center_daymet=array([0.5*(i.left+i.right) for i in daymet_seasonal['DTIME_bins'].data]) 


snap_seasonal=snap[['TBOT','PRECTmms','FSDS']].groupby_bins(snap['DTIME']%365,52).mean()
snap_seasonal_std=snap[['TBOT','PRECTmms','FSDS']].groupby_bins(snap['DTIME']%365,52).std()
t_center_snap=array([0.5*(i.left+i.right) for i in snap_seasonal['DTIME_bins'].data]) 

f,axs=subplots(2,3,num='Met comparison',clear=True)

axs[0,0].errorbar(t_center_daymet,daymet_seasonal['TBOT'],yerr=daymet_seasonal_std['TBOT'],label='Daymet')
axs[0,0].errorbar(t_center_snap,snap_seasonal['TBOT'],yerr=snap_seasonal_std['TBOT'],label='SNAP')

axs[0,0].set_title('Weekly temperature')
axs[0,0].set_xlabel('Day of year')
axs[0,0].set_ylabel('Temperature (K)')
axs[0,0].legend()

axs[0,1].errorbar(t_center_daymet,daymet_seasonal['PRECTmms'],yerr=daymet_seasonal_std['PRECTmms'])
axs[0,1].errorbar(t_center_snap,snap_seasonal['PRECTmms'],yerr=snap_seasonal_std['PRECTmms'])

axs[0,1].set_title('Weekly precip')
axs[0,1].set_xlabel('Day of year')
axs[0,1].set_ylabel('Precip (mm/s)')

axs[1,0].plot(t_center_daymet,snap_seasonal['PRECTmms'].values/daymet_seasonal['PRECTmms'].values)
axs[1,0].set_title('Ratio of SNAP to Daymet precip')
axs[1,0].set_xlabel('Day of year')
axs[1,0].set_ylabel('Precip ratio')

axs[1,1].plot(daymet_seasonal['PRECTmms'].values,snap_seasonal['PRECTmms'].values,'o')
axs[1,1].plot(linspace(daymet_seasonal['PRECTmms'].min(),daymet_seasonal['PRECTmms'].max(),10),linspace(daymet_seasonal['PRECTmms'].min(),daymet_seasonal['PRECTmms'].max(),10),'k:')
axs[1,1].set_title('SNAP vs Daymet precip')
axs[1,1].set_xlabel('Daymet precip (mm/s)')
axs[1,1].set_ylabel('SNAP precip (mm/s)')


precip_overlap=pandas.merge_asof(pandas.DataFrame(snap['PRECTmms'].data,index=snap['DTIME']/365+1970,columns=['SNAP_precip']) ,
                            pandas.DataFrame(daymet['PRECTmms'].data,index=daymet['DTIME']/365+1980,columns=['Daymet_precip']) ,
                            left_index=True,right_index=True,direction='nearest',tolerance=1e-5).dropna()
                            

snap_precip_corrected=snap['PRECTmms'].copy()
nonzero_snap=(snap_precip_corrected.to_masked_array()>5e-6).nonzero()[0]
nonzero_daymet=(daymet['PRECTmms'].to_masked_array()>5e-6).nonzero()[0]
# How many more nonzero precip points in SNAP than expected from Daymet ratio
extra_snap_nonzero = len(nonzero_snap)-len(snap['PRECTmms'])*len(nonzero_daymet)//len(daymet['PRECTmms']) 
# Random sample of that many points, to set to zero. Probably includes fewer than total points because of repetitions
x=randint(0,len(nonzero_snap),extra_snap_nonzero)
snap_precip_corrected[nonzero_snap[x]]=0.0
snap_precip_corrected[snap['TBOT']<273.15]=snap_precip_corrected[snap['TBOT']<273.15]/2

snap_precip_rcp85_corrected=snap_rcp85['PRECTmms'].copy()
nonzero_snap_rcp85=(snap_precip_rcp85_corrected.to_masked_array()>5e-6).nonzero()[0]
# How many more nonzero precip points in SNAP than expected from Daymet ratio
extra_snap_rcp85_nonzero = len(nonzero_snap_rcp85)-len(snap_rcp85['PRECTmms'])*len(nonzero_daymet)//len(daymet['PRECTmms']) 
# Random sample of that many points, to set to zero. Probably includes fewer than total points because of repetitions
x=randint(0,len(nonzero_snap_rcp85),extra_snap_rcp85_nonzero)
snap_precip_rcp85_corrected[nonzero_snap_rcp85[x]]=0.0
snap_precip_rcp85_corrected[snap_rcp85['TBOT']<273.15]=snap_precip_rcp85_corrected[snap_rcp85['TBOT']<273.15]/2



snap_corrected_seasonal=snap_precip_corrected.groupby_bins(snap_precip_corrected['DTIME']%365,52).mean()
snap_corrected_seasonal_std=snap_precip_corrected.groupby_bins(snap_precip_corrected['DTIME']%365,52).std()

axs[0,1].errorbar(t_center_snap,snap_corrected_seasonal,yerr=snap_corrected_seasonal_std)
axs[1,0].plot(t_center_daymet,snap_corrected_seasonal.values/daymet_seasonal['PRECTmms'].values)
axs[1,1].plot(daymet_seasonal['PRECTmms'].values,snap_corrected_seasonal.values,'o')

axs[0,2].errorbar(t_center_daymet,daymet_seasonal['FSDS'],yerr=daymet_seasonal_std['FSDS'],label='Daymet')
axs[0,2].errorbar(t_center_snap,snap_seasonal['FSDS'],yerr=snap_seasonal_std['FSDS'],label='SNAP')

axs[0,2].set_title('Weekly shortwave flux')
axs[0,2].set_xlabel('Day of year')
axs[0,2].set_ylabel('Shortwave flux (W/m2)')
axs[0,2].legend()