from plot_kougarok import *

import warnings

# Don't need to print this warning a billion times
warnings.filterwarnings(action='ignore',message='All-NaN slice encountered')  
warnings.filterwarnings(action='ignore',message='Mean of empty slice')  


e3sm_all=xarray.open_dataset('../../output_data/Kougarok_SNAP_rhizomefixed_E3SMpfts_20200518_h2_processed.nc',decode_times=True)
arcticpfts_all=xarray.open_dataset('../../output_data/Kougarok_SNAP_bzdormancy_Arcticpfts_20200716_h2_processed.nc',decode_times=True) # This is with dormancy set to -1 C and root_a params fixed
# soilthickness_all=xarray.open_dataset('../../output_data/Kougarok_SNAP_rhizomefixed_E3SMpfts_soilthickness_20200422_h2_all_processed.nc',decode_times=True)
communities_all=xarray.open_dataset('../../output_data/Kougarok_SNAP_rhizomefixed_E3SMpfts_communities_20200515_h2_processed.nc',decode_times=True) # 515 changed to RCP8.5 CO2 etc

e3sm_all['PFTcolors'][:]=pft_colors_default[:17]
communities_all['PFTcolors'][:]=pft_colors_default[:17]
arcticpfts_all['PFTcolors'][:]=pft_colors

start='1990-01-01'
end='2010-01-01'

data_global=e3sm_all.sel(time=slice(start,end))
# data_soildepth=soilthickness_all.sel(time=slice(start,end))
data_communities=communities_all.sel(time=slice(start,end))
data_Arcticpfts=arcticpfts_all.sel(time=slice(start,end))


meas_leaf_C.loc[:,:,'bryophyte'] = 0.0
meas_leaf_C.loc[:,:,'lichen'] = 0.0



f=figure(num='Ecosystem biomass',figsize=(14.76, 5.26),clear=True)
f.set_constrained_layout_pads(hspace=0,h_pad=0.02)
gs=f.add_gridspec(nrows=3,ncols=6,height_ratios=[1,0.5,1])
axs=zeros((2,6),dtype='object')
offsetw=2.0
for row in range(6):
    axs[0,row]=f.add_subplot(gs[0,row])
    axs[1,row]=f.add_subplot(gs[2,row])

for simnum,data_to_plot in enumerate([data_global,data_communities,data_Arcticpfts]):
    for econum in range(len(ecotype_names_list)):
        names=[]
        x=0.22*(simnum+1)
        ax=axs[0,econum]
        sca(ax)

        offset=0
        # axvline(offset-0.5,c='k',ls=':')
        
        leafc=get_var_PFTs('LEAFC',data_to_plot)
        nonvasc=(leafc.PFTnames=='arctic_lichen')|(leafc.PFTnames=='arctic_bryophyte')

        sca(axs[1,econum])
        h=plot_mod_bar_stack(x+offset,leafc[:,:,nonvasc],econum,width=0.2,hatch=None)
        sca(ax)

        plot_mod_bar_stack(x+offset,get_var_PFTs(['LIVESTEMC','DEADSTEMC'],data_to_plot)+leafc[:,:,~nonvasc],min(econum,len(data_to_plot.ecotype)-1),width=0.2,hatch=None)

        plot_mod_bar_stack(x+offset,-get_var_PFTs(['LIVECROOTC','DEADCROOTC'],data_to_plot),min(econum,len(data_to_plot.ecotype)-1),width=0.2,hatch=None,op='min')

        
        bar(x+offset,-get_var_PFTs('FROOTC',data_to_plot).sel(ecotype=min(econum,len(data_to_plot.ecotype)-1)).max(dim='time').sum(dim='PFT'),width=0.2,
                bottom=-get_var_PFTs(['LIVECROOTC','DEADCROOTC'],data_to_plot).sel(ecotype=min(econum,len(data_to_plot.ecotype)-1)).max(dim='time').sum(dim='PFT'),
                hatch=None,facecolor=[0.9,0.9,0.9],edgecolor='k',linewidth=0.5,linestyle='--')


ticlocs=[]
names=[]
pft_order=['dwarf shrub evergreen','dwarf shrub deciduous', 'low shrub deciduous','potential tall shrub deciduous non-alder','potential tall shrub deciduous alder',
       'forb', 'graminoid']
for econum in range(len(ecotype_names_list)):
    x=0
    ax=axs[0,econum]
    sca(ax)
    offset=0#offsetw*econum
    # axvline(offset-0.5,c='k',ls=':')

    sca(axs[1,econum])
    plot_obs_bar_stack(x+offset,meas_nonvasc_C ,econum,width=0.2,pfts=['lichen','bryophyte'])
    names.append('Nonvascular')
    ticlocs.append(x+offset)
    
    sca(ax)
    # x+=1.3
    
    plot_obs_bar_stack(x+offset,meas_leaf_C+meas_stem_C ,econum,width=0.2,pfts=pft_order)
    names.append('Vascular')

    plot_obs_bar_stack(x+offset,-meas_rhizome_C,econum,width=0.2)

    
    bar(x+offset,-meas_root_C[landscape_ecotypes[econum]].mean(),yerr=meas_root_C[landscape_ecotypes[econum]].std(),width=0.2,facecolor=[0.9,0.9,0.9],edgecolor='k',linewidth=0.5,linestyle='--',
            bottom=-meas_rhizome_C[landscape_ecotypes[econum]].sum(level='PlotID').mean() )
    
    ticlocs.append(x+offset)
    
    x+=1
    
    axhline(0.0,c='k',ls=':',lw=1.0)


pfts_global=data_global['PFT'][(data_global.weights>0).any(dim='ecotype')].data.tolist()
handles=[]
for pftnum in pfts_global[1:]:
    name = prettify_pft_name(data_global.PFTnames.values[pftnum])
    if name == 'Dry Graminoid':
        name = 'Graminoid'
    handles.append(Rectangle([0,0],0,0,facecolor=data_global.PFTcolors.values[pftnum],label=name ))
l=axs[0,0].legend(handles=handles,fontsize='medium',ncol=2,title='E3SM PFTs (Sim 1)',loc=(0,-0.6))
l.set_in_layout(False)

pfts_Arctic=data_Arcticpfts['PFT'][(data_Arcticpfts.weights>0).any(dim='ecotype')].data.tolist()
handles=[]
for pftnum in pfts_Arctic[1:]:
    name = prettify_pft_name(data_Arcticpfts.PFTnames.values[pftnum])
    if name == 'Dry Graminoid':
        name = 'Graminoid'
    handles.append(Rectangle([0,0],0,0,facecolor=data_Arcticpfts.PFTcolors.values[pftnum],label=name ))
handles.append(Rectangle([0,0],0,0,facecolor=[0.9,0.9,0.9],edgecolor='k',linewidth=0.5,linestyle='--',label='Pooled fine roots' ))
l=axs[0,2].legend(handles=handles,fontsize='medium',ncol=5,title='Arctic PFTs',loc=(0,-0.6))
l.set_in_layout(False)


from string import ascii_lowercase
for econum in range(len(ecotype_names_list)):
    for row in [0,1]:
        
        axs[row,econum].set_title(landscape_ecotypes[econum])
        axs[row,econum].set_xticks(arange(4)*0.22)
        
        axs[row,econum].set_xticklabels(['Meas.','Sim 1','Sim 2','Sim 3'])

        axs[row,econum].text(0.05,1.05,'('+ascii_lowercase[econum+row*len(ecotype_names_list)]+')',transform=axs[row,econum].transAxes,fontsize='large')

        axs[1,econum].set_ylim(bottom=0)

        
axs[0,0].set_ylabel('Vascular\nbiomass (gC m$^{-2}$)')
axs[1,0].set_ylabel('Nonvascular\nbiomass (gC m$^{-2}$)')
    
f.suptitle('Biomass by plant community')





f,axs=subplots(ncols=3,nrows=1,num='PFT areas',clear=True,figsize=(11.4,5.0))
plot_PFT_distributions(axs)
axs[0].set_title('(a) Simulation 1')
axs[1].set_title('(b) Simulation 2')
axs[2].set_title('(c) Simulation 3')


barfig_PFT,axs=subplots(nrows=2,ncols=5,num='Biomass comparison by PFT',figsize=(12,6),clear=True,squeeze=False);
plot_pft_biomass_bars(data_to_plot=data_Arcticpfts,axs=axs.ravel(),per_area=False,use_pooled_default_PFTs=True,include_froot=False,include_storage=False)

# Remove PFTs with zero area
pfts_inuse=data_Arcticpfts['PFT'][(data_Arcticpfts.weights>0).any(dim='ecotype')].data.tolist()

barfig_PFT_perarea,axs=subplots(nrows=2,ncols=5,num='Biomass comparison per area of each PFT',figsize=(12,6),clear=True,squeeze=False);
plot_pft_biomass_bars(data_to_plot=data_Arcticpfts,axs=axs.ravel(),per_area=True,use_pooled_default_PFTs=False,pfts_inuse=pfts_inuse,include_storage=False)
axs[1,4].set_visible(False) 


import cartopy.crs as ccrs
import cartopy


f=figure('Site locations',clear=True,constrained_layout=False,figsize=(8,6))
gs=f.add_gridspec(nrows=2,ncols=2,width_ratios=[0.3,0.7])

ax=f.add_subplot(gs[0,0],projection=ccrs.Mercator())
ax.set_extent((-175,-137.8,50,72)) 
ax.coastlines(resolution='50m')
ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'))
ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'))

x1,x2=-168.5,-162
y1,y2=64,67

gl=ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,xlocs=arange(-175,-128,10),ylocs=arange(50,80,5),linestyle='--')
gl.xformatter=cartopy.mpl.gridliner.LONGITUDE_FORMATTER
gl.yformatter=cartopy.mpl.gridliner.LATITUDE_FORMATTER
gl.ylabels_right=False


ax=f.add_subplot(gs[1,0],projection=ccrs.Mercator())

ax.set_extent((x1,x2,y1,y2)) 
ax.coastlines(resolution='50m')
ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'))
ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'))
gl=ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,xlocs=arange(-170,-160,2),ylocs=arange(63,70,1),linestyle='--')
ax.plot(-165.418,64.508,'y*',transform=ccrs.PlateCarree(),ms=15.0,markeredgecolor='k')
ax.text(-165.418,64.508+(y2-y1)/30,'Nome',transform=ccrs.PlateCarree(),fontsize='large',ha='center')
gl.xformatter=cartopy.mpl.gridliner.LONGITUDE_FORMATTER
gl.yformatter=cartopy.mpl.gridliner.LATITUDE_FORMATTER
gl.ylabels_right=False
x1,x2=-164.840,-164.8001
y1,y2=65.153,65.17422

f.axes[0].indicate_inset_zoom(f.axes[1],edgecolor='0.3')

plotdata=pandas.read_csv('../../Amy_data/ngee_arctic_kougarok_2016_veg_comp_env_table_v1_20180828.csv',header=0,skiprows=[1,2])
plot_names=plotdata['plot_number'].dropna()
plot_lat=plotdata['local_latitude'].dropna()
plot_lon=plotdata['local_longitude'].dropna()

crs_site=ccrs.epsg(3338)
ax=f.add_subplot(gs[:,1],projection=crs_site)
ax.set_extent((x1,x2,y1,y2))  


# Applied this command to subset the DEM data: 
# gdal_translate -strict -projwin -164.9 65.18 -164.75 65.14 -projwin_srs EPSG:4269 -epo dds4/ifsar/dtm/DTM_n6500w16500P/DTM_n6500w16500P.tif dem_subset.tif

def gdal_command(x1y1,x2y2,dataset,newname):
    import os
    x1,y1=x1y1
    x2,y2=x2y2
    cmd='gdal_translate -strict -projwin {ulx:1.4f} {uly:1.4f} {lrx:1.4f} {lry:1.4f} -projwin_srs EPSG:4269 -epo {dataset:s} {newname:s}'.format(
            ulx=min(x1,x2),lrx=max(x1,x2),uly=max(y1,y2),lry=min(y1,y2),dataset=dataset,newname=newname)
    return cmd
    
elevdata=xarray.open_rasterio('../../dem_subset.tif').isel(band=0,drop=True)
import pyproj
# Transform from 3338 (NAD83 / Alaska Albers) to 4269 (NAD83 geodetic)
t=pyproj.Transformer.from_crs(3338,4269)
# In m coordinates:
x1_m,y1_m=t.transform(y1,x1,direction='INVERSE')
x2_m,y2_m=t.transform(y2,x2,direction='INVERSE')


x1_merc,y1_merc=f.axes[1].projection.transform_point(x1,y1,ccrs.PlateCarree())
x2_merc,y2_merc=f.axes[1].projection.transform_point(x2,y2,ccrs.PlateCarree())
f.axes[1].indicate_inset([x1_merc,y1_merc,x2_merc-x1_merc,y2_merc-y1_merc],f.axes[2],edgecolor='0.3') 

xx,yy=meshgrid(elevdata.x,elevdata.y) 
lat,lon=t.transform(xx,yy)

cs=ax.contour(elevdata.x,elevdata.y,elevdata.data,levels=arange(0,120,20),colors='w')
clabels=cs.clabel(fmt='%1.0f',manual=((178.9+x1_m,1314.3+y1_m),(600+x1_m,1921+y1_m),(911+x1_m,2004+y1_m),(1111+x1_m,1421+y1_m),(2035+x1_m,1001+y1_m),(767+x1_m,1377+y1_m)))
markers=['o','s','D','^','>','H']
for plotnum in range(len(plot_lat)):
    old_new_ecos={'NAMC':'DLST','TT-WBT':'ASV','DSLT':'BEL'}
    eco=plot_names[plotnum].split('_')[1][:-1]
    n=plot_names[plotnum].split('_')[1][-1]
    c='C'+str(landscape_ecotypes.index(old_new_ecos.get(eco,eco)))
    if c=='C2':
        c='y'
    x,y=t.transform(plot_lat[plotnum],plot_lon[plotnum],direction='INVERSE')
    
    h=ax.scatter(plot_lon[plotnum],plot_lat[plotnum],transform=ccrs.PlateCarree(),marker=markers[landscape_ecotypes.index(old_new_ecos.get(eco,eco))],facecolor=c,edgecolor='w',linewidth=0.4,s=30.0,zorder=10)
    if n=='1':
        h.set_label(old_new_ecos.get(eco,eco))

# gdal_translate -strict -projwin -164.9 65.18 -164.75 65.14 -projwin_srs EPSG:4269 -epo ../../2016_oneyear_layer_subset_good_metrics_ver16m1_3.tif modis_subset.tif
# modis=xarray.open_rasterio('../../modis_subset.tif')
# t=pyproj.Transformer.from_proj(modis.crs,'epsg:3338')
# xx,yy=meshgrid(modis.x,modis.y)
# x,y=t.transform(xx,yy)
# ax.pcolormesh(x-x1_m,y-y1_m,modis.sel(band=11))
# ax.set_xlabel('Easting (m)')
# ax.set_ylabel('Northing (m)')
# ax.set_aspect(1.0)
# 
# landsat_r=xarray.open_rasterio('../../landsat_band2_red_subset.tif')
# landsat_g=xarray.open_rasterio('../../landsat_band3_green_subset.tif')
# landsat_b=xarray.open_rasterio('../../landsat_band4_blue_subset.tif')
# rgb=xarray.concat((landsat_r,landsat_g,landsat_b),dim='band').transpose('y','x','band')
# rgb_corr=(rgb.where(rgb!=rgb.nodatavals[0])*0.0001)
# 
# t=pyproj.Transformer.from_proj(rgb.crs.strip('+init='),elevdata.crs.strip('+init=')) 
# xx,yy=meshgrid(rgb_corr.x,rgb_corr.y)
# x,y=t.transform(xx,yy)
# # clipval=rgb_corr.quantile(0.98).item()
# clipval=0.1
# 
# sentinel_rgb=xarray.open_rasterio('../../sentinel_subset.tif') 
# t_sentinel=pyproj.Transformer.from_proj(sentinel_rgb.crs.strip('+init='),elevdata.crs.strip('+init=')) 
# xx,yy=meshgrid(sentinel_rgb.x,sentinel_rgb.y)
# x_sentinel,y_sentinel=t_sentinel.transform(xx,yy)
# clipval_sentinel=130

mosaic=xarray.open_rasterio('../../Kougarok_mosaic.tif')

im=ax.imshow(mosaic.transpose('y','x','band'),extent=(mosaic.x[0],mosaic.x[-1],mosaic.y[0],mosaic.y[-1]),origin='lower',transform=ccrs.PlateCarree(),interpolation='hermite',regrid_shape=2500)
ax.set_xlabel('Easting (m)')
ax.set_ylabel('Northing (m)')
ax.set_aspect(1.0)


ax.legend(markerscale=1.2,fontsize='medium')

# Scale bar
xmin,xmax=ax.get_xlim()
ymin,ymax=ax.get_ylim()

ax.annotate('',xy=(xmin+(xmax-xmin)*0.1,ymin+(ymax-ymin)*0.1),xycoords='data',xytext=(xmin+(xmax-xmin)*0.1+1000,ymin+(ymax-ymin)*0.1),arrowprops={'arrowstyle':'|-|','color':'k'})

# North arrow
x0,y0=(xmin+(xmax-xmin)*0.9,ymin+(ymax-ymin)*0.1)
lt0,ln0=t.transform(x0,y0)
x1,y1=(xmin+(xmax-xmin)*0.9,ymin+(ymax-ymin)*0.1+250)
lt1,ln1=t.transform(x1,y1)
dlt=sqrt((ln1-ln0)**2+(lt1-lt0)**2)/3
ln2,lt2=ln1,lt1+dlt
x2,y2=t.transform(lt2,ln2,direction='INVERSE')
ax.annotate('',xytext=(x0,y0),xycoords='data',xy=(x2,y2),transform=ccrs.PlateCarree(),
    arrowprops={'arrowstyle':'fancy','color':'w'},color='w')
ax.text(x2,y2,'N',color='w',weight='bold',fontsize='large',ha='center')

ax.text(xmin+(xmax-xmin)*0.1+500,ymin+(ymax-ymin)*0.1+20,'1 km',ha='center')

f.axes[0].text(0.05,0.94,'(a)',transform=f.axes[0].transAxes)
f.axes[1].text(0.05,0.94,'(b)',transform=f.axes[1].transAxes)
f.axes[2].text(0.04,0.96,'(c)',transform=f.axes[2].transAxes,c='k')

gs.wspace=0.3
gs.right=0.98
gs.top=0.95
gs.bottom=0.04
gs.left=0.05
gs.update()



def plot_timeseries(d,axs,plot_sum=False,plot_PFTs=True,resample='1Y',resample_op='mean',ymin=-150,ymax=4600,force_all_panels=False,leg_ax=0,leg_ncol=2,**kwargs):
    econums=arange(len(d.ecotype),dtype=int)
    if force_all_panels and len(econums)==1:
        econums=zeros(len(axs.ravel()),dtype=int)
    for plotnum in range(len(econums)):
        ax=axs.ravel()[plotnum]
        econum=econums[plotnum]
        if resample is not None:
            data=getattr(d.resample(time=resample),resample_op)()
        else:
            data=d
        handles=[]
        if plot_PFTs:
            if len(data.PFT)==12:
                pfts=range(1,11)
            else:
                pfts=data.PFT[data.max(dim='time').isel(ecotype=econum)>0].values
            for pft in pfts:
                if data.isel(PFT=pft).dropna(dim='ecotype').any():
                    handles.append(ax.plot(data['time'],data.isel(PFT=pft,ecotype=econum),c=data['PFTcolors'][pft].item(),label=prettify_pft_name(data['PFTnames'][pft].item()),**kwargs)[0])
        if plot_sum:
            h=ax.plot(data['time'],data.isel(ecotype=econum).sum(dim='PFT'),**kwargs)[0]
            if plot_PFTs:
                h.set(linestyle='--',c='k',label='Sum over PFTs')
            handles.append(h)
        
        ax.set_title(landscape_ecotypes[econum])
        ax.set_xlabel('Year')
        ax.set_ylabel('Biomass (g C m$^{-2}$)')
        ax.axvline(2016,ls=':',lw=0.5,c='k')
        ax.set_ylim(ymin,ymax)
        if plotnum==leg_ax:
            ax.legend(handles=handles,ncol=leg_ncol)


f,axs=subplots(nrows=2,ncols=3,num='Total biomass (Arctic PFTs)',clear=True,figsize=(12,6))
f.set_constrained_layout_pads(hspace=0.18)
plot_timeseries(get_var_PFTs('TOTVEGC', arcticpfts_all),axs,plot_sum=False,ymax=None,ymin=-10,leg_ax=None)
l=axs[0,0].legend(loc=(0,-0.4),ncol=5)
l.set_in_layout(False)

    
f,axs=subplots(nrows=2,ncols=3,num='Total biomass (communities)',clear=True,figsize=(12,6))
plot_timeseries(get_var_PFTs('TOTVEGC', communities_all),axs,True)
    

f,axs=subplots(1,1,num='E3SM PFTs',clear=True,squeeze=False)
plot_timeseries(get_var_PFTs('TOTVEGC', e3sm_all),axs,True)
axs[0,0].set_title('E3SM original PFTs')

# Total biomass from three simulations compared
f,axs=subplots(nrows=2,ncols=3,num='Total biomass (sim comp)',clear=True,figsize=(12,6))
plot_timeseries(get_var_PFTs('TOTVEGC', e3sm_all),axs,plot_sum=True,plot_PFTs=False,c='navy',ls=':',label='Simulation 1',ymax=6300,force_all_panels=True)
plot_timeseries(get_var_PFTs('TOTVEGC', communities_all),axs,plot_sum=True,plot_PFTs=False,c='purple',ls='--',label='Simulation 2',ymax=6300)
plot_timeseries(get_var_PFTs('TOTVEGC', arcticpfts_all),axs,plot_sum=True,plot_PFTs=False,c='k',ls='-',label='Simulation 3',ymax=6300)
axs[0,0].legend()

# Aboveground biomass
obs_ag_C = (meas_leaf_C + meas_stem_C).add(meas_nonvasc_C,fill_value=0.0)
obs_bg_C = meas_rhizome_C.groupby(level=('Ecotype','PlotID')).sum() + meas_root_C
ag_mod_C_arctic = get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC'],data_Arcticpfts)
bg_mod_C_arctic = get_var_PFTs(['LIVECROOTC','DEADCROOTC','FROOTC'],data_Arcticpfts)
ag_mod_C_global = get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC'],data_communities)
bg_mod_C_global = get_var_PFTs(['LIVECROOTC','DEADCROOTC','FROOTC'],data_communities)

obs=(meas_leaf_C + meas_stem_C + meas_rhizome_C).add(meas_nonvasc_C,fill_value=0.0)
mod_new=get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC','LIVECROOTC','DEADCROOTC'],data_Arcticpfts)
mod_soildepth=get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC','LIVECROOTC','DEADCROOTC'],data_communities)

markers=['o','s','x','^','>','*']

f,axs=subplots(ncols=3,nrows=2,num='Model-data comparison by community',clear=True,squeeze=False)
ecotypes_included=range(6)#[4,5]
for econum in ecotypes_included:#range(len(ecotype_names_list)):
    ax=axs.ravel()[econum]
    for pft in unique(obs.index.get_level_values('ELM_PFT')):
        if pft in obs[landscape_ecotypes[econum]].index.get_level_values('ELM_PFT') :
            x=obs[landscape_ecotypes[econum]][:,pft].mean()
            xerr=obs[landscape_ecotypes[econum]][:,pft].std()
            pftnum=list(data_Arcticpfts.PFTnames.to_masked_array()).index(obsdata_PFT_mappings[pft])
            y=mod_new.isel(ecotype=econum,PFT=pftnum).max(dim='time')
            if x>0 or y>0:
                ax.errorbar(x,y,c=data_Arcticpfts.PFTcolors[pftnum].item(),ls='None',marker='o',xerr=xerr,ms=5.0)
                
                x_old=obs_to_defaultPFTs(obs,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,obsdata_E3SMPFT_mappings[pft]].mean()
                if obsdata_E3SMPFT_mappings[pft] in pft_names_default:
                    y_old=mod_soildepth.isel(ecotype=min(econum,len(ag_mod_C_global.ecotype)-1),PFT=pft_names_default.index(obsdata_E3SMPFT_mappings[pft])).max(dim='time')
                else:
                    y_old=0.0
                ax.annotate('',(x,y),(x_old,y_old),arrowprops=dict(arrowstyle="->",color=data_Arcticpfts.PFTcolors[pftnum].item(),linestyle='--',alpha=0.5))
                
            
    for pft in [11,12]:
        x=obs_to_defaultPFTs(obs,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,data_global.PFTnames[pft].item()].mean()
        xerr=obs_to_defaultPFTs(obs,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,data_global.PFTnames[pft].item()].std()
        y=mod_soildepth.isel(ecotype=min(econum,len(mod_soildepth.ecotype)-1),PFT=pft).max(dim='time')
        if x>0 or y>0:
            ax.errorbar(x,y,c=data_global.PFTcolors[pft].item(),ls='None',marker='o',xerr=xerr,ms=8.0,mfc='None')
    x=obs_to_defaultPFTs(obs,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,'nonvascular'].mean()
    y=0.0
    if x>0:
        ax.errorbar(x,y,c=data_Arcticpfts.PFTcolors[1].item(),ls='None',marker='o',xerr=xerr,ms=8.0,mfc='None',alpha=0.5)

    ax.plot(linspace(0,800,10),linspace(0,800,10),'k:')
    ax.plot(linspace(0,800,10),linspace(0,800,10),'k:')
    ax.set_title(ecotype_names_list[econum],fontsize='large')
    ax.set_xlabel('Observed biomass (g m$^{-2}$)',fontsize='large')
    ax.set_ylabel('Modeled biomass (g m$^{-2}$)',fontsize='large')

axs[0,0].set_xlim(-20,305)
axs[0,0].set_ylim(-20,305)
axs[0,1].set_xlim(-20,405)
axs[0,1].set_ylim(-20,405)
handles=[Line2D([0,0],[0,0],ls='None',marker='o',ms=8.0,c=data_Arcticpfts.PFTcolors[pft].item()) for pft in range(1,11)]
labels=[prettify_pft_name(data_Arcticpfts.PFTnames[pft].item()) for pft in range(1,11)]

axs[0,0].legend(handles=handles,labels=labels,fontsize='medium')


f,axs=subplots(ncols=2,nrows=1,num='Model-data comparison by pft type',clear=True,squeeze=True)
ecotypes_included=range(6)#[4,5]
for econum in ecotypes_included:#range(len(ecotype_names_list)):
    for pft in unique(obs.index.get_level_values('ELM_PFT')):
        if pft in obs[landscape_ecotypes[econum]].index.get_level_values('ELM_PFT'):
            if 'shrub' in pft:
                ax=axs[1]
            else:
                ax=axs[0]
            x=obs[landscape_ecotypes[econum]][:,pft].mean()
            xerr=obs[landscape_ecotypes[econum]][:,pft].std()
            pftnum=list(data_Arcticpfts.PFTnames.to_masked_array()).index(obsdata_PFT_mappings[pft])
            y=mod_new.isel(ecotype=econum,PFT=pftnum).max(dim='time')
            if x>0 or y>0:
                ax.errorbar(x,y,c=data_Arcticpfts.PFTcolors[pftnum].item(),ls='None',marker=markers[econum],xerr=xerr,ms=5.0)
                
                x_old=obs_to_defaultPFTs(obs,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,obsdata_E3SMPFT_mappings[pft]].mean()
                if obsdata_E3SMPFT_mappings[pft] in pft_names_default:
                    y_old=mod_soildepth.isel(ecotype=min(econum,len(ag_mod_C_global.ecotype)-1),PFT=pft_names_default.index(obsdata_E3SMPFT_mappings[pft])).max(dim='time')
                else:
                    y_old=0.0
                ax.annotate('',(x,y),(x_old,y_old),arrowprops=dict(arrowstyle="->",color=data_Arcticpfts.PFTcolors[pftnum].item(),alpha=0.5,linestyle='--'))
                
            

    for pft in [11,12]:
        if 'shrub' in mod_soildepth['PFTnames'][pft].item():
            ax=axs[1]
        else:
            ax=axs[0]
        x=obs_to_defaultPFTs(obs,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,data_global.PFTnames[pft].item()].mean()
        xerr=obs_to_defaultPFTs(obs,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,data_global.PFTnames[pft].item()].std()
        y=mod_soildepth.isel(ecotype=min(econum,len(mod_soildepth.ecotype)-1),PFT=pft).max(dim='time')
        if x>0 or y>0:
            ax.errorbar(x,y,c=data_global.PFTcolors[pft].item(),ls='None',marker=markers[econum],xerr=xerr,ms=8.0,mfc='None')
    x=obs_to_defaultPFTs(obs,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,'nonvascular'].mean()
    y=0.0
    if x>0:
        axs[0].errorbar(x,y,c=data_Arcticpfts.PFTcolors[1].item(),ls='None',marker=markers[econum],xerr=xerr,ms=8.0,alpha=0.5,mfc='None')

axs[1].plot(linspace(0,5000,10),linspace(0,5000,10),'k:')
axs[0].plot(linspace(0,280,10),linspace(0,280,10),'k:')
axs[0].set_title('Graminoid, forb, and nonvascular',fontsize='large')
axs[0].set_xlabel('Observed biomass (g m$^{-2}$)',fontsize='large')
axs[0].set_ylabel('Modeled biomass (g m$^{-2}$)',fontsize='large')
axs[1].set_title('Shrub',fontsize='large')
axs[1].set_xlabel('Observed biomass (g m$^{-2}$)',fontsize='large')
axs[1].set_ylabel('Modeled biomass (g m$^{-2}$)',fontsize='large')

handles=[Line2D([0,0],[0,0],ls='None',marker='o',ms=8.0,c=data_Arcticpfts.PFTcolors[pft].item()) for pft in range(1,11)]
labels=[prettify_pft_name(data_Arcticpfts.PFTnames[pft].item()) for pft in range(1,11)]
handles=handles+[Line2D([0,0],[0,0],ls='None',marker=markers[m],c='k',ms=8.0,mfc='None') for m in ecotypes_included]
labels.extend(ecotype_names_list[e] for e in ecotypes_included)
axs[0].legend(handles=handles,labels=labels,fontsize='medium')


f=figure('PFT comp',clear=True,figsize=(15,6))
gs=f.add_gridspec(ncols=3,nrows=4,height_ratios=[0.1,0.8,0.2,1],width_ratios=[1.5,1,1])
axs={}
axs['broadleaf_deciduous_boreal_shrub']=f.add_subplot(gs[1,0])
axs['c3_arctic_grass']=f.add_subplot(gs[1,1])
axs['nonvascular']=f.add_subplot(gs[1,2])
for s in axs['nonvascular'].spines.values():
    s.set_linestyle(':')


gs_new=gs[3,:].subgridspec(ncols=6,nrows=1,width_ratios=[1.5,1.5,1,1,1,1])
gs_shrub=gs_new[0:2].subgridspec(ncols=4,nrows=2,width_ratios=[1,1,1,.1])
axs['arctic_evergreen_shrub_dwarf']=f.add_subplot(gs_shrub[0,0])
axs['arctic_deciduous_shrub_dwarf']=f.add_subplot(gs_shrub[0,1])
axs['arctic_deciduous_shrub_low']=f.add_subplot(gs_shrub[0,2])
axs['arctic_deciduous_shrub_tall']=f.add_subplot(gs_shrub[1,0])
axs['arctic_deciduous_shrub_alder']=f.add_subplot(gs_shrub[1,1])

axs['arctic_dry_graminoid']=f.add_subplot(gs_new[0,2])
axs['arctic_forb']=f.add_subplot(gs_new[0,3])

axs['arctic_bryophyte']=f.add_subplot(gs_new[0,4])
axs['arctic_lichen']=f.add_subplot(gs_new[0,5])

plot_pft_biomass_bars(axs=[axs['nonvascular'],axs[data_communities.PFTnames[11].item()],axs[data_communities.PFTnames[12].item()]],
                    data_to_plot=data_communities,per_area=False,include_storage=False,include_froot=False,use_pooled_default_PFTs=True,leg_axis=1,pfts_inuse=[0,11,12])
axs['broadleaf_deciduous_boreal_shrub'].set_title('Arctic Deciduous Shrubs')
axs['c3_arctic_grass'].set_title('Arctic grasses')
axs['nonvascular'].set_title('Nonvascular')

pfts_inuse=data_Arcticpfts['PFT'][(data_Arcticpfts.weights>0).any(dim='ecotype')].data.tolist()
plot_pft_biomass_bars(axs=[axs[data_Arcticpfts.PFTnames[pft].item()] for pft in pfts_inuse[1:]],data_to_plot=data_Arcticpfts,per_area=False,include_storage=False,include_froot=False,use_pooled_default_PFTs=True,leg_axis=None,pfts_inuse=pfts_inuse)


f.add_artist(Line2D([0.42,0.42],[0.0,0.95],ls='--',c='k')) 
f.add_artist(Line2D([0.71,0.71],[0.0,0.95],ls='--',c='k')) 
f.add_artist(matplotlib.patches.FancyArrow(0.23,0.56,0.0,-0.06,width=0.01,length_includes_head=True,head_width=0.025))
f.add_artist(matplotlib.patches.FancyArrow(0.56,0.56,0.0,-0.06,width=0.01,length_includes_head=True,head_width=0.025))
f.add_artist(matplotlib.patches.FancyArrow(0.86,0.56,0.0,-0.06,width=0.01,length_includes_head=True,head_width=0.025))



f=figure('PFT comp 2',clear=True,figsize=(10,9.3))
gs=f.add_gridspec(ncols=2,nrows=3,height_ratios=[0.7,0.1,1],left=0.03,right=0.97)
axs={}
axs['broadleaf_deciduous_boreal_shrub']=f.add_subplot(gs[0,0])
axs['c3_arctic_grass']=f.add_subplot(gs[0,1])



gs_new=gs[2,:].subgridspec(ncols=7,nrows=2,width_ratios=[1,10,10,1,10,10,1])

axs['arctic_evergreen_shrub_dwarf']=f.add_subplot(gs_new[0,1])

axs['arctic_deciduous_shrub_low']=f.add_subplot(gs_new[0,2])
axs['arctic_deciduous_shrub_tall']=f.add_subplot(gs_new[1,1])
axs['arctic_deciduous_shrub_alder']=f.add_subplot(gs_new[1,2])

axs['arctic_dry_graminoid']=f.add_subplot(gs_new[0,4])
axs['arctic_forb']=f.add_subplot(gs_new[0,5])

axs['arctic_bryophyte']=f.add_subplot(gs_new[1,4])
axs['arctic_lichen']=f.add_subplot(gs_new[1,5])

axs['broadleaf_deciduous_boreal_shrub'].set_title('Arctic Deciduous Shrubs')
axs['c3_arctic_grass'].set_title('Arctic grasses')

pfts_inuse=data_Arcticpfts['PFT'][(data_Arcticpfts.weights>0).any(dim='ecotype')].data.tolist()

a1=f.add_artist(matplotlib.patches.FancyArrow(0.25,0.595,0.0,-0.06,width=0.02,length_includes_head=True,head_width=0.05,head_length=0.03,facecolor='k'))
a2=f.add_artist(matplotlib.patches.FancyArrow(0.75,0.595,0.0,-0.06,width=0.02,length_includes_head=True,head_width=0.05,head_length=0.03,facecolor='k'))

r1=f.add_artist(Rectangle(xy=(0.045,0.59),width=0.4,height=0.38,zorder=-10,facecolor='0.9'))
r2=f.add_artist(Rectangle(xy=(0.55,0.59),width=0.4,height=0.38,zorder=-10,facecolor='0.9'))
r3=f.add_artist(Rectangle(xy=(0.02,0.01),width=0.50,height=0.53,zorder=-10,facecolor='0.9'))
r4=f.add_artist(Rectangle(xy=(0.53,0.01),width=.98-.53,height=0.53,zorder=-10,facecolor='0.9'))
r5=f.add_artist(Rectangle(xy=(0.03,0.58),width=0.93,height=0.999-0.58,facecolor='0.95',zorder=-11))
r6=f.add_artist(Rectangle(xy=(0.00,0.0),width=1,height=0.57,facecolor='0.95',zorder=-11))

t=f.text(0.5,0.55,'Simulation 3',fontsize='x-large',ha='center')    
t2=f.suptitle('Simulation 2',fontsize='x-large')

ecotypes_included=range(6)#[4,5]
for econum in ecotypes_included:#range(len(ecotype_names_list)):
    for pft in unique(obs.index.get_level_values('ELM_PFT')):
        if obsdata_PFT_mappings[pft] not in axs or 'birch' in pft or 'willow' in pft:
            continue
        if pft in obs[landscape_ecotypes[econum]].index.get_level_values('ELM_PFT'):
            ax=axs[obsdata_PFT_mappings[pft]]
            x=obs[landscape_ecotypes[econum]][:,pft].mean()
            xerr=obs[landscape_ecotypes[econum]][:,pft].std(ddof=0)
            pftnum=list(data_Arcticpfts.PFTnames.to_masked_array()).index(obsdata_PFT_mappings[pft])
            y=mod_new.isel(ecotype=econum,PFT=pftnum).max(dim='time')
            if x>0 or y>0:
                ax.errorbar(x,y,c=data_Arcticpfts.PFTcolors[pftnum].item(),ls='None',marker=markers[econum],xerr=xerr,ms=7.0,mfc='None')


    for pft in [11,12]:
        if 'shrub' in mod_soildepth['PFTnames'][pft].item():
            ax=axs['broadleaf_deciduous_boreal_shrub']
        else:
            ax=axs['c3_arctic_grass']
        x=obs_to_defaultPFTs(obs,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,data_global.PFTnames[pft].item()].mean()
        xerr=obs_to_defaultPFTs(obs,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,data_global.PFTnames[pft].item()].std(ddof=0)
        y=mod_soildepth.isel(ecotype=min(econum,len(mod_soildepth.ecotype)-1),PFT=pft).max(dim='time')
        if x>0 or y>0:
            ax.errorbar(x,y,c=data_global.PFTcolors[pft].item(),ls='None',marker=markers[econum],xerr=xerr,ms=8.0,mfc='None')
    x=obs_to_defaultPFTs(obs,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,'nonvascular'].mean()
    y=0.0
    if x>0:
        ax.errorbar(x,y,c=data_Arcticpfts.PFTcolors[1].item(),ls='None',marker=markers[econum],xerr=xerr,ms=8.0,alpha=0.5,mfc='None')

for num,pft in enumerate(axs):
    axs[pft].set(title=prettify_pft_name(pft),xlabel='Obs Biomass (gC m$^{-2}$)',ylabel='Mod Biomass (gC m$^{-2}$)',aspect=1.0)
    maxval=max(axs[pft].get_ylim()[1],axs[pft].get_xlim()[1])
    axs[pft].plot([0,maxval],[0,maxval],'k:',lw=0.5)
    axs[pft].text(0.02,0.92,'('+ascii_lowercase[num]+')',transform=axs[pft].transAxes)
    
# Calculate R2 for each plot. Reading the data from the plots themselves since the plotting order was not by PFT
# linestyle and alpha business is to exclude 1-1 line and the pooled nonvascular points in Sim2 grass panel
for pft in axs:
    xdata=ma.masked_invalid([l.get_xdata() for l in axs[pft].lines if l.get_linestyle() is 'None' and l.get_alpha() is None])
    ydata=ma.masked_invalid([l.get_ydata() for l in axs[pft].lines if l.get_linestyle() is 'None' and l.get_alpha() is None])
    data=row_stack((xdata.ravel(),ydata.ravel()))
    print('%s: R2 = %1.2g, N = %d'%(pft,corrcoef(data)[0,1]**2,data.shape[1]))

handles=[Line2D([0,0],[0,0],ls='None',marker=markers[m],c='k',ms=8.0,mfc='None') for m in ecotypes_included]
labels=[landscape_ecotypes[e] for e in ecotypes_included]
axs['broadleaf_deciduous_boreal_shrub'].legend(handles=handles,labels=labels,fontsize='medium',loc='upper right')



bins=arange(0,6400,300)

def plot_hist(data,ax,color=None,label=None,bottom=None,bins=bins,alpha=0.9,norm=1.0,smooth='spline',show_points=False,**kwargs):
    x=(bins[0:-1]+bins[1:])/2
    h,b=histogram(data[data>0],bins=bins)
    if color is None:
        color=data.PFTcolors.item()
    if label is None:
        label=prettify_pft_name(data.PFTnames.item())
    if bottom is None:
        bottom=zeros(len(x))
    c=matplotlib.colors.to_rgb(color)
    fc=(c[0],c[1],c[2],alpha)
    if smooth=='kde':
        
        from scipy.stats import kde
        x_smooth=linspace(0,bins.max(),100)
        
        h_smooth=kde.gaussian_kde(data.dropna('ecotype'))(x_smooth)
        
        ax.fill_between(x_smooth,h_smooth,edgecolor=color,facecolor=fc,label=label,**kwargs)
    elif smooth=='spline':
        import scipy.interpolate
        
        x_smooth=arange(0,x.max(),10)
        
        h_smooth=scipy.interpolate.PchipInterpolator(concatenate(([0],x)),concatenate(([0],h)),extrapolate=False)(x_smooth) 
        h_smooth[h_smooth<0]=0
        if len(bottom) == len(x_smooth):
            bottom_smooth=bottom
        else:
            bottom_smooth=scipy.interpolate.PchipInterpolator(concatenate(([0],x)),concatenate(([0],bottom)),extrapolate=False)(x_smooth) 
            bottom_smooth[bottom_smooth<0]=0
        ax.fill_between(x_smooth,bottom_smooth,h_smooth*norm+bottom_smooth,edgecolor=color,facecolor=fc,label=label,where=(h_smooth>0)|(bottom_smooth==0),**kwargs)
        if show_points:
            ax.plot(x,h*norm+bottom_smooth[[abs(x_smooth-xx).argmin() for xx in x]],'o',c=color,ms=2.5)
        return h_smooth*norm+bottom_smooth
    else:
        ax.fill_between(x,bottom,h*norm+bottom,edgecolor=color,facecolor=fc,label=label,**kwargs)
    return h*norm+bottom

overlap=False
f,axs=subplots(nrows=4,ncols=1,num='Histograms',clear=True,figsize=(3.8,9.8))
shrubs=data_global.PFTnames.data.tolist().index('broadleaf_deciduous_boreal_shrub')
gram=data_global.PFTnames.data.tolist().index('c3_arctic_grass')
norm=1/(get_var_PFTs('TOTVEGC',data_global).mean(dim='time').count().item()-2)
ax=axs[1]
bottom=plot_hist(get_var_PFTs('TOTVEGC',data_global).isel(PFT=shrubs).mean(dim='time'),ax,linestyle='-',norm=norm,zorder=1)
if overlap:
    bottom[:]=0
bottom=plot_hist(get_var_PFTs('TOTVEGC',data_global).isel(PFT=gram).mean(dim='time'),ax,linestyle='--',bottom=bottom,norm=norm)
ax.set_title('Simulation 1')
ax.legend(handles=ax.collections[::-1])
ax.set_ylabel('Fraction of patches')

shrubs=data_communities.PFTnames.data.tolist().index('broadleaf_deciduous_boreal_shrub')
gram=data_communities.PFTnames.data.tolist().index('c3_arctic_grass')
norm=1/(get_var_PFTs('TOTVEGC',data_communities).mean(dim='time').count().item()-2)
ax=axs[2]
bottom=plot_hist(get_var_PFTs('TOTVEGC',data_communities).isel(PFT=shrubs).mean(dim='time'),ax,linestyle='-',norm=norm,zorder=1)
if overlap:
    bottom[:]=0
bottom=plot_hist(get_var_PFTs('TOTVEGC',data_communities).isel(PFT=gram).mean(dim='time'),ax,linestyle='--',bottom=bottom,norm=norm)
ax.set_title('Simulation 2')

ax.set_ylabel('Fraction of patches')

pft_order_meas=['bryophyte', 'dwarf shrub deciduous', 'dwarf shrub evergreen',
       'forb', 'graminoid', 'lichen', 'low shrub deciduous', 'mixed',
       'potential tall shrub deciduous non-alder','potential tall shrub deciduous alder']

pft_order=[pft_names.index(obsdata_PFT_mappings[name]) for name in pft_order_meas]  
norm=1/(get_var_PFTs('TOTVEGC',data_Arcticpfts).mean(dim='time').count().item()-1)
bottom[:]=0
ax=axs[3]
for pft in pft_order:
    if 'shrub' in pft_names[pft] and pft_names.index(pft_names[pft]) in pfts_inuse:
        bottom=plot_hist(get_var_PFTs('TOTVEGC',data_Arcticpfts).isel(PFT=pft).mean(dim='time'),ax=ax,bottom=bottom,linestyle='-',norm=norm,zorder=2)

if overlap:
    bottom[:]=0
for pft in pft_order:
    if 'forb' in pft_names[pft] or 'graminoid' in pft_names[pft] and pft_names.index(pft_names[pft]) in pfts_inuse:
        bottom=plot_hist(get_var_PFTs('TOTVEGC',data_Arcticpfts).isel(PFT=pft).mean(dim='time'),ax=ax,bottom=bottom,linestyle='--',norm=norm)
ax.set_title('Arctic PFTs')

if overlap:
    bottom[:]=0
for pft in pft_order:
    if 'bryophyte' in pft_names[pft] or 'lichen' in pft_names[pft] and pft_names.index(pft_names[pft]) in pfts_inuse:
        bottom=plot_hist(get_var_PFTs('TOTVEGC',data_Arcticpfts).isel(PFT=pft).mean(dim='time'),ax=ax,bottom=bottom,linestyle=':',norm=norm,zorder=-1)
ax.set_title('Simulation 3')

ax.set_ylabel('Fraction of patches')

totalC=meas_leaf_C+meas_rhizome_C+meas_root_C+meas_stem_C
norm=1/len(totalC)
bottom[:]=0
ax=axs[0]
for pft in pft_order_meas:
    if 'shrub' in pft:
        bottom=plot_hist(totalC[:,:,pft],ax=ax,bottom=bottom,color=pft_colors[pft_names.index(obsdata_PFT_mappings[pft])],label='',linestyle='-',norm=norm,zorder=2)

if overlap:
    bottom[:]=0
for pft in pft_order_meas:
    if 'graminoid' in pft or 'forb' in pft:
        bottom=plot_hist(totalC[:,:,pft],ax=ax,bottom=bottom,color=pft_colors[pft_names.index(obsdata_PFT_mappings[pft])],label='',linestyle='--',norm=norm)

if overlap:
    bottom[:]=0
for pft in pft_order_meas:
    if 'bryophyte' in pft or 'lichen' in pft:
        bottom=plot_hist(meas_nonvasc_C[:,:,pft],ax=ax,bottom=bottom,color=pft_colors[pft_names.index(obsdata_PFT_mappings[pft])],label='',linestyle=':',norm=norm,zorder=-1)

ax.set_title('Measurements')
ax.set_ylabel('Fraction of plots')
ax.legend(handles=axs[3].collections[::-1])

for num,ax in enumerate(axs.ravel()):
    ax.set_xlabel('Biomass (gC m$^{-2}$)')
    ax.text(0.02,1.05,ascii_lowercase[num],transform=ax.transAxes,fontsize='large')

# Soil depth measurements
soildepth_mean=soildepth.replace({'Rock':0.0}).groupby('ecotype').mean().mean(axis=1)

# Model max rooting depths from a and b parameters
def zfunc(d,a,b):
    return 1-0.5*(exp(-a*d)+exp(-b*d))
    
def max_rooting_depth(a,b,thresh=0.99):
    import scipy.optimize as opt
    def zfunc0(d,a,b,thresh):
        z=zfunc(d,a,b)
        return log(z/thresh)
    return opt.fsolve(zfunc0,x0=1/min(a,b),args=(a,b,thresh))

# Active layer thickness
columndata_arcticpfts=xarray.open_dataset('../../output_data/Kougarok_SNAP_bzdormancy_Arcticpfts_20200715_h0.nc').sel(time=slice(start,end))
columndata_soildepth=xarray.open_dataset('../../output_data/E3SMpfts_soilthickness_20200316_h0.nc')
month=array([t.item().month for t in columndata_arcticpfts.time])
july_ALT=columndata_arcticpfts['ALT'].isel(time=(month==7)).mean(dim='time')

