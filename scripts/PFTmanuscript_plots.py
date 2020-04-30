from plot_kougarok import *

import warnings

# Don't need to print this warning a billion times
warnings.filterwarnings(action='ignore',message='All-NaN slice encountered')  
warnings.filterwarnings(action='ignore',message='Mean of empty slice')  

# data_global=read_pftfile('../../output_data/E3SMpfts_20200316_h2.nc',maxyear=None)
# data_soildepth=read_pftfile('../../output_data/E3SMpfts_soilthickness_20200316_h2.nc',maxyear=None)
# data_oldpfts=read_pftfile('../../output_data/E3SMpfts_communities_20200316_h2.nc',maxyear=None)
# data_newpfts=read_pftfile('../../output_data/Arcticpfts_20200323_h2.nc',maxyear=None)

# data_global=xarray.open_dataset('../../output_data/E3SMpfts_processed_20200316.nc')
# data_soildepth=xarray.open_dataset('../../output_data/E3SMpfts_soilthickness_processed_20200316.nc')
# data_communities=xarray.open_dataset('../../output_data/E3SMpfts_communities_processed_20200316.nc')
# data_Arcticpfts=xarray.open_dataset('../../output_data/Arcticpfts_processed_20200323.nc')

# data_global=xarray.open_dataset('../../output_data/E3SMpfts_processed_20200316.nc')
# data_soildepth=xarray.open_dataset('../../output_data/E3SMpfts_soilthickness_processed_20200316.nc')
# data_communities=xarray.open_dataset('../../output_data/E3SMpfts_communities_processed_20200316.nc')
# data_Arcticpfts=xarray.open_dataset('../../output_data/Arcticpfts_processed_20200323.nc')

e3sm_all=xarray.open_dataset('../../output_data/E3SMpfts_20200415_h2_all_processed.nc')
arcticpfts_all=xarray.open_dataset('../../output_data/Kougarok_SNAP_rhizomefixed_Arcticpfts_20200422_h2_all_processed.nc')
soilthickness_all=xarray.open_dataset('../../output_data/Kougarok_SNAP_rhizomefixed_E3SMpfts_soilthickness_20200422_h2_all_processed.nc')
communities_all=xarray.open_dataset('../../output_data/Kougarok_SNAP_rhizomefixed_E3SMpfts_communities_20200422_h2_all_processed.nc')

start='1990-01-01'
end='2010-01-01'

data_global=e3sm_all.sel(time=slice(start,end))
data_soildepth=soilthickness_all.sel(time=slice(start,end))
data_communities=communities_all.sel(time=slice(start,end))
data_Arcticpfts=arcticpfts_all.sel(time=slice(start,end))


meas_leaf_C.loc[:,:,'bryophyte'] = 0.0
meas_leaf_C.loc[:,:,'lichen'] = 0.0



f=figure(num='Ecosystem biomass',figsize=(12, 8),clear=True)
gs=f.add_gridspec(nrows=7,ncols=6,height_ratios=[1,.25,1,1,1,0.3,1])
axs=zeros((5,6),dtype='object')
for r,row in enumerate([0,2,3,4,6]):
    for col in range(6):
        axs[r,col]=f.add_subplot(gs[row,col])
for simnum,data_to_plot in enumerate([data_global,data_soildepth,data_communities,data_Arcticpfts]):
    for econum in range(len(ecotype_names_list)):
        names=[]
        x=0.0
        ax=axs[simnum,econum]
        sca(ax)
        
        leafc=get_var_PFTs('LEAFC',data_to_plot)
        nonvasc=(leafc.PFTnames=='arctic_lichen')|(leafc.PFTnames=='arctic_bryophyte')

        h=plot_mod_bar_stack(x,leafc[:,:,nonvasc],econum,width=0.8,hatch='//')
        # plot_obs_bar_stack(x,meas_nonvasc_C ,econum,width=0.8)
        names.append('Nonvasc')
        # text(x+0.4,h[-1][0].get_y()+h[-1][0].get_height()+50,'Modeled value',rotation=90,va='bottom',ha='center')
        # text(x,h[-1][0].get_y()+h[-1][0].get_height()+50,'Measured value',rotation=90,va='bottom',ha='center')
        x+=1
        
        # h=plot_mod_bar_stack(x,leafc[:,:,~nonvasc],min(econum,len(data_to_plot.ecotype)-1),width=0.8,hatch='//')
        # # plot_obs_bar_stack(x,meas_leaf_C ,econum,width=0.8)
        # names.append('Leaf')
        # # text(x+0.4,h[-1][0].get_y()+h[-1][0].get_height()+50,'Modeled value',rotation=90,va='bottom',ha='center')
        # # text(x,h[-1][0].get_y()+h[-1][0].get_height()+50,'Measured value',rotation=90,va='bottom',ha='center')
        # x+=1

        plot_mod_bar_stack(x,get_var_PFTs(['LIVESTEMC','DEADSTEMC'],data_to_plot)+leafc[:,:,~nonvasc],min(econum,len(data_to_plot.ecotype)-1),width=0.8,hatch='//')
        # plot_obs_bar_stack(x,meas_stem_C,econum,width=0.8)
        names.append('AG')
        x+=1

        plot_mod_bar_stack(x,get_var_PFTs(['LIVECROOTC','DEADCROOTC','FROOTC'],data_to_plot),min(econum,len(data_to_plot.ecotype)-1),width=0.8,hatch='//')
        # bar(x,meas_root_C[landscape_ecotypes[econum]].mean(),yerr=meas_root_C[landscape_ecotypes[econum]].std(),width=0.8,facecolor=[0.9,0.9,0.9],edgecolor='k')
        names.append('BG')
        x+=1

        # plot_mod_bar_stack(x,get_var_PFTs(['LIVECROOTC','DEADCROOTC'],data_to_plot),min(econum,len(data_to_plot.ecotype)-1),width=0.8,hatch='//')
        # # text(x+0.4,100,'? ? ? ? ? ?',fontsize='large',rotation=90,va='bottom',ha='center')
        # # plot_obs_bar_stack(x,meas_rhizome_C,econum,width=0.8)
        # names.append('Croot')
        # x+=1
        
        # plot_mod_bar_stack(x+0.4,get_var_PFTs(['STORVEGC'],data_to_plot),econum,width=0.4,hatch='//')
        # names.append('Storage')

        ylabel('Biomass (gC m$^{-2}$)')
        xticks(arange(len(names)),names)
        ylim(-50,1050)
        ax.set_xlim(left=-0.64)
        # title(ecotype_names_list[econum])
        # ylim(0,3000)
        # ax.set_xlim(right=x+0.6)

for econum in range(len(ecotype_names_list)):
    names=[]
    x=0.0
    ax=axs[4,econum]
    sca(ax)

    # h=plot_mod_bar_stack(x+0.4,get_var_PFTs('LEAFC',data_to_plot),econum,width=0.4,hatch='//')
    plot_obs_bar_stack(x,meas_nonvasc_C ,econum,width=0.8)
    names.append('Nonvasc')
    # text(x+0.4,h[-1][0].get_y()+h[-1][0].get_height()+50,'Modeled value',rotation=90,va='bottom',ha='center')
    # text(x,h[-1][0].get_y()+h[-1][0].get_height()+50,'Measured value',rotation=90,va='bottom',ha='center')
    x+=1
    
    plot_obs_bar_stack(x,meas_leaf_C+meas_stem_C ,econum,width=0.8)
    names.append('AG')
    # text(x+0.4,h[-1][0].get_y()+h[-1][0].get_height()+50,'Modeled value',rotation=90,va='bottom',ha='center')
    # text(x,h[-1][0].get_y()+h[-1][0].get_height()+50,'Measured value',rotation=90,va='bottom',ha='center')
    x+=1

    # # plot_mod_bar_stack(x+0.4,get_var_PFTs(['LIVESTEMC','DEADSTEMC'],data_to_plot),econum,width=0.4,hatch='//')
    # plot_obs_bar_stack(x,meas_stem_C,econum,width=0.8)
    # names.append('Stem')
    # x+=1

    # plot_mod_bar_stack(x+0.4,get_var_PFTs(['LIVECROOTC','DEADCROOTC'],data_to_plot),econum,width=0.4,hatch='//')
    # text(x+0.4,100,'? ? ? ? ? ?',fontsize='large',rotation=90,va='bottom',ha='center')
    plot_obs_bar_stack(x,meas_rhizome_C,econum,width=0.8)
    # names.append('Rhizome')
    # x+=1
    
    # plot_mod_bar_stack(x+0.4,get_var_PFTs('FROOTC',data_to_plot),econum,width=0.4,hatch='//')
    bar(x,meas_root_C[landscape_ecotypes[econum]].mean(),yerr=meas_root_C[landscape_ecotypes[econum]].std(),width=0.8,facecolor=[0.9,0.9,0.9],edgecolor='k',linewidth=0.5,linestyle='--',bottom=meas_rhizome_C[landscape_ecotypes[econum]].sum(level='PlotID').mean() )
    names.append('BG')
    x+=1
    
    # plot_mod_bar_stack(x+0.4,get_var_PFTs(['STORVEGC'],data_to_plot),econum,width=0.4,hatch='//')
    # names.append('Storage')

    ylabel('Biomass (gC m$^{-2}$)')
    xticks(arange(len(names)),names)
    ylim(-50,1250)
    # title(ecotype_names_list[econum])
    # ylim(0,3000)
    # ax.set_xlim(right=x+0.6)

axs[4,2].set_ylim(top=3050)

pfts_global=data_global['PFT'][(data_global.weights>0).any(dim='ecotype')].data.tolist()
handles=[]
for pftnum in pfts_global[1:]:
    name = prettify_pft_name(data_global.PFTnames.values[pftnum])
    if name == 'Dry Graminoid':
        name = 'Graminoid'
    handles.append(Rectangle([0,0],0,0,facecolor=data_global.PFTcolors.values[pftnum],label=name ))
l=axs[1,0].legend(handles=handles,fontsize='small',ncol=2,title='E3SM PFTs',loc=(0.0,1.1))
l.set_in_layout(False)

pfts_Arctic=data_Arcticpfts['PFT'][(data_Arcticpfts.weights>0).any(dim='ecotype')].data.tolist()
handles=[]
for pftnum in pfts_Arctic[1:]:
    name = prettify_pft_name(data_Arcticpfts.PFTnames.values[pftnum])
    if name == 'Dry Graminoid':
        name = 'Graminoid'
    handles.append(Rectangle([0,0],0,0,facecolor=data_Arcticpfts.PFTcolors.values[pftnum],label=name ))
l=axs[4,0].legend(handles=handles,fontsize='small',ncol=5,title='Arctic PFTs',loc=(0.0,1.1))
l.set_in_layout(False)

for econum in range(len(ecotype_names_list)):
    axs[0,econum].set_title(ecotype_names_list[econum])
    
axs[0,0].text(-0.45,0.5,'Level 1:\nE3SM grid cell',rotation=90,ha='center',va='center',fontsize='large',transform=axs[0,0].transAxes)
axs[1,0].text(-0.45,0.5,'Level 2:\nE3SM PFTs, site soil depths',rotation=90,ha='center',va='center',fontsize='large',transform=axs[1,0].transAxes)
axs[2,0].text(-0.45,0.5,'Level 3:\nE3SM PFTs, site areas',rotation=90,ha='center',va='center',fontsize='large',transform=axs[2,0].transAxes)
axs[3,0].text(-0.45,0.5,'Level 4:\nArctic PFTs',rotation=90,ha='center',va='center',fontsize='large',transform=axs[3,0].transAxes)
axs[4,0].text(-0.45,0.5,'Site measurements',rotation=90,ha='center',va='center',fontsize='large',transform=axs[4,0].transAxes)
    
f.suptitle('Biomass by plant community for different site and PFT definitions')




f,axs=subplots(nrows=5,ncols=6,num='Ecosystem NPP',figsize=(12, 8),clear=True)
for simnum,data_to_plot in enumerate([data_global,data_soildepth,data_communities,data_Arcticpfts]):
    for econum in range(len(ecotype_names_list)):
        names=[]
        x=0.0
        ax=axs[simnum,econum]
        sca(ax)
        

        plot_mod_bar_stack(x,get_var_PFTs(['AGNPP'],data_to_plot)*3600*24*365,min(econum,len(data_to_plot.ecotype)-1),width=0.8,hatch='//',op='mean')
        # plot_obs_bar_stack(x,meas_stem_C,econum,width=0.8)
        names.append('AGNPP')
        x+=1

        plot_mod_bar_stack(x,get_var_PFTs(['BGNPP'],data_to_plot)*3600*24*365,min(econum,len(data_to_plot.ecotype)-1),width=0.8,hatch='//',op='mean')
        # plot_obs_bar_stack(x,meas_stem_C,econum,width=0.8)
        names.append('BGNPP')
        x+=1


        ylabel('NPP (gC m$^{-2}$ year$^{-1}$)')
        xticks(arange(len(names)),names)
        # ylim(-50,1050)
        ax.set_xlim(left=-0.64)
        # title(ecotype_names_list[econum])
        # ylim(0,3000)
        # ax.set_xlim(right=x+0.6)

aboveground_NPP = (meas_leaf_NPP+meas_stem_NPP).groupby(level=('Ecotype','PlotID')).sum()
belowground_NPP = meas_rhizome_NPP.groupby(level=('Ecotype','PlotID')).sum() + meas_froot_NPP
for econum in range(len(ecotype_names_list)):
    names=[]
    x=0.0
    ax=axs[4,econum]
    sca(ax)

    
    bar(x,aboveground_NPP[landscape_ecotypes[econum]].mean(),yerr=meas_root_C[landscape_ecotypes[econum]].std(),width=0.8,facecolor=[0.9,0.9,0.9],edgecolor='k')
    names.append('AGNPP')
    x+=1


    # plot_mod_bar_stack(x+0.4,get_var_PFTs(['LIVESTEMC','DEADSTEMC'],data_to_plot),econum,width=0.4,hatch='//')
    bar(x,belowground_NPP[landscape_ecotypes[econum]].mean(),yerr=meas_root_C[landscape_ecotypes[econum]].std(),width=0.8,facecolor=[0.9,0.9,0.9],edgecolor='k')
    names.append('BGNPP')
    x+=1


    # plot_mod_bar_stack(x+0.4,get_var_PFTs(['STORVEGC'],data_to_plot),econum,width=0.4,hatch='//')
    # names.append('Storage')

    ylabel('NPP (gC m$^{-2}$ year$^{-1}$)')
    xticks(arange(len(names)),names)
    # ylim(-50,1050)
    # title(ecotype_names_list[econum])
    # ylim(0,3000)
    # ax.set_xlim(right=x+0.6)

# axs[4,2].set_ylim(top=3050)

pfts_global=data_global['PFT'][(data_global.weights>0).any(dim='ecotype')].data.tolist()
handles=[]
for pftnum in pfts_global[1:]:
    name = prettify_pft_name(data_global.PFTnames.values[pftnum])
    if name == 'Dry Graminoid':
        name = 'Graminoid'
    handles.append(Rectangle([0,0],0,0,facecolor=data_global.PFTcolors.values[pftnum],label=name ))
l=axs[0,0].legend(handles=handles,fontsize='small',ncol=1,title='E3SM PFTs')

pfts_Arctic=data_Arcticpfts['PFT'][(data_Arcticpfts.weights>0).any(dim='ecotype')].data.tolist()
handles=[]
for pftnum in pfts_Arctic[1:]:
    name = prettify_pft_name(data_Arcticpfts.PFTnames.values[pftnum])
    if name == 'Dry Graminoid':
        name = 'Graminoid'
    handles.append(Rectangle([0,0],0,0,facecolor=data_Arcticpfts.PFTcolors.values[pftnum],label=name ))
l=axs[3,0].legend(handles=handles,fontsize='small',ncol=1,title='Arctic PFTs',loc='upper right')

for econum in range(len(ecotype_names_list)):
    axs[0,econum].set_title(ecotype_names_list[econum])
    
axs[0,0].text(-0.38,0.5,'Level 1:\nE3SM grid cell',rotation=90,ha='center',va='center',fontsize='large',transform=axs[0,0].transAxes)
axs[1,0].text(-0.38,0.5,'Level 2:\nE3SM PFTs, site soil depths',rotation=90,ha='center',va='center',fontsize='large',transform=axs[1,0].transAxes)
axs[2,0].text(-0.38,0.5,'Level 3:\nE3SM PFTs, site areas',rotation=90,ha='center',va='center',fontsize='large',transform=axs[2,0].transAxes)
axs[3,0].text(-0.38,0.5,'Level 4:\nArctic PFTs',rotation=90,ha='center',va='center',fontsize='large',transform=axs[3,0].transAxes)
axs[4,0].text(-0.38,0.5,'Site measurements',rotation=90,ha='center',va='center',fontsize='large',transform=axs[4,0].transAxes)
    
f.suptitle('NPP by plant community for different site and PFT definitions')



f,axs=subplots(ncols=3,nrows=1,num='PFT areas',clear=True)
plot_PFT_distributions(axs)


barfig_PFT,axs=subplots(nrows=2,ncols=5,num='Biomass comparison by PFT',figsize=(12,6),clear=True,squeeze=False);
plot_pft_biomass_bars(data_to_plot=data_Arcticpfts,axs=axs.ravel(),per_area=False,use_pooled_default_PFTs=True,include_froot=False,include_storage=False)

# Remove some problematic PFTs. There isn't any tall evergreen shrub, and deciduous dwarf shrub area is too low
pfts_inuse=data_Arcticpfts['PFT'][(data_Arcticpfts.weights>0).any(dim='ecotype')].data.tolist()
# pfts_inuse.remove(data_Arcticpfts.PFTnames.data.tolist().index('arctic_evergreen_shrub_tall'))
pfts_inuse.remove(data_Arcticpfts.PFTnames.data.tolist().index('arctic_deciduous_shrub_dwarf'))
barfig_PFT_perarea,axs=subplots(nrows=2,ncols=4,num='Biomass comparison per area of each PFT',figsize=(12,6),clear=True,squeeze=False);
plot_pft_biomass_bars(data_to_plot=data_Arcticpfts,axs=axs.ravel(),per_area=True,use_pooled_default_PFTs=True,pfts_inuse=pfts_inuse)

f,axs=subplots(ncols=3,nrows=4,num='Default PFTs biomass per area',clear=True,figsize=(11.3,8.8))
pfts_inuse=data_communities['PFT'][(data_communities.weights>0).any(dim='ecotype')].data.tolist()
plot_pft_biomass_bars(axs=axs[0,:],data_to_plot=data_global,per_area=True,use_pooled_default_PFTs=True,leg_axis=1,pfts_inuse=pfts_inuse)
plot_pft_biomass_bars(axs=axs[1,:],data_to_plot=data_soildepth,per_area=True,use_pooled_default_PFTs=True,leg_axis=1,pfts_inuse=pfts_inuse)
plot_pft_biomass_bars(axs=axs[2,:],data_to_plot=data_communities,per_area=True,use_pooled_default_PFTs=True,leg_axis=1,pfts_inuse=pfts_inuse)
plot_pft_biomass_bars(axs=axs[3,:],data_to_plot=new_to_defaultPFTs(data_Arcticpfts,data_soildepth),per_area=True,use_pooled_default_PFTs=True,leg_axis=1,pfts_inuse=pfts_inuse)
# x=axs[0,1].get_xlim()
# y=axs[0,1].get_ylim()
# axs[1,1].set_xlim(*x)
# axs[1,1].set_ylim(*y)
# axs[2,1].set_xlim(*x)
# axs[2,1].set_ylim(*y)
axs[1,1].get_legend().set_visible(False)
axs[2,1].get_legend().set_visible(False)
axs[3,1].get_legend().set_visible(False)
axs[0,0].text(-0.3,0.5,'Level 1:\nE3SM grid cell',rotation=90,ha='center',va='center',fontsize='large',transform=axs[0,0].transAxes)
axs[1,0].text(-0.3,0.5,'Level 2:\nE3SM PFTs, site soil depths',rotation=90,ha='center',va='center',fontsize='large',transform=axs[1,0].transAxes)
axs[2,0].text(-0.3,0.5,'Level 3:\nE3SM PFTs, site areas',rotation=90,ha='center',va='center',fontsize='large',transform=axs[2,0].transAxes)
axs[3,0].text(-0.3,0.5,'Level 4:\nArctic PFTs, site areas',rotation=90,ha='center',va='center',fontsize='large',transform=axs[3,0].transAxes)
for num in range(4):
    axs[num,1].set_title('Shrubs')
    axs[num,2].set_title('Grasses')



import cartopy.crs as ccrs
import cartopy
figure('Site location',clear=True)
ax=subplot(111,projection=ccrs.Miller(central_longitude=-155))
ax.set_extent((-175,-138,50,75)) 
ax.coastlines(resolution='50m')
ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'))
ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'))
ax.scatter(-164.82,65.16,marker='*',s=50,c='g',transform=ccrs.Geodetic())



def plot_timeseries(d,axs,plot_sum=False,resample='1Y',resample_op='mean',ymin=-150,ymax=4600,**kwargs):
    for econum in range(0,len(d.ecotype)):
        ax=axs.ravel()[econum]
        if resample is not None:
            data=getattr(d.resample(time=resample),resample_op)()
        else:
            data=d
        handles=[]
        if len(data.PFT)==12:
            pfts=range(1,11)
        else:
            pfts=data.PFT[data.max(dim='time').isel(ecotype=econum)>0]
        for pft in pfts:
            handles.append(ax.plot(data['time'],data.isel(PFT=pft,ecotype=econum),c=data['PFTcolors'][pft].item(),label=prettify_pft_name(data['PFTnames'][pft].item()),**kwargs)[0])
        if plot_sum:
            handles.append(ax.plot(data['time'],data.isel(ecotype=econum).sum(dim='PFT'),c='k',ls='--',label='Sum over PFTs')[0])
        # plot_var_PFTs(var,hist_pfts,econum,ls='-',minyear=1850,maxyear=2100,ax=ax)
        ax.set_title(ecotype_names[landscape_ecotypes[econum]])
        ax.set_xlabel('Time (years)')
        ax.set_ylabel('Biomass (g m$^{-2}$)')
        ax.axvline(2016,ls=':',lw=0.5,c='k')
        ax.set_ylim(ymin,ymax)
    ax.legend(handles=handles,ncol=2)


f,axs=subplots(nrows=2,ncols=3,num='Stored fraction',clear=True,figsize=(12,6))
plot_timeseries(get_var_PFTs('STORVEGC', arcticpfts_all)/get_var_PFTs('TOTVEGC',arcticpfts_all),axs,resample_op='min',ymin=-.1,ymax=1.1)

f,axs=subplots(nrows=2,ncols=3,num='Total biomass',clear=True,figsize=(12,6))
plot_timeseries(get_var_PFTs('TOTVEGC', arcticpfts_all),axs,True)

f,axs=subplots(nrows=2,ncols=3,num='Total biomass (soilthickness)',clear=True,figsize=(12,6))
plot_timeseries(get_var_PFTs('TOTVEGC', soilthickness_all),axs,True)
    
f,axs=subplots(nrows=2,ncols=3,num='Total biomass (communities)',clear=True,figsize=(12,6))
plot_timeseries(get_var_PFTs('TOTVEGC', communities_all),axs,True)
    
    
f,axs=subplots(nrows=2,ncols=3,num='Total biomass per unit PFT area',clear=True,figsize=(12,6))
plot_timeseries(get_var_PFTs('TOTVEGC', arcticpfts_all,weight_area=False),axs,ymax=6000)
    
f,axs=subplots(nrows=2,ncols=3,num='Total Leaf C',clear=True,figsize=(12,6))
plot_timeseries(get_var_PFTs('LEAFC', arcticpfts_all,weight_area=True),axs,True,resample_op='max',ymax=400)

f,axs=subplots(1,1,num='E3SM PFTs',clear=True,squeeze=False)
plot_timeseries(get_var_PFTs('TOTVEGC', e3sm_all),axs,True)
axs[0,0].set_title('E3SM original PFTs')





# Aboveground biomass
obs_ag_C = (meas_leaf_C + meas_stem_C).add(meas_nonvasc_C,fill_value=0.0)
obs_bg_C = meas_rhizome_C.groupby(level=('Ecotype','PlotID')).sum() + meas_root_C
ag_mod_C_arctic = get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC'],data_Arcticpfts)
bg_mod_C_arctic = get_var_PFTs(['LIVECROOTC','DEADCROOTC','FROOTC'],data_Arcticpfts)
ag_mod_C_global = get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC'],data_soildepth)
bg_mod_C_global = get_var_PFTs(['LIVECROOTC','DEADCROOTC','FROOTC'],data_soildepth)

obs=(meas_leaf_C + meas_stem_C + meas_rhizome_C).add(meas_nonvasc_C,fill_value=0.0)
mod_new=get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC','LIVECROOTC','DEADCROOTC'],data_Arcticpfts)
mod_soildepth=get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC','LIVECROOTC','DEADCROOTC'],data_soildepth)

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
                
            
    # x=obs_bg_C[landscape_ecotypes[econum]].mean()
    # xerr=obs_bg_C[econum].std()
    # y=bg_mod_C_arctic.max(dim='time').sum(dim='PFT').isel(ecotype=econum)
    # axs[1,1].errorbar(x,y,c='k',ls='None',marker=markers[econum],xerr=xerr)
    # 
    # y=bg_mod_C_global.max(dim='time').sum(dim='PFT').isel(ecotype=0)
    # axs[0,1].errorbar(x,y,c='k',ls='None',marker=markers[econum],xerr=xerr)
    # 
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
# handles=handles+[Line2D([0,0],[0,0],ls='None',marker=markers[m],c='k',ms=8.0) for m in ecotypes_included]
# labels.extend(ecotype_names_list[e] for e in ecotypes_included)
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
                
            
    # x=obs_bg_C[landscape_ecotypes[econum]].mean()
    # xerr=obs_bg_C[econum].std()
    # y=bg_mod_C_arctic.max(dim='time').sum(dim='PFT').isel(ecotype=econum)
    # axs[1,1].errorbar(x,y,c='k',ls='None',marker=markers[econum],xerr=xerr)
    # 
    # y=bg_mod_C_global.max(dim='time').sum(dim='PFT').isel(ecotype=0)
    # axs[0,1].errorbar(x,y,c='k',ls='None',marker=markers[econum],xerr=xerr)
    # 
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

axs[1].plot(linspace(0,1800,10),linspace(0,1800,10),'k:')
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


f=figure('PFT comp',clear=True)
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

plot_pft_biomass_bars(axs=[axs['nonvascular'],axs[data_soildepth.PFTnames[11].item()],axs[data_soildepth.PFTnames[12].item()]],
                    data_to_plot=data_soildepth,per_area=False,include_storage=False,include_froot=False,use_pooled_default_PFTs=True,leg_axis=1,pfts_inuse=[0,11,12])
axs['broadleaf_deciduous_boreal_shrub'].set_title('Arctic Deciduous Shrubs')
axs['c3_arctic_grass'].set_title('Arctic grasses')
axs['nonvascular'].set_title('Nonvascular')

pfts_inuse=data_Arcticpfts['PFT'][(data_Arcticpfts.weights>0).any(dim='ecotype')].data.tolist()
plot_pft_biomass_bars(axs=[axs[data_Arcticpfts.PFTnames[pft].item()] for pft in pfts_inuse[1:]],data_to_plot=data_Arcticpfts,per_area=False,include_storage=False,include_froot=False,use_pooled_default_PFTs=True,leg_axis=None,pfts_inuse=pfts_inuse)
ev_dwarf_ax.set_title('Evergreen dwarf shrubs')
dec_dwarf_ax.set_title('Deciduous dwarf shrubs')
dec_low_ax.set_title('Deciduous low shrubs')
dec_tall_ax.set_title('Deciduous tall shrubs')
alder_ax.set_title('Alder shrubs')
gram_new_ax.set_title('Arctic graminoids')
forb_ax.set_title('Forbs')
bryophyte_ax.set_title('Bryophytes')
lichen_ax.set_title('Lichens')

# f.text(0.5,0.52,'New Arctic PFTs',fontsize='x-large',ha='center') 
# f.text(0.5,0.95,'Original PFTs',fontsize='x-large',ha='center') 
f.add_artist(Line2D([0.42,0.42],[0.0,0.95],ls='--',c='k')) 
f.add_artist(Line2D([0.71,0.71],[0.0,0.95],ls='--',c='k')) 
f.add_artist(matplotlib.patches.FancyArrow(0.23,0.56,0.0,-0.06,width=0.01,length_includes_head=True,head_width=0.025))
f.add_artist(matplotlib.patches.FancyArrow(0.56,0.56,0.0,-0.06,width=0.01,length_includes_head=True,head_width=0.025))
f.add_artist(matplotlib.patches.FancyArrow(0.86,0.56,0.0,-0.06,width=0.01,length_includes_head=True,head_width=0.025))




bins=arange(0,6400,300)
# bins=concatenate((arange(0,2500,300),[4000,5000,6000]))
def plot_hist(data,ax,color=None,label=None,bottom=None,bins=bins,alpha=0.5,norm=1.0,smooth='spline',show_points=False,**kwargs):
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
        # x_smooth=linspace(x.min(),x.max(),len(x)*10+1)
        x_smooth=arange(0,x.max(),10)
        # x_smooth=concatenate((arange(x.min(),x.max(),(x[1]-x[0])/10),[x.max()]))
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
f,axs=subplots(nrows=2,ncols=2,num='Histograms',clear=True)
shrubs=data_global.PFTnames.data.tolist().index('broadleaf_deciduous_boreal_shrub')
gram=data_global.PFTnames.data.tolist().index('c3_arctic_grass')
norm=1/(get_var_PFTs('TOTVEGC',data_global).mean(dim='time').count().item()-2)
ax=axs[0,1]
bottom=plot_hist(get_var_PFTs('TOTVEGC',data_global).isel(PFT=shrubs).mean(dim='time'),ax,linestyle='-',norm=norm,zorder=1)
if overlap:
    bottom[:]=0
bottom=plot_hist(get_var_PFTs('TOTVEGC',data_global).isel(PFT=gram).mean(dim='time'),ax,linestyle='--',bottom=bottom,norm=norm)
ax.set_title('E3SM grid cell')
ax.legend(handles=ax.collections[::-1])
ax.set_ylabel('Fraction of patches')

shrubs=data_soildepth.PFTnames.data.tolist().index('broadleaf_deciduous_boreal_shrub')
gram=data_soildepth.PFTnames.data.tolist().index('c3_arctic_grass')
norm=1/(get_var_PFTs('TOTVEGC',data_soildepth).mean(dim='time').count().item()-2)
ax=axs[1,1]
bottom=plot_hist(get_var_PFTs('TOTVEGC',data_soildepth).isel(PFT=shrubs).mean(dim='time'),ax,linestyle='-',norm=norm,zorder=1)
if overlap:
    bottom[:]=0
bottom=plot_hist(get_var_PFTs('TOTVEGC',data_soildepth).isel(PFT=gram).mean(dim='time'),ax,linestyle='--',bottom=bottom,norm=norm)
ax.set_title('E3SM PFTs, site soil depths')
# ax.legend(handles=ax.collections[::-1])
ax.set_ylabel('Fraction of patches')

pft_order_meas=['bryophyte', 'dwarf shrub deciduous', 'dwarf shrub evergreen',
       'forb', 'graminoid', 'lichen', 'low shrub deciduous', 'mixed',
       'potential tall shrub deciduous non-alder','potential tall shrub deciduous alder']

pft_order=[pft_names.index(obsdata_PFT_mappings[name]) for name in pft_order_meas]  
norm=1/(get_var_PFTs('TOTVEGC',data_Arcticpfts).mean(dim='time').count().item()-1)
bottom[:]=0
ax=axs[1,0]
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
ax.set_title('Arctic PFTs')

ax.set_ylabel('Fraction of patches')

totalC=meas_leaf_C+meas_rhizome_C+meas_root_C+meas_stem_C
norm=1/len(totalC)
bottom[:]=0
ax=axs[0,0]
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
ax.legend(handles=axs[1,0].collections[::-1])

axs.ravel()[0].set_xlabel('Biomass (gC m$^{-2}$)')
axs.ravel()[1].set_xlabel('Biomass (gC m$^{-2}$)')
axs.ravel()[2].set_xlabel('Biomass (gC m$^{-2}$)')
axs.ravel()[3].set_xlabel('Biomass (gC m$^{-2}$)')


# axs[0].set_xlim(0,2050)
# axs[1].set_xlim(0,2050)
# axs[2].set_xlim(0,2050)