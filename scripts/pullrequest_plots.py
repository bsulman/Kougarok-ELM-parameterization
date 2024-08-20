from plot_kougarok import *

import warnings

# Don't need to print this warning a billion times
warnings.filterwarnings(action='ignore',message='All-NaN slice encountered')  
warnings.filterwarnings(action='ignore',message='Mean of empty slice')  

matplotlib.rcParams['figure.constrained_layout.use']=True


e3sm_all=read_pftfile('/gpfs/wolf2/cades/cli185/proj-shared/f9y/archives/AK-K64_usrpft_arctic_ICB20TRCNPRDCTCBC/run/run_default/AK-K64_usrpft_arctic_ICB20TRCNPRDCTCBC.elm.h0.2008-01-01-00000.nc',decode_times=True)
# arcticpfts_all=xarray.open_dataset('../../output_data/Kougarok_SNAP_bzdormancy_Arcticpfts_20200716_h2_processed.nc',decode_times=True) # This is with dormancy set to -1 C and root_a params fixed
arcticpfts_all=read_pftfile('/gpfs/wolf2/cades/cli185/proj-shared/f9y/archives/AK-K64_usrpft_arctic_ICB20TRCNPRDCTCBC/run/run_arcticpft_vsoilthick/AK-K64_usrpft_arctic_ICB20TRCNPRDCTCBC.elm.h0.2008-01-01-00000.nc',decode_times=True)
# soilthickness_all=xarray.open_dataset('../../output_data/Kougarok_SNAP_rhizomefixed_E3SMpfts_soilthickness_20200422_h2_all_processed.nc',decode_times=True)
communities_all=read_pftfile('/gpfs/wolf2/cades/cli185/proj-shared/f9y/archives/AK-K64_usrpft_arctic_ICB20TRCNPRDCTCBC/run/run_default_vsoilthick/AK-K64_usrpft_arctic_ICB20TRCNPRDCTCBC.elm.h0.2008-01-01-00000.nc',decode_times=True) # 515 changed to RCP8.5 CO2 etc

newparams=xarray.open_dataset('/gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata/lnd/clm2/paramdata/clm_params_c180524-sub12_updated20240201.nc')

meas_leaf_C.loc[:,:,'bryophyte'] = 0.0
meas_leaf_C.loc[:,:,'lichen'] = 0.0


f=figure(num='Ecosystem biomass',figsize=(14.76, 5.26),clear=True,constrained_layout=True)
f.set_constrained_layout_pads(hspace=0,h_pad=0.02)
gs=f.add_gridspec(nrows=3,ncols=6,height_ratios=[1,0.5,1])
axs=zeros((2,6),dtype='object')
offsetw=2.0
for row in range(6):
    axs[0,row]=f.add_subplot(gs[0,row])
    axs[1,row]=f.add_subplot(gs[2,row])

for simnum,data_to_plot in enumerate([e3sm_all,arcticpfts_all]):
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
    plot_obs_bar_stack(x+offset,meas_nonvasc_C ,econum,width=0.2,pfts=['lichen','bryophyte'],ebar_args={'elinewidth':1.1,'capsize':3.0})
    names.append('Nonvascular')
    ticlocs.append(x+offset)
    
    sca(ax)
    # x+=1.3
    
    plot_obs_bar_stack(x+offset,meas_leaf_C+meas_stem_C ,econum,width=0.2,pfts=pft_order,ebar_args={'elinewidth':1.1,'capsize':3.0})
    names.append('Vascular')

    plot_obs_bar_stack(x+offset,-meas_rhizome_C,econum,width=0.2,ebar_args={'elinewidth':1.1,'capsize':3.0})

    
    bar(x+offset,-meas_root_C[landscape_ecotypes[econum]].mean(),yerr=meas_root_C[landscape_ecotypes[econum]].std(),width=0.2,facecolor=[0.9,0.9,0.9],edgecolor='k',linewidth=0.5,linestyle='--',
            bottom=-meas_rhizome_C[landscape_ecotypes[econum]].groupby('PlotID').sum().mean() ,capsize=3.0)
    
    ticlocs.append(x+offset)
    
    x+=1
    
    axhline(0.0,c='k',ls=':',lw=1.0)


pfts_global=e3sm_all['PFT'][(e3sm_all.weights>0).any(dim='ecotype')].data.tolist()
handles=[]
for pftnum in pfts_global[1:]:
    name = prettify_pft_name(e3sm_all.PFTnames.values[pftnum])
    if name == 'Dry Graminoid':
        name = 'Graminoid'
    handles.append(Rectangle([0,0],0,0,facecolor=e3sm_all.PFTcolors.values[pftnum],label=name ))
l=axs[0,0].legend(handles=handles,fontsize='medium',ncol=2,title='E3SM PFTs (Sim 1)',loc=(0,-0.6))
l.set_in_layout(False)

pfts_Arctic=arcticpfts_all['PFT'][(arcticpfts_all.weights>0).any(dim='ecotype')].data.tolist()
handles=[]
for pftnum in pfts_Arctic[1:]:
    name = prettify_pft_name(arcticpfts_all.PFTnames.values[pftnum])
    if name == 'Dry Graminoid':
        name = 'Graminoid'
    handles.append(Rectangle([0,0],0,0,facecolor=arcticpfts_all.PFTcolors.values[pftnum],label=name ))
handles.append(Rectangle([0,0],0,0,facecolor=[0.9,0.9,0.9],edgecolor='k',linewidth=0.5,linestyle='--',label='Pooled fine roots' ))
handles[0].set_hatch('..')
l=axs[0,2].legend(handles=handles,fontsize='medium',ncol=5,title='Arctic PFTs',loc=(0,-0.6))
l.set_in_layout(False)

# One reviewer had trouble distinguishing the lichen and bryophyte colors, so adding hatching to lichen
for ax in axs[1,:]:
    ax.patches[0].set_hatch('..')
    ax.patches[2].set_hatch('..')

from string import ascii_lowercase
for econum in range(len(ecotype_names_list)):
    for row in [0,1]:
        
        axs[row,econum].set_title(landscape_ecotypes[econum])
        axs[row,econum].set_xticks(arange(3)*0.22)
        
        axs[row,econum].set_xticklabels(['Meas.','E3SM PFTs','Arctic PFTs'])

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
axs[0].set_ylim(top=100)
axs[1].set_ylim(top=100)
axs[2].set_ylim(top=100)


barfig_PFT,axs=subplots(nrows=2,ncols=5,num='Biomass comparison by PFT',figsize=(12,6),clear=True,squeeze=False);
plot_pft_biomass_bars(data_to_plot=arcticpfts_all,axs=axs.ravel(),per_area=False,use_pooled_default_PFTs=True,include_froot=False,include_storage=False)

# Remove PFTs with zero area
pfts_inuse=arcticpfts_all['PFT'][(arcticpfts_all.weights>0).any(dim='ecotype')].data.tolist()

barfig_PFT_perarea,axs=subplots(nrows=2,ncols=5,num='Biomass comparison per area of each PFT',figsize=(12,6),clear=True,squeeze=False);
plot_pft_biomass_bars(data_to_plot=arcticpfts_all,axs=axs.ravel(),per_area=True,use_pooled_default_PFTs=False,pfts_inuse=pfts_inuse,include_storage=False)
axs[1,4].set_visible(False) 


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

import nc_time_axis

# Aboveground biomass
obs_ag_C = (meas_leaf_C + meas_stem_C).add(meas_nonvasc_C,fill_value=0.0)
obs_bg_C = meas_rhizome_C.groupby(level=('Ecotype','PlotID')).sum() + meas_root_C
ag_mod_C_arctic = get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC'],arcticpfts_all)
bg_mod_C_arctic = get_var_PFTs(['LIVECROOTC','DEADCROOTC','FROOTC'],arcticpfts_all)
ag_mod_C_global = get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC'],communities_all)
bg_mod_C_global = get_var_PFTs(['LIVECROOTC','DEADCROOTC','FROOTC'],communities_all)

obs=(meas_leaf_C + meas_stem_C + meas_rhizome_C).add(meas_nonvasc_C,fill_value=0.0)
mod_new=get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC','LIVECROOTC','DEADCROOTC'],arcticpfts_all)
mod_soildepth=get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC','LIVECROOTC','DEADCROOTC'],communities_all)

markers=['o','s','h','^','>','*']



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

pfts_inuse=arcticpfts_all['PFT'][(arcticpfts_all.weights>0).any(dim='ecotype')].data.tolist()

a1=f.add_artist(matplotlib.patches.FancyArrow(0.25,0.595,0.0,-0.06,width=0.02,length_includes_head=True,head_width=0.05,head_length=0.03,facecolor='k'))
a2=f.add_artist(matplotlib.patches.FancyArrow(0.75,0.595,0.0,-0.06,width=0.02,length_includes_head=True,head_width=0.05,head_length=0.03,facecolor='k'))

r1=f.add_artist(Rectangle(xy=(0.045,0.59),width=0.4,height=0.38,zorder=-10,facecolor='0.9'))
r2=f.add_artist(Rectangle(xy=(0.55,0.59),width=0.4,height=0.38,zorder=-10,facecolor='0.9'))
r3=f.add_artist(Rectangle(xy=(0.02,0.01),width=0.50,height=0.53,zorder=-10,facecolor='0.9'))
r4=f.add_artist(Rectangle(xy=(0.53,0.01),width=.98-.53,height=0.53,zorder=-10,facecolor='0.9'))
r5=f.add_artist(Rectangle(xy=(0.03,0.58),width=0.93,height=0.999-0.58,facecolor='0.95',zorder=-11))
r6=f.add_artist(Rectangle(xy=(0.00,0.0),width=1,height=0.57,facecolor='0.95',zorder=-11))

t=f.text(0.5,0.55,'Arctic PFTs',fontsize='x-large',ha='center')    
t2=f.suptitle('Default PFTs',fontsize='x-large')

ecotypes_included=range(6)#[4,5]
for econum in ecotypes_included:#range(len(ecotype_names_list)):
    for pft in unique(obs.index.get_level_values('ELM_PFT')):
        if obsdata_PFT_mappings[pft] not in axs or 'birch' in pft or 'willow' in pft:
            continue
        if pft in obs[landscape_ecotypes[econum]].index.get_level_values('ELM_PFT'):
            ax=axs[obsdata_PFT_mappings[pft]]
            x=obs[landscape_ecotypes[econum]][:,pft].mean()
            xerr=obs[landscape_ecotypes[econum]][:,pft].std(ddof=0)
            pftnum=list(arcticpfts_all.PFTnames.to_masked_array()).index(obsdata_PFT_mappings[pft])
            y=mod_new.isel(ecotype=econum,PFT=pftnum).max(dim='time')
            if x>0 or y>0:
                ax.errorbar(x,y,c=arcticpfts_all.PFTcolors[pftnum].item(),ls='None',marker=markers[econum],xerr=xerr,ms=7.0,mfc='None')


    for pft in [11,12]:
        if 'shrub' in mod_soildepth['PFTnames'][pft].item():
            ax=axs['broadleaf_deciduous_boreal_shrub']
            y2=mod_new.isel(ecotype=econum,PFT=mod_new.PFTnames.str.contains('shrub')).sum(dim='PFT').max(dim='time')
        else:
            ax=axs['c3_arctic_grass']

        x=obs_to_defaultPFTs(obs,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,e3sm_all.PFTnames[pft].item()].mean()
        xerr=obs_to_defaultPFTs(obs,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,e3sm_all.PFTnames[pft].item()].std(ddof=0)
        y=mod_soildepth.isel(ecotype=min(econum,len(mod_soildepth.ecotype)-1),PFT=pft).max(dim='time')

        if x>0 or y>0:
            ax.errorbar(x,y,c=e3sm_all.PFTcolors[pft].item(),ls='None',marker=markers[econum],xerr=xerr,ms=8.0,mfc='None')
        if x>0 and 'shrub' in mod_soildepth['PFTnames'][pft].item():
            ax.errorbar(x,y2,c=e3sm_all.PFTcolors[pft].item(),ls='None',marker=markers[econum],xerr=xerr,ms=8.0)
    x=obs_to_defaultPFTs(obs,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,'nonvascular'].mean()
    # y=0.0
    # if x>0:
    #     ax.errorbar(x,y,c=arcticpfts_all.PFTcolors[1].item(),ls='None',marker=markers[econum],xerr=xerr,ms=8.0,alpha=0.5,mfc='None')

for num,pft in enumerate(axs):
    axs[pft].set(title=prettify_pft_name(pft),xlabel='Obs Biomass (gC m$^{-2}$)',ylabel='Mod Biomass (gC m$^{-2}$)',aspect=1.0)
    maxval=max(axs[pft].get_ylim()[1],axs[pft].get_xlim()[1])
    axs[pft].plot([0,maxval],[0,maxval],'k:',lw=0.5)
    axs[pft].text(0.02,0.92,'('+ascii_lowercase[num]+')',transform=axs[pft].transAxes)
    
# Calculate R2 for each plot. Reading the data from the plots themselves since the plotting order was not by PFT
# linestyle and alpha business is to exclude 1-1 line and the pooled nonvascular points in Sim2 grass panel
for pft in axs:
    xdata=ma.masked_invalid([l.get_xdata().astype(float) for l in axs[pft].lines if l.get_linestyle() == 'None' and l.get_alpha() is None])
    ydata=ma.masked_invalid([l.get_ydata().astype(float) for l in axs[pft].lines if l.get_linestyle() == 'None' and l.get_alpha() is None])
    data=ma.row_stack((xdata.ravel(),ydata.ravel()))
    print('%s: R2 = %1.2g, N = %d'%(pft,corrcoef(ma.compress_cols(data))[0,1]**2,data.shape[1]))

handles=[Line2D([0,0],[0,0],ls='None',marker=markers[m],c='k',ms=8.0,mfc='None') for m in ecotypes_included]
labels=[landscape_ecotypes[e] for e in ecotypes_included]
axs['broadleaf_deciduous_boreal_shrub'].legend(handles=handles,labels=labels,fontsize='medium',loc='upper right',ncol=2)

axs['broadleaf_deciduous_boreal_shrub'].set_yticks([0,2000,4000,6000])
axs['c3_arctic_grass'].set_yticks([0,25,50])
axs['arctic_evergreen_shrub_dwarf'].set_xticks([0,150,300])
axs['arctic_evergreen_shrub_dwarf'].set_yticks([0,150,300])
axs['arctic_deciduous_shrub_low'].set_yticks([0,500,1000])
axs['arctic_deciduous_shrub_tall'].set_yticks([0,250,500])
axs['arctic_deciduous_shrub_alder'].set_yticks([0,2500,5000])
axs['arctic_dry_graminoid'].set_xticks([0,20,40,60])
axs['arctic_forb'].set_xticks([0,2,4,6])
axs['arctic_bryophyte'].set_xticks([0,50,100,150])


f,ax=subplots(num='Total biomass 1-1',clear=True)
x=obs.groupby(['Ecotype','PlotID']).sum().groupby('Ecotype').mean()
xerr=obs.groupby(['Ecotype','PlotID']).sum().groupby('Ecotype').std()
y=mod_soildepth.max(dim='time').sum(dim='PFT')
y2=mod_new.max(dim='time').sum(dim='PFT')
for n,e in enumerate(landscape_ecotypes):
    ax.errorbar(x[e],y[n],xerr=xerr[e],marker=markers[n],c='k',ms=8.0,mfc='None')
    ax.errorbar(x[e],y2[n],xerr=xerr[e],marker=markers[n],c='b',ms=8.0,label=e)
maxval=max([x.max(),y.max(),y2.max()]) 
ax.plot([0,maxval],[0,maxval],'k:',lw=0.5)
ax.set(aspect=1.0,xlabel='Obs Biomass (gC m$^{-2}$)',ylabel='Mod Biomass (gC m$^{-2}$)',title='Total biomass comparison')
ax.legend()


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

overlap=True
alpha=0.1
f,axs=subplots(nrows=4,ncols=1,num='Histograms',clear=True,figsize=(3.8,9.8))
shrubs=e3sm_all.PFTnames.data.tolist().index('broadleaf_deciduous_boreal_shrub')
gram=e3sm_all.PFTnames.data.tolist().index('c3_arctic_grass')
if not overlap:
    norm=1/(get_var_PFTs('TOTVEGC',e3sm_all).mean(dim='time').count().item()-2)
else:
    norm=1
ax=axs[1]
bottom=plot_hist(get_var_PFTs('TOTVEGC',e3sm_all).isel(PFT=shrubs).mean(dim='time'),ax,linestyle='-',norm=norm,zorder=1,alpha=alpha)
if overlap:
    bottom[:]=0
bottom=plot_hist(get_var_PFTs('TOTVEGC',e3sm_all).isel(PFT=gram).mean(dim='time'),ax,linestyle='--',bottom=bottom,norm=norm,alpha=alpha)
ax.set_title('Simulation 1')
ax.legend(handles=ax.collections[::-1])
ax.set_ylabel('Fraction of patches')

shrubs=communities_all.PFTnames.data.tolist().index('broadleaf_deciduous_boreal_shrub')
gram=communities_all.PFTnames.data.tolist().index('c3_arctic_grass')
if not overlap:
    norm=1/(get_var_PFTs('TOTVEGC',communities_all).mean(dim='time').count().item()-2)
else:
    norm=1/len(landscape_ecotypes)
ax=axs[2]
if overlap:
    bottom[:]=0
bottom=plot_hist(get_var_PFTs('TOTVEGC',communities_all).isel(PFT=shrubs).mean(dim='time'),ax,linestyle='-',norm=norm,zorder=1,alpha=alpha)
if overlap:
    bottom[:]=0
bottom=plot_hist(get_var_PFTs('TOTVEGC',communities_all).isel(PFT=gram).mean(dim='time'),ax,linestyle='--',bottom=bottom,norm=norm,alpha=alpha)
ax.set_title('Simulation 2')

ax.set_ylabel('Fraction of patches')

pft_order_meas=['bryophyte', 'dwarf shrub deciduous', 'dwarf shrub evergreen',
       'forb', 'graminoid', 'lichen', 'low shrub deciduous', 'mixed',
       'potential tall shrub deciduous non-alder','potential tall shrub deciduous alder']

pft_order=[pft_names.index(obsdata_PFT_mappings[name]) for name in pft_order_meas]  
if not overlap:
    norm=1/(get_var_PFTs('TOTVEGC',arcticpfts_all).mean(dim='time').count().item()-1)
else:
    norm=1/len(landscape_ecotypes)
bottom[:]=0
ax=axs[3]
for pft in pft_order:
    if 'shrub' in pft_names[pft] and pft_names.index(pft_names[pft]) in pfts_inuse:
        bottom=plot_hist(get_var_PFTs('TOTVEGC',arcticpfts_all).isel(PFT=pft).mean(dim='time'),ax=ax,bottom=bottom,linestyle='-',norm=norm,zorder=2,alpha=alpha)
        if overlap:
            bottom[:]=0

if overlap:
    bottom[:]=0
for pft in pft_order:
    if 'forb' in pft_names[pft] or 'graminoid' in pft_names[pft] and pft_names.index(pft_names[pft]) in pfts_inuse:
        bottom=plot_hist(get_var_PFTs('TOTVEGC',arcticpfts_all).isel(PFT=pft).mean(dim='time'),ax=ax,bottom=bottom,linestyle='--',norm=norm,alpha=alpha)
        if overlap:
            bottom[:]=0
ax.set_title('Arctic PFTs')

if overlap:
    bottom[:]=0
for pft in pft_order:
    if 'bryophyte' in pft_names[pft] or 'lichen' in pft_names[pft] and pft_names.index(pft_names[pft]) in pfts_inuse:
        bottom=plot_hist(get_var_PFTs('TOTVEGC',arcticpfts_all).isel(PFT=pft).mean(dim='time'),ax=ax,bottom=bottom,linestyle=':',norm=norm,zorder=-1,alpha=alpha)
        if overlap:
            bottom[:]=0
ax.set_title('Simulation 3')

ax.set_ylabel('Fraction of patches')

totalC=meas_leaf_C+meas_rhizome_C+meas_root_C+meas_stem_C
if not overlap:
    norm=1/len(totalC)
else:
    norm=1/len(unique(totalC.index.get_level_values('PlotID')))
bottom[:]=0
ax=axs[0]
for pft in pft_order_meas:
    if 'shrub' in pft:
        bottom=plot_hist(totalC[:,:,pft],ax=ax,bottom=bottom,color=pft_colors[pft_names.index(obsdata_PFT_mappings[pft])],label='',linestyle='-',norm=norm,zorder=2,alpha=alpha)
        if overlap:
            bottom[:]=0

if overlap:
    bottom[:]=0
for pft in pft_order_meas:
    if 'graminoid' in pft or 'forb' in pft:
        bottom=plot_hist(totalC[:,:,pft],ax=ax,bottom=bottom,color=pft_colors[pft_names.index(obsdata_PFT_mappings[pft])],label='',linestyle='--',norm=norm,alpha=alpha)
        if overlap:
            bottom[:]=0

if overlap:
    bottom[:]=0
for pft in pft_order_meas:
    if 'bryophyte' in pft or 'lichen' in pft:
        bottom=plot_hist(meas_nonvasc_C[:,:,pft],ax=ax,bottom=bottom,color=pft_colors[pft_names.index(obsdata_PFT_mappings[pft])],label='',linestyle=':',norm=norm,zorder=-1,alpha=alpha)
        if overlap:
            bottom[:]=0

ax.set_title('Measurements')
ax.set_ylabel('Fraction of plots')
ax.legend(handles=axs[3].collections[::-1])

for num,ax in enumerate(axs.ravel()):
    ax.set_xlabel('Biomass (gC m$^{-2}$)')
    ax.text(0.02,1.05,ascii_lowercase[num],transform=ax.transAxes,fontsize='large')


f,axs=subplots(nrows=2,ncols=3,num='Total biomass (Arctic PFTs)',clear=True,figsize=(12,6))
f.set_constrained_layout_pads(hspace=0.18)
for pft in mod_new.PFT.values[:-1]:
    for eco in range(6):
        mod_new.isel(PFT=pft,ecotype=eco).plot(ax=axs.ravel()[eco],c=pft_colors[pft],label=prettify_pft_name(pft_names[pft]))
l=axs[0,0].legend(loc=(0,-0.4),ncol=5)
l.set_in_layout(False)
for eco in range(6):
    ax=axs.ravel()[eco]
    ax.set_title(landscape_ecotypes[eco])
    ax.set_xlabel('Year')
    ax.set_ylabel('Biomass (g C m$^{-2}$)')

f,axs=subplots(nrows=2,ncols=3,num='Total biomass (Default PFTs)',clear=True,figsize=(12,6))
f.set_constrained_layout_pads(hspace=0.18)
for pft in mod_soildepth.PFT.values:
    for eco in range(6):
        if mod_soildepth.isel(PFT=pft,ecotype=eco).max()>0.0:
            mod_soildepth.isel(PFT=pft,ecotype=eco).plot(ax=axs.ravel()[eco],c=pft_colors_default[pft],label=prettify_pft_name(pft_names_default[pft]))
l=axs[0,0].legend(loc=(0,-0.4),ncol=5)
l.set_in_layout(False)
for eco in range(6):
    ax=axs.ravel()[eco]
    ax.set_title(landscape_ecotypes[eco])
    ax.set_xlabel('Year')
    ax.set_ylabel('Biomass (g C m$^{-2}$)')


f,a=subplots(num='N fixation',figsize=(6.4,6.7),clear=True)
x=arange(6)
nfix_obs={'AS':1.95,'ASV':0.53}
nfix_obs_error={'AS':0.68,'ASV':0.19}
nfix=e3sm_all['NFIX_TO_SMINN'].dropna(dim='column').mean(dim='time')*365*24*3600
bar(x,nfix,label='N fixation (default)',hatch='//',width=0.3,fc='C0')
bar(x,e3sm_all['NDEP_TO_SMINN'].dropna(dim='column').mean(dim='time')*365*24*3600,bottom=nfix,label='N deposition (default)',hatch='//',width=0.3,fc='C1')
bar(x+0.3,arcticpfts_all['NFIX_TO_SMINN'].dropna(dim='column').mean(dim='time')*365*24*3600,label='N fixation (updated)',hatch='o',width=0.3,fc='C0')
bar(x+0.3,arcticpfts_all['NDEP_TO_SMINN'].dropna(dim='column').mean(dim='time')*365*24*3600,bottom=arcticpfts_all['NFIX_TO_SMINN'].dropna(dim='column').mean(dim='time')*365*24*3600,label='N deposition (updated)',hatch='o',width=0.3,fc='C1')
bar(x+0.6,[nfix_obs.get(e,nan) for e in landscape_ecotypes], yerr=[nfix_obs_error.get(e,nan) for e in landscape_ecotypes]  ,width=0.3,label='Obs N fixation',fc='C0')
title('N deposition and inputs')
xticks(x,landscape_ecotypes)
ylabel('N input rate (gN m$^{-2}$ year$^{-1}$)')
legend(loc='upper left',fontsize='medium')