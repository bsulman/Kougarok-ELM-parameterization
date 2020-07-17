from plot_kougarok import *

import warnings

# Don't need to print this warning a billion times
warnings.filterwarnings(action='ignore',message='All-NaN slice encountered')  
warnings.filterwarnings(action='ignore',message='Mean of empty slice')  


data_soilthickness=xarray.open_dataset('../../output_data/E3SMpfts_soilthickness_processed_20200316.nc')
data_global=xarray.open_dataset('../../output_data/E3SMpfts_processed_20200316.nc')
data_Arcticpfts=xarray.open_dataset('../../output_data/Arcticpfts_processed_20200323.nc')

meas_leaf_C.loc[:,:,'bryophyte'] = 0.0
meas_leaf_C.loc[:,:,'lichen'] = 0.0

f,axs=subplots(nrows=3,ncols=6,num='Ecosystem biomass',figsize=(15, 10),clear=True)
for simnum,data_to_plot in enumerate([data_global,data_Arcticpfts]):
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
        
        h=plot_mod_bar_stack(x,leafc[:,:,~nonvasc],min(econum,len(data_to_plot.ecotype)-1),width=0.8,hatch='//')
        # plot_obs_bar_stack(x,meas_leaf_C ,econum,width=0.8)
        names.append('Leaf')
        # text(x+0.4,h[-1][0].get_y()+h[-1][0].get_height()+50,'Modeled value',rotation=90,va='bottom',ha='center')
        # text(x,h[-1][0].get_y()+h[-1][0].get_height()+50,'Measured value',rotation=90,va='bottom',ha='center')
        x+=1

        plot_mod_bar_stack(x,get_var_PFTs(['LIVESTEMC','DEADSTEMC'],data_to_plot),min(econum,len(data_to_plot.ecotype)-1),width=0.8,hatch='//')
        # plot_obs_bar_stack(x,meas_stem_C,econum,width=0.8)
        names.append('Stem')
        x+=1

        plot_mod_bar_stack(x,get_var_PFTs('FROOTC',data_to_plot),min(econum,len(data_to_plot.ecotype)-1),width=0.8,hatch='//')
        # bar(x,meas_root_C[landscape_ecotypes[econum]].mean(),yerr=meas_root_C[landscape_ecotypes[econum]].std(),width=0.8,facecolor=[0.9,0.9,0.9],edgecolor='k')
        names.append('FRoot')
        x+=1

        plot_mod_bar_stack(x,get_var_PFTs(['LIVECROOTC','DEADCROOTC'],data_to_plot),min(econum,len(data_to_plot.ecotype)-1),width=0.8,hatch='//')
        # text(x+0.4,100,'? ? ? ? ? ?',fontsize='large',rotation=90,va='bottom',ha='center')
        # plot_obs_bar_stack(x,meas_rhizome_C,econum,width=0.8)
        names.append('Croot')
        x+=1
        
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
    ax=axs[2,econum]
    sca(ax)

    # h=plot_mod_bar_stack(x+0.4,get_var_PFTs('LEAFC',data_to_plot),econum,width=0.4,hatch='//')
    plot_obs_bar_stack(x,meas_nonvasc_C ,econum,width=0.8)
    names.append('Nonvasc')
    # text(x+0.4,h[-1][0].get_y()+h[-1][0].get_height()+50,'Modeled value',rotation=90,va='bottom',ha='center')
    # text(x,h[-1][0].get_y()+h[-1][0].get_height()+50,'Measured value',rotation=90,va='bottom',ha='center')
    x+=1
    
    plot_obs_bar_stack(x,meas_leaf_C ,econum,width=0.8)
    names.append('Leaf')
    # text(x+0.4,h[-1][0].get_y()+h[-1][0].get_height()+50,'Modeled value',rotation=90,va='bottom',ha='center')
    # text(x,h[-1][0].get_y()+h[-1][0].get_height()+50,'Measured value',rotation=90,va='bottom',ha='center')
    x+=1

    # plot_mod_bar_stack(x+0.4,get_var_PFTs(['LIVESTEMC','DEADSTEMC'],data_to_plot),econum,width=0.4,hatch='//')
    plot_obs_bar_stack(x,meas_stem_C,econum,width=0.8)
    names.append('Stem')
    x+=1

    # plot_mod_bar_stack(x+0.4,get_var_PFTs('FROOTC',data_to_plot),econum,width=0.4,hatch='//')
    bar(x,meas_root_C[landscape_ecotypes[econum]].mean(),yerr=meas_root_C[landscape_ecotypes[econum]].std(),width=0.8,facecolor=[0.9,0.9,0.9],edgecolor='k')
    names.append('FRoot')
    x+=1

    # plot_mod_bar_stack(x+0.4,get_var_PFTs(['LIVECROOTC','DEADCROOTC'],data_to_plot),econum,width=0.4,hatch='//')
    # text(x+0.4,100,'? ? ? ? ? ?',fontsize='large',rotation=90,va='bottom',ha='center')
    plot_obs_bar_stack(x,meas_rhizome_C,econum,width=0.8)
    names.append('Rhizome')
    x+=1
    
    # plot_mod_bar_stack(x+0.4,get_var_PFTs(['STORVEGC'],data_to_plot),econum,width=0.4,hatch='//')
    # names.append('Storage')

    ylabel('Biomass (gC m$^{-2}$)')
    xticks(arange(len(names)),names)
    ylim(-50,1050)
    # title(ecotype_names_list[econum])
    # ylim(0,3000)
    # ax.set_xlim(right=x+0.6)

axs[2,2].set_ylim(top=3050)

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
l=axs[1,0].legend(handles=handles,fontsize='small',ncol=1,title='Arctic PFTs',loc='upper right')

for econum in range(len(ecotype_names_list)):
    axs[0,econum].set_title(ecotype_names_list[econum])
    
axs[0,0].text(-0.38,0.5,'E3SM grid cell',rotation=90,ha='center',va='center',fontsize='large',transform=axs[0,0].transAxes)
# axs[1,0].text(-0.38,0.5,'Level 2:\nE3SM PFTs, site soil depths',rotation=90,ha='center',va='center',fontsize='large',transform=axs[1,0].transAxes)
# axs[2,0].text(-0.38,0.5,'Level 3:\nE3SM PFTs, site areas',rotation=90,ha='center',va='center',fontsize='large',transform=axs[2,0].transAxes)
axs[1,0].text(-0.38,0.5,'Arctic PFTs',rotation=90,ha='center',va='center',fontsize='large',transform=axs[1,0].transAxes)
axs[2,0].text(-0.38,0.5,'Site measurements',rotation=90,ha='center',va='center',fontsize='large',transform=axs[2,0].transAxes)


f,axs=subplots(2,1,num='Model-data comparison',clear=True,squeeze=False)
# Aboveground biomass
obs_ag_C = (meas_leaf_C + meas_stem_C).add(meas_nonvasc_C,fill_value=0.0)
obs_bg_C = meas_rhizome_C.groupby(level=('Ecotype','PlotID')).sum() + meas_root_C
ag_mod_C_arctic = get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC'],data_Arcticpfts)
bg_mod_C_arctic = get_var_PFTs(['LIVECROOTC','DEADCROOTC','FROOTC'],data_Arcticpfts)
ag_mod_C_global = get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC'],data_global)
bg_mod_C_global = get_var_PFTs(['LIVECROOTC','DEADCROOTC','FROOTC'],data_global)

markers=['o','s','x','^','>','*']

for econum in range(len(ecotype_names_list)):
    for pft in unique(obs_ag_C.index.get_level_values('ELM_PFT')):
        if pft in obs_ag_C[landscape_ecotypes[econum]].index.get_level_values('ELM_PFT') :
            x=obs_ag_C[landscape_ecotypes[econum]][:,pft].mean()
            xerr=obs_ag_C[landscape_ecotypes[econum]][:,pft].std()
            pftnum=list(data_Arcticpfts.PFTnames.to_masked_array()).index(obsdata_PFT_mappings[pft])
            y=ag_mod_C_arctic.isel(ecotype=econum,PFT=pftnum).max(dim='time')
            axs[1,0].errorbar(x,y,c=data_Arcticpfts.PFTcolors[pftnum].item(),ls='None',marker=markers[econum],xerr=xerr,ms=5.0)
            
    # x=obs_bg_C[landscape_ecotypes[econum]].mean()
    # xerr=obs_bg_C[econum].std()
    # y=bg_mod_C_arctic.max(dim='time').sum(dim='PFT').isel(ecotype=econum)
    # axs[1,1].errorbar(x,y,c='k',ls='None',marker=markers[econum],xerr=xerr)
    # 
    # y=bg_mod_C_global.max(dim='time').sum(dim='PFT').isel(ecotype=0)
    # axs[0,1].errorbar(x,y,c='k',ls='None',marker=markers[econum],xerr=xerr)
    # 
    for pft in [11,12]:
        x=obs_to_defaultPFTs(obs_ag_C,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,data_global.PFTnames[pft].item()].mean()
        xerr=obs_to_defaultPFTs(obs_ag_C,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,data_global.PFTnames[pft].item()].std()
        y=ag_mod_C_global.isel(ecotype=min(econum,len(ag_mod_C_global.ecotype)-1),PFT=pft).max(dim='time')
        axs[0,0].errorbar(x,y,c=data_global.PFTcolors[pft].item(),ls='None',marker=markers[econum],xerr=xerr,ms=5.0)

axs[0,0].plot(linspace(0,800,10),linspace(0,800,10),'k:')
axs[1,0].plot(linspace(0,800,10),linspace(0,800,10),'k:')
axs[0,0].set_title('Default E3SM')
axs[1,0].set_title('Arctic parameterized PFTs')
axs[0,0].set_xlabel('Observed biomass')
axs[1,0].set_xlabel('Observed biomass')
axs[0,0].set_ylabel('Modeled biomass')
axs[1,0].set_ylabel('Modeled biomass')

handles=[Line2D([0,0],[0,0],ls='None',marker='o',ms=8.0,c=data_Arcticpfts.PFTcolors[pft].item()) for pft in range(1,11)]
labels=[prettify_pft_name(data_Arcticpfts.PFTnames[pft].item()) for pft in range(1,11)]
handles=handles+[Line2D([0,0],[0,0],ls='None',marker=m,c='k',ms=8.0) for m in markers]
labels=labels+ecotype_names_list
axs[1,0].legend(handles=handles,labels=labels)



f,axs=subplots(1,1,num='Model-data comparison 2',clear=True,squeeze=False)
# Aboveground biomass
obs_ag_C = (meas_leaf_C + meas_stem_C).add(meas_nonvasc_C,fill_value=0.0)
obs_bg_C = meas_rhizome_C.groupby(level=('Ecotype','PlotID')).sum() + meas_root_C
ag_mod_C_arctic = get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC'],data_Arcticpfts)
bg_mod_C_arctic = get_var_PFTs(['LIVECROOTC','DEADCROOTC','FROOTC'],data_Arcticpfts)
ag_mod_C_global = get_var_PFTs(['LIVESTEMC','DEADSTEMC','LEAFC'],data_global)
bg_mod_C_global = get_var_PFTs(['LIVECROOTC','DEADCROOTC','FROOTC'],data_global)

markers=['o','s','x','^','>','*']

ecotypes_included=[4,5]
for econum in ecotypes_included:#range(len(ecotype_names_list)):
    for pft in unique(obs_ag_C.index.get_level_values('ELM_PFT')):
        if pft in obs_ag_C[landscape_ecotypes[econum]].index.get_level_values('ELM_PFT') :
            x=obs_ag_C[landscape_ecotypes[econum]][:,pft].mean()
            xerr=obs_ag_C[landscape_ecotypes[econum]][:,pft].std()
            pftnum=list(data_Arcticpfts.PFTnames.to_masked_array()).index(obsdata_PFT_mappings[pft])
            y=ag_mod_C_arctic.isel(ecotype=econum,PFT=pftnum).max(dim='time')
            axs[0,0].errorbar(x,y,c=data_Arcticpfts.PFTcolors[pftnum].item(),ls='None',marker=markers[econum],xerr=xerr,ms=5.0)
            
            x_old=obs_to_defaultPFTs(obs_ag_C,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,obsdata_E3SMPFT_mappings[pft]].mean()
            if obsdata_E3SMPFT_mappings[pft] in pft_names_default:
                y_old=ag_mod_C_global.isel(ecotype=min(econum,len(ag_mod_C_global.ecotype)-1),PFT=pft_names_default.index(obsdata_E3SMPFT_mappings[pft])).max(dim='time')
            else:
                y_old=0.0
            annotate('',(x,y),(x_old,y_old),arrowprops=dict(arrowstyle="->",color=data_Arcticpfts.PFTcolors[pftnum].item()))
            
            
    # x=obs_bg_C[landscape_ecotypes[econum]].mean()
    # xerr=obs_bg_C[econum].std()
    # y=bg_mod_C_arctic.max(dim='time').sum(dim='PFT').isel(ecotype=econum)
    # axs[1,1].errorbar(x,y,c='k',ls='None',marker=markers[econum],xerr=xerr)
    # 
    # y=bg_mod_C_global.max(dim='time').sum(dim='PFT').isel(ecotype=0)
    # axs[0,1].errorbar(x,y,c='k',ls='None',marker=markers[econum],xerr=xerr)
    # 
    for pft in [11,12]:
        x=obs_to_defaultPFTs(obs_ag_C,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,data_global.PFTnames[pft].item()].mean()
        xerr=obs_to_defaultPFTs(obs_ag_C,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,data_global.PFTnames[pft].item()].std()
        y=ag_mod_C_global.isel(ecotype=min(econum,len(ag_mod_C_global.ecotype)-1),PFT=pft).max(dim='time')
        axs[0,0].errorbar(x,y,c=data_global.PFTcolors[pft].item(),ls='None',marker=markers[econum],xerr=xerr,ms=8.0)
    x=obs_to_defaultPFTs(obs_ag_C,obsdata_defaultPFT_mappings)[landscape_ecotypes[econum]][:,'nonvascular'].mean()
    y=0.0
    axs[0,0].errorbar(x,y,c=data_Arcticpfts.PFTcolors[1].item(),ls='None',marker=markers[econum],xerr=xerr,ms=8.0,alpha=0.5,mfc='w')

axs[0,0].plot(linspace(0,800,10),linspace(0,800,10),'k:')
axs[0,0].plot(linspace(0,800,10),linspace(0,800,10),'k:')
axs[0,0].set_title('Obs vs modeled biomass',fontsize='large')
axs[0,0].set_xlabel('Observed biomass (g m$^{-2}$)',fontsize='large')
axs[0,0].set_ylabel('Modeled biomass (g m$^{-2}$)',fontsize='large')

handles=[Line2D([0,0],[0,0],ls='None',marker='o',ms=8.0,c=data_Arcticpfts.PFTcolors[pft].item()) for pft in range(1,11)]
labels=[prettify_pft_name(data_Arcticpfts.PFTnames[pft].item()) for pft in range(1,11)]
handles=handles+[Line2D([0,0],[0,0],ls='None',marker=markers[m],c='k',ms=8.0) for m in ecotypes_included]
labels.extend(ecotype_names_list[e] for e in ecotypes_included)
axs[0,0].legend(handles=handles,labels=labels,fontsize='large')



show()