figure('Obs biomass');clf()
for econum in range(len(ecotype_names_list)):
    names=[]
    x=0.0
    ax=subplot(1,6,econum+1)

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
    ylim(-50,3050)
    # title(ecotype_names_list[econum])
    # ylim(0,3000)
    # ax.set_xlim(right=x+0.6)

    tight_layout()
    
data_global=read_pftfile('../../output_data/E3SMpfts_20191122_h2.nc',maxyear=None)
data_oldpfts=read_pftfile('../../output_data/E3SMpfts_communities_20191122_h2.nc',maxyear=None)
# data_newpfts=read_pftfile('../../output_data/Arcticpfts_20191122_h2.nc',maxyear=None)
data_newpfts=read_pftfile('../../output_data/soildepth_Nfix_hist_h1_20190911.nc',maxyear=None)

data_to_plot=data_global
figure('Mod biomass',figsize=(25.19, 3.76));clf()
for econum in range(len(ecotype_names_list)):
    names=[]
    x=0.0
    ax=subplot(1,6,econum+1)
    
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
    ylim(-50,3050)
    ax.set_xlim(left=-0.64)
    # title(ecotype_names_list[econum])
    # ylim(0,3000)
    # ax.set_xlim(right=x+0.6)

    tight_layout()
    
    
figure('N fixation',figsize=(6.4,6.7));clf()
subplot(111)
columndata=xarray.open_dataset('../../output_data/E3SMpfts_20191122_h0.nc')
columndata_new=xarray.open_dataset('../../output_data/Arcticpfts_20191122_h0.nc')
x=arange(6)
nfix_obs={'AS':1.95,'ASV':0.53}
nfix_obs_error={'AS':0.68,'ASV':0.19}
nfix=columndata['NFIX_TO_SMINN'].sel(time=t_col>(t_col.max()-100)).mean(dim='time')*365*24*3600
nfix_new=columndata_new['NFIX_TO_SMINN'].sel(time=t_col>(t_col.max()-100)).mean(dim='time')*365*24*3600
bar(x,nfix,label='N fixation (default)',hatch='//',width=0.3)
bar(x+0.3,nfix_new,label='N fixation (updated)',hatch='//',width=0.3)
# bar(x,columndata['NDEP_TO_SMINN'].sel(time=t_col>(t_col.max()-100)).mean(dim='time')*365*24*3600,bottom=nfix,label='N deposition',hatch='//',width=0.4)
bar(x+0.6,[nfix_obs.get(e,nan) for e in landscape_ecotypes], yerr=[nfix_obs_error.get(e,nan) for e in landscape_ecotypes]  ,width=0.3,label='Obs N fixation')
title('N deposition and inputs')
xticks(x,landscape_ecotypes)
ylabel('N input rate (gN m$^{-2}$ year$^{-1}$)')
legend()
tight_layout()

figure('Non-vascular');clf()
meas_lichen_BM=Koug_meas_biomass['Nonvascular_gperm2'][:,:,'lichen']
meas_moss_BM=Koug_meas_biomass['Nonvascular_gperm2'][:,:,'bryophyte']
# meas_other_nonvasc_BM=Koug_meas_biomass['Nonvascular_gperm2'][:,'other']
meas_vasc_AG_BM=(Koug_meas_biomass['Leaf_gperm2']+Koug_meas_biomass['Stem_gperm2']).sum(level='PlotID')

w=0.2
for econum in range(len(landscape_ecotypes)):
    plotIDs=meas_lichen_BM[landscape_ecotypes[econum]].index.get_level_values('PlotID')
    for plotnum,plotID in enumerate(plotIDs):
        bottom=0.0
        y=meas_lichen_BM[landscape_ecotypes[econum]][plotID]
        lichenbar=bar(econum+plotnum*w/len(plotIDs),y,bottom=bottom,width=w/len(plotIDs),label='Lichen',color=pft_colors[pft_names.index('arctic_lichen')],align='edge')
        bottom=bottom+y
        y=meas_moss_BM[landscape_ecotypes[econum]][plotID]
        mossbar=bar(econum+plotnum*w/len(plotIDs),y,bottom=bottom,width=w/len(plotIDs),label='Moss',color=pft_colors[pft_names.index('arctic_bryophyte')],align='edge')
        y=meas_vasc_AG_BM[plotID]
        vascbar=bar(econum+w+plotnum*w/len(plotIDs),y,facecolor='brown',width=w/len(plotIDs),label='Vascular',align='edge')


title('Aboveground vascular and non-vascular biomass')
ylabel('Total biomass (g DW m$^{-2}$)')
l=legend(handles=(mossbar,lichenbar,vascbar))
l.set_draggable(True)
xticks(arange(len(landscape_ecotypes))+w,[ecotype_names[name] for name in landscape_ecotypes] ,rotation=45,ha='right')
ylim(0,1600)
tight_layout()


