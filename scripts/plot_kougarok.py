from kougarok_plotting import *


if __name__=='__main__':
    
    import sys
    if len(sys.argv)==2:
        filename=sys.argv[1]
    else:
        raise RuntimeError('Must supply data file name')

    import os
    dataname=os.path.basename(filename)
    vegdata_PFTs=read_pftfile(filename,maxyear=400)
    
    columndata=xarray.open_dataset(filename.replace('h1','h0'),autoclose=True)
    Tsoil10cm=columndata['TSOI_10CM']

    minyear=0
    maxyear=300
    
    Koug_meas_biomass.rename(index={'NAMC':'DL','TTWBT':'ASV'},inplace=True)
    Koug_meas_chem.rename(index={'NAMC':'DL','TTWBT':'ASV'},inplace=True)

    meas_leaf_C=(Koug_meas_biomass['LeafBiomass_gperm2']*Koug_meas_chem['LeafC_percent']/100).groupby('Ecotype').sum()
    meas_nonvasc_C=(Koug_meas_biomass['NonvascularBiomass_gperm2']*0.5).groupby('Ecotype').sum()
    # Should this include nonvascular biomass?
    obsdata=meas_leaf_C+meas_nonvasc_C
    # for ecotype in landscape_ecotypes:
    #     obsdata=meas_leaf_C+meas_nonvasc_C
    #     defaultparams_ax.plot([mindate,maxdate],[obsdata[ecotype],obsdata[ecotype]],c='C%d'%(landscape_ecotypes.index(ecotype)),ls=':')
    #     newparams_ax.plot([mindate,maxdate],[obsdata[ecotype],obsdata[ecotype]],c='C%d'%(landscape_ecotypes.index(ecotype)),ls=':')

    meas_root_C=Koug_meas_biomass['FineRootBiomass_gperm2'][:,'mixed']*(Koug_meas_chem['FineRootC_percent']/100).groupby('Ecotype').mean()
    # Should rhizomes be treated as fine or coarse roots?

    meas_stem_C=(Koug_meas_biomass['StemBiomass_gperm2']*Koug_meas_chem['StemC_percent']/100).groupby('Ecotype').sum()
    # Check this. WOODC appears to be sum of LIVESTEMC+DEADSTEMC+LIVECROOTC+DEADCROOTC

    # Try comparing rhizomes to LIVECROOTC
    meas_rhizome_C=(Koug_meas_biomass['RhizomeBiomass_gperm2']*Koug_meas_chem['RhizomeC_percent']/100).groupby('Ecotype').sum()


    meas_leaf_C=(Koug_meas_biomass['LeafBiomass_gperm2']*Koug_meas_chem['LeafC_percent']/100).copy()
    meas_nonvasc_C=(Koug_meas_biomass['NonvascularBiomass_gperm2']*0.5)
    # Should this include nonvascular biomass?
    for ecotype in landscape_ecotypes:
        meas_leaf_C[ecotype,'moss']=meas_nonvasc_C.loc[ecotype,'moss']
        meas_leaf_C[ecotype,'lichen']=meas_nonvasc_C.loc[ecotype,'lichen']

    meas_root_C=Koug_meas_biomass['FineRootBiomass_gperm2'][:,'mixed']*(Koug_meas_chem['FineRootC_percent']/100)
    meas_stem_C=(Koug_meas_biomass['StemBiomass_gperm2']*Koug_meas_chem['StemC_percent']/100)
    meas_rhizome_C=(Koug_meas_biomass['RhizomeBiomass_gperm2']*Koug_meas_chem['RhizomeC_percent']/100)
    meas_NPP_C=(Koug_meas_biomass['LeafNPP_gperm2peryr']*Koug_meas_chem['LeafC_percent']/100)+\
        (Koug_meas_biomass['StemNPP_gperm2peryr']*Koug_meas_chem['StemC_percent']/100)#+\

    meas_froot_NPP =(Koug_meas_biomass['FineRootNPP_gperm2peryr']*Koug_meas_chem['FineRootC_percent']/100).groupby('Ecotype').sum()

    plotvars=['leaf','froot','croot','store','stem','cnpp','height','downreg']
    nplots=len(plotvars)
        
    year=xarray.DataArray(array([tt.year for tt in vegdata_PFTs.time.data]),dims='time',coords={'time':vegdata_PFTs.time},name='year')
    annualmax=vegdata_PFTs.groupby(year,squeeze=False).max(dim='time',keep_attrs=True)
    annualmax['time']=vegdata_PFTs['time'].groupby(year).first()
        
    fig=figure('Ecotype time series (%s)'%dataname,figsize=(15,10))
    fig.clf()
    gs=fig.add_gridspec(ncols=nplots,nrows=len(landscape_ecotypes))
    subplot_handles={}
    

    for econum in range(len(landscape_ecotypes)):
        ecotype=landscape_ecotypes[econum]
        subplot_handles[ecotype]={}
        for var in plotvars:
            subplot_handles[ecotype][var]=fig.add_subplot(gs[econum,plotvars.index(var)])

        plot_var_PFTs('LEAFC',annualmax,obsdata=meas_leaf_C,ecotype_num=econum,ax=subplot_handles[ecotype]['leaf'],maxyear=maxyear)
        subplot_handles[ecotype]['leaf'].text(-1.0,0.5,ecotype,rotation=90,va='center',ha='center',transform=subplot_handles[ecotype]['leaf'].transAxes)
        # Should rhizomes be treated as fine or coarse roots?
        plot_var_PFTs('FROOTC',annualmax,obsdata=meas_root_C,plotsum=True,ecotype_num=econum,ax=subplot_handles[ecotype]['froot'],maxyear=maxyear)

        plot_var_PFTs(['LIVECROOTC','DEADCROOTC'],vegdata_PFTs,longname='C Root',plotsum=True,ecotype_num=econum,ax=subplot_handles[ecotype]['croot'],maxyear=maxyear)

        plot_var_PFTs(['LIVESTEMC','DEADSTEMC'],vegdata_PFTs,obsdata=meas_stem_C,longname='Stem C',ecotype_num=econum,ax=subplot_handles[ecotype]['stem'],maxyear=maxyear)

        # plot_var_PFTs('NPP',vegdata_PFTs,longname='NPP',modfactor=3600*24,units='gC m$^{-2}$ day $^{-1}$',ecotype_num=econum,ax=subplot_handles[ecotype]['npp'],maxyear=maxyear)

        dat=plot_var_PFTs('AGNPP',vegdata_PFTs,longname='AGNPP',modfactor=3600*24*365,units='gC m$^{-2}$',cumulative=True,obsdata=meas_NPP_C,maxyear=maxyear,minyear=maxyear-1,ecotype_num=econum,ax=subplot_handles[ecotype]['cnpp'])
        t=array([tt.year + (tt.dayofyr-1)/365 for tt in dat['time'].data])
        datmax=dat.isel(time=nonzero(t<=maxyear)[0][-1],ecotype=econum).max()
        subplot_handles[ecotype]['cnpp'].set_ylim(bottom=-1,top=max(datmax*1.1,meas_NPP_C[ecotype].max()*1.1))

        plot_var_PFTs('HTOP',vegdata_PFTs,weight_area=False,ecotype_num=econum,ax=subplot_handles[ecotype]['height'],maxyear=maxyear)

        plot_var_PFTs('STORVEGC',vegdata_PFTs,ecotype_num=econum,longname='Stored C',obsdata=meas_rhizome_C,ax=subplot_handles[ecotype]['store'],maxyear=maxyear)

        plot_var_PFTs('DOWNREG',annualmax,weight_area=False,ecotype_num=econum,ax=subplot_handles[ecotype]['downreg'],maxyear=maxyear,longname='Nut. lim',units='frac')

        for ax in subplot_handles[ecotype].values():
            ax.title.set_fontsize('small')
            ax.xaxis.label.set_fontsize('small')
            ax.yaxis.label.set_fontsize('small')
            for label in ax.get_xticklabels():
                label.set_fontsize('small')
            for label in ax.get_yticklabels():
                label.set_fontsize('small')

    

    tight_layout(rect=(0.03,0.0,1.0,0.95))
    subplots_adjust(hspace=0.8)
    subplot_handles[landscape_ecotypes[0]]['froot'].legend(loc=(-1,1.2),ncol=7,fontsize='small')

    # Plot some biomass ratios

    meas_leaf_C=(Koug_meas_biomass['LeafBiomass_gperm2']*Koug_meas_chem['LeafC_percent']/100)
    meas_root_C=(Koug_meas_biomass['FineRootBiomass_gperm2'][:,'mixed']*(Koug_meas_chem['FineRootC_percent']/100))[:,'mixed']
    meas_stem_C=(Koug_meas_biomass['StemBiomass_gperm2']*(Koug_meas_chem['StemC_percent']/100))

    obs_leafCN = Koug_meas_chem['LeafC_percent']/Koug_meas_chem['LeafN_percent']
    obs_stemCN = Koug_meas_chem['StemC_percent']/Koug_meas_chem['StemN_percent']
    obs_frootCN = 54.6/1.3 # Obs uses a single value for fine roots
    # Measured tissue C:N ratios are in weight units, and model expects gC/gN so they should be comparable.


    # For now: let's assume that relative amount of leaf biomass is proportional to relative amount of root biomass
    # But we may want to change to a different approach like PFT % coverage
    # Or does it make more sense to calculate these at the ecotype scale?
    # Also: This param is really defined to be NPP ratio, not biomass ratio. Should we parameterize using NPP measurements instead?
    leafCfrac=meas_leaf_C/meas_leaf_C.groupby('Ecotype').sum()

    def froot_leaf(ecotype,pft):
        froot_leaf=meas_root_C[ecotype]*leafCfrac[ecotype][pft]/meas_leaf_C[ecotype][pft]
        #print('{ecotype:s}: {pft:s} leaf frac {leafCfrac:1.2f}, froot_leaf = {froot_leaf:1.2f}'.format(leafCfrac=leafCfrac[ecotype][pft],ecotype=ecotype,pft=pft,froot_leaf=froot_leaf))
        return froot_leaf

    def leaf_stem(ecotype,pft):
        return meas_stem_C[ecotype]/meas_leaf_C[ecotype][pft]

    leaf=vegdata_PFTs['LEAFC_unweighted']
    stemleafratio_newparams=xarray.Dataset({'stem_leaf_ratio_unweighted':(vegdata_PFTs['LIVESTEMC_unweighted']+vegdata_PFTs['DEADSTEMC_unweighted'])/leaf.where(leaf>0.01).rolling(time=365*2).max()})
    frootleafratio_newparams=xarray.Dataset({'froot_leaf_ratio_unweighted':vegdata_PFTs['FROOTC_unweighted']/vegdata_PFTs['LEAFC_unweighted']})

    froot_leaf_obs = meas_root_C/meas_leaf_C.sum(level='Ecotype')


    figure('Temperature and root respiration cumulative (%s)'%dataname);clf()
    
    t2=array([tt.year + (tt.month-.5)/12 for tt in Tsoil10cm.time.data])
    plot(t2,Tsoil10cm.isel(lndgrid=1)-273.15,'b-')
    plot([0,maxyear],[0.0,0.0],'k--')
    ylabel('Soil temperature (C)')

    ax2=twinx()
    t=array([tt.year + (tt.dayofyr-1)/365 for tt in vegdata_PFTs['time'].data])
    ecotype_num=1
    startyear=10
    endyear=25
    for yr in range(startyear,endyear):
        growingseason_start=nonzero((t>=yr)&(vegdata_PFTs.sel(PFT=3,ecotype=ecotype_num)['GPP_unweighted'].values>0))[0][0]
        growingseason_start_nextyear=nonzero((t>=yr+1)&(vegdata_PFTs.sel(PFT=3,ecotype=ecotype_num)['GPP_unweighted'].values>0))[0][0]
        xx=arange(growingseason_start,growingseason_start_nextyear)
        plot_var_PFTs('MR',vegdata_PFTs.isel(time=xx),ecotype_num,cumulative=True,ax=ax2,modfactor=3600*24,units='gC/m2',minyear=yr,ls='--')
        plot_var_PFTs('FROOT_MR',vegdata_PFTs.isel(time=xx),ecotype_num,cumulative=True,ax=ax2,modfactor=3600*24,units='gC/m2',minyear=yr,ls=':')
        plot_var_PFTs('GPP',vegdata_PFTs.isel(time=xx),ecotype_num,cumulative=True,ax=ax2,modfactor=3600*24,units='gC/m2')
    xlim(startyear,endyear)

    tight_layout()
    
    barplots=True
    figure('SOM stocks (%s)'%dataname,figsize=(10,5));clf()
    t_col=array([tt.year + (tt.dayofyr-1)/365 for tt in columndata['time'].data])
    if not barplots:
        subplot(131)
        plot(t_col,columndata.TOTSOMC/1000)
        title('Soil C')
        xlabel('Time (years)')
        ylabel('{name:s} ({units:s})'.format(name=columndata.TOTSOMC.long_name,units='kgC m$^{-2}$'))
        subplot(132)
        plot(t_col,columndata.TOTSOMN/1000)
        title('Soil N')
        xlabel('Time (years)')
        ylabel('{name:s} ({units:s})'.format(name=columndata.TOTSOMN.long_name,units='kgN m$^{-2}$'))
        subplot(133)
        plot(t_col,columndata.TOTSOMP)
        title('Soil P')
        xlabel('Time (years)')
        ylabel('{name:s} ({units:s})'.format(name=columndata.TOTSOMP.long_name,units=columndata.TOTSOMP.units))
        

    else:
        x=arange(6)
        subplot(131)
        bar(x,columndata.TOTSOMC.sel(time=t_col>(t_col.max()-100)).mean(dim='time')/1000,hatch='//')
        title('Soil C')
        xticks(x,landscape_ecotypes)
        ylabel('{name:s} ({units:s})'.format(name=columndata.TOTSOMC.long_name,units='kgC m$^{-2}$'))
        subplot(132)
        bar(x,columndata.TOTSOMN.sel(time=t_col>(t_col.max()-100)).mean(dim='time')/1000,hatch='//')
        title('Soil N')
        xticks(x,landscape_ecotypes)
        ylabel('{name:s} ({units:s})'.format(name=columndata.TOTSOMN.long_name,units='kgN m$^{-2}$'))
        subplot(133)
        bar(x,columndata.TOTSOMP.sel(time=t_col>(t_col.max()-100)).mean(dim='time'),hatch='//')
        title('Soil P')
        xticks(x,landscape_ecotypes)
        ylabel('{name:s} ({units:s})'.format(name=columndata.TOTSOMP.long_name,units=columndata.TOTSOMP.units))
        

    tight_layout()


    figure('N fixation(%s)'%dataname);clf()
    x=arange(6)
    nfix_obs={'AS':1.95,'ASV':0.53}
    nfix_obs_error={'AS':0.68,'ASV':0.19}
    nfix=columndata['NFIX_TO_SMINN'].sel(time=t_col>(t_col.max()-100)).mean(dim='time')*365*24*3600
    bar(x,nfix,label='N fixation',hatch='//',width=0.4)
    bar(x,columndata['NDEP_TO_SMINN'].sel(time=t_col>(t_col.max()-100)).mean(dim='time')*365*24*3600,bottom=nfix,label='N deposition',hatch='//',width=0.4)
    bar(x+0.4,[nfix_obs.get(e,nan) for e in landscape_ecotypes], yerr=[nfix_obs_error.get(e,nan) for e in landscape_ecotypes]  ,width=0.4,label='Obs N fixation')
    title('N deposition and inputs')
    xticks(x,landscape_ecotypes)
    ylabel('N input rate (gN m$^{-2}$ year$^{-1}$)')
    legend()
    


    Koug_meas_biomass.rename(index={'NAMC':'DL','TTWBT':'ASV'},inplace=True)
    Koug_meas_chem.rename(index={'NAMC':'DL','TTWBT':'ASV'},inplace=True)

    meas_leaf_C=(Koug_meas_biomass['LeafBiomass_gperm2']*Koug_meas_chem['LeafC_percent']/100)
    meas_nonvasc_C=(Koug_meas_biomass['NonvascularBiomass_gperm2']*0.5)
    meas_rhizome_C=(Koug_meas_biomass['RhizomeBiomass_gperm2']*Koug_meas_chem['RhizomeC_percent']/100).groupby('Ecotype').sum()
    meas_root_C=Koug_meas_biomass['FineRootBiomass_gperm2'][:,'mixed']*(Koug_meas_chem['FineRootC_percent']/100)
    meas_stem_C=(Koug_meas_biomass['StemBiomass_gperm2']*Koug_meas_chem['StemC_percent']/100)
    meas_rhizome_C=(Koug_meas_biomass['RhizomeBiomass_gperm2']*Koug_meas_chem['RhizomeC_percent']/100)
    meas_NPP_C=(Koug_meas_biomass['LeafNPP_gperm2peryr']*Koug_meas_chem['LeafC_percent']/100)+\
        (Koug_meas_biomass['StemNPP_gperm2peryr']*Koug_meas_chem['StemC_percent']/100)#+\

    meas_froot_NPP =(Koug_meas_biomass['FineRootNPP_gperm2peryr']*Koug_meas_chem['FineRootC_percent']/100).groupby('Ecotype').sum()

    years=array([xx.year for xx in vegdata_PFTs.time.values])
    last_100_years=vegdata_PFTs.sel(time=((years.max()-years)<100))

    data_to_plot=last_100_years

    figure('PFT distributions',figsize=(9,5));clf()
    plot_PFT_distributions()
    tight_layout()


    def plot_mod_bar_stack(x,dat,econum,do_legend=False,**kwargs):
        handles=[]
        bottom=0.0
        for pftnum in range(len(pft_names)):
            if pftnum==10:
                val=dat.sel(PFT=10,ecotype=econum).max(dim='time')+dat.sel(PFT=11,ecotype=econum).max(dim='time')
            elif pftnum==11:
                continue
            else:
                val=dat.sel(PFT=pftnum,ecotype=econum).max(dim='time')
            if ~isnan(val):
                if pft_names[pftnum].startswith('arctic'):
                    name=pft_names[pftnum][len('arctic_'):]
                else:
                    name=pft_names[pftnum]
                if name=='dry_graminoid':
                    name='graminoid'
                handles.append(bar(x,val,bottom=bottom,facecolor=pft_colors[pftnum],label=name,**kwargs))
                bottom=bottom+val
        return handles
        
        


    def plot_obs_bar_stack(x,obsdata,ecotype_num,**kwargs):
        names=[]
        bottom=0.0
        for pft in obsdata[landscape_ecotypes[ecotype_num]].index:
            val=obsdata[(landscape_ecotypes[ecotype_num],pft)]

            if ~isnan(val):
                bar(x,val,bottom=bottom,facecolor=pft_colors[(pft_names.index(obsdata_PFT_mappings[pft]))],**kwargs)
                bottom=bottom+val

    def prettify_pft_name(name):
        if name.startswith('arctic'):
            pretty_name=name[len('arctic_'):]
        else:
            pretty_name=name
        pretty_name = ' '.join(pretty_name.split('_')).title()
        return pretty_name

    barfig=figure('Biomass comparison by ecotype (%s)'%dataname,figsize=(15,8));clf()
    for econum in range(6):
        names=[]
        x=0.0
        ax=subplot(2,3,econum+1)

        h=plot_mod_bar_stack(x+0.4,get_var_PFTs('LEAFC',data_to_plot),econum,width=0.4,hatch='//')
        plot_obs_bar_stack(x,meas_leaf_C.add(meas_nonvasc_C,fill_value=0.0) ,econum,width=0.4)
        names.append('Leaf')
        text(x+0.4,h[-1][0].get_y()+h[-1][0].get_height()+50,'Modeled value',rotation=90,va='bottom',ha='center')
        text(x,h[-1][0].get_y()+h[-1][0].get_height()+50,'Measured value',rotation=90,va='bottom',ha='center')
        x+=1

        plot_mod_bar_stack(x+0.4,get_var_PFTs(['LIVESTEMC','DEADSTEMC'],data_to_plot),econum,width=0.4,hatch='//')
        plot_obs_bar_stack(x,meas_stem_C,econum,width=0.4)
        names.append('Stem')
        x+=1

        plot_mod_bar_stack(x+0.4,get_var_PFTs('FROOTC',data_to_plot),econum,width=0.4,hatch='//')
        bar(x,meas_root_C[landscape_ecotypes[econum]],width=0.4,facecolor=[0.9,0.9,0.9],edgecolor='k')
        names.append('Fine root')
        x+=1

        plot_mod_bar_stack(x+0.4,get_var_PFTs(['LIVECROOTC','DEADCROOTC'],data_to_plot),econum,width=0.4,hatch='//')
        # text(x+0.4,100,'? ? ? ? ? ?',fontsize='large',rotation=90,va='bottom',ha='center')
        plot_obs_bar_stack(x,meas_rhizome_C,econum,width=0.4)
        names.append('Rhizome/\nCoarse root')
        x+=1
        
        plot_mod_bar_stack(x+0.4,get_var_PFTs(['STORVEGC'],data_to_plot),econum,width=0.4,hatch='//')
        names.append('Storage')

        ylabel('Biomass (gC m$^{-2}$)')
        xticks(arange(len(names))+0.3,names)
        title('%s biomass'%ecotype_names_list[econum])
        # ylim(0,3000)
        # ax.set_xlim(right=x+0.6)


    handles=[]
    for pftnum in range(1,len(pft_names)-1):
        name = prettify_pft_name(pft_names[pftnum])
        if name == 'Dry Graminoid':
            name = 'Graminoid'
        handles.append(Rectangle([0,0],0,0,facecolor=pft_colors[pftnum],label=name ))
    l=barfig.axes[1].legend(handles=handles,fontsize='small',ncol=2)
    l.set_draggable(True)

    tight_layout()



    barfig_PFT=figure('Biomass comparison by PFT (%s)'%dataname,figsize=(15,8));clf()

    w=0.8/6/2
    for pftnum in range(1,len(pft_names)-1):
        subplot(2,5,pftnum)
        bottom=zeros(6)
        bottom_obs=zeros(6)
        handles=[]
        
        vals=get_var_PFTs('FROOTC',data_to_plot).max(dim='time').sel(PFT=pftnum)
        h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='brown',label='Fine root',hatch='//')
        bottom=bottom+vals.fillna(0)
        handles.append(h)
        
        vals=get_var_PFTs(['LIVECROOTC','DEADCROOTC'],data_to_plot).mean(dim='time').sel(PFT=pftnum)
        h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='orange',label='Coarse root',hatch='//')
        bottom=bottom+vals.fillna(0)
        handles.append(h)
        
        vals=get_var_PFTs('STORVEGC',data_to_plot).mean(dim='time').sel(PFT=pftnum)
        h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='yellow',label='Storage',hatch='//')
        bottom=bottom+vals.fillna(0)
        handles.append(h)
        # bar(arange(6)/6+w,meas_rhizome_C)
        
        vals=get_var_PFTs(['LIVESTEMC','DEADSTEMC'],data_to_plot).mean(dim='time').sel(PFT=pftnum)
        h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='blue',label='Stem',hatch='//')
        bottom=bottom+vals.fillna(0)
        handles.append(h)
        
        vals=get_var_PFTs('LEAFC',data_to_plot).max(dim='time').sel(PFT=pftnum)
        h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='green',label='Leaf',hatch='//')
        bottom=bottom+vals.fillna(0)
        handles.append(h)  

        ylabel('Biomass (gC m$^{-2}$)')
        # xticks(arange(1,len(pft_names)-1)+0.3,[prettify_pft_name(name) for name in pft_names])
        xticks(arange(6)/6,landscape_ecotypes)
        title(prettify_pft_name(pft_names[pftnum]))

    for pftname in obsdata_PFT_mappings.keys():
        if pftname in ['other','mixed']:
             continue
        ax=barfig_PFT.axes[pft_names.index(obsdata_PFT_mappings[pftname])-1]
        for econum in range(len(landscape_ecotypes)):
            bottom=0.0
            val=meas_rhizome_C[landscape_ecotypes[econum]].fillna(0.0).get(pftname,0.0)
            ax.bar(econum/6+w,val,bottom=bottom,facecolor='orange',width=w)
            bottom=bottom+val
            val=meas_stem_C[landscape_ecotypes[econum]].fillna(0.0).get(pftname,0.0)
            ax.bar(econum/6+w,val,bottom=bottom,facecolor='blue',width=w)
            bottom=bottom+val
            val=meas_leaf_C[landscape_ecotypes[econum]].fillna(0.0).get(pftname,0.0)
            ax.bar(econum/6+w,val,bottom=bottom,facecolor='green',width=w)
            bottom=bottom+val
            val=meas_nonvasc_C[landscape_ecotypes[econum]].fillna(0.0).get(pftname,0.0)
            h_nonvasc=ax.bar(econum/6+w,val,bottom=bottom,facecolor='cyan',width=w,label='Non-vascular')
            bottom=bottom+val

    handles.append(h_nonvasc)
    l=barfig_PFT.axes[3].legend(fontsize='small',ncol=2,handles=handles)
    l.set_draggable(True)

    tight_layout()



    barfig_PFT_perarea=figure('Biomass comparison per area of each PFT (%s)'%dataname,figsize=(15,8));clf()

    w=0.8/6/3
    for pftnum in range(1,len(pft_names)-1):
        subplot(2,5,pftnum)
        bottom=zeros(6)
        bottom_obs=zeros(6)
        handles=[]
        
        vals=get_var_PFTs('FROOTC',data_to_plot,weight_area=False).max(dim='time').sel(PFT=pftnum)
        h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='brown',label='Fine root',hatch='//')
        bottom=bottom+vals.fillna(0)
        handles.append(h)
        
        vals=get_var_PFTs(['LIVECROOTC','DEADCROOTC'],data_to_plot,weight_area=False).mean(dim='time').sel(PFT=pftnum)
        h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='orange',label='Coarse root',hatch='//')
        bottom=bottom+vals.fillna(0)
        handles.append(h)
        
        vals=get_var_PFTs('STORVEGC',data_to_plot,weight_area=False).mean(dim='time').sel(PFT=pftnum)
        h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='yellow',label='Storage',hatch='//')
        bottom=bottom+vals.fillna(0)
        handles.append(h)
        # bar(arange(6)/6+w,meas_rhizome_C)
        
        vals=get_var_PFTs(['LIVESTEMC','DEADSTEMC'],data_to_plot,weight_area=False).mean(dim='time').sel(PFT=pftnum)
        h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='blue',label='Stem',hatch='//')
        bottom=bottom+vals.fillna(0)
        handles.append(h)
        
        vals=get_var_PFTs('LEAFC',data_to_plot,weight_area=False).max(dim='time').sel(PFT=pftnum)
        h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='green',label='Leaf',hatch='//')
        bottom=bottom+vals.fillna(0)
        handles.append(h)  

        ylabel('Biomass (gC m$^{-2}$)')
        # xticks(arange(1,len(pft_names)-1)+0.3,[prettify_pft_name(name) for name in pft_names])
        xticks(arange(6)/6,landscape_ecotypes)
        title(prettify_pft_name(pft_names[pftnum]))

    for pftname in obsdata_PFT_mappings.keys():
        if pftname in ['other','mixed']:
             continue
        ax=barfig_PFT_perarea.axes[pft_names.index(obsdata_PFT_mappings[pftname])-1]
        for econum in range(len(landscape_ecotypes)):
            bottom=0.0
            areafrac=PFT_percents[landscape_ecotypes[econum]][obsdata_PFT_mappings[pftname]]/100
            if areafrac==0.0:
                invareafrac=0.0
            else:
                invareafrac=1.0/areafrac
                
            if pftname not in ['moss','lichen']:
                val=meas_root_C[landscape_ecotypes[econum]]
                h_obsroot=ax.bar(econum/6+w*2,val,bottom=bottom,facecolor='brown',width=w,hatch='.',label='Pooled obs fine roots')
                # bottom=bottom+val
            val=meas_rhizome_C[landscape_ecotypes[econum]].fillna(0.0).get(pftname,0.0)*invareafrac
            h_rhizome=ax.bar(econum/6+w,val,bottom=bottom,facecolor='orange',width=w,label='Obs rhizome')
            bottom=bottom+val
            val=meas_stem_C[landscape_ecotypes[econum]].fillna(0.0).get(pftname,0.0)*invareafrac
            ax.bar(econum/6+w,val,bottom=bottom,facecolor='blue',width=w)
            bottom=bottom+val
            val=meas_leaf_C[landscape_ecotypes[econum]].fillna(0.0).get(pftname,0.0)*invareafrac
            ax.bar(econum/6+w,val,bottom=bottom,facecolor='green',width=w)
            bottom=bottom+val
            val=meas_nonvasc_C[landscape_ecotypes[econum]].fillna(0.0).get(pftname,0.0)*invareafrac
            h_nonvasc=ax.bar(econum/6+w,val,bottom=bottom,facecolor='cyan',width=w,label='Obs non-vascular')
            bottom=bottom+val
            

    handles.append(h_nonvasc)
    handles.append(h_obsroot)
    handles.append(h_rhizome)
    l=barfig_PFT_perarea.axes[3].legend(fontsize='small',ncol=1,handles=handles)
    l.set_draggable(True)

    tight_layout()


    figure('Height (%s)'%dataname,figsize=(13,7));clf()


    for econum in range(6):
        names=[]
        x=0.0
        ax=subplot(2,3,econum+1)
        
        heights=get_var_PFTs('HTOP',data_to_plot,weight_area=False).sel(ecotype=econum).max(dim='time')   
        good=~isnan(heights.values)&(heights.values>0)
        good[11]=False # wet graminoids that we are ignoring
        x=arange(len(heights[good]))
        pftnames=[prettify_pft_name(name) for name in array(pft_names)[good]]
        pftnames[-1]='Graminoid'
        for num in range(len(x)):
            bar(x[num],heights[good][num],facecolor=array(pft_colors)[good][num],hatch='//',width=0.4)
            if pftnames[num] in ['Deciduous Shrub Tall','Deciduous Shrub Alder']:
                bar(x[num]+0.4,obs_heights['tall_shrub_height_mean'],facecolor=array(pft_colors)[good][num],width=0.4)
            elif 'Low' in pftnames[num]:
                bar(x[num]+0.4,obs_heights['low_shrub_height_mean'],facecolor=array(pft_colors)[good][num],width=0.4)
            elif 'Dwarf' in pftnames[num]:
                bar(x[num]+0.4,obs_heights['dwarf_shrub_height_mean'],facecolor=array(pft_colors)[good][num],width=0.4)
            elif pftnames[num] is 'Forb':
                bar(x[num]+0.4,obs_heights['forb_height_mean'],facecolor=array(pft_colors)[good][num],width=0.4)
        xticks(x,pftnames,rotation=45,ha='right')
        ylabel('Height (m)')
        title(ecotype_names_list[econum])
    tight_layout()

    figure('Aboveground-belowground (%s)'%dataname,figsize=(10.8,4.8));clf()
    meas_leaf_C=(Koug_meas_biomass['LeafBiomass_gperm2']*Koug_meas_chem['LeafC_percent']/100)
    meas_root_C=(Koug_meas_biomass['FineRootBiomass_gperm2'][:,'mixed']*(Koug_meas_chem['FineRootC_percent']/100))[:,'mixed']
    meas_rhizome_C=(Koug_meas_biomass['RhizomeBiomass_gperm2']*Koug_meas_chem['RhizomeC_percent']/100)
    meas_rhizome_NPP=(Koug_meas_biomass['RhizomeNPP_gperm2peryr']*Koug_meas_chem['RhizomeC_percent']/100)
    meas_leaf_NPP   =(Koug_meas_biomass['LeafNPP_gperm2peryr']*Koug_meas_chem['LeafC_percent']/100)
    meas_stem_NPP   =(Koug_meas_biomass['StemNPP_gperm2peryr']*Koug_meas_chem['StemC_percent']/100)
    meas_root_NPP   =(Koug_meas_biomass['FineRootNPP_gperm2peryr']*Koug_meas_chem['FineRootC_percent']/100)[:,'mixed']

    aboveground_C = (meas_leaf_C+meas_stem_C).sum(level='Ecotype')
    aboveground_NPP = (meas_leaf_NPP+meas_stem_NPP).sum(level='Ecotype')
    belowground_C = meas_rhizome_C.sum(level='Ecotype') + meas_root_C
    belowground_NPP = meas_rhizome_NPP.sum(level='Ecotype') + meas_root_NPP

    mod_aboveground_NPP = get_var_PFTs('AGNPP',data_to_plot).mean(dim='time').sum(dim='PFT')*3600*24*365
    mod_belowground_NPP = get_var_PFTs('BGNPP',data_to_plot).mean(dim='time').sum(dim='PFT')*3600*24*365
    mod_aboveground_C  = get_var_PFTs(['LEAFC','LIVESTEMC','DEADSTEMC'],data_to_plot).max(dim='time').sum(dim='PFT')
    mod_belowground_C  = get_var_PFTs(['FROOTC','LIVECROOTC','DEADCROOTC'],data_to_plot).max(dim='time').sum(dim='PFT')

    bar(arange(len(landscape_ecotypes)),100*(belowground_C/(aboveground_C+belowground_C))[landscape_ecotypes],width=0.2,label='Biomass')
    bar(arange(len(landscape_ecotypes))+0.2,100*(belowground_NPP/(aboveground_NPP+belowground_NPP))[landscape_ecotypes],width=0.2,label='NPP')
    bar(arange(len(landscape_ecotypes))+0.4,100*(mod_belowground_C/(mod_aboveground_C+mod_belowground_C)),width=0.2,label='Mod biomass',facecolor='C0',hatch='//')
    bar(arange(len(landscape_ecotypes))+0.6,100*(mod_belowground_NPP/(mod_aboveground_NPP+mod_belowground_NPP)),width=0.2,label='Mod NPP',facecolor='C1',hatch='//')
    names=[ecotype_names[name] for name in landscape_ecotypes]
    names[0]='Dryas-lichen\ndwarf shrub tundra'
    xticks(arange(len(landscape_ecotypes))+0.3, names,rotation=30,ha='right',fontsize='x-large')
    title('Measured belowground percentage',fontsize='x-large')
    ylabel('Belowground fraction (%)',fontsize='x-large')
    ylim(0,100)
    legend(fontsize='large')
    tight_layout()


    figure('Non-vascular (%s)'%dataname);clf()
    meas_lichen_BM=Koug_meas_biomass['NonvascularBiomass_gperm2'][:,'lichen']
    meas_moss_BM=Koug_meas_biomass['NonvascularBiomass_gperm2'][:,'moss']
    meas_other_nonvasc_BM=Koug_meas_biomass['NonvascularBiomass_gperm2'][:,'other']
    meas_vasc_AG_BM=(Koug_meas_biomass['LeafBiomass_gperm2']+Koug_meas_biomass['StemBiomass_gperm2']).sum(level='Ecotype')

    bottom=meas_lichen_BM.reindex(landscape_ecotypes,fill_value=0.0)*0

    y=meas_lichen_BM.reindex(landscape_ecotypes,fill_value=0.0)
    bar(arange(len(landscape_ecotypes)),y,width=0.3,label='Lichen',color=pft_colors[pft_names.index('arctic_lichen')],bottom=bottom)
    bottom=bottom+y
    y=meas_moss_BM.reindex(landscape_ecotypes,fill_value=0.0)
    bar(arange(len(landscape_ecotypes)),y,width=0.3,label='Moss',color=pft_colors[pft_names.index('arctic_bryophyte')],bottom=bottom)
    bottom=bottom+y
    y=meas_other_nonvasc_BM.reindex(landscape_ecotypes,fill_value=0.0)
    bar(arange(len(landscape_ecotypes)),y,width=0.3,label='Other nonvasc.',color='gray',bottom=bottom)
    bottom=bottom+y
    y=meas_vasc_AG_BM.reindex(landscape_ecotypes,fill_value=0.0)
    bar(arange(len(landscape_ecotypes))+0.3,y,width=0.3,label='Vascular',color='brown')
    # bottom=bottom.add(meas_vasc_AG_BM,fill_value=0.0)

    mod_moss=get_var_PFTs('TOTVEGC',data_to_plot,0).max(dim='time').sel(PFT=pft_names.index('arctic_bryophyte'))
    mod_lichen=get_var_PFTs('TOTVEGC',data_to_plot,0).max(dim='time').sel(PFT=pft_names.index('arctic_lichen'))
    bar(arange(len(landscape_ecotypes))+0.3+0.3,mod_lichen,label='Mod moss',color=[.3,1,.3],width=0.3)
    bar(arange(len(landscape_ecotypes))+0.3+0.3,mod_moss,bottom=mod_lichen,label='Mod lichen',color='orange',width=0.3)
    title('Aboveground vascular and non-vascular biomass')
    ylabel('Total biomass (g DW m$^{-2}$)')
    l=legend()
    l.set_draggable(True)
    xticks(arange(len(landscape_ecotypes))+0.3,[ecotype_names[name] for name in landscape_ecotypes] ,rotation=45,ha='right')
    tight_layout()
    
    figure('VCMAX (%s)'%dataname);clf()
    for n in range(6):
        ax=subplot(3,2,n+1)
        plot_var_PFTs('VCMAX25TOP',vegdata_PFTs,weight_area=False,ecotype_num=n,maxyear=maxyear,minyear=maxyear-1)
        title(ecotype_names_list[n])
        ylabel('Vcmax(25C)\n(umol m$^{-2}$ s$^{-1}$)')
    legend(fontsize='small',ncol=3)
    tight_layout()


    import cartopy.crs as ccrs
    import cartopy

    figure('Site location');clf()
    ax=subplot(111,projection=ccrs.Miller(central_longitude=-155))
    ax.set_extent((-175,-138,50,75)) 
    ax.coastlines(resolution='50m')
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'))
    ax.scatter(-164.82,65.16,marker='*',s=50,c='g',transform=ccrs.Geodetic())



    show()
